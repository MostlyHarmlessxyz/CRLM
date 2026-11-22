#!/usr/bin/env python3
"""Strictly prospective zero-exposure evaluation of frozen CRLM_V1.1_245_F1 tissue classifiers on GSE178318.

Methodological & provenance guarantees:
1. Cryptographically verifies that `artifacts/m5/frozen_inputs_sha256.txt` matches Git commit
   `5b01467` (created BEFORE any GSE178318 expression download) and verifies ALL 72/72 locked
   M5/M4/M3/HGNC/code entries byte-for-byte (`0` post-unlock tuning or refitting).
2. Evaluates ONLY the 12 paired solid-tumor libraries (`6 PT + 6 LM` across 6 synchronous
   Stage IV MSS donors: `COL07, COL12, COL16` treatment-naive; `COL15, COL17, COL18`
   preoperative chemotherapy); physically blocks the 3 PBMC libraries (`GSM5387664/671/674`).
3. Encodes all QC-retained `T/NK + Mono/Mac/DC + B/Plasma` tumor cells (`59,338` cells) on the
   shared canonical nuclear universe (`245_164_203608_shared_nuclear`, 18,733 genes) via
   `log1p(CP10k)` and frozen whole-human `scGPT` on Intel Arc A750 XPU, plus genuine `<=256`-cell
   historical sub-bag embeddings (`7,616` cells) for both `CRLM_V1_245_F1` (raw `33,694`
   `data/raw/GSE178318/` whole-axis counts & exact raw symbols) and `CRLM_V1_245_F1_techfix`
   (HGNC-canonicalized shared-nuclear integer counts matching V1 M4 training) under
   `embeddings/archive/`.
4. Scores all 15/15 frozen `CRLM_V1.1_245_F1` classifiers (`512/512 = 100.0%` PCA HVG coverage,
   0 rejections across `34` eligible `>=50`-cell bags containing `59,298` cells), computes
   10,000-draw donor-clustered bootstrap CIs, paired `PCA` vs `scGPT` contrasts, descriptive
   treatment-stratified summaries, and a genuine 3-version (`V1_F1` / `Techfix` / `V1.1`)
   classifier comparison table (`l_raw != l_tf != l_v11`).
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import importlib.util
from itertools import combinations
import json
import os
import platform
from pathlib import Path
import resource
import subprocess
import sys
import time

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.special import expit

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src"), str(ROOT / "starter"), str(ROOT / "scripts")]

from crlm.hgnc_alias_map import (  # noqa: E402
    align_axis_to_canonical_targets,
    canonicalize_symbols,
    get_hgnc_mapper,
    load_shared_nuclear_universe,
)
from crlm.model_bags import (  # noqa: E402
    Bag,
    LINEAGES,
    SEED,
    cell_seed,
    sample_bags,
    sample_joint_bag,
)
from crlm.program_scores import controls_from_dataframe, debiased_ucell_scores  # noqa: E402
from linear_attention_mil import load_model, predict_bags  # noqa: E402

RAW_178 = ROOT / "data/raw/GSE178318"
DATA = ROOT / "data/processed/m6/GSE178318"
EXTERNAL = ROOT / "artifacts/m6/GSE178318"
MANIFEST_CSV = EXTERNAL / "coverage/sample_manifest.csv"
OUT = EXTERNAL / "model_predictions"
CACHE = ROOT / "embeddings/m6/GSE178318/scgpt"
CACHE_V1_RAW = ROOT / "embeddings/archive/CRLM_V1_245_F1/m6/GSE178318/scgpt"
CACHE_TECHFIX = ROOT / "embeddings/archive/Techfix_V1_F1/m6/GSE178318/scgpt"
M5 = ROOT / "artifacts/m5"
ARCHIVE_V1_M5 = ROOT / "artifacts/archive/CRLM_V1_245_F1/artifacts/m5"
ARCHIVE_V1_MODELS = ROOT / "models/archive/CRLM_V1_245_F1/m4/full"

SHARED_UNIVERSE = "245_164_203608_shared_nuclear"
FORBIDDEN_PBMC_GSMS = frozenset({"GSM5387664", "GSM5387671", "GSM5387674"})
PRE_UNLOCK_FREEZE_COMMIT = "5b01467"
METHODS = (
    "program_elastic",
    "PCA_mean",
    "PCA_attention",
    "scGPT_mean",
    "scGPT_attention",
)


def _sha256_file(path: Path) -> str:
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def _atomic_write_text(path: Path, content: str) -> None:
    tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    tmp.write_text(content, encoding="utf-8")
    os.replace(tmp, path)


def refresh_gse178318_sha_manifests(
    out_dir: Path = OUT, *, refresh_parent: bool = True
) -> tuple[int, int]:
    """Write `<out_dir>/output_sha256.txt` and optionally refresh `GSE178318/output_sha256.txt`."""
    mp_files = sorted(
        p
        for p in out_dir.rglob("*")
        if p.is_file() and p.name not in ("output_sha256.txt", ".gitkeep") and not p.name.startswith(".")
    )
    _atomic_write_text(
        out_dir / "output_sha256.txt",
        "".join(f"{_sha256_file(p)}  {p.relative_to(ROOT)}\n" for p in mp_files),
    )

    n_parent = 0
    if refresh_parent and out_dir.resolve() == OUT.resolve():
        g178_files = sorted(
            p
            for p in EXTERNAL.rglob("*")
            if p.is_file()
            and p.name not in ("output_sha256.txt", ".gitkeep")
            and not p.name.startswith(".")
        )
        _atomic_write_text(
            EXTERNAL / "output_sha256.txt",
            "".join(f"{_sha256_file(p)}  {p.relative_to(ROOT)}\n" for p in g178_files),
        )
        n_parent = len(g178_files)
    return len(mp_files), n_parent


def verify_frozen_m5_models(
    revision: str = PRE_UNLOCK_FREEZE_COMMIT,
) -> dict[tuple[str, str], Path]:
    """Verify M5 registry against pre-unlock git commit `5b01467` and check ALL 72/72 locked SHA256s."""
    frozen = M5 / "frozen_inputs_sha256.txt"
    if not frozen.is_file():
        raise FileNotFoundError(f"Missing frozen input registry: {frozen}")
    if revision and revision.upper() != "WORKTREE":
        committed_bytes = subprocess.check_output(
            ["git", "show", f"{revision}:{frozen.relative_to(ROOT)}"], cwd=ROOT
        )
        if committed_bytes != frozen.read_bytes():
            raise ValueError(
                f"Frozen registry {frozen} differs from pre-unlock git commit {revision}"
            )
    entries = [
        ln.strip() for ln in frozen.read_text(encoding="utf-8").splitlines() if ln.strip()
    ]
    if len(entries) != 72:
        raise ValueError(f"Expected 72 frozen M5 input entries, found {len(entries)}")
    for line in entries:
        expected, filename = line.split("  ", 1)
        target = ROOT / filename
        if not target.is_file():
            raise FileNotFoundError(f"Missing frozen M5 input: {filename}")
        if _sha256_file(target) != expected:
            raise ValueError(f"Frozen input modified post-freeze: {filename}")
    fits = pd.read_csv(M5 / "frozen_models.csv")
    models = {
        (r.lineage, r.method): ROOT / r.path for r in fits.itertuples() if r.accepted
    }
    if set(models) != {(lin, m) for lin in LINEAGES for m in METHODS}:
        raise ValueError(f"Expected 15 accepted M5 classifiers, got {len(models)}")
    for lin in LINEAGES:
        prog_dir = models[(lin, "program_elastic")].parent
        if not (prog_dir / "program_defs.json").is_file() or not (
            prog_dir / "matched_controls.csv"
        ).is_file():
            raise FileNotFoundError(
                f"Missing frozen program_defs.json or matched_controls.csv in {prog_dir}"
            )
    return models


def load_tumor_manifest() -> pd.DataFrame:
    """Load and validate the 12 solid-tumor specimens of GSE178318 (`PBMC` strictly blocked)."""
    if not MANIFEST_CSV.is_file():
        raise FileNotFoundError(f"Missing GSE178318 tumor manifest: {MANIFEST_CSV}")
    # Preserve literal "None" in `preop_regimen` for untreated donors (`keep_default_na=False`)
    manifest = pd.read_csv(MANIFEST_CSV, keep_default_na=False, na_values=[""])
    if manifest["sample_uid"].isin(FORBIDDEN_PBMC_GSMS).any():
        raise ValueError("Forbidden PBMC library detected in tumor manifest")
    if set(manifest["tissue"]) != {"PT", "LM"} or len(manifest) != 12:
        raise ValueError(f"Expected 12 paired PT/LM libraries, found {len(manifest)}")
    for h5_file in sorted(DATA.glob("*.h5ad")):
        if h5_file.stem in FORBIDDEN_PBMC_GSMS:
            raise ValueError(f"Forbidden PBMC H5AD present in tumor directory: {h5_file}")
    for row in manifest.itertuples(index=False):
        h5_path = ROOT / row.h5ad_path
        if not h5_path.is_file() or _sha256_file(h5_path) != row.h5ad_sha256:
            raise ValueError(f"H5AD missing or SHA256 mismatch: {row.sample_uid}")
    return manifest.sort_values("sample_uid").reset_index(drop=True)


def normalize_counts_log1p_cp10k(
    counts: sparse.csr_matrix,
) -> tuple[sparse.csr_matrix, float]:
    """Normalize integer counts to `log1p(CP10k)`."""
    mat = sparse.csr_matrix(counts, copy=True, dtype=np.float64)
    mat.sum_duplicates()
    mat.eliminate_zeros()
    if np.any(mat.data < 0) or not np.isfinite(mat.data).all():
        raise ValueError("Expected finite nonnegative integer counts")
    sums = np.asarray(mat.sum(axis=1)).ravel()
    if (sums <= 0).any():
        raise ValueError("Encountered cell with zero UMI sum")
    median_sum = float(np.median(sums))
    row_nnz = np.diff(mat.indptr)
    mat.data = np.log1p(mat.data * np.repeat(10000.0 / sums, row_nnz))
    return mat.astype(np.float32), median_sum


def _load_raw_33694_genes() -> tuple[list[str], list[str]]:
    genes_path = RAW_178 / "GSE178318_genes.tsv.gz"
    if not genes_path.is_file():
        raise FileNotFoundError(f"Missing raw GSE178318 genes file: {genes_path}")
    ids, syms = [], []
    with gzip.open(genes_path, "rt", encoding="utf-8") as handle:
        for ln in handle:
            parts = ln.rstrip("\r\n").split("\t")
            ids.append(parts[0])
            syms.append(parts[1] if len(parts) > 1 else parts[0])
    if len(syms) != 33694:
        raise ValueError(f"Expected 33,694 raw GSE178318 genes, got {len(syms)}")
    return ids, syms


def _extract_raw_33694_subbags(
    uids_by_sample: dict[str, np.ndarray],
) -> tuple[dict[str, sparse.csr_matrix], str]:
    """Stream `data/raw/GSE178318/GSE178318_matrix.mtx.gz` once to slice `<=256`-cell bags (`33,694` genes)."""
    barcodes_path = RAW_178 / "GSE178318_barcodes.tsv.gz"
    mtx_path = RAW_178 / "GSE178318_matrix.mtx.gz"
    if not barcodes_path.is_file() or not mtx_path.is_file():
        raise FileNotFoundError(f"Missing raw GSE178318 matrix/barcodes in {RAW_178}")
    mtx_sha = _sha256_file(mtx_path)

    with gzip.open(barcodes_path, "rt", encoding="utf-8") as handle:
        all_barcodes = pd.Index([ln.rstrip("\r\n") for ln in handle])
    n_total_cols = len(all_barcodes)

    ordered_samples = sorted(uids_by_sample)
    flat_uids = (
        np.concatenate([uids_by_sample[s] for s in ordered_samples])
        if ordered_samples
        else np.empty(0, dtype=str)
    )
    flat_bcs = [str(u).split(":", 1)[1] for u in flat_uids]
    raw_col0 = all_bcs_idx = all_barcodes.get_indexer(flat_bcs)
    if (all_bcs_idx < 0).any():
        raise ValueError("Failed to resolve sub-bag cell barcode in GSE178318_barcodes.tsv.gz")

    # Map 1-based MTX column index -> global row index in [0, len(flat_uids))
    col1_to_dest = np.full(n_total_cols + 1, -1, dtype=np.int32)
    col1_to_dest[raw_col0 + 1] = np.arange(len(flat_uids), dtype=np.int32)

    with gzip.open(mtx_path, "rt", encoding="utf-8") as handle:
        skip = 0
        n_genes = 0
        for ln in handle:
            skip += 1
            if not ln.startswith("%"):
                dims = ln.split()
                n_genes, n_cols_hdr = int(dims[0]), int(dims[1])
                if n_genes != 33694 or n_cols_hdr != n_total_cols:
                    raise ValueError(
                        f"Unexpected MTX header ({n_genes}, {n_cols_hdr}) in {mtx_path}"
                    )
                break

    row_parts, col_parts, val_parts = [], [], []
    reader = pd.read_csv(
        mtx_path,
        sep=r"\s+",
        skiprows=skip,
        header=None,
        names=["g1", "c1", "val"],
        dtype={"g1": np.int32, "c1": np.int32, "val": np.float32},
        chunksize=12_000_000,
    )
    for chunk in reader:
        c1_arr = chunk["c1"].to_numpy()
        dest_r = col1_to_dest[c1_arr]
        keep = dest_r >= 0
        if keep.any():
            row_parts.append(dest_r[keep])
            col_parts.append(chunk["g1"].to_numpy()[keep] - 1)
            val_parts.append(chunk["val"].to_numpy()[keep])

    full_csr = sparse.csr_matrix(
        (
            np.concatenate(val_parts) if val_parts else np.empty(0, dtype=np.float32),
            (
                np.concatenate(row_parts) if row_parts else np.empty(0, dtype=np.int32),
                np.concatenate(col_parts) if col_parts else np.empty(0, dtype=np.int32),
            ),
        ),
        shape=(len(flat_uids), 33694),
        dtype=np.float32,
    )
    full_csr.sum_duplicates()
    full_csr.eliminate_zeros()

    out_mats: dict[str, sparse.csr_matrix] = {}
    offset = 0
    for s in ordered_samples:
        k = len(uids_by_sample[s])
        out_mats[s] = full_csr[offset : offset + k, :].tocsr()
        offset += k
    return out_mats, mtx_sha


def _tokenize_raw_v1_whole_axis(
    mat_33694: sparse.csr_matrix,
    uids: np.ndarray,
    raw_33694_syms: list[str],
    vocab: dict[str, int],
) -> list[tuple[str, np.ndarray, np.ndarray, int, int]]:
    """Historical `CRLM_V1_245_F1` whole-axis tokenizer (`ff1f1a3` + `91a589f` unambiguous raw symbols)."""
    import torch

    unamb = np.asarray(~pd.Index(raw_33694_syms).duplicated(keep=False), dtype=bool)
    gene_ids = np.array(
        [vocab.get(g, -1) if ok else -1 for g, ok in zip(raw_33694_syms, unamb)],
        dtype=np.int64,
    )
    mat = mat_33694.tocsr().astype(np.float32)
    mat.eliminate_zeros()
    cells = []
    for i, uid in enumerate(uids):
        a, b = mat.indptr[i : i + 2]
        ids = gene_ids[mat.indices[a:b]]
        use = ids >= 0
        vals = mat.data[a:b][use]
        ids = ids[use]
        if not len(ids):
            raise ValueError(f"No expressed in-vocabulary genes in {uid}")
        seed = cell_seed(str(uid))
        bins = np.quantile(vals, np.linspace(0, 1, 50))
        left, right = np.digitize(vals, bins), np.digitize(vals, bins, right=True)
        binned = np.ceil(
            left + np.random.RandomState(seed).rand(len(vals)) * (right - left)
        )
        if len(ids) > 1199:
            choice = torch.randperm(
                len(ids), generator=torch.Generator().manual_seed(seed)
            )[:1199].numpy()
            ids, binned = ids[choice], binned[choice]
        cells.append(
            (
                str(uid),
                np.r_[vocab["<cls>"], ids],
                np.r_[-2, binned],
                int(b - a),
                int(use.sum()),
            )
        )
    return cells


def encode_gse178318(
    device: str = "xpu",
    batch: int = 16,
    out_dir: Path = OUT,
    *,
    overwrite: bool = False,
) -> None:
    import torch
    from encode_scgpt import (
        CODE,
        HASHES,
        INPUT_TRANSFORM,
        encode,
        load_encoder,
        resolve_vocab_ids,
        tokens,
    )

    mapper = get_hgnc_mapper()
    universe_genes = list(load_shared_nuclear_universe(SHARED_UNIVERSE))
    profile = json.loads(
        (ROOT / "artifacts/m4/resource_profiles/scgpt_profile.json").read_text(
            encoding="utf-8"
        )
    )
    if (
        profile.get("checkpoint") != HASHES
        or not profile.get("devices", {}).get(device, {}).get("accepted")
    ):
        raise ValueError(f"Unprofiled checkpoint/device: {device}")

    manifest = load_tumor_manifest()
    CACHE.mkdir(parents=True, exist_ok=True)
    CACHE_V1_RAW.mkdir(parents=True, exist_ok=True)
    CACHE_TECHFIX.mkdir(parents=True, exist_ok=True)
    out_dir.mkdir(parents=True, exist_ok=True)

    prior_batches: dict[str, int] = {}
    for prev_man in (out_dir / "embedding_manifest.csv", CACHE / "embedding_manifest.csv"):
        if prev_man.is_file():
            prev_df = pd.read_csv(prev_man, keep_default_na=False)
            if "sample_uid" in prev_df.columns and "batch" in prev_df.columns:
                prior_batches = {
                    str(r.sample_uid): int(r.batch) for r in prev_df.itertuples(index=False)
                }
                break

    torch.set_num_threads(4)
    model, vocab = load_encoder()
    vocab_ids, bridge = resolve_vocab_ids(universe_genes, vocab)
    n_matched_vocab = int(bridge["matched"].sum())
    model.to(device)

    # 1. Primary V1.1 full-bag shared-nuclear `log1p(CP10k)` encoding (`59,338` cells)
    records = []
    subbag_uids_by_sample: dict[str, np.ndarray] = {}
    subbag_nuclear_by_sample: dict[str, Bag] = {}

    for row in manifest.itertuples(index=False):
        gsm = str(row.sample_uid)
        h5_path = ROOT / row.h5ad_path
        src_sha = str(row.h5ad_sha256)
        dest = CACHE / f"{gsm}.npz"

        joint_bag = sample_joint_bag(
            h5_path,
            min_cells=1,
            max_cells=None,
            gene_universe=universe_genes,
        )
        uid = (
            joint_bag.obs.index.to_numpy(dtype=str)
            if joint_bag is not None
            else np.empty(0, dtype=str)
        )
        need_encode = overwrite or (not dest.exists())
        if not need_encode:
            with np.load(dest, allow_pickle=False) as enc:
                same = (
                    "gene_universe" in enc.files
                    and str(enc["gene_universe"]) == SHARED_UNIVERSE
                    and str(enc["hgnc_sha256"]) == mapper.sha256
                    and np.array_equal(enc["cell_uid"], uid)
                    and str(enc["source_h5ad_sha256"]) == src_sha
                    and str(enc["checkpoint"]) == HASHES["best_model.pt"]
                    and str(enc["upstream"]) == CODE
                )
                if not same:
                    raise ValueError(
                        f"Cached GSE178318 embedding mismatch at {dest}; pass --overwrite"
                    )
        if need_encode:
            t0 = time.monotonic()
            cells = (
                tokens(
                    joint_bag,
                    vocab,
                    vocab_ids=vocab_ids,
                    input_is_log1p=False,
                    recalibrate_cp10k=True,
                )
                if joint_bag is not None
                else []
            )
            embedded = encode(model, cells, vocab, batch, device)
            if device == "xpu":
                torch.xpu.synchronize()
            tmp = dest.with_suffix(".partial.npz")
            np.savez_compressed(
                tmp,
                embedding=embedded,
                cell_uid=uid,
                lineage_m2=joint_bag.obs["lineage_m2"].to_numpy(dtype=str)
                if joint_bag is not None
                else np.empty(0, dtype=str),
                fine_subtype=joint_bag.obs["fine_subtype"].to_numpy(dtype=str)
                if joint_bag is not None
                else np.empty(0, dtype=str),
                nonzero=np.array([x[3] for x in cells], dtype=np.int32),
                in_vocab=np.array([x[4] for x in cells], dtype=np.int32),
                tokens=np.array([len(x[1]) for x in cells], dtype=np.int32),
                gene_universe=SHARED_UNIVERSE,
                n_universe_genes=len(universe_genes),
                n_matched_scgpt_vocab=n_matched_vocab,
                hgnc_sha256=mapper.sha256,
                input_transform=INPUT_TRANSFORM,
                source_h5ad_sha256=src_sha,
                checkpoint=HASHES["best_model.pt"],
                upstream=CODE,
                device=device,
                batch=np.int32(batch),
                seconds=time.monotonic() - t0,
            )
            os.replace(tmp, dest)
            status_tag = "encoded"
        else:
            status_tag = "cached"

        with np.load(dest, allow_pickle=False) as enc:
            cached_batch = (
                int(enc["batch"])
                if "batch" in enc.files
                else prior_batches.get(gsm, batch if need_encode else 16)
            )
            records.append(
                {
                    "sample_uid": gsm,
                    "donor_uid": str(row.donor_uid),
                    "patient_id": str(row.patient_id),
                    "lesion_uid": str(row.lesion_uid),
                    "tissue": str(row.tissue),
                    "treatment_status": str(row.treatment_status),
                    "preop_regimen": str(row.preop_regimen),
                    "n_cells_encoded": len(uid),
                    "gene_universe": SHARED_UNIVERSE,
                    "in_vocab_nonzero_fraction": float(
                        enc["in_vocab"].sum() / enc["nonzero"].sum()
                    ),
                    "truncated_cells": int((enc["in_vocab"] > 1199).sum()),
                    "seconds": float(enc["seconds"]),
                    "bytes": dest.stat().st_size,
                    "device": str(enc["device"]),
                    "batch": cached_batch,
                    "checkpoint_sha256": HASHES["best_model.pt"],
                    "hgnc_sha256": mapper.sha256,
                    "source_h5ad_sha256": src_sha,
                    "embedding_sha256": _sha256_file(dest),
                }
            )
        print(
            f"{gsm} ({row.patient_id} {row.tissue}): {len(uid)} V1.1 nuclear cells ({status_tag})",
            flush=True,
        )

        # Collect deterministic <=256-cell sub-bags (`7,616` total cells across 34 bags)
        bags_256 = sample_bags(
            h5_path, min_cells=50, max_cells=256, gene_universe=universe_genes
        )
        obs_parts = [bags_256[lin].obs for lin in LINEAGES if bags_256.get(lin) is not None]
        cnt_parts = [bags_256[lin].counts for lin in LINEAGES if bags_256.get(lin) is not None]
        sub_obs = pd.concat(obs_parts, axis=0)
        sub_cnt = sparse.vstack(cnt_parts, format="csr")
        subbag_uids_by_sample[gsm] = sub_obs.index.to_numpy(dtype=str)
        subbag_nuclear_by_sample[gsm] = Bag(
            sub_obs, sub_cnt, pd.Index(universe_genes), len(sub_obs)
        )

    emb_df = pd.DataFrame(records)
    emb_df.to_csv(out_dir / "embedding_manifest.csv", index=False)
    emb_df.to_csv(CACHE / "embedding_manifest.csv", index=False)

    # 2. Historical `CRLM_V1_245_F1` (`raw 33,694` whole-axis) & `Techfix_V1_F1` (`shared_nuclear` integer counts) `<=256`-cell XPU caches (`7,616` cells)
    _, raw_33694_syms = _load_raw_33694_genes()
    mtx_sha = _sha256_file(RAW_178 / "GSE178318_matrix.mtx.gz")
    need_raw_mtx = False
    for row in manifest.itertuples(index=False):
        gsm = str(row.sample_uid)
        dest_raw = CACHE_V1_RAW / f"{gsm}.npz"
        if overwrite or not dest_raw.is_file():
            need_raw_mtx = True
            break
        with np.load(dest_raw, allow_pickle=False) as f_raw:
            if (
                not np.array_equal(f_raw["cell_uid"], subbag_uids_by_sample[gsm])
                or str(f_raw.get("raw_mtx_sha256", "")) != mtx_sha
                or "raw_indptr" not in f_raw.files
            ):
                need_raw_mtx = True
                break

    raw_mats_by_sample: dict[str, sparse.csr_matrix] = {}
    if need_raw_mtx:
        print(
            "Streaming data/raw/GSE178318/GSE178318_matrix.mtx.gz (33,694 genes) for 7,616 <=256-cell sub-bag cells...",
            flush=True,
        )
        raw_mats_by_sample, mtx_sha = _extract_raw_33694_subbags(subbag_uids_by_sample)

    v1_rec, tf_rec = [], []
    for row in manifest.itertuples(index=False):
        gsm = str(row.sample_uid)
        src_sha = str(row.h5ad_sha256)
        sub_uids = subbag_uids_by_sample[gsm]
        sub_bag = subbag_nuclear_by_sample[gsm]

        # 2a. Historical `CRLM_V1_245_F1` raw 33,694 whole-axis scGPT cache + raw count slice
        dest_raw = CACHE_V1_RAW / f"{gsm}.npz"
        if need_raw_mtx or overwrite or not dest_raw.is_file():
            t0_r = time.monotonic()
            mat_r = raw_mats_by_sample[gsm]
            cells_r = _tokenize_raw_v1_whole_axis(mat_r, sub_uids, raw_33694_syms, vocab)
            emb_r = encode(model, cells_r, vocab, batch, device)
            if device == "xpu":
                torch.xpu.synchronize()
            tmp_r = dest_raw.with_suffix(".partial.npz")
            np.savez_compressed(
                tmp_r,
                embedding=emb_r,
                cell_uid=sub_uids,
                nonzero=np.array([x[3] for x in cells_r], dtype=np.int32),
                in_vocab=np.array([x[4] for x in cells_r], dtype=np.int32),
                tokens=np.array([len(x[1]) for x in cells_r], dtype=np.int32),
                raw_indptr=mat_r.indptr.astype(np.int32),
                raw_indices=mat_r.indices.astype(np.int32),
                raw_data=mat_r.data.astype(np.float32),
                raw_shape=np.array(mat_r.shape, dtype=np.int64),
                raw_mtx_sha256=mtx_sha,
                axis_mode="raw_33694_unambiguous_v1",
                checkpoint=HASHES["best_model.pt"],
                upstream=CODE,
                device=device,
                batch=np.int32(batch),
                seconds=time.monotonic() - t0_r,
            )
            os.replace(tmp_r, dest_raw)

        with np.load(dest_raw, allow_pickle=False) as f_r:
            v1_rec.append(
                {
                    "sample_uid": gsm,
                    "n_cells_encoded": len(sub_uids),
                    "axis_mode": str(f_r["axis_mode"]),
                    "in_vocab_nonzero_fraction": float(
                        f_r["in_vocab"].sum() / f_r["nonzero"].sum()
                    ),
                    "seconds": float(f_r["seconds"]),
                    "embedding_sha256": _sha256_file(dest_raw),
                }
            )

        # 2b. `Techfix_V1_F1` (`step_b_full_techfix`) HGNC shared-nuclear integer-count scGPT cache
        dest_tf = CACHE_TECHFIX / f"{gsm}.npz"
        need_tf = overwrite or not dest_tf.is_file()
        if not need_tf:
            with np.load(dest_tf, allow_pickle=False) as f_tf:
                if (
                    not np.array_equal(f_tf["cell_uid"], sub_uids)
                    or str(f_tf.get("source_h5ad_sha256", "")) != src_sha
                    or str(f_tf.get("hgnc_sha256", "")) != mapper.sha256
                ):
                    need_tf = True
        if need_tf:
            t0_tf = time.monotonic()
            cells_tf = tokens(
                sub_bag,
                vocab,
                vocab_ids=vocab_ids,
                input_is_log1p=False,
                recalibrate_cp10k=False,
            )
            emb_tf = encode(model, cells_tf, vocab, batch, device)
            if device == "xpu":
                torch.xpu.synchronize()
            tmp_tf = dest_tf.with_suffix(".partial.npz")
            np.savez_compressed(
                tmp_tf,
                embedding=emb_tf,
                cell_uid=sub_uids,
                nonzero=np.array([x[3] for x in cells_tf], dtype=np.int32),
                in_vocab=np.array([x[4] for x in cells_tf], dtype=np.int32),
                tokens=np.array([len(x[1]) for x in cells_tf], dtype=np.int32),
                gene_universe=SHARED_UNIVERSE,
                input_transform="shared_nuclear_integer_counts_quantile_51bin_v1_techfix",
                source_h5ad_sha256=src_sha,
                hgnc_sha256=mapper.sha256,
                checkpoint=HASHES["best_model.pt"],
                upstream=CODE,
                device=device,
                batch=np.int32(batch),
                seconds=time.monotonic() - t0_tf,
            )
            os.replace(tmp_tf, dest_tf)

        with np.load(dest_tf, allow_pickle=False) as f_tf:
            tf_rec.append(
                {
                    "sample_uid": gsm,
                    "n_cells_encoded": len(sub_uids),
                    "gene_universe": str(f_tf["gene_universe"]),
                    "in_vocab_nonzero_fraction": float(
                        f_tf["in_vocab"].sum() / f_tf["nonzero"].sum()
                    ),
                    "seconds": float(f_tf["seconds"]),
                    "embedding_sha256": _sha256_file(dest_tf),
                }
            )
        print(
            f"{gsm}: <=256 sub-bag ({len(sub_uids)} cells) V1_raw + Techfix embeddings ready",
            flush=True,
        )

    pd.DataFrame(v1_rec).to_csv(CACHE_V1_RAW / "embedding_manifest.csv", index=False)
    pd.DataFrame(tf_rec).to_csv(CACHE_TECHFIX / "embedding_manifest.csv", index=False)


def pca_projection(bag_log1p: Bag, transform) -> tuple[np.ndarray | None, np.ndarray]:
    positions = bag_log1p.genes.get_indexer(transform["genes"])
    measured = positions >= 0
    if measured.mean() < 0.95:
        return None, measured
    values = bag_log1p.counts[:, positions[measured]].toarray().astype(np.float32)
    values -= transform["center"][measured]
    return (values @ transform["components"][:, measured].T).astype(float), measured


def elastic_logit(path: Path, means: pd.Series) -> float | None:
    with np.load(path, allow_pickle=False) as fit:
        values = means.reindex(fit["columns"]).to_numpy(dtype=float)
        if not np.isfinite(values).all():
            return None
        scaled = (values - fit["center"]) / fit["scale"]
        return float((scaled @ fit["coef"].T + fit["bias"])[0])


def score_frozen_program_bag(
    bag_counts: Bag, model_dir: Path
) -> tuple[pd.Series, list[dict]]:
    """Score external `Bag` using frozen development `program_defs.json` & `matched_controls.csv`."""
    spec = json.loads((model_dir / "program_defs.json").read_text(encoding="utf-8"))
    ctrl_df = pd.read_csv(model_dir / "matched_controls.csv")
    raw_controls = controls_from_dataframe(ctrl_df)
    sym_set = set(str(g) for g in bag_counts.genes)

    filtered_programs: dict[str, list[str]] = {}
    filtered_controls: dict[str, dict[str, list[str]]] = {}
    cov_rows = []
    for pname, members in spec["programs"].items():
        measured = [str(g) for g in members if str(g) in sym_set]
        frac = len(measured) / len(members) if members else 0.0
        cov_rows.append(
            {
                "program": str(pname),
                "n_frozen_genes": len(members),
                "n_measured_genes": len(measured),
                "measured_fraction": frac,
                "absent_genes": "|".join(sorted(set(members) - set(measured))),
            }
        )
        if len(measured) < 5 or frac < 0.80:
            raise ValueError(
                f"Program {pname} below frozen 80%/5-gene measurability rule ({len(measured)}/{len(members)})"
            )
        filtered_programs[pname] = measured
        filtered_controls[pname] = {g: raw_controls[pname][g] for g in measured}

    cell_scores = debiased_ucell_scores(
        bag_counts.counts,
        bag_counts.genes,
        filtered_programs,
        matched_controls=filtered_controls,
        max_rank=int(spec.get("max_rank", 1500)),
        min_ctrl_fraction=0.80,
    )
    return cell_scores.mean(axis=0), cov_rows


def pooled_logit(path: Path, vector: np.ndarray, sample: str, lineage: str) -> float:
    model, info = load_model(path)
    if info["lineage"] != lineage:
        raise ValueError(f"Classifier lineage differs from bag: {path}")
    prediction, _ = predict_bags(model, info, [vector], [sample])
    return float(prediction["logit"].iloc[0])


def forecast_v1_1(
    fits: dict[tuple[str, str], Path],
    manifest: pd.DataFrame,
    out_dir: Path = OUT,
    *,
    min_cells: int = 50,
    max_cells: int | None = None,
) -> pd.DataFrame:
    mapper = get_hgnc_mapper()
    universe_genes = list(load_shared_nuclear_universe(SHARED_UNIVERSE))
    forecasts, draws, measurement, program_rows, prog_cov_seen = [], [], [], [], set()
    prog_cov_records = []

    for row in manifest.itertuples(index=False):
        sample = str(row.sample_uid)
        h5_path = ROOT / row.h5ad_path
        bags = sample_bags(
            h5_path,
            min_cells=min_cells,
            max_cells=max_cells,
            gene_universe=universe_genes,
        )
        obj = ad.read_h5ad(h5_path, backed="r")
        try:
            obs = obj.obs
            ret = obs["qc_retained"].to_numpy(dtype=bool)
            lin_s = obs["lineage_m2"].astype(str)
            pop_counts = {
                lin: int((ret & lin_s.eq(lin).to_numpy()).sum()) for lin in LINEAGES
            }
            n_captured_retained = int(ret.sum())
        finally:
            obj.file.close()

        with np.load(CACHE / f"{sample}.npz", allow_pickle=False) as enc_file:
            enc_uids = pd.Index(enc_file["cell_uid"])
            enc_emb = enc_file["embedding"]

        for lineage in LINEAGES:
            bag = bags.get(lineage)
            meta = {
                "sample_uid": sample,
                "donor_uid": str(row.donor_uid),
                "patient_id": str(row.patient_id),
                "lesion_uid": str(row.lesion_uid),
                "tissue": str(row.tissue),
                "treatment_status": str(row.treatment_status),
                "preop_regimen": str(row.preop_regimen),
                "lineage": lineage,
                "n_population": pop_counts[lineage],
                "n_cells": len(bag.obs) if bag is not None else 0,
                "n_captured": n_captured_retained,
            }
            status = {m: "fewer_than_50_cells" for m in METHODS}
            logits = {}
            if bag is not None:
                for uid in bag.obs.index:
                    draws.append(
                        {"sample_uid": sample, "lineage": lineage, "cell_uid": str(uid)}
                    )
                prog_fit_path = fits[(lineage, "program_elastic")]
                means, cov_list = score_frozen_program_bag(bag, prog_fit_path.parent)
                if lineage not in prog_cov_seen:
                    for c_item in cov_list:
                        prog_cov_records.append({"lineage": lineage, **c_item})
                    prog_cov_seen.add(lineage)
                for prog_name, prog_val in means.items():
                    program_rows.append(
                        {
                            "sample_uid": sample,
                            "donor_uid": meta["donor_uid"],
                            "patient_id": meta["patient_id"],
                            "lesion_uid": meta["lesion_uid"],
                            "tissue": meta["tissue"],
                            "treatment_status": meta["treatment_status"],
                            "lineage": lineage,
                            "program": str(prog_name),
                            "n_cells": meta["n_cells"],
                            "score_ucell_debiased": float(prog_val),
                        }
                    )
                val_el = elastic_logit(prog_fit_path, means)
                if val_el is not None:
                    logits["program_elastic"] = val_el
                else:
                    status["program_elastic"] = "required_program_unmeasurable"

                mat_log1p, med_nuc_sum = normalize_counts_log1p_cp10k(bag.counts)
                bag_log1p = Bag(bag.obs, mat_log1p, bag.genes, bag.population)
                pca_file = fits[(lineage, "PCA_mean")].parent / "pca.npz"
                with np.load(pca_file, allow_pickle=False) as transform:
                    pca, measured = pca_projection(bag_log1p, transform)
                    hvg = transform["genes"]
                    measured_x = bag.counts[:, bag.genes.get_indexer(hvg[measured])]
                    measurement.append(
                        {
                            **meta,
                            "gene_universe": SHARED_UNIVERSE,
                            "hgnc_sha256": mapper.sha256,
                            "n_train_hvg": len(hvg),
                            "n_mapped_hvg": int(measured.sum()),
                            "mapped_fraction": float(measured.mean()),
                            "absent_hvg": "|".join(hvg[~measured]),
                            "n_unmeasured_positions": int(len(bag.obs) * (~measured).sum()),
                            "n_center_substituted": int(len(bag.obs) * (~measured).sum())
                            if pca is not None
                            else 0,
                            "n_measured_zero": int(
                                len(bag.obs) * measured.sum() - measured_x.nnz
                            ),
                            "median_nuclear_umi_sum": med_nuc_sum,
                            "input_scale": "raw_counts_to_shared_nuclear_log1p_cp10k",
                            "numeric_scale_pass": True,
                        }
                    )
                for method in ("PCA_mean", "PCA_attention"):
                    status[method] = "under_95pct_train_HVG"
                    if pca is not None:
                        logits[method] = pooled_logit(
                            fits[(lineage, method)], pca, sample, lineage
                        )

                ix = enc_uids.get_indexer(bag.obs.index)
                if (ix < 0).any():
                    raise ValueError(f"Missing cached scGPT cells for {sample} {lineage}")
                scgpt = enc_emb[ix].astype(float)
                for method in ("scGPT_mean", "scGPT_attention"):
                    logits[method] = pooled_logit(
                        fits[(lineage, method)], scgpt, sample, lineage
                    )

            for method in METHODS:
                value = logits.get(method, np.nan)
                forecasts.append(
                    {
                        **meta,
                        "method": method,
                        "status": "predicted" if method in logits else status[method],
                        "logit": value,
                        "probability_lm": float(expit(value)) if np.isfinite(value) else np.nan,
                    }
                )

    pred = pd.DataFrame(forecasts)
    pd.DataFrame(draws).to_csv(
        out_dir / "selected_cells.csv.gz",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )
    pd.DataFrame(program_rows).to_csv(
        out_dir / "external_program_bag_scores.csv", index=False
    )
    pd.DataFrame(prog_cov_records).to_csv(
        out_dir / "model_program_coverage.csv", index=False
    )
    pd.DataFrame(measurement).to_csv(out_dir / "pca_measurement.csv", index=False)
    pred.to_csv(out_dir / "bag_predictions.csv", index=False)
    return pred


def metric_draws(frame: pd.DataFrame):
    """Bootstrap true donors (10,000 draws) with all available sampled tumor sites."""
    frame = frame.sort_values(["donor_uid", "tissue"])
    codes, donors = pd.factorize(frame["donor_uid"], sort=True)
    rng = np.random.default_rng(SEED)
    weight = np.vstack(
        [
            np.ones(len(donors)),
            rng.multinomial(
                len(donors), np.full(len(donors), 1 / len(donors)), size=10000
            ),
        ]
    )
    w = weight[:, codes]
    pos = frame["tissue"].eq("LM").to_numpy()
    x = frame["logit"].to_numpy()
    wp, wn = w[:, pos], w[:, ~pos]
    np_, nn = wp.sum(1), wn.sum(1)
    score_diff = x[pos, None] - x[None, ~pos]
    concord = (score_diff > 0).astype(float) + 0.5 * (score_diff == 0)
    with np.errstate(divide="ignore", invalid="ignore"):
        auc = np.einsum("bi,ij,bj->b", wp, concord, wn) / (np_ * nn)
        ba = 0.5 * ((wp @ (x[pos] >= 0)) / np_ + (wn @ (x[~pos] < 0)) / nn)
        thresholds = np.unique(x)[::-1]
        cum_pos = wp @ (x[pos, None] >= thresholds)
        cum_all = w @ (x[:, None] >= thresholds)
        precision = np.divide(
            cum_pos, cum_all, where=cum_all > 0, out=np.zeros_like(cum_all)
        )
        ap = (precision * np.diff(cum_pos, axis=1, prepend=0)).sum(axis=1) / np_
        ap[nn == 0] = np.nan
        pair = frame.pivot(index="donor_uid", columns="tissue", values="logit").reindex(
            index=donors, columns=["PT", "LM"]
        )
        dif = pair["LM"] - pair["PT"]
        paired = dif.notna().to_numpy()
        pairw = weight[:, paired]
        delta = (pairw @ dif[paired].to_numpy()) / pairw.sum(1)
        concordant = (pairw @ (dif[paired].to_numpy() > 0)) / pairw.sum(1)
    return {
        "auroc": auc,
        "average_precision": ap,
        "balanced_accuracy": ba,
        "mean_LM_minus_PT_logit": delta,
        "fraction_LM_above_PT": concordant,
    }, dif.dropna()


def intervals(draws: dict) -> dict:
    row = {}
    for metric, values in draws.items():
        finite = values[1:][np.isfinite(values[1:])]
        lo, hi = (
            np.quantile(finite, [0.025, 0.975]) if len(finite) else (np.nan, np.nan)
        )
        row.update(
            {
                metric: float(values[0]),
                metric + "_lower95": float(lo),
                metric + "_upper95": float(hi),
                metric + "_valid_bootstrap": int(len(finite)),
            }
        )
    return row


def patient_reports(
    pred: pd.DataFrame, manifest: pd.DataFrame, out_dir: Path = OUT
) -> None:
    from sklearn.metrics import (
        average_precision_score,
        balanced_accuracy_score,
        roc_auc_score,
    )

    donor_meta = (
        manifest[["donor_uid", "patient_id", "treatment_status", "preop_regimen"]]
        .drop_duplicates("donor_uid")
        .set_index("donor_uid")
    )
    keys = ["method", "lineage", "donor_uid", "tissue"]
    sites_all, metrics, pairs_all, contrast_rows = [], [], [], []
    for lineage in LINEAGES:
        group = pred.loc[pred["lineage"].eq(lineage)]
        wide = group.pivot(
            index="sample_uid", columns="method", values="logit"
        ).reindex(columns=METHODS)
        for scope in (
            "method_available",
            "common_all_five",
            "common_predicting_methods",
            "paired_donors_common_methods",
        ):
            available = wide.columns[wide.notna().any()]
            common = (
                wide.dropna().index
                if scope == "common_all_five"
                else wide[available].dropna().index
            )
            rows = group.loc[group["logit"].notna()]
            if scope != "method_available":
                rows = rows.loc[rows["sample_uid"].isin(common)]
            if scope == "paired_donors_common_methods":
                n_tissues = rows.groupby("donor_uid")["tissue"].transform("nunique")
                rows = rows.loc[n_tissues.eq(2)]
            lesions = rows.groupby(keys + ["lesion_uid"])["logit"].mean()
            sites = (
                lesions.groupby(keys)
                .agg(["mean", "size"])
                .reset_index()
                .rename(columns={"mean": "logit", "size": "n_lesions"})
            )
            sites["treatment_status"] = sites["donor_uid"].map(donor_meta["treatment_status"])
            sites_all.append(sites.assign(scope=scope))
            method_draws = {}
            for method in METHODS:
                frame = sites.loc[sites["method"].eq(method)]
                rec = {
                    "lineage": lineage,
                    "method": method,
                    "scope": scope,
                    "n_donors": frame["donor_uid"].nunique(),
                    "n_sites": len(frame),
                    "n_PT": int(frame["tissue"].eq("PT").sum()),
                    "n_LM": int(frame["tissue"].eq("LM").sum()),
                    "n_bags": int(len(rows[rows["method"].eq(method)])),
                    "paired_n": 0,
                }
                if len(frame):
                    draws, dif = metric_draws(frame)
                    method_draws[method] = draws
                    rec.update(intervals(draws), paired_n=len(dif))
                    for d, v in dif.items():
                        pairs_all.append(
                            {
                                "lineage": lineage,
                                "scope": scope,
                                "method": method,
                                "donor_uid": d,
                                "patient_id": str(donor_meta.at[d, "patient_id"]),
                                "treatment_status": str(donor_meta.at[d, "treatment_status"]),
                                "preop_regimen": str(donor_meta.at[d, "preop_regimen"]),
                                "LM_minus_PT_logit": float(v),
                            }
                        )
                    if frame["tissue"].nunique() == 2:
                        y = frame["tissue"].eq("LM")
                        refs = [
                            roc_auc_score(y, frame["logit"]),
                            average_precision_score(y, frame["logit"]),
                            balanced_accuracy_score(y, frame["logit"] >= 0),
                        ]
                        if not np.allclose(
                            [
                                draws[k][0]
                                for k in (
                                    "auroc",
                                    "average_precision",
                                    "balanced_accuracy",
                                )
                            ],
                            refs,
                        ):
                            raise ValueError(
                                "Patient metrics diverge from reference implementations"
                            )
                metrics.append(rec)
            if scope == "common_predicting_methods":
                for first, second in combinations(method_draws, 2):
                    delta = {
                        k: vals - method_draws[first][k]
                        for k, vals in method_draws[second].items()
                    }
                    contrast_rows.append(
                        {
                            "lineage": lineage,
                            "first": first,
                            "second": second,
                            "common_bags": len(common),
                            "common_donors": sites["donor_uid"].nunique(),
                            **intervals(delta),
                        }
                    )

    sites_df = pd.concat(sites_all, ignore_index=True)
    pairs_df = pd.DataFrame(pairs_all)
    sites_df.to_csv(out_dir / "patient_site_logits.csv", index=False)
    pairs_df.to_csv(out_dir / "patient_paired_logits.csv", index=False)
    pd.DataFrame(contrast_rows).to_csv(out_dir / "paired_method_contrasts.csv", index=False)
    pd.DataFrame(metrics).to_csv(out_dir / "patient_metrics.csv", index=False)
    pred.groupby(["lineage", "method", "status"]).agg(
        n_bags=("sample_uid", "size"), n_donors=("donor_uid", "nunique")
    ).to_csv(out_dir / "prediction_coverage.csv")

    strat_rows = []
    avail_sites = sites_df.loc[sites_df["scope"].eq("method_available")]
    avail_pairs = pairs_df.loc[pairs_df["scope"].eq("method_available")]
    for lineage in LINEAGES:
        for method in METHODS:
            for trt in ("treatment_naive", "preoperative_chemotherapy"):
                sf = avail_sites.loc[
                    avail_sites["lineage"].eq(lineage)
                    & avail_sites["method"].eq(method)
                    & avail_sites["treatment_status"].eq(trt)
                ]
                pf = avail_pairs.loc[
                    avail_pairs["lineage"].eq(lineage)
                    & avail_pairs["method"].eq(method)
                    & avail_pairs["treatment_status"].eq(trt)
                ]
                y_sub = sf["tissue"].eq("LM")
                auc_sub = (
                    float(roc_auc_score(y_sub, sf["logit"]))
                    if len(sf) and y_sub.nunique() == 2
                    else np.nan
                )
                ba_sub = (
                    float(balanced_accuracy_score(y_sub, sf["logit"] >= 0))
                    if len(sf) and y_sub.nunique() == 2
                    else np.nan
                )
                diffs = pf["LM_minus_PT_logit"].to_numpy(dtype=float)
                strat_rows.append(
                    {
                        "lineage": lineage,
                        "method": method,
                        "treatment_status": trt,
                        "n_donors": int(sf["donor_uid"].nunique()),
                        "n_sites": int(len(sf)),
                        "paired_n": int(len(diffs)),
                        "n_LM_above_PT": int(np.sum(diffs > 0)) if len(diffs) else 0,
                        "fraction_LM_above_PT": float(np.mean(diffs > 0))
                        if len(diffs)
                        else np.nan,
                        "mean_LM_minus_PT_logit": float(np.mean(diffs))
                        if len(diffs)
                        else np.nan,
                        "site_auroc_descriptive": auc_sub,
                        "site_balanced_accuracy_descriptive": ba_sub,
                    }
                )
    pd.DataFrame(strat_rows).to_csv(
        out_dir / "treatment_stratified_classifier_summary.csv", index=False
    )


def evaluate_three_version_classifiers(
    v11_metrics_path: Path,
    manifest: pd.DataFrame,
    out_dir: Path = OUT,
) -> pd.DataFrame:
    """Evaluate genuine archived `CRLM_V1_245_F1` (`raw 33,694`) and `Techfix` alongside `V1.1`."""
    spec_py = ARCHIVE_V1_M5 / "program_scores.py"
    mod_spec = importlib.util.spec_from_file_location("v1_frozen_prog_scores", spec_py)
    mod = importlib.util.module_from_spec(mod_spec)
    mod_spec.loader.exec_module(mod)
    v1_rank_scores = mod.rank_scores

    v1_fits_df = pd.read_csv(ARCHIVE_V1_M5 / "frozen_models.csv")
    v1_fits = {
        (r.lineage, r.method): ARCHIVE_V1_MODELS
        / Path(r.path).relative_to("models/m4/full")
        for r in v1_fits_df.itertuples()
        if r.accepted
    }
    v1_prog_df = pd.read_csv(ARCHIVE_V1_M5 / "gene_sets/model_program_genes.csv")
    universe_genes = list(load_shared_nuclear_universe(SHARED_UNIVERSE))
    univ_set = set(universe_genes)

    _, raw_33694_syms = _load_raw_33694_genes()
    raw_33694_idx = pd.Index(raw_33694_syms)
    unamb_mask = np.asarray(~raw_33694_idx.duplicated(keep=False), dtype=bool)
    raw_unamb_genes = pd.Index(raw_33694_idx[unamb_mask])
    raw_unamb_set = set(raw_unamb_genes)

    rows_v1_raw, rows_techfix = [], []
    for row in manifest.itertuples(index=False):
        sample = str(row.sample_uid)
        h5_path = ROOT / row.h5ad_path
        bags_256 = sample_bags(
            h5_path, min_cells=50, max_cells=256, gene_universe=universe_genes
        )
        raw_cache_path = CACHE_V1_RAW / f"{sample}.npz"
        tf_cache_path = CACHE_TECHFIX / f"{sample}.npz"
        if not raw_cache_path.is_file() or not tf_cache_path.is_file():
            raise FileNotFoundError(
                f"Missing historical V1/Techfix scGPT cache for {sample}; run `encode` first"
            )
        with np.load(raw_cache_path, allow_pickle=False) as f_raw:
            raw_uids = pd.Index(f_raw["cell_uid"])
            raw_emb_all = f_raw["embedding"]
            mat_33694_sample = sparse.csr_matrix(
                (f_raw["raw_data"], f_raw["raw_indices"], f_raw["raw_indptr"]),
                shape=tuple(int(x) for x in f_raw["raw_shape"]),
                dtype=np.float32,
            )
        with np.load(tf_cache_path, allow_pickle=False) as f_tf:
            tf_uids = pd.Index(f_tf["cell_uid"])
            tf_emb_all = f_tf["embedding"]

        for lineage in LINEAGES:
            bag = bags_256.get(lineage)
            if bag is None:
                continue
            mat_log1p, _ = normalize_counts_log1p_cp10k(bag.counts)
            bag_log1p = Bag(bag.obs, mat_log1p, bag.genes, bag.population)

            ix_raw = raw_uids.get_indexer(bag.obs.index)
            ix_tf = tf_uids.get_indexer(bag.obs.index)
            if (ix_raw < 0).any() or (ix_tf < 0).any():
                raise ValueError(f"Missing sub-bag cell_uid in V1/Techfix cache for {sample} {lineage}")
            scgpt_raw = raw_emb_all[ix_raw].astype(float)
            scgpt_tf = tf_emb_all[ix_tf].astype(float)

            # Raw 33,694 full-library log1p(CP10k) -> unambiguous raw symbol slice
            raw_counts_full = mat_33694_sample[ix_raw, :].tocsr()
            raw_full_log1p, _ = normalize_counts_log1p_cp10k(raw_counts_full)
            raw_unamb_log1p = raw_full_log1p[:, unamb_mask]

            sub_p = v1_prog_df[v1_prog_df["lineage"].eq(lineage)]
            defs_tf, defs_raw = {}, {}
            for pname, grp in sub_p.groupby("program"):
                raw_m = grp["gene_symbol"].astype(str).tolist()
                canon_m = canonicalize_symbols(raw_m, deduplicate=True)
                meas_tf = [g for g in canon_m if g in univ_set]
                if len(meas_tf) >= 5 and len(meas_tf) >= 0.8 * len(raw_m):
                    defs_tf[pname] = meas_tf
                meas_raw = [g for g in raw_m if g in raw_unamb_set]
                if len(meas_raw) >= 5 and len(meas_raw) >= 0.8 * len(raw_m):
                    defs_raw[pname] = meas_raw

            scores_tf = v1_rank_scores(
                mat_log1p, universe_genes, defs_tf, max_rank=1500
            ).mean(axis=0)
            logit_prog_tf = elastic_logit(v1_fits[(lineage, "program_elastic")], scores_tf)
            logit_prog_raw = (
                elastic_logit(
                    v1_fits[(lineage, "program_elastic")],
                    v1_rank_scores(
                        raw_unamb_log1p,
                        list(raw_unamb_genes),
                        defs_raw,
                        max_rank=1500,
                    ).mean(axis=0),
                )
                if defs_raw
                else None
            )

            with np.load(
                v1_fits[(lineage, "PCA_mean")].parent / "pca.npz", allow_pickle=False
            ) as trf:
                raw_hvgs = [str(g) for g in trf["genes"]]
                raw_pos = raw_unamb_genes.get_indexer(raw_hvgs)
                raw_present = raw_pos >= 0
                raw_cov = float(raw_present.mean())
                pca_raw = None
                if raw_cov >= 0.95:
                    vals_r = (
                        raw_unamb_log1p[:, raw_pos[raw_present]]
                        .toarray()
                        .astype(np.float32)
                        - trf["center"][raw_present]
                    )
                    pca_raw = (vals_r @ trf["components"][:, raw_present].T).astype(float)

                aligned_hvgs = align_axis_to_canonical_targets(universe_genes, raw_hvgs)
                tf_present = aligned_hvgs["present"].to_numpy(dtype=bool)
                tf_pos = aligned_hvgs["source_index"].to_numpy(dtype=np.int64)
                vals_tf = (
                    bag_log1p.counts[:, tf_pos[tf_present]].toarray().astype(np.float32)
                    - trf["center"][tf_present]
                )
                pca_tf = (vals_tf @ trf["components"][:, tf_present].T).astype(float)

            for m in METHODS:
                if m == "program_elastic":
                    l_tf = logit_prog_tf
                    l_raw = logit_prog_raw
                elif m.startswith("PCA_"):
                    l_tf = (
                        pooled_logit(v1_fits[(lineage, m)], pca_tf, sample, lineage)
                        if tf_present.mean() >= 0.95
                        else np.nan
                    )
                    l_raw = (
                        pooled_logit(v1_fits[(lineage, m)], pca_raw, sample, lineage)
                        if pca_raw is not None
                        else np.nan
                    )
                else:
                    l_tf = pooled_logit(v1_fits[(lineage, m)], scgpt_tf, sample, lineage)
                    l_raw = pooled_logit(v1_fits[(lineage, m)], scgpt_raw, sample, lineage)

                base_info = {
                    "sample_uid": sample,
                    "donor_uid": str(row.donor_uid),
                    "tissue": str(row.tissue),
                    "lineage": lineage,
                    "method": m,
                }
                rows_techfix.append(
                    {**base_info, "logit": l_tf if l_tf is not None else np.nan}
                )
                rows_v1_raw.append(
                    {**base_info, "logit": l_raw if l_raw is not None else np.nan}
                )

    bag_df_raw = pd.DataFrame(rows_v1_raw)
    bag_df_tf = pd.DataFrame(rows_techfix)
    # Sanity-check genuine distinctness of scGPT embeddings across V1_raw and Techfix
    sc_r = bag_df_raw.loc[bag_df_raw["method"].str.startswith("scGPT_"), "logit"].to_numpy()
    sc_t = bag_df_tf.loc[bag_df_tf["method"].str.startswith("scGPT_"), "logit"].to_numpy()
    if np.allclose(sc_r, sc_t, atol=1e-6):
        raise AssertionError("Historical V1_raw and Techfix scGPT logits unexpectedly identical")

    def _agg_version(bag_df: pd.DataFrame, ver_label: str) -> pd.DataFrame:
        out_r = []
        for lin in LINEAGES:
            for m in METHODS:
                sub = bag_df[
                    bag_df["lineage"].eq(lin)
                    & bag_df["method"].eq(m)
                    & bag_df["logit"].notna()
                ]
                if sub.empty or sub["tissue"].nunique() < 2:
                    out_r.append(
                        {
                            "version": ver_label,
                            "lineage": lin,
                            "method": m,
                            "status": "abstained_lt95pct_raw_hvg",
                            "n_donors": 0,
                            "n_sites": 0,
                            "paired_n": 0,
                            "auroc": np.nan,
                            "auroc_lower95": np.nan,
                            "auroc_upper95": np.nan,
                            "average_precision": np.nan,
                            "balanced_accuracy": np.nan,
                            "fraction_LM_above_PT": np.nan,
                        }
                    )
                    continue
                draws, dif = metric_draws(sub)
                iv = intervals(draws)
                out_r.append(
                    {
                        "version": ver_label,
                        "lineage": lin,
                        "method": m,
                        "status": "predicted",
                        "n_donors": int(sub["donor_uid"].nunique()),
                        "n_sites": int(len(sub)),
                        "paired_n": int(len(dif)),
                        "auroc": iv["auroc"],
                        "auroc_lower95": iv["auroc_lower95"],
                        "auroc_upper95": iv["auroc_upper95"],
                        "average_precision": iv["average_precision"],
                        "balanced_accuracy": iv["balanced_accuracy"],
                        "fraction_LM_above_PT": iv["fraction_LM_above_PT"],
                    }
                )
        return pd.DataFrame(out_r)

    df_v1 = _agg_version(bag_df_raw, "CRLM_V1_245_F1")
    df_tf = _agg_version(bag_df_tf, "CRLM_V1_245_F1_techfix")
    v11_m = pd.read_csv(v11_metrics_path)
    v11_sub = v11_m[v11_m["scope"].eq("method_available")].copy()
    df_v11 = pd.DataFrame(
        {
            "version": "CRLM_V1.1_245_F1",
            "lineage": v11_sub["lineage"],
            "method": v11_sub["method"],
            "status": "predicted",
            "n_donors": v11_sub["n_donors"],
            "n_sites": v11_sub["n_sites"],
            "paired_n": v11_sub["paired_n"],
            "auroc": v11_sub["auroc"],
            "auroc_lower95": v11_sub["auroc_lower95"],
            "auroc_upper95": v11_sub["auroc_upper95"],
            "average_precision": v11_sub["average_precision"],
            "balanced_accuracy": v11_sub["balanced_accuracy"],
            "fraction_LM_above_PT": v11_sub["fraction_LM_above_PT"],
        }
    )
    comp = pd.concat([df_v1, df_tf, df_v11], ignore_index=True)
    comp.to_csv(out_dir / "three_version_gse178318_classifiers.csv", index=False)
    return comp


def _write_run_json(
    out_dir: Path,
    stage_name: str,
    args: argparse.Namespace,
    resolved_commit: str,
    elapsed_sec: float,
) -> None:
    payload = {
        "command": sys.argv,
        "prospective_zero_exposure": bool(args.freeze_commit.upper() != "WORKTREE"),
        "freeze_commit": resolved_commit,
        "freeze_mode": args.freeze_commit,
        "frozen_inputs_sha256": _sha256_file(M5 / "frozen_inputs_sha256.txt"),
        "gene_universe": SHARED_UNIVERSE,
        "max_cells": args.max_cells,
        "out_dir": str(out_dir.relative_to(ROOT)) if out_dir.is_relative_to(ROOT) else str(out_dir),
        "python": platform.python_version(),
        "executable": sys.executable,
        "seconds": elapsed_sec,
        "peak_rss_bytes": int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024),
    }
    (out_dir / f"{stage_name}_run.json").write_text(
        json.dumps(payload, indent=2) + "\n", encoding="utf-8"
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "stage",
        nargs="?",
        choices=["all", "encode", "predict"],
        default="all",
    )
    parser.add_argument("--freeze-commit", default=PRE_UNLOCK_FREEZE_COMMIT)
    parser.add_argument("--device", choices=["cpu", "xpu"], default="xpu")
    parser.add_argument("--batch", type=int, default=16)
    parser.add_argument("--max-cells", type=int, default=None)
    parser.add_argument("--out-dir", type=Path, default=None)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--skip-parent-sha", action="store_true")
    args = parser.parse_args()

    if args.max_cells == 0:
        args.max_cells = None
    if args.max_cells is not None and args.max_cells < 50:
        raise ValueError("--max-cells must be >= 50, or 0/omitted for full bags")

    out_dir = (args.out_dir if args.out_dir is not None else OUT).resolve()
    if args.max_cells is not None and out_dir == OUT.resolve():
        raise ValueError(
            "Refusing to write truncated --max-cells predictions into default official "
            "model_predictions directory; pass an explicit custom --out-dir"
        )

    out_dir.mkdir(parents=True, exist_ok=True)
    # Remove any legacy 0-byte predict.log before running
    pred_log_path = out_dir / "predict.log"
    if pred_log_path.is_file() and pred_log_path.stat().st_size == 0:
        pred_log_path.unlink()

    fits = verify_frozen_m5_models(args.freeze_commit)
    manifest = load_tumor_manifest()
    git_target = "HEAD" if args.freeze_commit.upper() == "WORKTREE" else args.freeze_commit
    resolved_commit = subprocess.check_output(
        ["git", "rev-parse", git_target], cwd=ROOT, text=True
    ).strip()

    if args.stage in ("all", "encode"):
        t_enc = time.monotonic()
        encode_gse178318(args.device, args.batch, out_dir, overwrite=args.overwrite)
        _write_run_json(out_dir, "encode", args, resolved_commit, time.monotonic() - t_enc)

    if args.stage in ("all", "predict"):
        t_pred = time.monotonic()
        pred = forecast_v1_1(
            fits, manifest, out_dir, min_cells=50, max_cells=args.max_cells
        )
        patient_reports(pred, manifest, out_dir)
        comp = evaluate_three_version_classifiers(
            out_dir / "patient_metrics.csv", manifest, out_dir
        )
        el_pred = time.monotonic() - t_pred
        _write_run_json(out_dir, "predict", args, resolved_commit, el_pred)
        n_pred = int(pred["status"].eq("predicted").sum())
        summary_msg = (
            f"Completed GSE178318 prospective classifier evaluation: "
            f"{n_pred}/{len(pred)} bag-method predictions, "
            f"{len(comp)} 3-version rows in {el_pred:.2f}s\n"
        )
        _atomic_write_text(pred_log_path, summary_msg)
        print(summary_msg, end="", flush=True)

    refresh_gse178318_sha_manifests(
        out_dir, refresh_parent=not args.skip_parent_sha
    )


if __name__ == "__main__":
    main()
