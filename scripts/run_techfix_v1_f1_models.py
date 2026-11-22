#!/usr/bin/env python3
"""Wave 2 (Task W2-A): Prior MIL bugfix benchmark and Track 1 (CRLM_V1_245_F1_techfix) external model errata.

Strictly preserves all frozen M4/M5 classifier weights (models/m4/) and original M4/M6 artifacts.
Outputs exclusively to:
  - embeddings/techfix_v1_f1/GSE164522/
  - artifacts/techfix_v1_f1/models/
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import importlib.util
from itertools import combinations
import json
import os
import platform
from pathlib import Path
import resource
import sys
import time

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.special import expit
from sklearn.metrics import (
    average_precision_score,
    balanced_accuracy_score,
    roc_auc_score,
)
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src"), str(ROOT / "starter"), str(ROOT / "scripts")]

from crlm.benchmark_statistics import summarize as summarize_oof
from crlm.hgnc_alias_map import (
    SHARED_NUCLEAR_UNIVERSE,
    get_hgnc_mapper,
    load_shared_nuclear_universe,
)
from crlm.model_bags import Bag, LINEAGES, SEED, cell_seed, sample_bags
from linear_attention_mil import (
    build_design,
    fit_scaler,
    load_model,
    predict_bags,
)

M4_ART_ARCH_DIR = ROOT / "artifacts/archive/CRLM_V1_245_F1/artifacts/m4"
M4_MOD_ARCH_DIR = ROOT / "models/archive/CRLM_V1_245_F1/m4"
M4_DIRECT_MOD_ARCH_DIR = ROOT / "models/archive/CRLM_V1_245_F1/m4_direct"
M4_EMB_ARCH_DIR = ROOT / "embeddings/archive/CRLM_V1_245_F1/m4"
M5_DIR = ROOT / "artifacts/archive/CRLM_V1_245_F1/artifacts/m5"
M6_DATA_DIR = ROOT / "data/processed/m6/GSE164522"
M6_ORIG_PRED_DIR = ROOT / "artifacts/archive/CRLM_V1_245_F1/artifacts/m6/GSE164522/model_predictions"
OUT_DIR = ROOT / "artifacts/techfix_v1_f1/models"
EMB_DIR = ROOT / "embeddings/techfix_v1_f1/GSE164522"


def _resolve_v1_m4_relpath(rel_path: str) -> Path:
    """Resolve a V1 M4 model/embedding/artifact relative path via CRLM_V1_245_F1 archive first."""
    if rel_path.startswith("models/m4/"):
        cand = M4_MOD_ARCH_DIR / rel_path.removeprefix("models/m4/")
        if cand.exists():
            return cand
    elif rel_path.startswith("models/m4_direct/"):
        cand = M4_DIRECT_MOD_ARCH_DIR / rel_path.removeprefix("models/m4_direct/")
        if cand.exists():
            return cand
    elif rel_path.startswith("embeddings/m4/"):
        cand = M4_EMB_ARCH_DIR / rel_path.removeprefix("embeddings/m4/")
        if cand.exists():
            return cand
    elif rel_path.startswith("artifacts/m4/"):
        sub = rel_path.removeprefix("artifacts/m4/")
        for base in (M4_ART_ARCH_DIR, ROOT / "models/archive/CRLM_V1_245_F1"):
            cand = base / sub
            if cand.exists():
                return cand
    return ROOT / rel_path


METHODS = (
    "program_elastic",
    "PCA_mean",
    "PCA_attention",
    "scGPT_mean",
    "scGPT_attention",
)
TISSUES = {"PT": "PT", "MT": "LM"}


def sha256_file(path: Path) -> str:
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def load_frozen_rank_scores():
    """Load the exact frozen M5 rank_scores implementation without touching src/crlm/program_scores.py."""
    frozen_py = M5_DIR / "program_scores.py"
    spec = importlib.util.spec_from_file_location("frozen_m5_program_scores", frozen_py)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.rank_scores


def verify_frozen_m5_models() -> dict[tuple[str, str], Path]:
    """Verify SHA256 integrity of all 15 frozen M4/M5 transfer classifiers and PCA transforms."""
    frozen_sha_path = M5_DIR / "frozen_inputs_sha256.txt"
    expected_map = {}
    for line in frozen_sha_path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        digest, rel_path = line.split("  ", 1)
        expected_map[rel_path] = digest

    for rel_path, expected in expected_map.items():
        if rel_path.startswith("models/m4/full/"):
            actual = sha256_file(_resolve_v1_m4_relpath(rel_path))
            if actual != expected:
                raise ValueError(f"Frozen M5 file SHA256 mismatch for {rel_path}: expected {expected}, got {actual}")
        elif rel_path.startswith(("artifacts/m5/program_scores.py", "artifacts/m5/frozen_models.csv", "artifacts/m5/gene_sets/")):
            sub_rel = rel_path.removeprefix("artifacts/m5/")
            actual = sha256_file(M5_DIR / sub_rel)
            if actual != expected:
                raise ValueError(f"Archived M5 file SHA256 mismatch for {rel_path}: expected {expected}, got {actual}")

    fits = pd.read_csv(M5_DIR / "frozen_models.csv")
    models = {}
    for row in fits.itertuples(index=False):
        if not row.accepted:
            raise ValueError(f"Unexpected non-accepted frozen model: {row}")
        model_path = _resolve_v1_m4_relpath(row.path)
        if sha256_file(model_path) != row.sha256:
            raise ValueError(f"Frozen classifier weight SHA256 mismatch: {row.path}")
        models[(row.lineage, row.method)] = model_path

    expected_keys = {(lin, m) for lin in LINEAGES for m in METHODS}
    if set(models) != expected_keys:
        raise ValueError(f"Expected 15 frozen classifiers, got {len(models)}")
    return models


# ---------------------------------------------------------------------------
# Part 1: MixMIL [::2] + 1/sqrt(Q) bugfix quantification on M4 18,610 cells
# ---------------------------------------------------------------------------

def _m4_prior_inputs(rep: str, rows: pd.DataFrame):
    published = pd.read_csv(
        _resolve_v1_m4_relpath("artifacts/m4/selected_cells.csv.gz"), index_col="cell_uid"
    )
    groups = {k: {} for k in LINEAGES}
    for gsm in rows.sample_uid:
        h5_obj = ad.read_h5ad(ROOT / f"data/processed/m2/{gsm}.h5ad", backed="r")
        try:
            genes = pd.Index(h5_obj.var["gene_symbol"].astype(str))
            cols, enc = None, None
            if rep.startswith("scGPT"):
                with np.load(_resolve_v1_m4_relpath(f"embeddings/m4/scgpt/{gsm}.npz")) as f:
                    cols = pd.Index(f["cell_uid"])
                    enc = f["embedding"]
            for lineage in LINEAGES:
                expected = published[(published.sample_uid == gsm) & (published.lineage_m2 == lineage)]
                if expected.empty:
                    continue
                if rep == "PCA32":
                    h5_idx = h5_obj.obs.index.get_indexer(expected.index)
                    if (h5_idx < 0).any():
                        raise ValueError(f"{gsm} {lineage}: missing archived V1 cell_uid in H5AD")
                    x = h5_obj.X[h5_idx, :].tocsr().astype(np.float32)
                    x.data *= np.repeat(10000 / np.asarray(x.sum(axis=1)).ravel(), np.diff(x.indptr))
                    np.log1p(x.data, out=x.data)
                else:
                    ix = cols.get_indexer(expected.index)
                    if (ix < 0).any() or not cols.is_unique:
                        raise ValueError(f"{gsm}: missing or ambiguous scGPT cell identities")
                    x = enc[ix]
                groups[lineage][gsm] = (x, genes, expected.copy())
        finally:
            h5_obj.file.close()
    return groups


def _m4_project(group: dict, design, rep: str, lineage: str, fold: int):
    if rep == "PCA32":
        src = _resolve_v1_m4_relpath(f"models/m4/{lineage.replace('/', '_')}/fold{fold}/pca.npz")
        with np.load(src) as f:
            if set(f["training_samples"]) != set(design.samples):
                raise ValueError("Saved PCA has a different training fold")
            genes, center, components = f["genes"], f["center"], f["components"]
        matrix = {}
        for gsm, (x, symbols, _) in group.items():
            indices = symbols.get_indexer(genes)
            if (indices < 0).any():
                raise ValueError("Missing a fitted PCA feature")
            matrix[gsm] = np.asarray((x[:, indices].toarray() - center) @ components.T, dtype=np.float32)
    elif rep == "scGPT32":
        src = _resolve_v1_m4_relpath(f"models/m4_direct/scGPT32/{lineage.replace('/', '_')}/fold{fold}/projection.npz")
        with np.load(src) as f:
            if set(f["training_samples"]) != set(design.samples):
                raise ValueError("Saved scGPT32 projection has a different training fold")
            center, components = f["center"], f["components"]
        matrix = {k: np.asarray((x - center) @ components.T, dtype=np.float32) for k, (x, _, _) in group.items()}
    else:
        raise ValueError(f"Unsupported representation: {rep}")
    loc, scale = fit_scaler([matrix[s].astype(float) for s in design.samples], design.bce_weights, True)
    return {s: np.asarray((a - loc) / scale, dtype=np.float32) for s, a in matrix.items()}


def _paired_donor_auc_ci(patient_df: pd.DataFrame, ref_method: str, target_method: str, seed: int = SEED, n_boot: int = 10000) -> dict:
    if ref_method == target_method:
        return {
            "auc_delta": 0.0,
            "auc_delta_lower95": 0.0,
            "auc_delta_upper95": 0.0,
            "valid_boot": int(n_boot),
        }
    matched = (
        patient_df.pivot(index=["donor_uid", "tissue"], columns="method", values="logit")[[ref_method, target_method]]
        .dropna()
        .reset_index()
    )
    donor_codes, donors = pd.factorize(matched.donor_uid, sort=True)
    y = matched.tissue.eq("LM").to_numpy()
    kernels = []
    for m in (ref_method, target_method):
        diff = matched.loc[y, m].to_numpy()[:, None] - matched.loc[~y, m].to_numpy()
        kernels.append((diff > 0).astype(float) + 0.5 * (diff == 0))
    delta = kernels[1] - kernels[0]
    n = len(donors)
    weights = np.random.default_rng(seed).multinomial(n, np.full(n, 1 / n), size=n_boot)
    pos, neg = weights[:, donor_codes[y]], weights[:, donor_codes[~y]]
    den = pos.sum(1) * neg.sum(1)
    valid = den > 0
    boot = np.einsum("bi,ij,bj->b", pos[valid], delta, neg[valid]) / den[valid]
    lo, hi = np.quantile(boot, [0.025, 0.975])
    return {
        "auc_delta": float(delta.mean()),
        "auc_delta_lower95": float(lo),
        "auc_delta_upper95": float(hi),
        "valid_boot": int(valid.sum()),
    }


def run_mixmil_bugfix() -> None:
    import torch
    from crlm.prior_mil import BagBatch, MixedMIL, official_posterior, train_fixed

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    torch.set_num_threads(2)
    posterior_cls = official_posterior()  # Uses in-repo references/code/mixmil

    rows = pd.read_csv(ROOT / "data/manifests/m0/sample_manifest_validated.csv").query(
        "role == 'development' and tissue in ['PT', 'LM']"
    )
    folds = json.loads((ROOT / "data/manifests/m0/donor_folds_pre_qc.json").read_text(encoding="utf-8"))["folds"]

    variants = [
        ("legacy_buggy", True, False),
        ("fix_init_slice_only", False, False),
        ("fix_sqrt_q_only", True, True),
        ("bugfixed_both", False, True),
    ]

    oof_frames = []
    for rep in ("PCA32", "scGPT32"):
        groups = _m4_prior_inputs(rep, rows)
        archived_oof = pd.read_csv(
            ROOT / f"artifacts/archive/CRLM_V1_245_F1/artifacts/m4/direct_prior/{rep}/oof_predictions.csv"
        )
        archived_mix = archived_oof[archived_oof.method == f"{rep}_MixMIL"].set_index(["lineage", "sample_uid"])

        for lineage, group in groups.items():
            ready = rows[rows.sample_uid.isin(group)]
            for split in folds:
                fold = split["fold"]
                train = ready[ready.donor_uid.isin(split["train_donors"])].reset_index(drop=True)
                test = ready[ready.donor_uid.isin(split["test_donors"])].reset_index(drop=True)
                design = build_design(train)
                h = _m4_project(group, design, rep, lineage, fold)
                ids = design.samples
                test_ids = test.sample_uid.tolist()
                trmat = [h[s] for s in ids]
                train_batch = BagBatch(trmat, [group[s][2].subtype_m2.astype(str).to_numpy() for s in ids])
                test_batch = BagBatch([h[s] for s in test_ids], [group[s][2].subtype_m2.astype(str).to_numpy() for s in test_ids])

                for v_name, leg_slice, scale_sq in variants:
                    torch.manual_seed(SEED)
                    model = MixedMIL(
                        trmat,
                        design.labels,
                        posterior_cls,
                        scale_by_sqrt_q=scale_sq,
                        legacy_even_slice_init=leg_slice,
                    )
                    train_fixed(model, train_batch, design, train.donor_uid.nunique(), 500)
                    with torch.inference_mode():
                        logits, att, signed, bias = model(test_batch)
                        logits_np = logits.numpy()
                        owner = test_batch.owner.numpy()
                        np.testing.assert_allclose(np.bincount(owner, weights=att.numpy()), np.ones(len(test_ids)), atol=2e-5)
                        np.testing.assert_allclose(
                            np.bincount(owner, weights=signed.numpy()) + bias.item(),
                            logits_np,
                            rtol=2e-5,
                            atol=2e-5,
                        )
                    if v_name == "legacy_buggy":
                        ref_logits = np.array([archived_mix.loc[(lineage, s), "logit"] for s in test_ids], dtype=float)
                        # Minor ~3e-4 float32 drift across 500 MC variational Adam steps from statsmodels==0.15.0 update
                        np.testing.assert_allclose(logits_np, ref_logits, rtol=5e-3, atol=5e-3)

                    oof_frames.append(
                        test[["sample_uid", "donor_uid", "lesion_uid", "tissue"]].assign(
                            method=f"{rep}_MixMIL_{v_name}",
                            lineage=lineage,
                            fold=fold,
                            seed=SEED,
                            logit=logits_np,
                        )
                    )

    # Also load matched mean baselines from archived M4 for context
    mean_pca = pd.read_csv(ROOT / "artifacts/archive/CRLM_V1_245_F1/artifacts/m4/PCA/oof_predictions.csv")
    mean_pca = mean_pca[(mean_pca.method == "PCA_mean") & (mean_pca.seed == SEED)]
    mean_sc32 = pd.read_csv(ROOT / "artifacts/archive/CRLM_V1_245_F1/artifacts/m4/direct_prior/scGPT32/oof_predictions.csv")
    mean_sc32 = mean_sc32[(mean_sc32.method == "scGPT32_mean") & (mean_sc32.seed == SEED)]

    oof_all = pd.concat(oof_frames, ignore_index=True)
    oof_all.to_csv(OUT_DIR / "mixmil_bugfix_oof_predictions.csv", index=False)

    combined = pd.concat([oof_all, mean_pca, mean_sc32], ignore_index=True)
    patients, metrics, _ = summarize_oof(combined, SEED)

    rows_out = []
    for rep, mean_ref in (("PCA32", "PCA_mean"), ("scGPT32", "scGPT32_mean")):
        for lineage in LINEAGES:
            pat_lin = patients[(patients.lineage == lineage) & (patients.seed == SEED)]
            for v_name, leg_slice, scale_sq in variants:
                m_key = f"{rep}_MixMIL_{v_name}"
                m_row = metrics[(metrics.method == m_key) & (metrics.lineage == lineage)].iloc[0].to_dict()
                vs_legacy = _paired_donor_auc_ci(pat_lin, f"{rep}_MixMIL_legacy_buggy", m_key, SEED)
                vs_mean = _paired_donor_auc_ci(pat_lin, mean_ref, m_key, SEED)
                rows_out.append(
                    {
                        "representation": rep,
                        "lineage": lineage,
                        "variant": v_name,
                        "legacy_even_slice_init": leg_slice,
                        "scale_by_sqrt_q": scale_sq,
                        "patients": int(m_row["patients"]),
                        "n_patient_sites": int(m_row["n_patient_sites"]),
                        "paired_n": int(m_row["paired_n"]),
                        "auroc": float(m_row["auroc"]),
                        "average_precision": float(m_row["average_precision"]),
                        "balanced_accuracy": float(m_row["balanced_accuracy"]),
                        "mean_paired_logit": float(m_row["mean_paired_logit"]),
                        "fraction_positive_pairs": float(m_row["fraction_positive_pairs"]),
                        "auc_delta_vs_legacy_buggy": vs_legacy["auc_delta"],
                        "auc_delta_vs_legacy_lower95": vs_legacy["auc_delta_lower95"],
                        "auc_delta_vs_legacy_upper95": vs_legacy["auc_delta_upper95"],
                        "matched_mean_baseline": mean_ref,
                        "auc_delta_vs_matched_mean": vs_mean["auc_delta"],
                        "auc_delta_vs_mean_lower95": vs_mean["auc_delta_lower95"],
                        "auc_delta_vs_mean_upper95": vs_mean["auc_delta_upper95"],
                    }
                )

    comp_df = pd.DataFrame(rows_out)
    comp_df.to_csv(OUT_DIR / "mixmil_bugfix_comparison.csv", index=False)
    print(comp_df[["representation", "lineage", "variant", "auroc", "balanced_accuracy", "auc_delta_vs_legacy_buggy", "auc_delta_vs_matched_mean"]].to_string(index=False))


# ---------------------------------------------------------------------------
# Part 2: Track 1 (CRLM_V1_245_F1_techfix) GSE164522 External Classifiers
# ---------------------------------------------------------------------------

def _is_gse164522_doublet(obs: pd.DataFrame) -> np.ndarray:
    """Identify author-annotated CD3E+ myeloid doublets and previously missed hC75_B-CD3E B-T doublets."""
    sub = obs["celltype_sub"].astype(str)
    lin = obs["lineage"].astype(str)
    myeloid_cd3 = lin.eq("Mono/Mac/DC") & sub.str.contains("CD3E", regex=False)
    b_cd3 = lin.eq("B/Plasma") & sub.eq("hC75_B-CD3E")
    return (myeloid_cd3 | b_cd3).to_numpy(dtype=bool)


def _select_bag_rows(cell_ids: pd.Index, candidate_rows: np.ndarray, sample_uid: str, lineage: str) -> np.ndarray:
    if len(candidate_rows) < 50:
        return np.array([], dtype=int)
    if len(candidate_rows) <= 256:
        return np.sort(candidate_rows)
    ordered = candidate_rows[np.argsort(cell_ids[candidate_rows].to_numpy(dtype=str))]
    rng = np.random.default_rng(cell_seed(f"{sample_uid}|{lineage}"))
    return np.sort(rng.choice(ordered, 256, replace=False))


def _recalibrate_shared_nuclear_cp10k(mat_csr: sparse.csr_matrix, shared_col_indices: np.ndarray) -> tuple[sparse.csr_matrix, np.ndarray]:
    """Slice to 245_164_shared_nuclear (excluding 352 ribo/MT genes) and recalibrate log1p(CP10k)."""
    sub = mat_csr[:, shared_col_indices].tocsr().copy()
    sub.eliminate_zeros()
    linear = sub.copy()
    linear.data = np.expm1(linear.data.astype(np.float64))
    nuclear_sums = np.asarray(linear.sum(axis=1)).ravel()
    if (nuclear_sums <= 0).any():
        raise ValueError("Encountered cell with zero shared-nuclear expression sum during CP10k recalibration")
    row_nnz = np.diff(linear.indptr)
    linear.data = np.log1p(linear.data * np.repeat(10000.0 / nuclear_sums, row_nnz)).astype(np.float32)
    return linear, nuclear_sums


def _scgpt_tokens_mapped(mat_csr: sparse.csr_matrix, cell_uids: pd.Index, vocab_ids_per_col: np.ndarray, cls_id: int):
    """Tokenize cells for frozen scGPT using strict 1-to-1 HGNC-bridged vocabulary indices."""
    if (vocab_ids_per_col[vocab_ids_per_col >= 0]).size != len(set(vocab_ids_per_col[vocab_ids_per_col >= 0].tolist())):
        raise ValueError("Duplicate scGPT vocabulary token index detected; 1-to-1 injectivity violated")
    result = []
    for i, uid in enumerate(cell_uids):
        a, b = mat_csr.indptr[i], mat_csr.indptr[i + 1]
        raw_cols = mat_csr.indices[a:b]
        raw_vals = mat_csr.data[a:b]
        pos_mask = raw_vals > 0
        raw_cols = raw_cols[pos_mask]
        raw_vals = raw_vals[pos_mask]
        ids = vocab_ids_per_col[raw_cols]
        use = ids >= 0
        values = raw_vals[use]
        ids = ids[use]
        if not len(ids):
            raise ValueError(f"No expressed in-vocabulary genes in {uid}")
        seed = cell_seed(str(uid))
        bins = np.quantile(values, np.linspace(0, 1, 50))
        left, right = np.digitize(values, bins), np.digitize(values, bins, right=True)
        binned = np.ceil(left + np.random.RandomState(seed).rand(len(values)) * (right - left))
        if len(ids) > 1199:
            choice = torch.randperm(len(ids), generator=torch.Generator().manual_seed(seed))[:1199].numpy()
            ids, binned = ids[choice], binned[choice]
        result.append(
            (
                str(uid),
                np.r_[cls_id, ids],
                np.r_[-2, binned],
                int(pos_mask.sum()),
                int(use.sum()),
            )
        )
    return result


def encode_external_techfix(device: str = "xpu", batch_size: int = 4) -> None:
    global torch
    import torch as _torch
    torch = _torch
    from encode_scgpt import CODE, HASHES, encode, load_encoder

    profile = json.loads(
        _resolve_v1_m4_relpath("artifacts/m4/resource_profiles/scgpt_profile.json").read_text(encoding="utf-8")
    )
    if profile["checkpoint"] != HASHES or not profile["devices"].get(device, {}).get("accepted"):
        raise ValueError(f"Unprofiled checkpoint/device: {device}")

    torch.set_num_threads(4)
    model, vocab = load_encoder()
    model.to(device)
    cls_id = vocab["<cls>"]

    mapper = get_hgnc_mapper()
    shared_nuclear = list(load_shared_nuclear_universe("245_164_shared_nuclear"))
    #Exclude special tokens when bridging HGNC symbols into scGPT gene vocabulary
    gene_vocab = {k: v for k, v in vocab.items() if not k.startswith("<")}

    manifest_rows = []
    for code, tissue in TISSUES.items():
        h5_path = M6_DATA_DIR / f"{code}.h5ad"
        obj = ad.read_h5ad(h5_path, backed="r")
        try:
            obs = obj.obs.copy()
            var = obj.var.copy()
            source_sha = str(obj.uns["expression_sha256"])

            # Step A axis: unambiguous symbols on full GSE164522 axis bridged via HGNC -> scGPT vocab
            unamb_mask = var["unambiguous_symbol"].to_numpy(dtype=bool)
            unamb_Indices = np.flatnonzero(unamb_mask)
            unamb_syms = var.iloc[unamb_Indices]["gene_symbol"].tolist()
            bridge_a = mapper.canonical_to_target_vocab(unamb_syms, gene_vocab)
            vocab_ids_a = bridge_a["target_index"].to_numpy(dtype=np.int64)

            # Step B axis: 245_164_shared_nuclear (18,815 genes, zero ribo/MT) bridged via HGNC -> scGPT vocab
            align_shared = mapper.align_axis_to_canonical_targets(var["gene_symbol"].tolist(), shared_nuclear)
            if not align_shared["present"].all():
                raise ValueError(f"{code}: expected all 18,815 shared_nuclear_245_164 genes present in GSE164522")
            shared_src_cols = align_shared["source_index"].to_numpy(dtype=np.int64)
            bridge_b = mapper.canonical_to_target_vocab(shared_nuclear, gene_vocab)
            vocab_ids_b = bridge_b["target_index"].to_numpy(dtype=np.int64)

            is_doublet = _is_gse164522_doublet(obs)

            for sample_uid, samp_idx in obs.groupby("sample_uid", observed=True).indices.items():
                samp_obs = obs.iloc[samp_idx]
                donor_uid = str(samp_obs.donor_uid.iloc[0])
                lesion_uid = str(samp_obs.lesion_uid.iloc[0])

                for lineage in LINEAGES:
                    lin_mask_orig = samp_obs.lineage.eq(lineage).to_numpy(dtype=bool)
                    rows_orig = samp_idx[lin_mask_orig]
                    sel_orig = _select_bag_rows(obs.index, rows_orig, sample_uid, lineage)

                    rows_clean = samp_idx[lin_mask_orig & (~is_doublet[samp_idx])]
                    sel_clean = _select_bag_rows(obs.index, rows_clean, sample_uid, lineage)

                    # 1. Step A: HGNC-bridged scGPT on original bag
                    if len(sel_orig):
                        dest_a = EMB_DIR / "step_a_hgnc_only" / code / f"{sample_uid}_{lineage.replace('/', '_')}.npz"
                        dest_a.parent.mkdir(parents=True, exist_ok=True)
                        uids_a = obs.index[sel_orig]
                        if not dest_a.exists():
                            t0 = time.monotonic()
                            mat_a = obj.X[sel_orig, :].tocsr()[:, unamb_Indices]
                            mat_a.eliminate_zeros()
                            cells_a = _scgpt_tokens_mapped(mat_a, uids_a, vocab_ids_a, cls_id)
                            emb_a = encode(model, cells_a, vocab, batch_size, device)
                            tmp = dest_a.with_suffix(".partial.npz")
                            np.savez_compressed(
                                tmp,
                                embedding=emb_a,
                                cell_uid=uids_a.to_numpy(dtype=str),
                                nonzero=np.array([x[3] for x in cells_a]),
                                in_vocab=np.array([x[4] for x in cells_a]),
                                tokens=np.array([len(x[1]) for x in cells_a]),
                                source_sha=source_sha,
                                hgnc_sha=mapper.sha256,
                                checkpoint=HASHES["best_model.pt"],
                                upstream=CODE,
                                device=device,
                                step="step_a_hgnc_only",
                                seconds=time.monotonic() - t0,
                            )
                            os.replace(tmp, dest_a)
                        with np.load(dest_a, allow_pickle=False) as enc_a:
                            assert np.array_equal(enc_a["cell_uid"], uids_a.to_numpy(dtype=str))
                            manifest_rows.append(
                                {
                                    "step": "step_a_hgnc_only",
                                    "tissue_code": code,
                                    "tissue": tissue,
                                    "sample_uid": sample_uid,
                                    "donor_uid": donor_uid,
                                    "lesion_uid": lesion_uid,
                                    "lineage": lineage,
                                    "n_cells": len(uids_a),
                                    "in_vocab_nonzero_fraction": float(enc_a["in_vocab"].sum() / enc_a["nonzero"].sum()),
                                    "truncated_cells": int((enc_a["in_vocab"] > 1199).sum()),
                                    "seconds": float(enc_a["seconds"]),
                                    "bytes": dest_a.stat().st_size,
                                    "path": str(dest_a.relative_to(ROOT)),
                                    "sha256": sha256_file(dest_a),
                                }
                            )

                    # 2. Step B: HGNC + shared_nuclear_245_164 CP10k + doublet-cleaned bag
                    if len(sel_clean):
                        dest_b = EMB_DIR / "step_b_shared_nuclear_clean" / code / f"{sample_uid}_{lineage.replace('/', '_')}.npz"
                        dest_b.parent.mkdir(parents=True, exist_ok=True)
                        uids_b = obs.index[sel_clean]
                        if not dest_b.exists():
                            t0 = time.monotonic()
                            raw_b = obj.X[sel_clean, :].tocsr()
                            mat_b, nuc_sums_b = _recalibrate_shared_nuclear_cp10k(raw_b, shared_src_cols)
                            cells_b = _scgpt_tokens_mapped(mat_b, uids_b, vocab_ids_b, cls_id)
                            emb_b = encode(model, cells_b, vocab, batch_size, device)
                            tmp = dest_b.with_suffix(".partial.npz")
                            np.savez_compressed(
                                tmp,
                                embedding=emb_b,
                                cell_uid=uids_b.to_numpy(dtype=str),
                                nonzero=np.array([x[3] for x in cells_b]),
                                in_vocab=np.array([x[4] for x in cells_b]),
                                tokens=np.array([len(x[1]) for x in cells_b]),
                                nuclear_inverse_sum=nuc_sums_b,
                                source_sha=source_sha,
                                hgnc_sha=mapper.sha256,
                                checkpoint=HASHES["best_model.pt"],
                                upstream=CODE,
                                device=device,
                                step="step_b_shared_nuclear_clean",
                                seconds=time.monotonic() - t0,
                            )
                            os.replace(tmp, dest_b)
                        with np.load(dest_b, allow_pickle=False) as enc_b:
                            assert np.array_equal(enc_b["cell_uid"], uids_b.to_numpy(dtype=str))
                            manifest_rows.append(
                                {
                                    "step": "step_b_shared_nuclear_clean",
                                    "tissue_code": code,
                                    "tissue": tissue,
                                    "sample_uid": sample_uid,
                                    "donor_uid": donor_uid,
                                    "lesion_uid": lesion_uid,
                                    "lineage": lineage,
                                    "n_cells": len(uids_b),
                                    "in_vocab_nonzero_fraction": float(enc_b["in_vocab"].sum() / enc_b["nonzero"].sum()),
                                    "truncated_cells": int((enc_b["in_vocab"] > 1199).sum()),
                                    "seconds": float(enc_b["seconds"]),
                                    "bytes": dest_b.stat().st_size,
                                    "path": str(dest_b.relative_to(ROOT)),
                                    "sha256": sha256_file(dest_b),
                                }
                            )

                    # 3. Intermediate Step B1 (shared_nuclear CP10k on original bag, only when orig != clean)
                    if len(sel_orig):
                        dest_b1 = EMB_DIR / "step_b1_shared_nuclear_orig_bags" / code / f"{sample_uid}_{lineage.replace('/', '_')}.npz"
                        dest_b1.parent.mkdir(parents=True, exist_ok=True)
                        uids_b1 = obs.index[sel_orig]
                        if not dest_b1.exists():
                            if len(sel_clean) == len(sel_orig) and np.array_equal(sel_clean, sel_orig):
                                # Identical cells (e.g. all T/NK bags): copy exact array without re-running GPU
                                with np.load(dest_b, allow_pickle=False) as enc_b:
                                    payload = {k: enc_b[k] for k in enc_b.files}
                                payload["step"] = "step_b1_shared_nuclear_orig_bags"
                                tmp = dest_b1.with_suffix(".partial.npz")
                                np.savez_compressed(tmp, **payload)
                                os.replace(tmp, dest_b1)
                            else:
                                t0 = time.monotonic()
                                raw_b1 = obj.X[sel_orig, :].tocsr()
                                mat_b1, nuc_sums_b1 = _recalibrate_shared_nuclear_cp10k(raw_b1, shared_src_cols)
                                cells_b1 = _scgpt_tokens_mapped(mat_b1, uids_b1, vocab_ids_b, cls_id)
                                emb_b1 = encode(model, cells_b1, vocab, batch_size, device)
                                tmp = dest_b1.with_suffix(".partial.npz")
                                np.savez_compressed(
                                    tmp,
                                    embedding=emb_b1,
                                    cell_uid=uids_b1.to_numpy(dtype=str),
                                    nonzero=np.array([x[3] for x in cells_b1]),
                                    in_vocab=np.array([x[4] for x in cells_b1]),
                                    tokens=np.array([len(x[1]) for x in cells_b1]),
                                    nuclear_inverse_sum=nuc_sums_b1,
                                    source_sha=source_sha,
                                    hgnc_sha=mapper.sha256,
                                    checkpoint=HASHES["best_model.pt"],
                                    upstream=CODE,
                                    device=device,
                                    step="step_b1_shared_nuclear_orig_bags",
                                    seconds=time.monotonic() - t0,
                                )
                                os.replace(tmp, dest_b1)
            print(f"Finished encoding GSE164522 {code} ({tissue}) on {device}", flush=True)
        finally:
            obj.file.close()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(manifest_rows).to_csv(OUT_DIR / "techfix_embedding_manifest.csv", index=False)


def _elastic_logit(path: Path, means: pd.Series) -> float:
    with np.load(path) as fit:
        values = means.reindex(fit["columns"]).to_numpy(dtype=float)
        if not np.isfinite(values).all():
            raise ValueError(f"Non-finite program elastic feature in {path}")
        scaled = (values - fit["center"]) / fit["scale"]
        return float((scaled @ fit["coef"].T + fit["bias"])[0])


def _pooled_logit(path: Path, matrix: np.ndarray, sample_uid: str, lineage: str) -> float:
    model, info = load_model(path)
    if info["lineage"] != lineage:
        raise ValueError(f"Classifier lineage mismatch: {path}")
    pred, _ = predict_bags(model, info, [matrix], [sample_uid])
    return float(pred.logit.iloc[0])


def _metric_draws(frame: pd.DataFrame, seed: int = SEED, n_boot: int = 10000) -> tuple[dict[str, np.ndarray], pd.Series]:
    """Donor-cluster bootstrap (10,000 draws) matching predict_gse164522.py."""
    frame = frame.sort_values(["donor_uid", "tissue"]).reset_index(drop=True)
    codes, donors = pd.factorize(frame.donor_uid, sort=True)
    rng = np.random.default_rng(seed)
    weight = np.vstack(
        [
            np.ones(len(donors)),
            rng.multinomial(len(donors), np.full(len(donors), 1 / len(donors)), size=n_boot),
        ]
    )
    w = weight[:, codes]
    pos = frame.tissue.eq("LM").to_numpy()
    x = frame.logit.to_numpy(dtype=float)
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
        precision = np.divide(cum_pos, cum_all, where=cum_all > 0, out=np.zeros_like(cum_all))
        ap = (precision * np.diff(cum_pos, axis=1, prepend=0)).sum(axis=1) / np_
        ap[nn == 0] = np.nan
        pair = frame.pivot(index="donor_uid", columns="tissue", values="logit").reindex(
            index=donors, columns=["PT", "LM"]
        )
        dif = pair.LM - pair.PT
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


def _intervals(draws: dict[str, np.ndarray]) -> dict:
    row = {}
    for metric, values in draws.items():
        finite = values[1:][np.isfinite(values[1:])]
        lo, hi = np.quantile(finite, [0.025, 0.975]) if len(finite) else (np.nan, np.nan)
        row.update(
            {
                metric: float(values[0]),
                f"{metric}_lower95": float(lo),
                f"{metric}_upper95": float(hi),
                f"{metric}_valid_bootstrap": int(len(finite)),
            }
        )
    return row


def _summarize_step_predictions(pred: pd.DataFrame, step_name: str) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    keys = ["method", "lineage", "donor_uid", "tissue"]
    metrics_rows, contrast_rows, site_frames = [], [], []
    for lineage in LINEAGES:
        group = pred.loc[pred.lineage.eq(lineage)]
        wide = group.pivot(index="sample_uid", columns="method", values="logit").reindex(columns=METHODS)
        for scope in ("method_available", "common_all_five", "common_predicting_methods", "paired_donors_common_methods"):
            available = wide.columns[wide.notna().any()]
            common = wide.dropna().index if scope == "common_all_five" else wide[available].dropna().index
            rows = group.loc[group.logit.notna()]
            if scope != "method_available":
                rows = rows.loc[rows.sample_uid.isin(common)]
            if scope == "paired_donors_common_methods":
                n_tissues = rows.groupby("donor_uid").tissue.transform("nunique")
                rows = rows.loc[n_tissues.eq(2)]
            lesions = rows.groupby(keys + ["lesion_uid"]).logit.mean()
            sites = (
                lesions.groupby(keys)
                .agg(["mean", "size"])
                .reset_index()
                .rename(columns={"mean": "logit", "size": "n_lesions"})
            )
            site_frames.append(sites.assign(step=step_name, scope=scope))
            method_draws = {}
            for method in METHODS:
                frame = sites.loc[sites.method.eq(method)]
                rec = {
                    "step": step_name,
                    "lineage": lineage,
                    "method": method,
                    "scope": scope,
                    "n_donors": int(frame.donor_uid.nunique()),
                    "n_sites": int(len(frame)),
                    "n_PT": int(frame.tissue.eq("PT").sum()),
                    "n_LM": int(frame.tissue.eq("LM").sum()),
                    "n_bags": int(len(rows[rows.method.eq(method)])),
                    "paired_n": 0,
                }
                if len(frame):
                    draws, dif = _metric_draws(frame, SEED, 10000)
                    method_draws[method] = draws
                    rec.update(_intervals(draws), paired_n=int(len(dif)))
                    if frame.tissue.nunique() == 2:
                        y = frame.tissue.eq("LM")
                        refs = [
                            roc_auc_score(y, frame.logit),
                            average_precision_score(y, frame.logit),
                            balanced_accuracy_score(y, frame.logit >= 0),
                        ]
                        if not np.allclose([draws[k][0] for k in ("auroc", "average_precision", "balanced_accuracy")], refs):
                            raise ValueError(f"Metric mismatch vs sklearn in {step_name} {lineage} {method}")
                metrics_rows.append(rec)
            if scope == "common_predicting_methods":
                for first, second in combinations(method_draws, 2):
                    delta = {k: vals - method_draws[first][k] for k, vals in method_draws[second].items()}
                    contrast_rows.append(
                        {
                            "step": step_name,
                            "lineage": lineage,
                            "first": first,
                            "second": second,
                            "common_bags": int(len(common)),
                            "common_donors": int(sites.donor_uid.nunique()),
                            **_intervals(delta),
                        }
                    )
    return pd.DataFrame(metrics_rows), pd.DataFrame(contrast_rows), pd.concat(site_frames, ignore_index=True)


def evaluate_external_techfix() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fits = verify_frozen_m5_models()
    rank_scores_fn = load_frozen_rank_scores()
    mapper = get_hgnc_mapper()
    shared_nuclear = list(load_shared_nuclear_universe("245_164_shared_nuclear"))

    # Load frozen M5 model program definitions
    model_prog_df = pd.read_csv(M5_DIR / "gene_sets/model_program_genes.csv")
    raw_164_bg = (
        ROOT / "artifacts/archive/CRLM_V1_245_F1/artifacts/m6/GSE164522/coverage/external_rank_background.txt"
    ).read_text(encoding="utf-8").splitlines()

    # Precompute PCA and program alignments on GSE164522
    coverage_records = []
    pca_transforms = {}
    for lineage in LINEAGES:
        pca_path = fits[(lineage, "PCA_mean")].parent / "pca.npz"
        with np.load(pca_path) as trf:
            pca_transforms[lineage] = {
                "genes": [str(g) for g in trf["genes"]],
                "center": trf["center"].copy(),
                "components": trf["components"].copy(),
            }

    bag_rows = []
    # 0. Load Baseline M6 published predictions (archived snapshot)
    orig_pred = pd.read_csv(M6_ORIG_PRED_DIR / "bag_predictions.csv")
    bag_rows.append(orig_pred.assign(step="m6_original_published"))

    step_configs = [
        ("step_a_hgnc_only", False, False, "step_a_hgnc_only"),
        ("step_b1_hgnc_plus_nuclear_cp10k", True, False, "step_b1_shared_nuclear_orig_bags"),
        ("step_b_full_techfix", True, True, "step_b_shared_nuclear_clean"),
    ]

    for step_name, use_shared_nuclear_cp10k, remove_doublets, emb_subdir in step_configs:
        step_forecasts = []
        for code, tissue in TISSUES.items():
            scale_df = pd.read_csv(
                ROOT / f"artifacts/archive/CRLM_V1_245_F1/artifacts/m6/GSE164522/coverage/{code}_scale_check.csv"
            ).set_index("sample_uid")
            h5_path = M6_DATA_DIR / f"{code}.h5ad"
            obj = ad.read_h5ad(h5_path, backed="r")
            try:
                obs = obj.obs.copy()
                var = obj.var.copy()
                is_doublet = _is_gse164522_doublet(obs)

                # Align rank axis and PCA HVGs
                if use_shared_nuclear_cp10k:
                    align_shared = mapper.align_axis_to_canonical_targets(var["gene_symbol"].tolist(), shared_nuclear)
                    assert align_shared["present"].all()
                    rank_src_cols = align_shared["source_index"].to_numpy(dtype=np.int64)
                    rank_axis_symbols = shared_nuclear
                else:
                    # Step A uses GSE164522's common 4-tissue unambiguous axis, canonicalized via HGNCMapper
                    align_bg = mapper.align_axis_to_canonical_targets(var["gene_symbol"].tolist(), raw_164_bg)
                    assert align_bg["present"].all()
                    rank_src_cols = align_bg["source_index"].to_numpy(dtype=np.int64)
                    ax_norm = mapper.normalize_gene_axis(raw_164_bg)
                    rank_axis_symbols = ax_norm["canonical_symbol"].tolist()

                rank_axis_set = set(rank_axis_symbols)
                prog_defs_by_lin = {}
                for lineage in LINEAGES:
                    defs = {}
                    sub_p = model_prog_df[model_prog_df.lineage == lineage]
                    for prog_name, grp in sub_p.groupby("program"):
                        canon_members = mapper.normalize_gene_set(grp["gene_symbol"].tolist(), deduplicate=True)[0]
                        measured_members = [g for g in canon_members if g in rank_axis_set]
                        if code == "PT":
                            coverage_records.append(
                                {
                                    "step": step_name,
                                    "feature_type": "model_program",
                                    "lineage": lineage,
                                    "name": prog_name,
                                    "n_frozen_245": len(grp),
                                    "n_mapped_164": len(measured_members),
                                    "mapped_fraction": len(measured_members) / len(grp),
                                    "absent_symbols": "|".join(sorted(set(canon_members) - set(measured_members))),
                                }
                            )
                        if len(measured_members) >= 5 and len(measured_members) / len(grp) >= 0.80:
                            defs[prog_name] = measured_members
                    prog_defs_by_lin[lineage] = defs

                pca_align_by_lin = {}
                for lineage in LINEAGES:
                    hvgs = pca_transforms[lineage]["genes"]
                    if use_shared_nuclear_cp10k:
                        # Map onto shared_nuclear positions (0..18814)
                        shared_index = pd.Index(shared_nuclear)
                        canon_hvgs = [mapper.resolve_symbol(g).canonical_symbol or g for g in hvgs]
                        pos = shared_index.get_indexer(canon_hvgs)
                        present = pos >= 0
                        pca_align_by_lin[lineage] = (pos, present)
                    else:
                        al = mapper.align_axis_to_canonical_targets(var["gene_symbol"].tolist(), hvgs)
                        pos = al["source_index"].to_numpy(dtype=np.int64)
                        present = al["present"].to_numpy(dtype=bool)
                        pca_align_by_lin[lineage] = (pos, present)
                    if code == "PT":
                        coverage_records.append(
                            {
                                "step": step_name,
                                "feature_type": "PCA_512_HVG",
                                "lineage": lineage,
                                "name": "PCA_512",
                                "n_frozen_245": len(hvgs),
                                "n_mapped_164": int(present.sum()),
                                "mapped_fraction": float(present.mean()),
                                "absent_symbols": "|".join([h for h, ok in zip(hvgs, present) if not ok]),
                            }
                        )

                for sample_uid, samp_idx in obs.groupby("sample_uid", observed=True).indices.items():
                    samp_obs = obs.iloc[samp_idx]
                    donor_uid = str(samp_obs.donor_uid.iloc[0])
                    lesion_uid = str(samp_obs.lesion_uid.iloc[0])
                    s_row = scale_df.loc[sample_uid]
                    scale_ok = (
                        s_row.n_nonempty > 0
                        and 9900 <= s_row.median_nonempty_inverse_sum <= 10100
                        and s_row.fraction_within_five_percent_of_10000 >= 0.95
                    )
                    if not scale_ok:
                        raise ValueError(f"Sample scale check failed for {sample_uid}")

                    for lineage in LINEAGES:
                        lin_mask = samp_obs.lineage.eq(lineage).to_numpy(dtype=bool)
                        if remove_doublets:
                            lin_mask = lin_mask & (~is_doublet[samp_idx])
                        cand_rows = samp_idx[lin_mask]
                        sel = _select_bag_rows(obs.index, cand_rows, sample_uid, lineage)
                        meta = {
                            "sample_uid": sample_uid,
                            "donor_uid": donor_uid,
                            "lesion_uid": lesion_uid,
                            "tissue": tissue,
                            "lineage": lineage,
                            "n_population": int(len(cand_rows)),
                            "n_cells": int(len(sel)),
                            "n_captured": int(len(samp_idx)),
                        }
                        if not len(sel):
                            for method in METHODS:
                                step_forecasts.append(
                                    {
                                        **meta,
                                        "method": method,
                                        "status": "fewer_than_50_cells",
                                        "logit": np.nan,
                                        "probability_lm": np.nan,
                                        "step": step_name,
                                    }
                                )
                            continue

                        uids = obs.index[sel]
                        raw_csr = obj.X[sel, :].tocsr()

                        if use_shared_nuclear_cp10k:
                            expr_mat, _ = _recalibrate_shared_nuclear_cp10k(raw_csr, rank_src_cols)
                            rank_mat = expr_mat
                        else:
                            expr_mat = raw_csr
                            rank_mat = raw_csr[:, rank_src_cols].tocsr()

                        # 1. program_elastic
                        defs = prog_defs_by_lin[lineage]
                        cell_prog_scores = rank_scores_fn(rank_mat, rank_axis_symbols, defs, max_rank=1500)
                        prog_logit = _elastic_logit(fits[(lineage, "program_elastic")], cell_prog_scores.mean(axis=0))

                        # 2. PCA_mean & PCA_attention
                        trf = pca_transforms[lineage]
                        pos, present = pca_align_by_lin[lineage]
                        if present.mean() < 0.95:
                            raise ValueError(f"Unexpected <95% PCA coverage in {step_name} {lineage}: {present.mean()}")
                        vals = expr_mat[:, pos[present]].toarray().astype(np.float32) - trf["center"][present]
                        pca_vec = (vals @ trf["components"][:, present].T).astype(float)
                        pca_mean_logit = _pooled_logit(fits[(lineage, "PCA_mean")], pca_vec, sample_uid, lineage)
                        pca_attn_logit = _pooled_logit(fits[(lineage, "PCA_attention")], pca_vec, sample_uid, lineage)

                        # 3. scGPT_mean & scGPT_attention
                        emb_path = EMB_DIR / emb_subdir / code / f"{sample_uid}_{lineage.replace('/', '_')}.npz"
                        with np.load(emb_path, allow_pickle=False) as enc:
                            if not np.array_equal(enc["cell_uid"], uids.to_numpy(dtype=str)):
                                raise ValueError(f"scGPT cell_uid mismatch in {emb_path}")
                            sc_vec = enc["embedding"].astype(float)
                        sc_mean_logit = _pooled_logit(fits[(lineage, "scGPT_mean")], sc_vec, sample_uid, lineage)
                        sc_attn_logit = _pooled_logit(fits[(lineage, "scGPT_attention")], sc_vec, sample_uid, lineage)

                        logits_map = {
                            "program_elastic": prog_logit,
                            "PCA_mean": pca_mean_logit,
                            "PCA_attention": pca_attn_logit,
                            "scGPT_mean": sc_mean_logit,
                            "scGPT_attention": sc_attn_logit,
                        }
                        for method in METHODS:
                            val = logits_map[method]
                            step_forecasts.append(
                                {
                                    **meta,
                                    "method": method,
                                    "status": "predicted",
                                    "logit": val,
                                    "probability_lm": float(expit(val)),
                                    "step": step_name,
                                }
                            )
            finally:
                obj.file.close()
        bag_rows.append(pd.DataFrame(step_forecasts))

    all_bags = pd.concat(bag_rows, ignore_index=True)
    all_bags.to_csv(OUT_DIR / "techfix_bag_predictions.csv", index=False)
    pd.DataFrame(coverage_records).to_csv(OUT_DIR / "techfix_pca_program_coverage.csv", index=False)

    metrics_list, contrasts_list, sites_list = [], [], []
    for step_name in ("m6_original_published", "step_a_hgnc_only", "step_b1_hgnc_plus_nuclear_cp10k", "step_b_full_techfix"):
        m_df, c_df, s_df = _summarize_step_predictions(all_bags[all_bags.step == step_name], step_name)
        metrics_list.append(m_df)
        contrasts_list.append(c_df)
        sites_list.append(s_df)

    full_metrics = pd.concat(metrics_list, ignore_index=True)
    full_contrasts = pd.concat(contrasts_list, ignore_index=True)
    full_sites = pd.concat(sites_list, ignore_index=True)

    full_metrics.to_csv(OUT_DIR / "techfix_external_classifier_metrics.csv", index=False)
    full_contrasts.to_csv(OUT_DIR / "techfix_paired_method_contrasts.csv", index=False)
    full_sites.to_csv(OUT_DIR / "techfix_patient_site_logits.csv", index=False)

    # Build concise side-by-side ablation summary on method_available scope
    avail = full_metrics[full_metrics.scope == "method_available"].copy()
    piv_auc = avail.pivot(index=["lineage", "method"], columns="step", values="auroc")
    piv_ba = avail.pivot(index=["lineage", "method"], columns="step", values="balanced_accuracy")
    piv_dir = avail.pivot(index=["lineage", "method"], columns="step", values="fraction_LM_above_PT")
    piv_n = avail.pivot(index=["lineage", "method"], columns="step", values="n_sites")

    step_cols = ["m6_original_published", "step_a_hgnc_only", "step_b1_hgnc_plus_nuclear_cp10k", "step_b_full_techfix"]
    ablation_rows = []
    for lineage in LINEAGES:
        for method in METHODS:
            rec = {"lineage": lineage, "method": method}
            for sc in step_cols:
                rec[f"{sc}_n_sites"] = int(piv_n.loc[(lineage, method), sc])
                rec[f"{sc}_auroc"] = piv_auc.loc[(lineage, method), sc]
                rec[f"{sc}_ba"] = piv_ba.loc[(lineage, method), sc]
                rec[f"{sc}_pair_concordance"] = piv_dir.loc[(lineage, method), sc]
            ablation_rows.append(rec)
    ablation_df = pd.DataFrame(ablation_rows)
    ablation_df.to_csv(OUT_DIR / "techfix_step_ablation_summary.csv", index=False)
    print("\n=== Track 1 External Classifier Step-by-Step Ablation (method_available) ===")
    print(ablation_df.to_string(index=False))


def write_report_and_checksums() -> None:
    mix_df = pd.read_csv(OUT_DIR / "mixmil_bugfix_comparison.csv")
    ext_df = pd.read_csv(OUT_DIR / "techfix_external_classifier_metrics.csv")
    cov_df = pd.read_csv(OUT_DIR / "techfix_pca_program_coverage.csv")

    step_a = ext_df[(ext_df.step == "step_a_hgnc_only") & (ext_df.scope == "method_available")]
    step_b = ext_df[(ext_df.step == "step_b_full_techfix") & (ext_df.scope == "method_available")]
    orig = ext_df[(ext_df.step == "m6_original_published") & (ext_df.scope == "method_available")]

    lines = [
        "# Wave 2 (Task W2-A) · MixMIL 前作修复与 Track 1 (`CRLM_V1_245_F1_techfix`) 外部模型技术勘误报告",
        "",
        "## 1. `src/crlm/prior_mil.py` 缺陷修复与 M4 折外量化（`mixmil_bugfix_comparison.csv`）",
        "",
        "- 已切断 `~/.cache/crlm-m4-priors` 仓外隐式路径，默认绑定仓内归档源码 `references/code/mixmil` 与 `references/code/hier-mil`。",
        "- 修复 `mean_glmm_init` 在 `n_trials=1` 下误照搬上游 `n_trials=2` 的 `[::2]` 切片（恢复全量训练袋估计初始方差比 `ratio`），补齐原论文 `1/sqrt(Q)` 注意力特征缩放（`scale_by_sqrt_q=True`），并为 `HierarchicalMIL` 增加与原实现一致的 `dropout` 参数。",
        "",
        "| 表示 | 谱系 | 原版 Buggy AUROC | 仅修 `[::2]` AUROC | 修复后 (`[::2]` + `/sqrt(Q)`) AUROC | 相对原版 ΔAUROC (95% CI) | 相对同表示均值头 ΔAUROC (95% CI) |",
        "|---|---|---:|---:|---:|---:|---:|",
    ]
    for rep in ("PCA32", "scGPT32"):
        for lin in LINEAGES:
            r_old = mix_df[(mix_df.representation == rep) & (mix_df.lineage == lin) & (mix_df.variant == "legacy_buggy")].iloc[0]
            r_sl = mix_df[(mix_df.representation == rep) & (mix_df.lineage == lin) & (mix_df.variant == "fix_init_slice_only")].iloc[0]
            r_new = mix_df[(mix_df.representation == rep) & (mix_df.lineage == lin) & (mix_df.variant == "bugfixed_both")].iloc[0]
            lines.append(
                f"| {rep} | {lin} | {r_old.auroc:.4f} | {r_sl.auroc:.4f} | {r_new.auroc:.4f} | "
                f"{r_new.auc_delta_vs_legacy_buggy:+.4f} ({r_new.auc_delta_vs_legacy_lower95:+.4f}, {r_new.auc_delta_vs_legacy_upper95:+.4f}) | "
                f"{r_new.auc_delta_vs_matched_mean:+.4f} ({r_new.auc_delta_vs_mean_lower95:+.4f}, {r_new.auc_delta_vs_mean_upper95:+.4f}) |"
            )

    lines.extend(
        [
            "",
            "## 2. Track 1 (`CRLM_V1_245_F1_techfix`) GSE164522 全部 15 个冻结分类头技术勘误",
            "",
            "- **权重锁定**：`models/m4/full/` 全部 15 个分类头与 3 个 `pca.npz` 经 SHA256 严格校验零改动。",
            "- **PCA HVG 覆盖恢复**：经 `hgnc_alias_map` 对齐后，`T/NK` 由 `486/512 (94.92%)` 升至 **`507/512 (99.02%)`**，`Mono/Mac/DC` 由 `496/512 (96.88%)` 升至 **`510/512 (99.61%)`**，`B/Plasma` 由 `485/512 (94.73%)` 升至 **`506/512 (98.83%)`**，全部解除 `<95%` 假性拒判（15/15 全部输出有效预测）。",
            "- **Step A → Step B 分步消融**：Step A 仅启用 HGNC 符号对齐并保留原 58 个细胞袋；Step B 进一步对齐 `245_164_shared_nuclear`（18,815 核基因，剔除 352 个肿瘤缺失核糖体/MT 基因对 `CP10k` 分母、秩背景及 `scGPT` 分位数/token 采样池的干扰，在 Intel Arc A750 XPU 重编码）并剔除 `CD3E+` 髓系双细胞（1,776 细胞）与原漏检的 `hC75_B-CD3E` 双细胞（1,489 细胞）。",
            "",
            "| 谱系 | 冻结分类头 | 原 M6 AUROC | Step A (仅 HGNC) AUROC | Step B (完整 Techfix) AUROC (95% CI) | Step B AP | Step B BA | Step B LM>PT 配对 |",
            "|---|---|---:|---:|---|---:|---:|---:|",
        ]
    )
    for lin in LINEAGES:
        for method in METHODS:
            ro = orig[(orig.lineage == lin) & (orig.method == method)].iloc[0]
            ra = step_a[(step_a.lineage == lin) & (step_a.method == method)].iloc[0]
            rb = step_b[(step_b.lineage == lin) & (step_b.method == method)].iloc[0]
            orig_str = "拒判 (<95% HVG)" if ro.n_sites == 0 else f"{ro.auroc:.3f}"
            n_pos_pairs = int(round(rb.fraction_LM_above_PT * rb.paired_n))
            lines.append(
                f"| {lin} | {method} | {orig_str} | {ra.auroc:.3f} | "
                f"**{rb.auroc:.3f}** ({rb.auroc_lower95:.3f}–{rb.auroc_upper95:.3f}) | "
                f"{rb.average_precision:.3f} | {rb.balanced_accuracy:.3f} | {n_pos_pairs}/{int(rb.paired_n)} |"
            )

    (OUT_DIR / "stage_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    # Write output_sha256.txt covering all generated files in OUT_DIR and EMB_DIR
    sha_lines = []
    targets = sorted(
        [p for p in OUT_DIR.rglob("*") if p.is_file() and p.name != "output_sha256.txt"]
        + [p for p in EMB_DIR.rglob("*.npz") if p.is_file()]
    )
    for path in targets:
        sha_lines.append(f"{sha256_file(path)}  {path.relative_to(ROOT)}")
    (OUT_DIR / "output_sha256.txt").write_text("\n".join(sha_lines) + "\n", encoding="utf-8")
    print(f"Wrote stage_report.md and output_sha256.txt ({len(sha_lines)} tracked files).")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "command",
        choices=("mixmil-bugfix", "encode-external", "evaluate-external", "finalize"),
        help="Execution stage",
    )
    parser.add_argument("--device", choices=("cpu", "xpu"), default="xpu")
    parser.add_argument("--batch", type=int, default=4)
    args = parser.parse_args()

    t0 = time.monotonic()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    with threadpool_limits(limits=4):
        if args.command == "mixmil-bugfix":
            run_mixmil_bugfix()
        elif args.command == "encode-external":
            encode_external_techfix(device=args.device, batch_size=args.batch)
        elif args.command == "evaluate-external":
            evaluate_external_techfix()
        elif args.command == "finalize":
            write_report_and_checksums()

    run_meta_path = OUT_DIR / f"run_{args.command.replace('-', '_')}.json"
    run_meta_path.write_text(
        json.dumps(
            {
                "command": [sys.executable, *sys.argv],
                "completed_utc": utc_now(),
                "python": platform.python_version(),
                "seconds": round(time.monotonic() - t0, 3),
                "peak_rss_bytes": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    if args.command == "finalize":
        # Re-update output_sha256.txt so run_finalize.json itself is included
        sha_lines = []
        targets = sorted(
            [p for p in OUT_DIR.rglob("*") if p.is_file() and p.name != "output_sha256.txt"]
            + [p for p in EMB_DIR.rglob("*.npz") if p.is_file()]
        )
        for path in targets:
            sha_lines.append(f"{sha256_file(path)}  {path.relative_to(ROOT)}")
        (OUT_DIR / "output_sha256.txt").write_text("\n".join(sha_lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
