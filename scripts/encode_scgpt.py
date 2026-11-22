"""Profile and encode full QC-retained tumor bags on the shared nuclear universe with frozen whole-human scGPT."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import resource
import sys
import time

import numpy as np
import pandas as pd
from scipy import sparse
import torch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from crlm.hgnc_alias_map import (
    canonical_to_target_vocab,
    get_hgnc_mapper,
    load_shared_nuclear_universe,
)
from crlm.model_bags import Bag, cell_seed, sample_bags, sample_joint_bag

CODE = "cebd6fae655b9c585a4807daa3ac31bb764f06b4"
MODEL = ROOT / "models/scGPT-whole-human"
HASHES = {
    "args.json": "c18e075e018140cb8b2d9029387b9de26607a5ce6a8ccabd6ead70cd76b95d60",
    "vocab.json": "acca93d114ca62c3f0f50debbd23e8c87f0714f4737764454f6b2b13f2e8580f",
    "best_model.pt": "6cb5d451ab5c4b33eb673adbe4fddc61d2389df1b89b7651a9fe2e557572b922",
}
DEFAULT_UNIVERSE = "245_164_shared_nuclear"
INPUT_TRANSFORM = "log1p_cp10k_shared_nuclear_51bin"
PROFILE_CELLS_PER_LINEAGE = 24
OUT = ROOT / "artifacts/m4/resource_profiles"
CACHE = ROOT / "embeddings/m4/scgpt"


def file_digest(path: Path) -> str:
    with Path(path).open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def source_sha(gsm: str) -> str:
    return file_digest(ROOT / f"data/processed/m2/{gsm}.h5ad")


def implementation() -> str:
    return ":".join(
        [
            CODE,
            file_digest(Path(__file__)),
            file_digest(ROOT / "src/crlm/model_bags.py"),
            file_digest(ROOT / "src/crlm/hgnc_alias_map.py"),
        ]
    )


def load_encoder():
    for name, expected in HASHES.items():
        if file_digest(MODEL / name) != expected:
            raise ValueError(f"Unexpected whole-human checkpoint member {name}")
    source = MODEL / f"scGPT-{CODE}/scgpt/model"
    spec = importlib.util.spec_from_file_location(
        "pinned_scgpt_encoder",
        source / "model.py",
        submodule_search_locations=[str(source)],
    )
    upstream = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = upstream
    spec.loader.exec_module(upstream)
    vocab = json.loads((MODEL / "vocab.json").read_text(encoding="utf-8"))
    args = json.loads((MODEL / "args.json").read_text(encoding="utf-8"))
    for symbol in ("<pad>", "<cls>", "<eoc>"):
        if symbol not in vocab:
            vocab[symbol] = len(vocab)
    model = upstream.TransformerModel(
        ntoken=len(vocab),
        d_model=args["embsize"],
        nhead=args["nheads"],
        d_hid=args["d_hid"],
        nlayers=args["nlayers"],
        nlayers_cls=args["n_layers_cls"],
        n_cls=1,
        vocab=vocab,
        dropout=args["dropout"],
        pad_token=args["pad_token"],
        pad_value=args["pad_value"],
        do_mvc=True,
        use_fast_transformer=False,
        pre_norm=False,
    )
    saved = torch.load(MODEL / "best_model.pt", map_location="cpu", weights_only=True)
    converted = {
        k.replace("self_attn._impl.Wqkv.", "self_attn.in_proj_")
        .replace("self_attn.Wqkv.", "self_attn.in_proj_")
        .replace("self_attn._impl.out_proj.", "self_attn.out_proj."): v
        for k, v in saved.items()
    }
    matches = {
        k: v
        for k, v in converted.items()
        if k in model.state_dict() and model.state_dict()[k].shape == v.shape
    }
    missing = model.load_state_dict(matches, strict=False).missing_keys
    if any(
        k.startswith(("encoder.", "value_encoder.", "transformer_encoder."))
        for k in missing
    ):
        raise ValueError(f"Missing tensors used by CLS encoding: {missing}")
    model.requires_grad_(False).eval()
    model.transformer_encoder.use_nested_tensor = False
    return model, vocab


def resolve_vocab_ids(genes, vocab: dict[str, int]) -> tuple[np.ndarray, pd.DataFrame]:
    """Deterministically bridge canonical/raw gene symbols onto scGPT vocabulary with 1-to-1 injectivity."""
    gene_vocab = {k: int(v) for k, v in vocab.items() if not str(k).startswith("<")}
    bridge = canonical_to_target_vocab([str(g) for g in genes], gene_vocab)
    vocab_ids = bridge["target_index"].to_numpy(dtype=np.int64)
    matched = vocab_ids[vocab_ids >= 0]
    if len(matched) != len(set(matched.tolist())):
        raise ValueError("Duplicate scGPT vocabulary token index detected; 1-to-1 injectivity violated")
    return vocab_ids, bridge


def normalize_bag_log1p_cp10k(
    counts,
    *,
    input_is_log1p: bool = False,
    scale_factor: float = 10000.0,
) -> sparse.csr_matrix:
    """Normalize or recalibrate a sparse `(cells x genes)` matrix onto `log1p(CP10k)`."""
    mat = sparse.csr_matrix(counts, copy=True, dtype=np.float64)
    mat.sum_duplicates()
    mat.eliminate_zeros()
    if not np.isfinite(mat.data).all() or (mat.data < 0).any():
        raise ValueError("scGPT expression input must be finite and nonnegative")
    if input_is_log1p:
        mat.data = np.expm1(mat.data)
    row_sums = np.asarray(mat.sum(axis=1)).ravel()
    if (row_sums <= 0).any():
        raise ValueError("Encountered cell with zero total expression on target nuclear universe")
    row_nnz = np.diff(mat.indptr)
    mat.data = np.log1p(mat.data * np.repeat(scale_factor / row_sums, row_nnz))
    return mat.astype(np.float32)


def tokens(
    bag: Bag,
    vocab: dict[str, int],
    *,
    vocab_ids: np.ndarray | None = None,
    input_is_log1p: bool = False,
    recalibrate_cp10k: bool = True,
):
    """Tokenize cells for frozen scGPT using HGNC-bridged injective vocab IDs and 51-bin quantile binning."""
    if vocab_ids is None:
        gene_ids, _ = resolve_vocab_ids(bag.genes, vocab)
    else:
        gene_ids = np.asarray(vocab_ids, dtype=np.int64)
        if len(gene_ids) != len(bag.genes):
            raise ValueError("vocab_ids length does not match bag.genes")
        matched = gene_ids[gene_ids >= 0]
        if len(matched) != len(set(matched.tolist())):
            raise ValueError("Duplicate scGPT token index in supplied vocab_ids")

    mat = (
        normalize_bag_log1p_cp10k(bag.counts, input_is_log1p=input_is_log1p)
        if recalibrate_cp10k
        else sparse.csr_matrix(bag.counts, copy=False)
    )
    cls_id = int(vocab["<cls>"])
    result = []
    for i, uid in enumerate(bag.obs.index):
        a, b = mat.indptr[i : i + 2]
        raw_cols = mat.indices[a:b]
        raw_vals = mat.data[a:b]
        pos = raw_vals > 0
        raw_cols = raw_cols[pos]
        raw_vals = raw_vals[pos]
        ids = gene_ids[raw_cols]
        use = ids >= 0
        values = raw_vals[use]
        ids = ids[use]
        if not len(ids):
            raise ValueError(f"No expressed in-vocabulary genes in {uid}")
        seed = cell_seed(str(uid))
        bins = np.quantile(values, np.linspace(0, 1, 50))
        left, right = np.digitize(values, bins), np.digitize(values, bins, right=True)
        binned = np.ceil(
            left + np.random.RandomState(seed).rand(len(values)) * (right - left)
        )
        if len(ids) > 1199:
            choice = torch.randperm(
                len(ids), generator=torch.Generator().manual_seed(seed)
            )[:1199].numpy()
            ids, binned = ids[choice], binned[choice]
        result.append(
            (
                str(uid),
                np.r_[cls_id, ids],
                np.r_[-2, binned],
                int(pos.sum()),
                int(use.sum()),
            )
        )
    return result


@torch.inference_mode()
def encode(model, cells, vocab, batch: int, device: str) -> np.ndarray:
    if not len(cells):
        return np.empty((0, 512), dtype=np.float32)
    output = []
    pad_id = int(vocab["<pad>"])
    for start in range(0, len(cells), batch):
        block = cells[start : start + batch]
        length = max(len(c[1]) for c in block)
        ids = np.full((len(block), length), pad_id, dtype=np.int64)
        values = np.full(ids.shape, -2.0, dtype=np.float32)
        for j, (_, token, val, _, _) in enumerate(block):
            ids[j, : len(token)] = token
            values[j, : len(val)] = val
        x = torch.from_numpy(ids).to(device)
        h = model._encode(x, torch.from_numpy(values).to(device), x.eq(pad_id))[:, 0, :]
        output.append(torch.nn.functional.normalize(h, dim=1).cpu().numpy().astype(np.float32))
    result = np.vstack(output)
    if not np.isfinite(result).all() or not np.allclose(
        np.linalg.norm(result, axis=1), 1.0, atol=2e-5
    ):
        raise ValueError("Invalid or non-unit-norm CLS embeddings")
    return result


def main(args):
    torch.set_num_threads(4)
    mapper = get_hgnc_mapper()
    universe_genes = list(load_shared_nuclear_universe(args.gene_universe))
    model, vocab = load_encoder()
    vocab_ids, bridge = resolve_vocab_ids(universe_genes, vocab)
    n_matched_vocab = int(bridge["matched"].sum())

    OUT.mkdir(parents=True, exist_ok=True)
    if args.profile:
        bags = sample_bags(
            ROOT / "data/processed/m2/GSM7844847.h5ad",
            min_cells=50,
            max_cells=None,
            gene_universe=universe_genes,
        )
        cells = []
        for lin_bag in bags.values():
            sub_n = min(PROFILE_CELLS_PER_LINEAGE, len(lin_bag.obs))
            sub_bag = Bag(
                lin_bag.obs.iloc[:sub_n].copy(),
                lin_bag.counts[:sub_n, :].tocsr(),
                lin_bag.genes,
                lin_bag.population,
            )
            cells.extend(tokens(sub_bag, vocab, vocab_ids=vocab_ids))
        if len(cells) < 64:
            raise ValueError(f"Expected >= 64 real profile cells, got {len(cells)}")

        rebatch_k = min(16, len(cells))
        report = {
            "checkpoint": HASHES,
            "upstream_commit": CODE,
            "torch": torch.__version__,
            "gene_universe": args.gene_universe,
            "n_universe_genes": len(universe_genes),
            "n_matched_scgpt_vocab": n_matched_vocab,
            "hgnc_sha256": mapper.sha256,
            "input_transform": INPUT_TRANSFORM,
            "cells": len(cells),
            "rebatch_check_cells": rebatch_k,
            "ids": [c[0] for c in cells],
            "max_length": 1200,
            "batch_size": args.batch,
            "xpu_available": bool(torch.xpu.is_available()),
            "cuda_available": bool(torch.cuda.is_available()),
            "devices": {},
        }
        if report["xpu_available"]:
            report["xpu_name"] = torch.xpu.get_device_name(0)
        reference = None
        for device in ["cpu"] + (["xpu"] if report["xpu_available"] else []):
            try:
                model.to(device)
                encode(model, cells[:1], vocab, 1, device)
                if device == "xpu":
                    torch.xpu.synchronize()
                    torch.xpu.reset_peak_memory_stats()
                start = time.perf_counter()
                result = encode(model, cells, vocab, args.batch, device)
                if device == "xpu":
                    torch.xpu.synchronize()
                seconds = time.perf_counter() - start
                reference = result if reference is None else reference
                changed_batch = encode(
                    model, list(reversed(cells[:rebatch_k])), vocab, 2, device
                )[::-1]
                cos_vec = (result * reference).sum(axis=1)
                abs_diff = np.abs(result - reference)
                rebatch_cos = (result[:rebatch_k] * changed_batch).sum(axis=1)
                rebatch_abs = np.abs(result[:rebatch_k] - changed_batch)

                min_cos_cpu = float(cos_vec.min())
                mean_cos_cpu = float(cos_vec.mean())
                max_abs_cpu = float(abs_diff.max())
                mean_abs_cpu = float(abs_diff.mean())
                min_cos_reb = float(rebatch_cos.min())
                max_abs_reb = float(rebatch_abs.max())

                accepted = bool(
                    min(min_cos_cpu, min_cos_reb) > 0.9999
                    and max(max_abs_cpu, max_abs_reb) < 1e-3
                )
                report["devices"][device] = dict(
                    seconds=seconds,
                    cells_per_second=len(cells) / seconds,
                    min_cosine_vs_cpu=min_cos_cpu,
                    mean_cosine_vs_cpu=mean_cos_cpu,
                    max_abs_diff_vs_cpu=max_abs_cpu,
                    mean_abs_diff_vs_cpu=mean_abs_cpu,
                    min_cosine_rebatched=min_cos_reb,
                    max_abs_diff_rebatched=max_abs_reb,
                    accepted=accepted,
                    peak_allocated_bytes=int(torch.xpu.max_memory_allocated())
                    if device == "xpu"
                    else None,
                )
            except (RuntimeError, NotImplementedError, ValueError) as exc:
                report["devices"][device] = dict(accepted=False, error=str(exc))
        report["process_peak_rss_bytes"] = int(
            resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024
        )
        (OUT / "scgpt_profile.json").write_text(
            json.dumps(report, indent=2) + "\n", encoding="utf-8"
        )
        print(json.dumps(report["devices"], indent=2), flush=True)
        return

    profile = json.loads((OUT / "scgpt_profile.json").read_text(encoding="utf-8"))
    if (
        profile.get("checkpoint") != HASHES
        or profile.get("gene_universe") != args.gene_universe
        or not profile.get("devices", {}).get(args.device, {}).get("accepted")
    ):
        raise ValueError(
            f"Selected backend ({args.device}) or universe ({args.gene_universe}) has not passed scgpt_profile.json"
        )

    CACHE.mkdir(parents=True, exist_ok=True)
    model.to(args.device)
    impl_id = implementation()
    rows = pd.read_csv(ROOT / "data/manifests/m0/sample_manifest_validated.csv")
    tumor_rows = rows.loc[
        rows["role"].eq("development") & rows["tissue"].isin(["PT", "LM"])
    ].reset_index(drop=True)
    if set(tumor_rows["tissue"]) != {"PT", "LM"}:
        raise ValueError("Expected verified development PT and LM rows only")

    for gsm in tumor_rows["sample_uid"]:
        dest = CACHE / f"{gsm}.npz"
        start = time.perf_counter()
        joint_bag = sample_joint_bag(
            ROOT / f"data/processed/m2/{gsm}.h5ad",
            min_cells=1,
            max_cells=None,
            gene_universe=universe_genes,
        )
        all_cells = (
            tokens(joint_bag, vocab, vocab_ids=vocab_ids)
            if joint_bag is not None
            else []
        )
        ids = np.asarray([c[0] for c in all_cells], dtype=str)
        src_digest = source_sha(gsm)

        if dest.exists() and not args.overwrite:
            with np.load(dest, allow_pickle=False) as existing:
                same = (
                    "gene_universe" in existing.files
                    and str(existing["gene_universe"]) == args.gene_universe
                    and str(existing["hgnc_sha256"]) == mapper.sha256
                    and np.array_equal(existing["cell_uid"], ids)
                    and str(existing["source_h5ad_sha256"]) == src_digest
                    and str(existing["implementation"]) == impl_id
                    and str(existing["device"]) == args.device
                )
                if not same:
                    raise ValueError(
                        f"Cache identity mismatch at {dest}; pass --overwrite to regenerate for {args.gene_universe}"
                    )
            continue

        embedding = encode(model, all_cells, vocab, args.batch, args.device)
        if args.device == "xpu":
            torch.xpu.synchronize()
        lin_arr = (
            joint_bag.obs["lineage_m2"].to_numpy(dtype=str)
            if joint_bag is not None
            else np.empty(0, dtype=str)
        )
        sub_arr = (
            joint_bag.obs["subtype_m2"].to_numpy(dtype=str)
            if joint_bag is not None
            else np.empty(0, dtype=str)
        )
        fine_arr = (
            joint_bag.obs["fine_subtype"].to_numpy(dtype=str)
            if joint_bag is not None
            else np.empty(0, dtype=str)
        )
        payload = dict(
            embedding=embedding,
            cell_uid=ids,
            lineage_m2=lin_arr,
            subtype_m2=sub_arr,
            fine_subtype=fine_arr,
            expressed_counts=np.array([c[3] for c in all_cells], dtype=np.int32),
            in_vocab=np.array([c[4] for c in all_cells], dtype=np.int32),
            tokens=np.array([len(c[1]) for c in all_cells], dtype=np.int32),
            gene_universe=args.gene_universe,
            n_universe_genes=len(universe_genes),
            n_matched_scgpt_vocab=n_matched_vocab,
            hgnc_sha256=mapper.sha256,
            input_transform=INPUT_TRANSFORM,
            weights_sha256=HASHES["best_model.pt"],
            source_h5ad_sha256=src_digest,
            implementation=impl_id,
            device=args.device,
            seconds=time.perf_counter() - start,
        )
        tmp = dest.with_suffix(".tmp.npz")
        np.savez_compressed(tmp, **payload)
        os.replace(tmp, dest)
        print(gsm, len(ids), round(payload["seconds"], 2), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", action="store_true")
    parser.add_argument("--batch", type=int, default=16)
    parser.add_argument("--device", choices=["cpu", "xpu"], default="cpu")
    parser.add_argument("--gene-universe", default=DEFAULT_UNIVERSE)
    parser.add_argument("--overwrite", action="store_true")
    main(parser.parse_args())
