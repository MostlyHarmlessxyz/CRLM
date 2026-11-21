"""Profile and encode the prespecified tumor bags with the frozen whole-human scGPT encoder."""

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
import torch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from crlm.model_bags import cell_seed, sample_bags

CODE = "cebd6fae655b9c585a4807daa3ac31bb764f06b4"
MODEL = ROOT / "models/scGPT-whole-human"
HASHES = {
    "args.json": "c18e075e018140cb8b2d9029387b9de26607a5ce6a8ccabd6ead70cd76b95d60",
    "vocab.json": "acca93d114ca62c3f0f50debbd23e8c87f0714f4737764454f6b2b13f2e8580f",
    "best_model.pt": "6cb5d451ab5c4b33eb673adbe4fddc61d2389df1b89b7651a9fe2e557572b922",
}
OUT = ROOT / "artifacts/m4/resource_profiles"
CACHE = ROOT / "embeddings/m4/scgpt"


def load_encoder():
    for name, expected in HASHES.items():
        with (MODEL / name).open("rb") as f:
            if hashlib.file_digest(f, "sha256").hexdigest() != expected:
                raise ValueError(f"Unexpected whole-human checkpoint member {name}")
    source = MODEL / f"scGPT-{CODE}/scgpt/model"
    # Load the upstream encoder and relative modules directly, without its training utilities.
    spec = importlib.util.spec_from_file_location(
        "pinned_scgpt_encoder",
        source / "model.py",
        submodule_search_locations=[str(source)],
    )
    upstream = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = upstream
    spec.loader.exec_module(upstream)
    vocab = json.loads((MODEL / "vocab.json").read_text())
    args = json.loads((MODEL / "args.json").read_text())
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
    # The released flash-attention checkpoint uses packed QKV weights. These are the upstream CPU aliases.
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


def tokens(bag, vocab):
    gene_ids = np.array([vocab.get(g, -1) for g in bag.genes])
    result = []
    for i, uid in enumerate(bag.obs.index):
        a, b = bag.counts.indptr[i : i + 2]
        ids = gene_ids[bag.counts.indices[a:b]]
        use = ids >= 0
        values = bag.counts.data[a:b][use]
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
                np.r_[vocab["<cls>"], ids],
                np.r_[-2, binned],
                b - a,
                int(use.sum()),
            )
        )
    return result


@torch.inference_mode()
def encode(model, cells, vocab, batch, device):
    output = []
    for start in range(0, len(cells), batch):
        block = cells[start : start + batch]
        length = max(len(c[1]) for c in block)
        ids = np.full((len(block), length), vocab["<pad>"], dtype=np.int64)
        values = np.full(ids.shape, -2, dtype=np.float32)
        for j, (_, token, val, _, _) in enumerate(block):
            ids[j, : len(token)] = token
            values[j, : len(val)] = val
        x = torch.from_numpy(ids).to(device)
        h = model._encode(x, torch.from_numpy(values).to(device), x.eq(vocab["<pad>"]))[
            :, 0, :
        ]
        output.append(torch.nn.functional.normalize(h, dim=1).cpu().numpy())
    result = np.vstack(output)
    if not np.isfinite(result).all() or not np.allclose(
        np.linalg.norm(result, axis=1), 1, atol=2e-5
    ):
        raise ValueError("Invalid CLS embeddings")
    return result


def main(args):
    torch.set_num_threads(4)
    model, vocab = load_encoder()
    OUT.mkdir(parents=True, exist_ok=True)
    if args.profile:
        bags = sample_bags(ROOT / "data/processed/m2/GSM7844847.h5ad")
        cells = [c for bag in bags.values() for c in tokens(bag, vocab)[:4]]
        report = {
            "checkpoint": HASHES,
            "upstream_commit": CODE,
            "torch": torch.__version__,
            "cells": len(cells),
            "ids": [c[0] for c in cells],
            "max_length": 1200,
            "batch_size": args.batch,
            "xpu_available": torch.xpu.is_available(),
            "cuda_available": torch.cuda.is_available(),
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
                    torch.xpu.reset_peak_memory_stats()
                start = time.perf_counter()
                result = encode(model, cells, vocab, args.batch, device)
                seconds = time.perf_counter() - start
                reference = result if reference is None else reference
                changed_batch = encode(
                    model, list(reversed(cells[:4])), vocab, 2, device
                )[::-1]
                cos = float((result * reference).sum(1).min())
                repeated = float((result[:4] * changed_batch).sum(1).min())
                report["devices"][device] = dict(
                    seconds=seconds,
                    cells_per_second=len(cells) / seconds,
                    min_cosine_vs_cpu=cos,
                    min_cosine_rebatched=repeated,
                    accepted=bool(min(cos, repeated) > 0.9999),
                    peak_allocated_bytes=torch.xpu.max_memory_allocated()
                    if device == "xpu"
                    else None,
                )
            except (RuntimeError, NotImplementedError, ValueError) as exc:
                report["devices"][device] = dict(accepted=False, error=str(exc))
        report["process_peak_rss_bytes"] = (
            resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024
        )
        (OUT / "scgpt_profile.json").write_text(json.dumps(report, indent=2) + "\n")
        print(json.dumps(report["devices"], indent=2), flush=True)
        return
    profile = json.loads((OUT / "scgpt_profile.json").read_text())
    if profile["checkpoint"] != HASHES or not profile["devices"].get(
        args.device, {}
    ).get("accepted"):
        raise ValueError("Selected backend has not passed the checkpoint profile")
    CACHE.mkdir(parents=True, exist_ok=True)
    model.to(args.device)
    rows = pd.read_csv(ROOT / "data/manifests/m0/sample_manifest_validated.csv")
    for gsm in rows.loc[rows.tissue.isin(["PT", "LM"]), "sample_uid"]:
        dest = CACHE / f"{gsm}.npz"
        start = time.perf_counter()
        bags = sample_bags(ROOT / f"data/processed/m2/{gsm}.h5ad")
        result = {}
        all_cells = []
        for bag in bags.values():
            all_cells.extend(tokens(bag, vocab))
        ids = np.asarray([c[0] for c in all_cells], dtype=str)
        if dest.exists():
            with np.load(dest, allow_pickle=False) as existing:
                if (
                    not np.array_equal(existing["cell_uid"], ids)
                    or str(existing["source_h5ad_sha256"]) != source_sha(gsm)
                    or str(existing["implementation"]) != implementation()
                    or str(existing["device"]) != args.device
                ):
                    raise ValueError(f"Cache identity mismatch: {dest}")
            continue
        result["embedding"] = (
            encode(model, all_cells, vocab, args.batch, args.device)
            if len(ids)
            else np.empty((0, 512), np.float32)
        )
        result.update(
            cell_uid=ids,
            expressed_counts=np.array([c[3] for c in all_cells]),
            in_vocab=np.array([c[4] for c in all_cells]),
            tokens=np.array([len(c[1]) for c in all_cells]),
            weights_sha256=HASHES["best_model.pt"],
            source_h5ad_sha256=source_sha(gsm),
            implementation=implementation(),
            device=args.device,
            seconds=time.perf_counter() - start,
        )
        tmp = dest.with_suffix(".tmp.npz")
        np.savez_compressed(tmp, **result)
        os.replace(tmp, dest)
        print(gsm, len(ids), round(time.perf_counter() - start, 1), flush=True)


def file_digest(path):
    with path.open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def source_sha(gsm):
    return file_digest(ROOT / f"data/processed/m2/{gsm}.h5ad")


def implementation():
    return ":".join(
        [
            CODE,
            file_digest(Path(__file__)),
            file_digest(ROOT / "src/crlm/model_bags.py"),
        ]
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", action="store_true")
    parser.add_argument("--batch", type=int, default=4)
    parser.add_argument("--device", choices=["cpu", "xpu"], default="cpu")
    main(parser.parse_args())
