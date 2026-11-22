"""Numerically compare adapted HA/CTA/MixMIL pooling to the vendored author implementations (`references/code/`)."""

from __future__ import annotations

import ast
import json
import math
from pathlib import Path
import sys

import numpy as np
import torch

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))
from crlm.prior_mil import (
    BagBatch,
    HierarchicalMIL,
    MixedMIL,
    official_posterior,
    resolve_prior_repo,
)


def _load_hiermil_without_pyg(model_py: Path):
    """Execute upstream `hier-mil/model.py` `Model` class AST using pure-PyTorch `global_add_pool`/`softmax` shims."""
    tree = ast.parse(model_py.read_text(encoding="utf-8"), filename=str(model_py))
    model_cls = next(
        c for c in tree.body if isinstance(c, ast.ClassDef) and c.name == "Model"
    )

    def global_add_pool(x: torch.Tensor, batch: torch.Tensor, size: int | None = None) -> torch.Tensor:
        n = int(batch.max().item()) + 1 if size is None else int(size)
        return x.new_zeros((n, x.shape[1])).index_add(0, batch, x)

    def global_mean_pool(x: torch.Tensor, batch: torch.Tensor, size: int | None = None) -> torch.Tensor:
        n = int(batch.max().item()) + 1 if size is None else int(size)
        sums = global_add_pool(x, batch, size=n)
        counts = torch.bincount(batch, minlength=n).clamp_min(1)[:, None].to(x.dtype)
        return sums / counts

    def pyg_softmax(src: torch.Tensor, index: torch.Tensor, num_nodes: int | None = None) -> torch.Tensor:
        was_1d = src.ndim == 1
        if was_1d:
            src = src[:, None]
        n = int(index.max().item()) + 1 if num_nodes is None else int(num_nodes)
        idx = index[:, None].expand_as(src)
        mx = src.new_full((n, src.shape[1]), -torch.inf)
        mx.scatter_reduce_(0, idx, src.detach(), reduce="amax")
        ex = (src - mx[index]).exp()
        den = global_add_pool(ex, index, size=n)[index]
        out = ex / den
        return out.squeeze(-1) if was_1d else out

    ns = {
        "torch": torch,
        "nn": torch.nn,
        "F": torch.nn.functional,
        "global_add_pool": global_add_pool,
        "global_mean_pool": global_mean_pool,
        "softmax": pyg_softmax,
    }
    exec(compile(ast.Module(body=[model_cls], type_ignores=[]), str(model_py), "exec"), ns)
    return ns["Model"]


def verify() -> None:
    prov = json.loads((ROOT / "references/code/PROVENANCE.json").read_text(encoding="utf-8"))
    revisions = {entry["id"]: entry["commit"] for entry in prov["repositories"]}
    hier_dir = resolve_prior_repo("hier-mil")
    mix_dir = resolve_prior_repo("mixmil")

    upstream_hier_cls = _load_hiermil_without_pyg(hier_dir / "model.py")
    torch.manual_seed(123)
    raw = np.random.default_rng(23).normal(size=(8, 8, 4)).astype("f")
    labels = [np.array(["a"] * 4 + ["b"] * 4) for _ in raw]
    batch = BagBatch(raw, labels)
    differences = {}

    for attention in (False, True):
        ref = upstream_hier_cls(4, attn1=attention, attn2=True, dropout=0.0, n_layers_lin=1, n_hid=32)
        ref.eval()
        model = HierarchicalMIL(4, attention, dropout=0.0)
        model.eval()
        model.transform.load_state_dict(ref.lin[0].state_dict())
        model.type_attention.load_state_dict(ref.w_ct[0].state_dict())
        model.head.load_state_dict(ref.lin_out.state_dict())
        if attention:
            model.cell_attention.load_state_dict(ref.w_c[0].state_dict())
        expected = ref(batch.x, batch.group, 16, 2).flatten()
        differences["HA" if attention else "CTA"] = float((expected - model(batch)[0]).abs().max().item())

    mix_path = mix_dir / "mixmil/model.py"
    module = ast.parse(mix_path.read_text(encoding="utf-8"), filename=str(mix_path))
    cls = next(c for c in module.body if isinstance(c, ast.ClassDef) and c.name == "MixMIL")
    method = next(f for f in cls.body if isinstance(f, ast.FunctionDef) and f.name == "_calc_bag_emb_effect_tensor")
    symbols = {"torch": torch}
    exec(compile(ast.Module(body=[method], type_ignores=[]), str(mix_path), "exec"), symbols)

    model = MixedMIL(
        raw,
        np.array([0, 1] * 4),
        official_posterior(mix_dir),
        scale_by_sqrt_q=True,
        legacy_even_slice_init=False,
    )
    beta = torch.randn(8, 8)
    pooled, _, _, _ = model.raw(batch, beta)
    u, z = beta.split(4)
    eta = z / z.square().mean(0).sqrt()
    raw_scaled = torch.tensor(raw, dtype=torch.float32) / math.sqrt(4.0)
    expected = symbols[method.name](None, u[:, None, :], eta[:, None, :], raw_scaled).squeeze(1)
    differences["MixMIL"] = float((expected - pooled).abs().max().item())

    if max(differences.values()) > 2e-6:
        raise ValueError(f"Upstream equivalence tolerance exceeded: {differences}")
    dest = Path(__file__).parent / "upstream_equivalence.json"
    dest.write_text(
        json.dumps(
            {
                "revisions": revisions,
                "max_absolute_difference": differences,
                "scope": "vendored references/code/; all subtypes observed for HA/CTA; 1/sqrt(Q)-scaled raw bag effects for MixMIL",
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    print(json.dumps(differences, indent=2))


if __name__ == "__main__":
    verify()
