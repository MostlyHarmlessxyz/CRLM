"""Direct HA, CTA and full-covariance Gaussian variational MixMIL comparators on full M4 bags.

Track 2 (`V1.1`) upgrades:
1. Eliminates implicit `~/.cache/crlm-m4-priors` dependency; verifies SHA256-locked
   in-repository upstream archives via `references/code/PROVENANCE.json`.
2. Evaluates both single-lineage full bags (`T/NK`, `Mono/Mac/DC`, `B/Plasma` with clean
   fine subtype hierarchy) AND whole-specimen multi-lineage bags (`All3Lineages` across
   all 11 clean immune cell types), faithful to Do & Lähdesmäki (2025) specimen-level HA/CTA.
3. Performs strict fold-internal 3-fold donor-grouped inner CV tuning for `HA` and `CTA`
   over `dropout in (0.0, 0.3)`, `weight_decay in (1e-3, 1e-2)`, `steps in (100, 250)`
   with inner-train-only feature projection/scaling (zero outer-test leakage).
4. Uses remediated `MixedMIL(..., scale_by_sqrt_q=True, legacy_even_slice_init=False)`.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import resource
import sys
import time

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.special import expit
from sklearn.decomposition import PCA
from sklearn.metrics import roc_auc_score
import torch
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src"), str(ROOT / "starter")]

from crlm.hgnc_alias_map import load_shared_nuclear_universe
from crlm.model_bags import (
    JOINT_LINEAGE,
    LINEAGES,
    SEED,
    cell_seed,
    sample_bags,
    sample_joint_bag,
)
from crlm.prior_mil import (
    DEFAULT_MIXMIL_DIR,
    BagBatch,
    HierarchicalMIL,
    MixedMIL,
    official_posterior,
    resolve_prior_repo,
    train_fixed,
)
from linear_attention_mil import (
    build_design,
    fit_model,
    fit_scaler,
    predict_bags,
    save_model,
)

OUT = ROOT / "artifacts/m4/direct_prior"
MODELS = ROOT / "models/m4_direct"
PROVENANCE_JSON = ROOT / "references/code/PROVENANCE.json"
DEFAULT_UNIVERSE = "245_164_shared_nuclear"
MIXMIL_REV = "16d92349ccf19b73c1dccd167934fcc3202b6be4"
HIERMIL_REV = "2fb155840394961595755bede506c4b54c0942ae"
TASKS = (*LINEAGES, JOINT_LINEAGE)

HA_CTA_DROPOUTS = (0.0, 0.3)
HA_CTA_WEIGHT_DECAYS = (1e-3, 1e-2)
HA_CTA_STEPS = (100, 250)


def _sha256_file(path: Path) -> str:
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def verify_prior_provenance(mixmil_repo: Path) -> dict[str, str]:
    """Verify vendored MixMIL and hier-mil sources against `references/code/PROVENANCE.json`."""
    mix_dir = resolve_prior_repo("mixmil", mixmil_repo)
    hier_dir = resolve_prior_repo("hier-mil")
    if not PROVENANCE_JSON.is_file():
        raise FileNotFoundError(f"Missing vendor provenance manifest: {PROVENANCE_JSON}")
    prov = json.loads(PROVENANCE_JSON.read_text(encoding="utf-8"))
    by_id = {entry["id"]: entry for entry in prov.get("repositories", [])}

    mix_meta = by_id.get("mixmil", {})
    hier_meta = by_id.get("hier-mil", {})
    if mix_meta.get("commit") != MIXMIL_REV:
        raise ValueError(f"Expected MixMIL commit {MIXMIL_REV}, got {mix_meta.get('commit')}")
    if hier_meta.get("commit") != HIERMIL_REV:
        raise ValueError(f"Expected hier-mil commit {HIERMIL_REV}, got {hier_meta.get('commit')}")

    post_sha = _sha256_file(mix_dir / "mixmil/posterior.py")
    if post_sha != mix_meta["key_file_sha256"]["mixmil/posterior.py"]:
        raise ValueError("SHA256 mismatch on vendored mixmil/posterior.py")
    hmod_sha = _sha256_file(hier_dir / "model.py")
    if hmod_sha != hier_meta["key_file_sha256"]["model.py"]:
        raise ValueError("SHA256 mismatch on vendored hier-mil/model.py")
    return {
        "mixmil_commit": MIXMIL_REV,
        "hiermil_commit": HIERMIL_REV,
        "mixmil_posterior_sha256": post_sha,
        "hiermil_model_sha256": hmod_sha,
    }


def _normalize_csr_log1p_cp10k(counts: sparse.csr_matrix) -> sparse.csr_matrix:
    x = counts.astype(np.float32, copy=True)
    x.sum_duplicates()
    x.eliminate_zeros()
    sums = np.asarray(x.sum(axis=1)).ravel()
    if (sums <= 0).any():
        raise ValueError("Encountered cell with zero total counts on shared nuclear universe")
    x.data *= np.repeat(10000.0 / sums, np.diff(x.indptr)).astype(np.float32)
    np.log1p(x.data, out=x.data)
    return x


def load_inputs(rep: str, rows: pd.DataFrame, *, gene_universe: str = DEFAULT_UNIVERSE):
    """Load single-lineage and `All3Lineages` full bags on `gene_universe`."""
    universe_genes = list(load_shared_nuclear_universe(gene_universe))
    sel_path = ROOT / "artifacts/m4/selected_cells.csv.gz"
    published = pd.read_csv(sel_path, index_col="cell_uid") if sel_path.is_file() else None

    groups: dict[str, dict[str, tuple]] = {k: {} for k in TASKS}
    for gsm in rows["sample_uid"]:
        h5_path = ROOT / f"data/processed/m2/{gsm}.h5ad"
        single_bags = sample_bags(
            h5_path, min_cells=50, max_cells=None, gene_universe=universe_genes
        )
        joint_bag = sample_joint_bag(
            h5_path, min_cells=50, max_cells=None, gene_universe=universe_genes
        )
        cols, enc = None, None
        if rep.startswith("scGPT"):
            with np.load(ROOT / f"embeddings/m4/scgpt/{gsm}.npz", allow_pickle=False) as f:
                cols = pd.Index(f["cell_uid"])
                enc = f["embedding"]
                if not cols.is_unique:
                    raise ValueError(f"{gsm}: duplicate cell_uid in scGPT cache")

        for lineage, bag in single_bags.items():
            if published is not None:
                expected = published[
                    (published["sample_uid"] == gsm) & (published["lineage_m2"] == lineage)
                ]
                if not bag.obs.index.equals(expected.index):
                    raise ValueError(f"{gsm} {lineage}: cell selection differs from M4 selected_cells.csv.gz")
            if rep == "PCA32":
                x = _normalize_csr_log1p_cp10k(bag.counts)
            else:
                ix = cols.get_indexer(bag.obs.index)
                if (ix < 0).any():
                    raise ValueError(f"{gsm} {lineage}: missing scGPT cell identities")
                x = enc[ix]
            obs_sub = bag.obs[["sample_uid", "subtype_m2", "fine_subtype", "fine_cell_type"]].copy()
            obs_sub["hierarchy_type"] = obs_sub["fine_subtype"].astype(str)
            groups[lineage][gsm] = (x, bag.genes, obs_sub)

        if joint_bag is not None:
            if rep == "PCA32":
                xj = _normalize_csr_log1p_cp10k(joint_bag.counts)
            else:
                ixj = cols.get_indexer(joint_bag.obs.index)
                if (ixj < 0).any():
                    raise ValueError(f"{gsm} {JOINT_LINEAGE}: missing scGPT cell identities")
                xj = enc[ixj]
            obs_joint = joint_bag.obs[
                ["sample_uid", "lineage_m2", "subtype_m2", "fine_subtype", "fine_cell_type"]
            ].copy()
            obs_joint["hierarchy_type"] = obs_joint["fine_cell_type"].astype(str)
            groups[JOINT_LINEAGE][gsm] = (xj, joint_bag.genes, obs_joint)

    return groups


def project_subset(
    group: dict,
    design,
    target_samples: list[str],
    rep: str,
    lineage: str,
    *,
    outer_fold: int | None = None,
    destination: Path | None = None,
) -> dict[str, np.ndarray]:
    """Fit projection & patient-weighted scaling strictly on `design.samples` and transform `target_samples`."""
    eval_samples = list(dict.fromkeys([*design.samples, *target_samples]))
    if rep == "PCA32":
        saved_m4 = (
            ROOT / f"models/m4/{lineage.replace('/', '_')}/fold{outer_fold}/pca.npz"
            if (outer_fold is not None and lineage in LINEAGES)
            else None
        )
        if saved_m4 is not None and saved_m4.is_file():
            with np.load(saved_m4, allow_pickle=False) as f:
                if set(f["training_samples"]) != set(design.samples):
                    raise ValueError("Saved M4 PCA has a different training fold")
                genes, center, components = f["genes"], f["center"], f["components"]
            matrix = {}
            for gsm in eval_samples:
                x, symbols, _ = group[gsm]
                indices = symbols.get_indexer(genes)
                if (indices < 0).any():
                    raise ValueError("Missing a fitted PCA feature")
                matrix[gsm] = np.asarray(
                    (x[:, indices].toarray() - center) @ components.T, dtype=np.float32
                )
        else:
            all_genes = group[design.samples[0]][1]
            x_tr = sparse.vstack([group[s][0] for s in design.samples], format="csr")
            var = (
                np.asarray(x_tr.power(2).mean(axis=0)).ravel()
                - np.asarray(x_tr.mean(axis=0)).ravel() ** 2
            )
            idx = np.sort(np.argsort(-var, kind="stable")[:512])
            pca = PCA(n_components=32, svd_solver="randomized", random_state=SEED)
            pca.fit(x_tr[:, idx].toarray())
            if destination is not None:
                np.savez_compressed(
                    destination / "pca.npz",
                    genes=np.asarray(all_genes[idx], dtype=str),
                    center=pca.mean_,
                    components=pca.components_,
                    training_samples=np.asarray(design.samples, dtype=str),
                    gene_universe=DEFAULT_UNIVERSE,
                )
            matrix = {
                gsm: np.asarray(pca.transform(group[gsm][0][:, idx].toarray()), dtype=np.float32)
                for gsm in eval_samples
            }
    elif rep == "scGPT32":
        xs = np.concatenate([group[s][0] for s in design.samples]).astype(np.float64)
        center = xs.mean(axis=0)
        xs -= center
        _, vectors = np.linalg.eigh(xs.T @ xs)
        components = vectors[:, -32:][:, ::-1].T
        signs = np.sign(
            components[np.arange(len(components)), np.abs(components).argmax(axis=1)]
        )
        components *= signs[:, None]
        matrix = {
            gsm: np.asarray((group[gsm][0] - center) @ components.T, dtype=np.float32)
            for gsm in eval_samples
        }
        if destination is not None:
            np.savez_compressed(
                destination / "projection.npz",
                center=center,
                components=components,
                training_samples=np.asarray(design.samples, dtype=str),
            )
    else:
        matrix = {gsm: np.asarray(group[gsm][0], dtype=np.float32) for gsm in eval_samples}

    loc, scale = fit_scaler(
        [matrix[s].astype(float) for s in design.samples], design.bce_weights, True
    )
    if destination is not None:
        np.savez_compressed(
            destination / "scaler.npz",
            center=loc,
            scale=scale,
            training_samples=np.asarray(design.samples, dtype=str),
        )
    return {s: np.asarray((a - loc) / scale, dtype=np.float32) for s, a in matrix.items()}


def make_inner_donor_splits(
    train_rows: pd.DataFrame,
    *,
    n_splits: int = 3,
    seed_key: str = "",
) -> list[tuple[pd.DataFrame, pd.DataFrame]]:
    """Construct `n_splits` deterministic donor-disjoint inner CV splits stratified by site pairing."""
    donor_sites = train_rows.groupby("donor_uid")["tissue"].agg(lambda s: set(s.tolist())).to_dict()
    paired = sorted(d for d, ts in donor_sites.items() if ts == {"PT", "LM"})
    pt_only = sorted(d for d, ts in donor_sites.items() if ts == {"PT"})
    lm_only = sorted(d for d, ts in donor_sites.items() if ts == {"LM"})

    rng = np.random.default_rng(cell_seed(f"inner_splits|{seed_key}"))
    buckets: list[set[str]] = [set() for _ in range(n_splits)]
    cursor = 0
    for stratum in (paired, pt_only, lm_only):
        order = rng.permutation(stratum).tolist() if len(stratum) > 1 else list(stratum)
        for d in order:
            buckets[cursor % n_splits].add(str(d))
            cursor += 1

    splits = []
    all_donors = set(donor_sites)
    for val_donors in buckets:
        if not val_donors:
            continue
        tr_donors = all_donors - val_donors
        if tr_donors & val_donors:
            raise ValueError("Inner CV donor leakage")
        in_tr = train_rows[train_rows["donor_uid"].isin(tr_donors)].reset_index(drop=True)
        in_va = train_rows[train_rows["donor_uid"].isin(val_donors)].reset_index(drop=True)
        if in_tr["tissue"].nunique() == 2 and len(in_va) > 0:
            splits.append((in_tr, in_va))
    if len(splits) < 2:
        raise ValueError("Insufficient donor support for inner CV tuning")
    return splits


def _batch_to_device(batch: BagBatch, device: str) -> BagBatch:
    if device == "cpu":
        return batch
    moved = object.__new__(BagBatch)
    moved.x = batch.x.to(device)
    moved.n = batch.n
    moved.owner = batch.owner.to(device)
    moved.group = batch.group.to(device)
    moved.group_owner = batch.group_owner.to(device)
    return moved


def _seed_all(seed: int) -> None:
    torch.manual_seed(seed)
    if torch.xpu.is_available():
        torch.xpu.manual_seed_all(seed)


def _train_hier_on_device(
    model: HierarchicalMIL,
    bags: BagBatch,
    design,
    steps: int,
    *,
    optim: torch.optim.Optimizer,
    lambda_pair: float = 0.1,
) -> float:
    dev = bags.x.device
    weights = torch.tensor(design.bce_weights, dtype=torch.float32, device=dev)
    labels = torch.tensor(design.labels, dtype=torch.float32, device=dev)
    contrasts = torch.tensor(design.pair_matrix, dtype=torch.float32, device=dev)
    model.train()
    loss = torch.tensor(0.0, device=dev)
    for _ in range(steps):
        logits = model(bags)[0][:, None]
        bce = torch.nn.functional.binary_cross_entropy_with_logits(
            logits, labels[:, None].expand_as(logits), reduction="none"
        )
        pair_term = (
            float(lambda_pair) * torch.nn.functional.softplus(-(contrasts @ logits)).mean()
            if len(contrasts) and lambda_pair > 0
            else 0.0
        )
        loss = (weights @ bce).mean() + pair_term
        if not torch.isfinite(loss):
            raise FloatingPointError("Non-finite training loss")
        optim.zero_grad()
        loss.backward()
        optim.step()
    model.eval()
    return float(loss.detach().cpu().item())


def prepare_inner_cv_batches(
    group: dict,
    train_rows: pd.DataFrame,
    rep: str,
    lineage: str,
    outer_fold: int,
    *,
    device: str = "cpu",
) -> list[tuple[object, pd.DataFrame, BagBatch, BagBatch, int]]:
    """Precompute strictly fold-isolated inner-CV projections and `BagBatch` objects."""
    raw_splits = make_inner_donor_splits(
        train_rows, n_splits=3, seed_key=f"{lineage}|fold{outer_fold}"
    )
    prepared = []
    for in_tr, in_va in raw_splits:
        in_design = build_design(in_tr)
        val_ids = in_va["sample_uid"].tolist()
        h_in = project_subset(
            group, in_design, val_ids, rep, lineage, outer_fold=None, destination=None
        )
        tr_batch = _batch_to_device(
            BagBatch(
                [h_in[s] for s in in_design.samples],
                [group[s][2]["hierarchy_type"].to_numpy(dtype=str) for s in in_design.samples],
            ),
            device,
        )
        va_batch = _batch_to_device(
            BagBatch(
                [h_in[s] for s in val_ids],
                [group[s][2]["hierarchy_type"].to_numpy(dtype=str) for s in val_ids],
            ),
            device,
        )
        prepared.append((in_design, in_va, tr_batch, va_batch, in_tr["donor_uid"].nunique()))
    return prepared


def tune_hierarchical_mil(
    width: int,
    cell_attention: bool,
    inner_batches: list[tuple[object, pd.DataFrame, BagBatch, BagBatch, int]],
    *,
    method: str,
    lineage: str,
    outer_fold: int,
    device: str = "cpu",
) -> tuple[dict[str, object], list[dict[str, object]]]:
    """Select `(dropout, weight_decay, steps)` strictly by inner-CV patient-site OOF AUROC."""
    step1, step2 = HA_CTA_STEPS
    records = []
    for dropout in HA_CTA_DROPOUTS:
        for wd in HA_CTA_WEIGHT_DECAYS:
            oof_step1, oof_step2 = [], []
            for k_idx, (in_design, in_va, tr_b, va_b, _) in enumerate(inner_batches):
                _seed_all(SEED + k_idx)
                model = HierarchicalMIL(width, cell_attention, dropout=dropout).to(device)
                optim = torch.optim.Adam(model.parameters(), lr=1e-3, weight_decay=float(wd))
                lam_p = 0.1 if in_design.pair_donors else 0.0
                _train_hier_on_device(
                    model,
                    tr_b,
                    in_design,
                    step1,
                    optim=optim,
                    lambda_pair=lam_p,
                )
                with torch.inference_mode():
                    logits_1 = model(va_b)[0].detach().cpu().numpy()
                oof_step1.append(
                    in_va[["sample_uid", "donor_uid", "lesion_uid", "tissue"]].assign(
                        logit=logits_1
                    )
                )

                _train_hier_on_device(
                    model,
                    tr_b,
                    in_design,
                    step2 - step1,
                    optim=optim,
                    lambda_pair=lam_p,
                )
                with torch.inference_mode():
                    logits_2 = model(va_b)[0].detach().cpu().numpy()
                oof_step2.append(
                    in_va[["sample_uid", "donor_uid", "lesion_uid", "tissue"]].assign(
                        logit=logits_2
                    )
                )

            for steps_val, frames in ((step1, oof_step1), (step2, oof_step2)):
                df_oof = pd.concat(frames, ignore_index=True)
                lesions = df_oof.groupby(
                    ["donor_uid", "tissue", "lesion_uid"], as_index=False
                )["logit"].mean()
                sites = lesions.groupby(["donor_uid", "tissue"], as_index=False)[
                    "logit"
                ].mean()
                y = sites["tissue"].eq("LM").to_numpy(dtype=float)
                scores = sites["logit"].to_numpy(dtype=float)
                auroc = float(roc_auc_score(y, scores)) if len(np.unique(y)) == 2 else -1.0
                bce = float(np.mean(np.logaddexp(0.0, scores) - y * scores))
                records.append(
                    {
                        "method": method,
                        "lineage": lineage,
                        "fold": outer_fold,
                        "dropout": float(dropout),
                        "weight_decay": float(wd),
                        "steps": int(steps_val),
                        "inner_oof_auroc": auroc,
                        "inner_oof_bce": bce,
                    }
                )

    # Tie-break toward higher AUROC, lower BCE, stronger dropout, stronger weight_decay, fewer steps
    best = max(
        records,
        key=lambda r: (
            round(r["inner_oof_auroc"], 8),
            -round(r["inner_oof_bce"], 8),
            r["dropout"],
            r["weight_decay"],
            -r["steps"],
        ),
    )
    return best, records


def run(rep: str, repo: Path, *, tune_inner: bool = True, device: str | None = None):
    start = time.perf_counter()
    hier_dev = (
        ("xpu" if torch.xpu.is_available() else "cpu")
        if device is None
        else str(device)
    )
    prov_info = verify_prior_provenance(repo)
    posterior = official_posterior(repo)
    out = OUT / rep
    out.mkdir(parents=True, exist_ok=True)

    rows = (
        pd.read_csv(ROOT / "data/manifests/m0/sample_manifest_validated.csv")
        .query("role == 'development' and tissue in ['PT', 'LM']")
        .reset_index(drop=True)
    )
    folds = json.loads((ROOT / "data/manifests/m0/donor_folds_pre_qc.json").read_text(encoding="utf-8"))[
        "folds"
    ]
    groups = load_inputs(rep, rows)
    torch.set_num_threads(4)

    runs, predictions, cell_parts, tuning_rows = [], [], [], []
    for lineage, group in groups.items():
        ready = rows[rows["sample_uid"].isin(group)].reset_index(drop=True)
        for split in folds:
            fold = split["fold"]
            train = ready[ready["donor_uid"].isin(split["train_donors"])].reset_index(drop=True)
            test = ready[ready["donor_uid"].isin(split["test_donors"])].reset_index(drop=True)
            if set(train["donor_uid"]) & set(test["donor_uid"]) or set(train["sample_uid"]) & set(
                test["sample_uid"]
            ):
                raise ValueError("Patient or sample leakage")
            design = build_design(train)
            if not design.pair_donors or test.empty:
                raise ValueError("A study fold lacks training pairs or test tissue")

            dest = MODELS / rep / lineage.replace("/", "_") / f"fold{fold}"
            dest.mkdir(parents=True, exist_ok=True)
            ids = design.samples
            test_ids = test["sample_uid"].tolist()
            h = project_subset(
                group, design, test_ids, rep, lineage, outer_fold=fold, destination=dest
            )
            trmat = [h[s] for s in ids]
            train_batch = BagBatch(
                trmat, [group[s][2]["hierarchy_type"].to_numpy(dtype=str) for s in ids]
            )
            train_batch_hier = _batch_to_device(train_batch, hier_dev)
            test_batch = BagBatch(
                [h[s] for s in test_ids],
                [group[s][2]["hierarchy_type"].to_numpy(dtype=str) for s in test_ids],
            )
            inner_batches = (
                prepare_inner_cv_batches(
                    group, train, rep, lineage, fold, device=hier_dev
                )
                if tune_inner
                else None
            )

            kinds = ["CTA", "HA"]
            if rep != "scGPT512":
                kinds.append("MixMIL")
            if rep == "scGPT32" or lineage == JOINT_LINEAGE:
                kinds += ["mean", "attention"]

            for kind in kinds:
                began = time.perf_counter()
                if kind in ("mean", "attention"):
                    if rep == "PCA32" and lineage == JOINT_LINEAGE:
                        method = f"PCA_{kind}"
                    elif rep == "scGPT512" and lineage == JOINT_LINEAGE:
                        method = f"scGPT_{kind}"
                    else:
                        method = f"{rep}_{kind}"
                else:
                    method = f"{rep}_{kind}"

                common = {
                    "method": method,
                    "lineage": lineage,
                    "fold": fold,
                    "seed": SEED,
                    "n_train_bags": len(train),
                    "n_train_donors": train["donor_uid"].nunique(),
                    "n_test_bags": len(test),
                }

                if kind in ("mean", "attention"):
                    model, fit = fit_model(
                        [x.astype(float) for x in trmat],
                        design,
                        pooling=kind,
                        lambda_pair=0.1,
                        lambda_l2=0.01,
                        seed=SEED,
                        standardize=False,
                        maxiter=500,
                    )
                    meta = {
                        "fit": fit,
                        "lineage": lineage,
                        "training_samples": ids,
                        "training_donors": sorted(train["donor_uid"].unique()),
                    }
                    save_model(dest / f"{kind}.npz", model, meta)
                    pred, _ = predict_bags(
                        model, meta, [h[s].astype(float) for s in test_ids], test_ids
                    )
                    logits = pred["logit"].to_numpy()
                    n_params = fit["n_parameters"]
                    steps = fit["iterations"]
                    status = "converged" if fit["converged"] else fit["message"]
                    obj = fit.get("final_training_objective")
                    configs_tried = 1
                    sel_drop, sel_wd, inner_auc = np.nan, np.nan, np.nan
                elif kind == "MixMIL":
                    torch.manual_seed(SEED)
                    width = trmat[0].shape[1]
                    model = MixedMIL(
                        trmat,
                        design.labels,
                        posterior,
                        scale_by_sqrt_q=True,
                        legacy_even_slice_init=False,
                    )
                    steps = 500
                    obj = train_fixed(
                        model, train_batch, design, train["donor_uid"].nunique(), steps
                    )
                    n_params = sum(p.numel() for p in model.parameters() if p.requires_grad)
                    torch.save(
                        {
                            "state_dict": model.state_dict(),
                            "width": width,
                            "method": method,
                            "seed": SEED,
                            "scale_by_sqrt_q": True,
                            "legacy_even_slice_init": False,
                            "training_samples": ids,
                            "training_donors": sorted(train["donor_uid"].unique()),
                        },
                        dest / f"{kind}.pt",
                    )
                    with torch.inference_mode():
                        lin, att, signed, bias = model(test_batch)
                        logits = lin.numpy()
                        first = test_ids[0]
                        single = BagBatch(
                            [h[first]],
                            [group[first][2]["hierarchy_type"].to_numpy(dtype=str)],
                        )
                        np.testing.assert_allclose(
                            model(single)[0].numpy(), logits[:1], rtol=2e-5, atol=2e-5
                        )
                        at = att.numpy()
                        si = signed.numpy()
                        owner = test_batch.owner.numpy()
                        cells = pd.concat(
                            [
                                group[s][2][["sample_uid", "subtype_m2", "hierarchy_type"]]
                                for s in test_ids
                            ]
                        ).copy()
                        np.testing.assert_allclose(
                            np.bincount(owner, weights=at), np.ones(len(test_ids)), atol=2e-5
                        )
                        np.testing.assert_allclose(
                            np.bincount(owner, weights=si) + bias.item(),
                            logits,
                            rtol=2e-5,
                            atol=2e-5,
                        )
                        cell_parts.append(
                            cells.assign(
                                method=method,
                                lineage=lineage,
                                fold=fold,
                                attention=at,
                                signed_logit_contribution=si,
                                bias=bias.item(),
                            )
                        )
                    status = "completed_fixed_steps"
                    configs_tried = 1
                    sel_drop, sel_wd, inner_auc = 0.0, 0.0, np.nan
                else:
                    width = trmat[0].shape[1]
                    if tune_inner and inner_batches is not None:
                        best_cfg, cfg_records = tune_hierarchical_mil(
                            width,
                            kind == "HA",
                            inner_batches,
                            method=method,
                            lineage=lineage,
                            outer_fold=fold,
                            device=hier_dev,
                        )
                        tuning_rows.extend(cfg_records)
                        sel_drop = float(best_cfg["dropout"])
                        sel_wd = float(best_cfg["weight_decay"])
                        steps = int(best_cfg["steps"])
                        inner_auc = float(best_cfg["inner_oof_auroc"])
                        configs_tried = len(cfg_records)
                    else:
                        sel_drop, sel_wd, steps, inner_auc, configs_tried = (
                            0.0,
                            1e-3,
                            100,
                            np.nan,
                            1,
                        )

                    _seed_all(SEED)
                    model = HierarchicalMIL(width, kind == "HA", dropout=sel_drop).to(
                        hier_dev
                    )
                    optim_outer = torch.optim.Adam(
                        model.parameters(), lr=1e-3, weight_decay=float(sel_wd)
                    )
                    obj = _train_hier_on_device(
                        model,
                        train_batch_hier,
                        design,
                        steps,
                        optim=optim_outer,
                        lambda_pair=0.1,
                    )
                    model = model.to("cpu").eval()
                    n_params = sum(p.numel() for p in model.parameters() if p.requires_grad)
                    torch.save(
                        {
                            "state_dict": model.state_dict(),
                            "width": width,
                            "method": method,
                            "seed": SEED,
                            "dropout": sel_drop,
                            "weight_decay": sel_wd,
                            "steps": steps,
                            "inner_oof_auroc": inner_auc,
                            "training_samples": ids,
                            "training_donors": sorted(train["donor_uid"].unique()),
                        },
                        dest / f"{kind}.pt",
                    )
                    with torch.inference_mode():
                        lin, att, signed, bias = model(test_batch)
                        logits = lin.numpy()
                        first = test_ids[0]
                        single = BagBatch(
                            [h[first]],
                            [group[first][2]["hierarchy_type"].to_numpy(dtype=str)],
                        )
                        np.testing.assert_allclose(
                            model(single)[0].numpy(), logits[:1], rtol=2e-5, atol=2e-5
                        )
                        at = att.numpy()
                        si = signed.numpy()
                        owner = test_batch.owner.numpy()
                        cells = pd.concat(
                            [
                                group[s][2][["sample_uid", "subtype_m2", "hierarchy_type"]]
                                for s in test_ids
                            ]
                        ).copy()
                        np.testing.assert_allclose(
                            np.bincount(owner, weights=at), np.ones(len(test_ids)), atol=2e-5
                        )
                        np.testing.assert_allclose(
                            np.bincount(owner, weights=si) + bias.item(),
                            logits,
                            rtol=2e-5,
                            atol=2e-5,
                        )
                        cell_parts.append(
                            cells.assign(
                                method=method,
                                lineage=lineage,
                                fold=fold,
                                attention=at,
                                signed_logit_contribution=si,
                                bias=bias.item(),
                            )
                        )
                    status = "tuned_inner_cv" if tune_inner else "completed_fixed_steps"

                if not np.isfinite(logits).all():
                    raise FloatingPointError(f"{method}: non-finite OOF logit")
                runs.append(
                    dict(
                        common,
                        configs_tried=configs_tried,
                        selected_dropout=sel_drop,
                        selected_weight_decay=sel_wd,
                        steps=steps,
                        inner_oof_auroc=inner_auc,
                        n_parameters=n_params,
                        seconds=time.perf_counter() - began,
                        status=status,
                        objective=obj,
                    )
                )
                if status in ("converged", "completed_fixed_steps", "tuned_inner_cv"):
                    predictions.append(
                        test[["sample_uid", "donor_uid", "lesion_uid", "tissue"]].assign(
                            method=method,
                            lineage=lineage,
                            fold=fold,
                            seed=SEED,
                            logit=logits,
                        )
                    )
                print(
                    f"{rep} {lineage} fold={fold} {kind}: {status} ({time.perf_counter()-began:.2f}s)",
                    flush=True,
                )
            pd.DataFrame(runs).to_csv(out / "benchmark_runs.csv", index=False)
            pd.concat(predictions, ignore_index=True).to_csv(
                out / "oof_predictions.csv", index=False
            )

    if tuning_rows:
        pd.DataFrame(tuning_rows).to_csv(out / "inner_cv_tuning.csv", index=False)
    cell_df = pd.concat(cell_parts)
    cell_df.index.name = "cell_uid"
    cell_df.to_csv(
        out / "cell_contributions.csv.gz", compression={"method": "gzip", "mtime": 0}
    )
    (out / "profile.json").write_text(
        json.dumps(
            {
                **prov_info,
                "gene_universe": DEFAULT_UNIVERSE,
                "tasks": list(TASKS),
                "inner_cv_tuned": bool(tune_inner),
                "hier_device": hier_dev,
                "torch": torch.__version__,
                "torch_threads": torch.get_num_threads(),
                "seconds": time.perf_counter() - start,
                "peak_rss_bytes": int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024),
                "command": " ".join([sys.executable, *sys.argv]),
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--representation", choices=["PCA32", "scGPT32", "scGPT512"], required=True
    )
    parser.add_argument("--mixmil-repo", type=Path, default=DEFAULT_MIXMIL_DIR)
    parser.add_argument("--device", choices=["cpu", "xpu"], default=None)
    parser.add_argument("--no-inner-tune", action="store_true")
    args = parser.parse_args()
    with threadpool_limits(limits=4):
        run(
            args.representation,
            args.mixmil_repo,
            tune_inner=not args.no_inner_tune,
            device=args.device,
        )
