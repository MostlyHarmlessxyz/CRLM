"""Patient-held-out classification of tumor specimens by lineage on full QC-retained bags.

Track 2 (`V1.1`) methodological guarantees:
1. Full QC-retained tumor bags (`min_cells=50, max_cells=None`) on the 18,815-gene
   `245_164_shared_nuclear` canonical axis (guaranteeing 100% external PCA HVG coverage),
   paired with a 10-draw `256`-cell sub-bag variance audit (`subbag_256x10_sensitivity.csv`).
2. Strict outer `train_donors` isolation for HVG selection, 32-D PCA, fold-internal
   `NL` vs `NC` Nonimmune parenchymal marker masking (`PCA_remove_organ_markers`),
   `debiased_ucell_scores` matched-control construction (`program_elastic`), and scaling.
3. Clean hierarchical subtype deconfounding (`PCA_subtype_balanced` with exact 50/50
   dominant-strata centroids and zero cell starvation, plus `PCA_pure_primary` and
   `PCA_pure_secondary` pure single-subtype controls).
4. Class-balanced donor-site weights (`sum(w_PT) == sum(w_LM) == 0.5`) + paired softplus
   ranking loss for `program_elastic`, `composition_elastic`, and `QC_elastic`, eliminating
   unpaired-PT double-weighting sign inversion.
5. 200-permutation donor-level Rademacher randomization null (`permutation_null_summary.csv`).
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
import time

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.optimize import minimize
from scipy.special import expit
from scipy.stats import spearmanr
from sklearn.decomposition import PCA
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src"), str(ROOT / "starter")]

from crlm.benchmark_statistics import summarize
from crlm.hgnc_alias_map import (
    canonicalize_symbols,
    get_hgnc_mapper,
    load_shared_nuclear_universe,
)
from crlm.model_bags import (
    LINEAGES,
    SEED,
    cell_seed,
    sample_bags,
    specimen_immune_composition,
    subsample_bag_indices,
)
from crlm.program_scores import (
    build_detection_matched_controls,
    compute_reference_gene_stats,
    controls_to_dataframe,
    debiased_ucell_scores,
)
from linear_attention_mil import (
    bag_forward,
    build_design,
    fit_model,
    predict_bags,
    save_model,
    unpack,
)

OUT = ROOT / "artifacts/m4"
MODELS = ROOT / "models/m4"
SPLIT = ROOT / "data/manifests/m0/donor_folds_pre_qc.json"
MANIFEST = ROOT / "data/manifests/m0/sample_manifest_validated.csv"
PROGRAMS = ROOT / "artifacts/m3/program_definitions_draft/program_genes.csv"
DEFAULT_UNIVERSE = "245_164_shared_nuclear"
N_PERMUTATIONS = 200
N_SUBBAG_REPEATS = 10

PRESPECIFIED_ORGAN_MARKERS = frozenset(
    {
        "ALB",
        "APOA1",
        "APOA2",
        "APOC1",
        "APOC3",
        "APOE",
        "APOH",
        "TTR",
        "FGA",
        "FGB",
        "FGG",
        "HP",
        "HPR",
        "HNF4A",
        "ASGR1",
        "ASGR2",
        "CPS1",
        "FTCD",
        "CYP2E1",
        "SERPINA1",
        "CRP",
        "SAA1",
        "SAA2",
        "EPCAM",
        "KRT8",
        "KRT18",
        "KRT19",
        "PIGR",
        "FABP1",
        "MUC2",
        "TFF3",
        "CEACAM5",
        "CEACAM6",
        "IGHA1",
        "IGHA2",
        "JCHAIN",
    }
)

DOMINANT_SUBTYPE_PAIRS: dict[str, tuple[str, str]] = {
    "T/NK": ("CD8 T cell", "CD4 T cell"),
    "Mono/Mac/DC": ("Macrophage", "Monocyte"),
    "B/Plasma": ("plasma cell", "B cell"),
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


def load_lineage_programs(
    universe_genes: list[str] | tuple[str, ...] | None = None,
    program_csv: Path = PROGRAMS,
) -> dict[str, dict[str, list[str]]]:
    """Load complete lineage programs canonicalized onto `universe_genes` (>=5 genes & >=80% coverage)."""
    if universe_genes is None:
        universe_genes = list(load_shared_nuclear_universe(DEFAULT_UNIVERSE))
    universe_set = set(str(g) for g in universe_genes)
    table = pd.read_csv(program_csv)
    complete = table.loc[table["role"].eq("complete_program")].copy()
    result: dict[str, dict[str, list[str]]] = {lin: {} for lin in LINEAGES}
    for (lin, prog), grp in complete.groupby(["lineage", "program"], sort=True):
        if lin not in result:
            continue
        canon_all = canonicalize_symbols(grp["gene_symbol"].astype(str).tolist(), deduplicate=True)
        measured = [g for g in canon_all if g in universe_set]
        if len(measured) >= 5 and len(measured) >= 0.8 * len(canon_all):
            result[lin][str(prog)] = measured
    for lin in LINEAGES:
        if not result[lin]:
            raise ValueError(f"Zero measurable complete programs for {lin} in {program_csv}")
    return result


def library_inputs(
    write_selection: bool = True,
    *,
    gene_universe: str = DEFAULT_UNIVERSE,
    max_cells: int | None = None,
):
    """Load full QC-retained development tumor bags (`PT`/`LM` only) on `gene_universe`."""
    universe_genes = list(load_shared_nuclear_universe(gene_universe))
    rows = pd.read_csv(MANIFEST).query(
        "role == 'development' and tissue in ['PT','LM']"
    ).reset_index(drop=True)
    if set(rows["tissue"]) != {"PT", "LM"}:
        raise ValueError("Expected verified development PT and LM rows only")

    groups: dict[str, dict[str, tuple]] = {lineage: {} for lineage in LINEAGES}
    samples = []
    cells = []
    for row in rows.itertuples(index=False):
        gsm = row.sample_uid
        h5_path = ROOT / f"data/processed/m2/{gsm}.h5ad"
        bags = sample_bags(
            h5_path,
            min_cells=50,
            max_cells=max_cells,
            gene_universe=universe_genes,
        )
        comp_full = specimen_immune_composition(h5_path)
        for lineage in LINEAGES:
            bag = bags.get(lineage)
            samples.append(
                dict(
                    sample_uid=gsm,
                    donor_uid=row.donor_uid,
                    tissue=row.tissue,
                    lineage=lineage,
                    n_population=bag.population if bag is not None else 0,
                    n_selected=len(bag.obs) if bag is not None else 0,
                    eligible=bag is not None,
                    reason="" if bag else "fewer_than_50_retained",
                )
            )
            if bag is None:
                continue
            x_log1p = _normalize_csr_log1p_cp10k(bag.counts)
            lin_comp = (
                comp_full.loc[comp_full["lineage_m2"].eq(lineage)]
                .set_index("fine_subtype")["prop_within_lineage"]
                .copy()
            )
            cell = bag.obs.copy()
            cell["population"] = bag.population
            cells.append(cell)
            groups[lineage][gsm] = (bag, x_log1p, lin_comp)

    if write_selection:
        OUT.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(samples).to_csv(OUT / "eligible_bags.csv", index=False)
        selected_df = pd.concat(cells)
        selected_df.index.name = "cell_uid"
        selected_df.to_csv(
            OUT / "selected_cells.csv.gz", compression={"method": "gzip", "mtime": 0}
        )
    return rows, groups


def load_normal_nonimmune_profiles(
    universe_genes: list[str] | tuple[str, ...],
) -> list[dict[str, object]]:
    """Load `NC`/`NL` normal `Nonimmune` QC-retained mean `log1p(CP10k)` profiles per library.

    Normal tissues are NEVER used for tumor classification labels or model fitting;
    they are exclusively sliced per training fold (`donor_uid in train_donors`) to
    identify empirical hepatic/colonic parenchymal background genes for `PCA_remove_organ_markers`.
    """
    rows = pd.read_csv(MANIFEST).query(
        "role == 'development' and tissue in ['NC','NL']"
    )
    universe_idx = pd.Index([str(g) for g in universe_genes])
    profiles: list[dict[str, object]] = []
    for row in rows.itertuples(index=False):
        path = ROOT / f"data/processed/m2/{row.sample_uid}.h5ad"
        obj = ad.read_h5ad(path, backed="r")
        try:
            obs = obj.obs
            var_syms = pd.Index(obj.var["canonical_symbol"].astype(str))
            col_idx = var_syms.get_indexer(universe_idx)
            if (col_idx < 0).any():
                raise ValueError(f"Missing shared nuclear genes in normal library {path}")
            keep = np.flatnonzero(
                obs["qc_retained"].to_numpy(dtype=bool)
                & obs["lineage_m2"].astype(str).eq("Nonimmune").to_numpy()
            )
            if len(keep) == 0:
                continue
            mat = obj.X[keep, :].tocsr()[:, col_idx].tocsr()
            x_log = _normalize_csr_log1p_cp10k(mat)
            mean_vec = np.asarray(x_log.mean(axis=0)).ravel().astype(np.float64)
            profiles.append(
                {
                    "sample_uid": str(row.sample_uid),
                    "donor_uid": str(row.donor_uid),
                    "tissue": str(row.tissue),
                    "n_cells": int(len(keep)),
                    "mean_log1p": mean_vec,
                }
            )
        finally:
            obj.file.close()
    return profiles


def build_fold_organ_mask(
    normal_profiles: list[dict[str, object]],
    train_donors: set[str],
    test_donors: set[str],
    universe_genes: list[str],
) -> tuple[set[str], dict[str, object]]:
    """Compute fold-internal `NL` vs `NC` `Nonimmune` parenchymal gene mask strictly from `train_donors`."""
    if train_donors & test_donors:
        raise ValueError("Train/test donor overlap in build_fold_organ_mask")
    tr_prof = [p for p in normal_profiles if str(p["donor_uid"]) in train_donors]
    if any(str(p["donor_uid"]) in test_donors for p in tr_prof):
        raise ValueError("Normal profile test-donor leakage detected")

    by_tissue_donor: dict[str, dict[str, list[np.ndarray]]] = {"NC": {}, "NL": {}}
    for p in tr_prof:
        by_tissue_donor[str(p["tissue"])].setdefault(str(p["donor_uid"]), []).append(
            np.asarray(p["mean_log1p"], dtype=np.float64)
        )
    if not by_tissue_donor["NC"] or not by_tissue_donor["NL"]:
        raise ValueError("Training fold lacks NC or NL Nonimmune normal profiles")

    nc_mean = np.mean(
        [np.mean(vecs, axis=0) for vecs in by_tissue_donor["NC"].values()], axis=0
    )
    nl_mean = np.mean(
        [np.mean(vecs, axis=0) for vecs in by_tissue_donor["NL"].values()], axis=0
    )
    emp_mask = (np.abs(nl_mean - nc_mean) >= 0.5) | (np.maximum(nl_mean, nc_mean) >= 1.5)
    emp_genes = {universe_genes[i] for i in np.flatnonzero(emp_mask)}
    prespec_in_univ = set(PRESPECIFIED_ORGAN_MARKERS) & set(universe_genes)
    union_genes = emp_genes | prespec_in_univ
    audit = {
        "n_train_nc_donors": len(by_tissue_donor["NC"]),
        "n_train_nl_donors": len(by_tissue_donor["NL"]),
        "n_empirical_nonimmune_genes": len(emp_genes),
        "n_prespecified_genes": len(prespec_in_univ),
        "n_total_masked_universe": len(union_genes),
    }
    return union_genes, audit


def pca_rep(
    group: dict,
    train: list[str],
    dest: Path | None = None,
    *,
    mask_genes: set[str] | None = None,
) -> tuple[dict[str, np.ndarray], np.ndarray]:
    """Fit top-512 HVGs and 32-D PCA strictly on `train` samples and project all samples in `group`."""
    all_genes = next(iter(group.values()))[0].genes
    if not all(all_genes.equals(v[0].genes) for v in group.values()):
        raise ValueError("Tumor gene axes differ")
    x_train = sparse.vstack([group[s][1] for s in train], format="csr")
    var = (
        np.asarray(x_train.power(2).mean(axis=0)).ravel()
        - np.asarray(x_train.mean(axis=0)).ravel() ** 2
    )
    if mask_genes:
        var = var.copy()
        var[all_genes.isin(mask_genes)] = -np.inf
    idx = np.sort(np.argsort(-var, kind="stable")[:512])
    selected_genes = np.asarray(all_genes[idx], dtype=str)
    model = PCA(n_components=32, svd_solver="randomized", random_state=SEED)
    model.fit(x_train[:, idx].toarray())
    if dest is not None:
        np.savez_compressed(
            dest,
            genes=selected_genes,
            center=model.mean_,
            components=model.components_,
            training_samples=np.asarray(train, dtype=str),
            gene_universe=DEFAULT_UNIVERSE,
        )
    projected = {
        s: model.transform(v[1][:, idx].toarray()).astype(float)
        for s, v in group.items()
    }
    return projected, selected_genes


def score_programs_train_locked(
    group: dict,
    train_samples: list[str],
    eval_samples: list[str],
    program_defs: dict[str, list[str]],
    universe_genes: list[str],
) -> tuple[dict[str, pd.DataFrame], dict[str, dict[str, list[str]]], pd.DataFrame]:
    """Build detection-matched controls strictly on `train_samples` and score `debiased_ucell_scores`."""
    ref_stats = compute_reference_gene_stats(
        [group[s][0].counts for s in train_samples],
        universe_genes,
        input_scale="counts",
        n_bins=20,
    )
    matched_controls = build_detection_matched_controls(
        program_defs,
        ref_stats,
        candidate_universe=universe_genes,
        n_bins=20,
        n_ctrl=10,
        seed=SEED,
    )
    scored: dict[str, pd.DataFrame] = {}
    for s in eval_samples:
        df = debiased_ucell_scores(
            group[s][0].counts,
            universe_genes,
            program_defs,
            matched_controls=matched_controls,
            max_rank=1500,
        )
        df.index = group[s][0].obs.index
        scored[s] = df
    return scored, matched_controls, ref_stats


def balanced_donor_site_weights(rows: pd.DataFrame) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute class-balanced donor-site sample weights (`sum(w_PT) == sum(w_LM) == 0.5`).

    Eliminates the 2x unpaired-PT weight distortion of `1 / (n_donors * n_sites)` while
    retaining exact lesion -> sample equal weighting within each `(donor_uid, tissue)`
    and returning the paired-donor contrast matrix `pair_matrix` (`LM - PT`).
    """
    rows = rows.reset_index(drop=True)
    design = build_design(rows)
    n_pt_donors = rows.loc[rows["tissue"].eq("PT"), "donor_uid"].nunique()
    n_lm_donors = rows.loc[rows["tissue"].eq("LM"), "donor_uid"].nunique()
    if n_pt_donors < 1 or n_lm_donors < 1:
        raise ValueError("Require at least one PT donor and one LM donor")

    weights = np.zeros(len(rows), dtype=np.float64)
    for (_, site), site_rows in rows.groupby(["donor_uid", "tissue"], sort=True):
        donor_w = 0.5 / n_lm_donors if site == "LM" else 0.5 / n_pt_donors
        n_lesions = site_rows["lesion_uid"].nunique()
        for _, lesion_rows in site_rows.groupby("lesion_uid", sort=True):
            idx = lesion_rows.index.to_numpy()
            weights[idx] = donor_w / (n_lesions * len(idx))

    if not np.isclose(weights.sum(), 1.0) or not np.isclose(
        weights[design.labels == 1].sum(), 0.5
    ):
        raise ValueError("Class-balanced donor-site weight normalization error")
    return weights, design.labels, design.pair_matrix


def fit_balanced_paired_elastic(
    x_train: np.ndarray,
    rows_train: pd.DataFrame,
    *,
    alpha: float = 0.05,
    l1_ratio: float = 0.5,
    lambda_pair: float = 0.25,
    maxiter: int = 2000,
) -> dict[str, object]:
    """Fit class-balanced + paired-ranking Elastic-Net logistic regression via smooth bound-constrained L-BFGS-B."""
    w_bal, y, pair_matrix = balanced_donor_site_weights(rows_train)
    center = w_bal @ x_train
    variance = w_bal @ ((x_train - center) ** 2)
    scale = np.sqrt(np.maximum(variance, 0.0))
    scale = np.where(scale >= 1e-12, scale, 1.0)
    z = (x_train - center) / scale

    p = z.shape[1]
    l1_pen = float(alpha) * float(l1_ratio)
    l2_pen = float(alpha) * (1.0 - float(l1_ratio))

    def obj_and_grad(u: np.ndarray) -> tuple[float, np.ndarray]:
        w_pos = u[:p]
        w_neg = u[p : 2 * p]
        bias = float(u[2 * p])
        w = w_pos - w_neg
        logits = z @ w + bias
        bce = float(np.dot(w_bal, np.logaddexp(0.0, logits) - y * logits))
        d_logits = w_bal * (expit(logits) - y)
        pair_loss = 0.0
        if len(pair_matrix) > 0 and lambda_pair > 0:
            delta = pair_matrix @ logits
            pair_loss = float(np.logaddexp(0.0, -delta).mean())
            d_logits = d_logits + float(lambda_pair) * (
                pair_matrix.T @ (-expit(-delta) / len(delta))
            )
        reg = l1_pen * float(w_pos.sum() + w_neg.sum()) + 0.5 * l2_pen * float(w @ w)
        grad_w_smooth = z.T @ d_logits + l2_pen * w
        grad = np.empty(2 * p + 1, dtype=np.float64)
        grad[:p] = grad_w_smooth + l1_pen
        grad[p : 2 * p] = -grad_w_smooth + l1_pen
        grad[2 * p] = float(d_logits.sum())
        return bce + float(lambda_pair) * pair_loss + reg, grad

    u0 = np.zeros(2 * p + 1, dtype=np.float64)
    bounds = [(0.0, None)] * (2 * p) + [(None, None)]
    res = minimize(
        obj_and_grad,
        u0,
        method="L-BFGS-B",
        jac=True,
        bounds=bounds,
        options={"maxiter": int(maxiter), "gtol": 1e-8, "ftol": 1e-13, "maxls": 50},
    )
    w_raw = res.x[:p] - res.x[p : 2 * p]
    w_clean = np.where(np.abs(w_raw) < 1e-9, 0.0, w_raw)
    bias = float(res.x[2 * p])
    return {
        "coef": w_clean[None, :],
        "bias": np.array([bias], dtype=np.float64),
        "center": center,
        "scale": scale,
        "bce_weights_balanced": w_bal,
        "converged": bool(res.success),
        "iterations": int(res.nit),
        "objective": float(res.fun),
        "n_parameters": int(p + 1),
    }


def elastic(
    group: dict,
    rows: pd.DataFrame,
    test: list[str],
    kind: str,
    dest: Path,
    *,
    lineage: str | None = None,
    program_defs: dict[str, list[str]] | None = None,
    universe_genes: list[str] | None = None,
    precomputed_scores: tuple[dict[str, pd.DataFrame], dict, pd.DataFrame] | None = None,
):
    """Fit fold-isolated `program_elastic`, `composition_elastic`, or `QC_elastic`."""
    train_samples = rows["sample_uid"].tolist()
    sample_names = train_samples + list(test)

    if kind == "program_elastic":
        if universe_genes is None:
            universe_genes = list(load_shared_nuclear_universe(DEFAULT_UNIVERSE))
        if program_defs is None:
            if lineage is None:
                lineage = str(next(iter(group.values()))[0].obs["lineage_m2"].iloc[0])
            program_defs = load_lineage_programs(universe_genes)[lineage]
        if precomputed_scores is None:
            scored_map, matched_controls, ref_stats = score_programs_train_locked(
                group, train_samples, sample_names, program_defs, universe_genes
            )
        else:
            scored_map, matched_controls, ref_stats = precomputed_scores
        col_names = list(program_defs.keys())
        vectors = [scored_map[s][col_names].mean(axis=0).to_numpy(dtype=float) for s in sample_names]
    elif kind == "composition_elastic":
        strata = sorted(
            {str(sub) for s in train_samples for sub in group[s][2].index}
        )
        col_names = strata
        vectors = [
            group[s][2].reindex(strata, fill_value=0.0).to_numpy(dtype=float)
            for s in sample_names
        ]
    elif kind == "QC_elastic":
        col_names = [
            f"{c}_q{q}"
            for q in (0.25, 0.5, 0.75)
            for c in ("n_counts_visible", "n_genes_visible")
        ]
        vectors = []
        for s in sample_names:
            q_df = group[s][0].obs[["n_counts_visible", "n_genes_visible"]].quantile(
                [0.25, 0.5, 0.75]
            )
            vectors.append(np.log1p(q_df.to_numpy(dtype=float)).ravel())
    else:
        raise ValueError(f"Unsupported elastic baseline kind: {kind}")

    x = np.stack(vectors, axis=0)
    n_train = len(train_samples)
    fitted = fit_balanced_paired_elastic(
        x[:n_train],
        rows,
        alpha=0.05,
        l1_ratio=0.5,
        lambda_pair=0.25,
    )
    dest.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        dest,
        coef=fitted["coef"],
        bias=fitted["bias"],
        center=fitted["center"],
        scale=fitted["scale"],
        columns=np.asarray(col_names, dtype=str),
        training_samples=np.asarray(train_samples, dtype=str),
        bce_weights_balanced=fitted["bce_weights_balanced"],
        scoring_method="ucell_debiased" if kind == "program_elastic" else kind,
        gene_universe=DEFAULT_UNIVERSE,
    )
    if kind == "program_elastic":
        mapper = get_hgnc_mapper()
        (dest.parent / "program_defs.json").write_text(
            json.dumps(
                {
                    "lineage": lineage,
                    "gene_universe": DEFAULT_UNIVERSE,
                    "hgnc_sha256": mapper.sha256,
                    "scoring_method": "ucell_debiased",
                    "max_rank": 1500,
                    "n_bins": 20,
                    "n_ctrl": 10,
                    "seed": SEED,
                    "programs": program_defs,
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
        controls_to_dataframe(matched_controls, ref_stats).to_csv(
            dest.parent / "matched_controls.csv", index=False
        )

    if test:
        z_test = (x[n_train:] - fitted["center"]) / fitted["scale"]
        test_logits = (z_test @ fitted["coef"].T + fitted["bias"]).ravel()
    else:
        test_logits = np.array([], dtype=float)
    return (
        test_logits,
        bool(fitted["converged"]),
        int(fitted["iterations"]),
        int(fitted["n_parameters"]),
    )


def subtype_balanced_centroids(
    bag,
    h: np.ndarray,
    lineage: str,
    *,
    min_per_stratum: int = 5,
) -> np.ndarray:
    """Exact 50/50 equal-stratum centroid representation on the two dominant subtypes.

    Uses ALL cells in `s1` and `s2` (zero subsampling starvation) while enforcing
    identical 50/50 subtype weighting across every PT and LM bag.
    """
    s1, s2 = DOMINANT_SUBTYPE_PAIRS[lineage]
    fine = bag.obs["fine_subtype"].astype(str).to_numpy()
    idx1 = np.flatnonzero(fine == s1)
    idx2 = np.flatnonzero(fine == s2)
    if len(idx1) < min_per_stratum or len(idx2) < min_per_stratum:
        return np.empty((0, h.shape[1]), dtype=float)
    return np.vstack([h[idx1].mean(axis=0), h[idx2].mean(axis=0)])


def pure_subtype_slice(
    bag,
    h: np.ndarray,
    target_subtype: str,
    *,
    min_cells: int = 20,
) -> np.ndarray:
    """Extract cells belonging exclusively to a single fine subtype (`>= min_cells`)."""
    fine = bag.obs["fine_subtype"].astype(str).to_numpy()
    idx = np.flatnonzero(fine == target_subtype)
    if len(idx) < min_cells:
        return np.empty((0, h.shape[1]), dtype=float)
    return h[idx]


def permute_training_donors(
    tr: pd.DataFrame,
    rng: np.random.Generator,
    *,
    balanced_half: bool = False,
) -> pd.DataFrame:
    """Permute PT/LM labels at the `donor_uid` level strictly within training rows."""
    donor_sites = tr.groupby("donor_uid")["tissue"].agg(lambda s: set(s.tolist())).to_dict()
    paired = sorted(d for d, ts in donor_sites.items() if ts == {"PT", "LM"})
    unpaired = sorted(d for d, ts in donor_sites.items() if ts != {"PT", "LM"})

    for _ in range(64):
        out = tr.copy()
        flip_donors: set[str] = set()
        if balanced_half:
            if paired:
                k_p = len(paired) // 2 + (
                    int(rng.integers(0, 2)) if len(paired) % 2 == 1 else 0
                )
                if k_p > 0:
                    flip_donors.update(rng.choice(paired, size=k_p, replace=False).tolist())
            if unpaired:
                k_u = len(unpaired) // 2 + (
                    int(rng.integers(0, 2)) if len(unpaired) % 2 == 1 else 0
                )
                if k_u > 0:
                    flip_donors.update(rng.choice(unpaired, size=k_u, replace=False).tolist())
        else:
            for d in paired + unpaired:
                if bool(rng.integers(0, 2)):
                    flip_donors.add(d)

        mask = out["donor_uid"].isin(flip_donors)
        out.loc[mask, "tissue"] = out.loc[mask, "tissue"].map({"PT": "LM", "LM": "PT"})
        if out["tissue"].nunique() == 2:
            return out
    raise RuntimeError("Failed to generate non-degenerate donor permutation")


def scgpt_rep(group: dict) -> dict[str, np.ndarray]:
    out = {}
    for gsm, (bag, _, _) in group.items():
        path = ROOT / f"embeddings/m4/scgpt/{gsm}.npz"
        with np.load(path, allow_pickle=False) as f:
            subset = pd.Index(f["cell_uid"]).get_indexer(bag.obs.index)
            if min(subset) < 0 or len(set(subset)) != len(subset):
                raise ValueError(f"{gsm} scGPT cell identity mismatch")
            out[gsm] = f["embedding"][subset].astype(float)
    return out


def cell_attribution(model, bag, values, scores, pooling, gsm, seed):
    z = (values - model["center"]) / model["scale"]
    w, bias, v = unpack(model["theta"], z.shape[1], pooling)
    logit, _, _, a = bag_forward(z, w, bias, v)
    signed = a * (z @ w)
    if not np.isclose(bias + signed.sum(), logit):
        raise ValueError("Cell contributions do not reconstruct the bag logit")
    k = max(1, len(a) // 10)
    top = np.argsort(-np.abs(signed), kind="stable")[:k]
    rng = np.random.default_rng(cell_seed(f"{gsm}|{seed}|remove"))
    subtype = bag.obs["fine_subtype"].astype(str).to_numpy()
    random = np.concatenate(
        [
            rng.choice(
                np.flatnonzero(subtype == t), np.sum(subtype[top] == t), replace=False
            )
            for t in np.unique(subtype[top])
        ]
    )

    def drop(indices):
        return logit - bag_forward(np.delete(z, indices, axis=0), w, bias, v)[0]

    frame = pd.DataFrame(
        {"cell_uid": bag.obs.index, "attention": a, "signed_contribution": signed}
    )
    corr = {
        f"rho_{c}": spearmanr(signed, scores[c]).statistic
        if scores[c].nunique() > 1
        else np.nan
        for c in scores
    }
    return frame, dict(
        delta_remove_top_abs=drop(top),
        delta_remove_random_subtype_matched=drop(random),
        n_removed=k,
        **corr,
    )


def run_permutation_null(
    fold_cache: dict[tuple[str, int], tuple[pd.DataFrame, pd.DataFrame, dict[str, np.ndarray]]],
    observed_metrics: pd.DataFrame,
    outdir: Path,
    *,
    n_permutations: int = N_PERMUTATIONS,
) -> None:
    """Run `n_permutations` unconditional donor-level sign-flip permutations on `PCA_mean`."""
    perm_frames = []
    for perm_id in range(n_permutations):
        for (lineage, f), (training, testing, vectors) in fold_cache.items():
            rng = np.random.default_rng(cell_seed(f"perm_rademacher|fold{f}|{perm_id}"))
            tr_perm = permute_training_donors(training, rng, balanced_half=False)
            design = build_design(tr_perm)
            pair_w = 0.1 if design.pair_donors else 0.0
            model, fit = fit_model(
                [vectors[s] for s in design.samples],
                design,
                pooling="mean",
                lambda_pair=pair_w,
                lambda_l2=0.01,
                seed=SEED,
                standardize=True,
                maxiter=500,
            )
            if not fit["converged"]:
                continue
            ids = testing["sample_uid"].tolist()
            meta = {"lineage": lineage, "fit": fit, "training_samples": design.samples}
            pred, _ = predict_bags(model, meta, [vectors[s] for s in ids], ids)
            perm_frames.append(
                testing[["sample_uid", "donor_uid", "lesion_uid", "tissue"]].assign(
                    method="PCA_permuted_draw",
                    lineage=lineage,
                    fold=f,
                    seed=perm_id,
                    logit=pred["logit"].to_numpy(),
                )
            )

    all_perm_oof = pd.concat(perm_frames, ignore_index=True)
    _, perm_metrics, _ = summarize(all_perm_oof, -1)
    perm_metrics = perm_metrics.rename(columns={"seed": "perm_id"}).sort_values(
        ["lineage", "perm_id"]
    )
    perm_metrics.to_csv(outdir / "permutation_null_distribution.csv", index=False)

    summaries = []
    for lineage in LINEAGES:
        obs_row = observed_metrics.loc[
            observed_metrics["method"].eq("PCA_mean")
            & observed_metrics["lineage"].eq(lineage)
            & observed_metrics["seed"].eq(SEED)
        ]
        obs_auc = float(obs_row["auroc"].iloc[0]) if len(obs_row) else np.nan
        sub = perm_metrics.loc[perm_metrics["lineage"].eq(lineage), "auroc"].dropna().to_numpy()
        p_one = float((1 + np.sum(sub >= obs_auc)) / (len(sub) + 1)) if len(sub) else np.nan
        p_two = (
            float((1 + np.sum(np.abs(sub - 0.5) >= abs(obs_auc - 0.5))) / (len(sub) + 1))
            if len(sub)
            else np.nan
        )
        summaries.append(
            {
                "lineage": lineage,
                "observed_pca_mean_auroc": obs_auc,
                "n_permutations": int(len(sub)),
                "null_auroc_mean": float(np.mean(sub)),
                "null_auroc_sd": float(np.std(sub, ddof=1)),
                "null_auroc_median": float(np.median(sub)),
                "null_auroc_q025": float(np.quantile(sub, 0.025)),
                "null_auroc_q975": float(np.quantile(sub, 0.975)),
                "empirical_p_one_sided": p_one,
                "empirical_p_two_sided": p_two,
            }
        )
    pd.DataFrame(summaries).to_csv(outdir / "permutation_null_summary.csv", index=False)


def run_subbag_sensitivity(
    representation: str,
    rows: pd.DataFrame,
    groups: dict,
    folds: list[dict],
    observed_metrics: pd.DataFrame,
    outdir: Path,
    *,
    max_cells: int = 256,
    n_repeats: int = N_SUBBAG_REPEATS,
) -> None:
    """Quantify full-bag vs `256 x 10` sub-bag OOF AUROC variance (`mean` and `attention`)."""
    sub_indices: dict[tuple[str, str], list[np.ndarray]] = {}
    sc_full: dict[str, dict[str, np.ndarray]] = {}
    for lineage, group in groups.items():
        if representation == "scGPT":
            sc_full[lineage] = scgpt_rep(group)
        for gsm, (bag, _, _) in group.items():
            sub_indices[(lineage, gsm)] = subsample_bag_indices(
                bag.obs.index, gsm, lineage, max_cells=max_cells, n_repeats=n_repeats, seed=SEED
            )

    draw_oof = []
    for r in range(n_repeats):
        for lineage, group in groups.items():
            eligible = rows[rows["sample_uid"].isin(group)]
            if representation == "PCA":
                sub_group = {}
                for gsm, (bag, x_log1p, comp) in group.items():
                    ix = sub_indices[(lineage, gsm)][r]
                    sub_group[gsm] = (bag, x_log1p[ix, :], comp)
            for fold in folds:
                f = fold["fold"]
                tr = eligible[eligible["donor_uid"].isin(fold["train_donors"])].reset_index(drop=True)
                te = eligible[eligible["donor_uid"].isin(fold["test_donors"])].reset_index(drop=True)
                if not len(te) or tr["tissue"].nunique() != 2:
                    continue
                design = build_design(tr)
                if not design.pair_donors:
                    continue
                if representation == "PCA":
                    h_sub, _ = pca_rep(sub_group, design.samples, dest=None)
                else:
                    h_sub = {
                        gsm: sc_full[lineage][gsm][sub_indices[(lineage, gsm)][r]]
                        for gsm in group
                    }
                for pooling in ("mean", "attention"):
                    model, fit = fit_model(
                        [h_sub[s] for s in design.samples],
                        design,
                        pooling=pooling,
                        lambda_pair=0.1,
                        lambda_l2=0.01,
                        seed=SEED,
                        standardize=True,
                        maxiter=500,
                    )
                    if not fit["converged"]:
                        continue
                    ids = te["sample_uid"].tolist()
                    meta = {"lineage": lineage, "fit": fit, "training_samples": design.samples}
                    pred, _ = predict_bags(model, meta, [h_sub[s] for s in ids], ids)
                    draw_oof.append(
                        te[["sample_uid", "donor_uid", "lesion_uid", "tissue"]].assign(
                            method=f"{representation}_{pooling}",
                            lineage=lineage,
                            fold=f,
                            seed=r,
                            logit=pred["logit"].to_numpy(),
                        )
                    )

    oof_df = pd.concat(draw_oof, ignore_index=True)
    _, draw_metrics, _ = summarize(oof_df, -1)
    draw_metrics = draw_metrics.rename(columns={"seed": "repeat_index"}).sort_values(
        ["lineage", "method", "repeat_index"]
    )
    draw_metrics.to_csv(outdir / "subbag_256x10_draws.csv", index=False)

    summary_rows = []
    for (lineage, method), grp in draw_metrics.groupby(["lineage", "method"], sort=True):
        full_row = observed_metrics.loc[
            observed_metrics["method"].eq(method)
            & observed_metrics["lineage"].eq(lineage)
            & observed_metrics["seed"].eq(SEED)
        ]
        full_auc = float(full_row["auroc"].iloc[0]) if len(full_row) else np.nan
        aucs = grp["auroc"].to_numpy(dtype=float)
        summary_rows.append(
            {
                "representation": representation,
                "lineage": lineage,
                "method": method,
                "full_bag_auroc": full_auc,
                "n_repeats": int(len(aucs)),
                "max_cells_subbag": int(max_cells),
                "sub256_mean_auroc": float(np.mean(aucs)),
                "sub256_sd_auroc": float(np.std(aucs, ddof=1)),
                "sub256_min_auroc": float(np.min(aucs)),
                "sub256_max_auroc": float(np.max(aucs)),
                "sub256_range_auroc": float(np.max(aucs) - np.min(aucs)),
                "delta_full_minus_sub256_mean": float(full_auc - np.mean(aucs)),
            }
        )
    pd.DataFrame(summary_rows).to_csv(outdir / "subbag_256x10_sensitivity.csv", index=False)


def run(args):
    OUT.mkdir(parents=True, exist_ok=True)
    MODELS.mkdir(parents=True, exist_ok=True)
    outdir = OUT / args.representation
    outdir.mkdir(exist_ok=True)

    universe_genes = list(load_shared_nuclear_universe(DEFAULT_UNIVERSE))
    all_program_defs = load_lineage_programs(universe_genes)
    normal_profiles = (
        load_normal_nonimmune_profiles(universe_genes)
        if args.representation == "PCA"
        else []
    )
    folds = json.loads(SPLIT.read_text(encoding="utf-8"))["folds"]
    rows, groups = library_inputs(write_selection=True)

    predictions, runs, perturb, attribution = [], [], [], []
    organ_audit_rows = []
    pca_fold_cache: dict[tuple[str, int], tuple[pd.DataFrame, pd.DataFrame, dict[str, np.ndarray]]] = {}

    for lineage, group in groups.items():
        eligible = rows[rows["sample_uid"].isin(group)]
        prog_defs = all_program_defs[lineage]
        s1_name, s2_name = DOMINANT_SUBTYPE_PAIRS[lineage]

        for fold in folds:
            f = fold["fold"]
            train_donor_set = set(fold["train_donors"])
            test_donor_set = set(fold["test_donors"])
            if train_donor_set & test_donor_set:
                raise ValueError("Donor leakage")
            training = eligible[
                eligible["donor_uid"].isin(train_donor_set)
            ].reset_index(drop=True)
            testing = eligible[
                eligible["donor_uid"].isin(test_donor_set)
            ].reset_index(drop=True)
            dest = MODELS / lineage.replace("/", "_") / f"fold{f}"
            dest.mkdir(parents=True, exist_ok=True)

            def record(method, seed, subset, logits):
                frame = subset[
                    ["sample_uid", "donor_uid", "lesion_uid", "tissue"]
                ].copy()
                frame = frame.assign(
                    method=method, lineage=lineage, fold=f, seed=seed, logit=logits
                )
                predictions.append(frame)

            if not len(testing) or training["tissue"].nunique() != 2:
                runs.append(
                    dict(
                        lineage=lineage,
                        fold=f,
                        method=args.representation,
                        status="insufficient support",
                    )
                )
                continue

            train = training["sample_uid"].tolist()
            test_samples = testing["sample_uid"].tolist()
            scored_fold = score_programs_train_locked(
                group,
                train,
                train + test_samples,
                prog_defs,
                universe_genes,
            )
            cell_scores_fold = scored_fold[0]

            if args.representation == "PCA":
                for kind in ("program_elastic", "composition_elastic", "QC_elastic"):
                    logits, ok, iterations, nparams = elastic(
                        group,
                        training,
                        test_samples,
                        kind,
                        dest / f"{kind}.npz",
                        lineage=lineage,
                        program_defs=prog_defs,
                        universe_genes=universe_genes,
                        precomputed_scores=scored_fold if kind == "program_elastic" else None,
                    )
                    runs.append(
                        dict(
                            method=kind,
                            lineage=lineage,
                            fold=f,
                            seed=SEED,
                            n_parameters=nparams,
                            iterations=iterations,
                            status="converged" if ok else "iteration limit",
                        )
                    )
                    if ok:
                        record(kind, SEED, testing, logits)

                vectors, naive_hvgs = pca_rep(group, train, dest / "pca.npz")
                pca_fold_cache[(lineage, f)] = (training, testing, vectors)

                organ_mask, mask_stats = build_fold_organ_mask(
                    normal_profiles, train_donor_set, test_donor_set, universe_genes
                )
                vectors_clean_organ, clean_hvgs = pca_rep(
                    group, train, dest / "pca_remove_organ.npz", mask_genes=organ_mask
                )
                displaced = sorted(set(naive_hvgs.tolist()) - set(clean_hvgs.tolist()))
                organ_audit_rows.append(
                    {
                        "lineage": lineage,
                        "fold": f,
                        **mask_stats,
                        "n_displaced_from_naive_top512": len(displaced),
                        "displaced_genes": "|".join(displaced),
                    }
                )
            else:
                vectors = scgpt_rep(group)
                vectors_clean_organ = None

            settings = [
                ("mean", 0.1, SEED, "mean"),
                ("attention", 0.1, SEED, "attention"),
                ("attention", 0.0, SEED, "no_pair"),
                ("attention", 0.1, SEED + 1, "attention"),
                ("attention", 0.1, SEED + 2, "attention"),
            ]
            if args.representation == "PCA":
                settings += [
                    ("mean", 0.1, SEED, "subtype_balanced"),
                    ("mean", 0.1, SEED, "pure_primary"),
                    ("mean", 0.1, SEED, "pure_secondary"),
                    ("mean", 0.1, SEED, "permuted"),
                    ("mean", 0.1, SEED, "remove_organ_markers"),
                ]

            for pooling, pair, seed, variant in settings:
                method = f"{args.representation}_{variant}"
                tr, te = training.copy(), testing.copy()
                h = vectors
                if variant == "subtype_balanced":
                    h = {
                        s: subtype_balanced_centroids(group[s][0], h[s], lineage)
                        for s in h
                    }
                    supported = [s for s, v in h.items() if len(v)]
                    tr = tr[tr["sample_uid"].isin(supported)].reset_index(drop=True)
                    te = te[te["sample_uid"].isin(supported)].reset_index(drop=True)
                elif variant in ("pure_primary", "pure_secondary"):
                    target_sub = s1_name if variant == "pure_primary" else s2_name
                    h = {
                        s: pure_subtype_slice(group[s][0], h[s], target_sub, min_cells=20)
                        for s in h
                    }
                    supported = [s for s, v in h.items() if len(v)]
                    tr = tr[tr["sample_uid"].isin(supported)].reset_index(drop=True)
                    te = te[te["sample_uid"].isin(supported)].reset_index(drop=True)
                elif variant == "remove_organ_markers":
                    h = vectors_clean_organ
                elif variant == "permuted":
                    rng_perm = np.random.default_rng(cell_seed(f"perm_balanced|fold{f}|{SEED}"))
                    tr = permute_training_donors(tr, rng_perm, balanced_half=True)

                info = dict(
                    method=method,
                    lineage=lineage,
                    fold=f,
                    seed=seed,
                    lambda_pair=pair,
                    n_train_bags=len(tr),
                    n_test_bags=len(te),
                    n_train_donors=tr["donor_uid"].nunique(),
                )
                if tr["tissue"].nunique() != 2 or not len(te):
                    runs.append(dict(info, status="insufficient support"))
                    continue
                design = build_design(tr)
                eff_pair = pair
                if eff_pair and not design.pair_donors:
                    if variant in ("pure_primary", "pure_secondary", "permuted"):
                        eff_pair = 0.0
                        info["lambda_pair"] = 0.0
                    else:
                        runs.append(dict(info, status="no training pair"))
                        continue

                start = time.perf_counter()
                model, fit = fit_model(
                    [h[s] for s in design.samples],
                    design,
                    pooling=pooling,
                    lambda_pair=eff_pair,
                    lambda_l2=0.01,
                    seed=seed,
                    standardize=True,
                    maxiter=500,
                )
                runs.append(
                    dict(
                        info,
                        n_parameters=fit["n_parameters"],
                        iterations=fit["iterations"],
                        seconds=time.perf_counter() - start,
                        status=fit["message"] if not fit["converged"] else "converged",
                    )
                )
                if not fit["converged"]:
                    continue
                meta = dict(
                    lineage=lineage,
                    fit=fit,
                    training_samples=design.samples,
                    training_donors=sorted(tr["donor_uid"].unique()),
                    test_donors=sorted(te["donor_uid"].unique()),
                    representation=args.representation,
                    variant=variant,
                    gene_universe=DEFAULT_UNIVERSE,
                )
                save_model(dest / f"{method}_{seed}.npz", model, meta)
                ids = te["sample_uid"].tolist()
                pred, _ = predict_bags(model, meta, [h[s] for s in ids], ids)
                record(method, seed, te, pred["logit"].to_numpy())
                if variant in ("mean", "attention") and seed == SEED:
                    for s in ids:
                        frame, stress = cell_attribution(
                            model,
                            group[s][0],
                            h[s],
                            cell_scores_fold[s],
                            pooling,
                            s,
                            seed,
                        )
                        attribution.append(
                            frame.assign(
                                method=method, fold=f, lineage=lineage, sample_uid=s
                            )
                        )
                        perturb.append(
                            dict(
                                method=method,
                                fold=f,
                                lineage=lineage,
                                sample_uid=s,
                                **stress,
                            )
                        )
            print(lineage, f, training["donor_uid"].nunique(), len(testing), flush=True)

    oof_this = pd.concat(predictions, ignore_index=True)
    pd.DataFrame(runs).to_csv(outdir / "benchmark_runs.csv", index=False)
    oof_this.to_csv(outdir / "oof_predictions.csv", index=False)
    pd.DataFrame(perturb).to_csv(outdir / "attribution_stress.csv", index=False)
    pd.concat(attribution, ignore_index=True).to_csv(
        outdir / "cell_contributions.csv.gz",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )

    _, local_metrics, _ = summarize(oof_this, SEED)
    if args.representation == "PCA":
        pd.DataFrame(organ_audit_rows).to_csv(
            outdir / "organ_marker_mask_audit.csv", index=False
        )
        run_permutation_null(
            pca_fold_cache,
            local_metrics,
            outdir,
            n_permutations=args.n_permutations,
        )

    run_subbag_sensitivity(
        args.representation,
        rows,
        groups,
        folds,
        local_metrics,
        outdir,
        max_cells=256,
        n_repeats=args.subbag_repeats,
    )

    available = [
        p for p in (OUT / "PCA", OUT / "scGPT") if (p / "oof_predictions.csv").exists()
    ]
    oof = pd.concat(
        [pd.read_csv(p / "oof_predictions.csv") for p in available], ignore_index=True
    )
    oof.to_csv(OUT / "oof_predictions.csv", index=False)
    pd.concat(
        [pd.read_csv(p / "benchmark_runs.csv") for p in available], ignore_index=True
    ).to_csv(OUT / "benchmark_runs.csv", index=False)
    patient, metrics, contrasts = summarize(oof, SEED)
    patient.to_csv(OUT / "patient_site_logits.csv", index=False)
    metrics.to_csv(OUT / "patient_metrics.csv", index=False)
    compare = OUT / "paired_model_comparisons"
    compare.mkdir(exist_ok=True)
    contrasts.to_csv(compare / "auc_differences.csv", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--representation", choices=["PCA", "scGPT"], required=True)
    parser.add_argument("--n-permutations", type=int, default=N_PERMUTATIONS)
    parser.add_argument("--subbag-repeats", type=int, default=N_SUBBAG_REPEATS)
    with threadpool_limits(limits=4):
        run(parser.parse_args())
