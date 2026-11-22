"""Small-sample robust paired pseudobulk negative-binomial and exact sign-flip contrasts in GSE245552.

Upgrades over V1:
1. Evaluates both broad immune lineages (`T/NK`, `Mono/Mac/DC`, `B/Plasma`) and hierarchical
   fine immune populations (`CD8 T cell`, `CD4 T cell`, `Macrophage`, `plasma cell`) from
   the remediated M2 pseudobulk manifest (`level in {'lineage', 'subtype', 'myeloid_fine'}`).
2. Attaches deterministic HGNC `canonical_symbol` alongside deposited `gene_symbol`.
3. Explicitly masks PyDESeq2 unconverged genes (`_genewise_converged == False`,
   `_MAP_converged == False`, or `_LFC_converged == False`) to `NaN` in p-value/FDR columns
   and logs unconverged counts plus dispersion-trend status (`parametric` vs `mean`).
4. Disables vacuous `cooks_filter` in paired `n=2`-per-donor designs and replaces
   anti-conservative asymptotic normal Wald p-values with residual-df Student-`t` Wald
   p-values (`df = N_pairs - 1`, retaining `pvalue_norm_wald` for methodological comparison).
5. Cross-checks every gene with a distribution-free full-enumeration (`2**N_pairs`) two-sided
   exact paired sign-flip permutation test on within-donor nuclear `log2(CP10k + 0.5)` differences.
"""

from __future__ import annotations

import argparse
import json
import resource
import sys
import time
from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from scipy import sparse
from scipy.stats import norm, t as student_t
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from crlm.hgnc_alias_map import get_hgnc_mapper

M2 = ROOT / "artifacts/m2/pseudobulk_manifest.csv"
AXIS = ROOT / "data/processed/m2/pseudobulk/gene_symbols.txt"
OUT = ROOT / "artifacts/m3/paired_counts"
MIN_CELLS_DEFAULT = 50
MIN_PAIRS_FIT = 3

# (level, parent_lineage, group_name, output_dir_slug)
POPULATION_SPECS: tuple[tuple[str, str, str, str], ...] = (
    ("lineage", "T/NK", "T/NK", "T_NK"),
    ("lineage", "Mono/Mac/DC", "Mono/Mac/DC", "Mono_Mac_DC"),
    ("lineage", "B/Plasma", "B/Plasma", "B_Plasma"),
    ("subtype", "T/NK", "CD8 T cell", "CD8_T_cell"),
    ("subtype", "T/NK", "CD4 T cell", "CD4_T_cell"),
    ("myeloid_fine", "Mono/Mac/DC", "Macrophage", "Macrophage"),
    ("subtype", "B/Plasma", "plasma cell", "plasma_cell"),
)

REFERENCE_READOUTS: dict[str, list[str]] = {
    "T/NK": [
        "CXCL13", "GZMB", "GZMK", "IFNG", "PDCD1", "TOX",
        "TCF7", "IL7R", "FOXP3", "CTLA4", "IL17A", "RORC",
        "TNFSF12", "STAT1", "IFIT1", "NKG7",
    ],
    "Mono/Mac/DC": [
        "SPP1", "TREM2", "GPNMB", "APOE", "APOC1", "FABP5",
        "LGALS3", "MMP9", "MARCO", "VSIG4", "TIMD4", "CD5L",
        "FOLR2", "CCL18", "CXCL9", "CXCL10", "STAT1", "CD74", "HLA-DRA",
    ],
    "B/Plasma": [
        "XBP1", "HSPA5", "HSP90B1", "MANF", "DERL3", "MZB1",
        "JCHAIN", "SDC1", "IGHA1", "IGHG1", "MS4A1", "CD79A", "HLA-DRA",
    ],
}


def bh_adjust(pvalues: np.ndarray) -> np.ndarray:
    """Exact Benjamini-Hochberg step-up FDR adjustment preserving NaNs."""
    p = np.asarray(pvalues, dtype=np.float64)
    out = np.full(p.shape, np.nan, dtype=np.float64)
    valid = np.flatnonzero(np.isfinite(p))
    m = len(valid)
    if m == 0:
        return out
    pv = np.clip(p[valid], 0.0, 1.0)
    order = np.argsort(pv, kind="mergesort")
    ranked = pv[order]
    ranks = np.arange(1, m + 1, dtype=np.float64)
    scaled = ranked * (float(m) / ranks)
    monotone = np.minimum.accumulate(scaled[::-1])[::-1]
    adjusted = np.clip(monotone, 0.0, 1.0)
    out_valid = np.empty(m, dtype=np.float64)
    out_valid[order] = adjusted
    out[valid] = out_valid
    return out


def exact_paired_signflip_test(
    counts: pd.DataFrame,
    meta: pd.DataFrame,
    *,
    pseudocount: float = 0.5,
    scale_factor: float = 10000.0,
    gene_chunk: int = 4096,
) -> pd.DataFrame:
    """Full-enumeration (`2**N_pairs`) two-sided exact paired sign-flip test on `log2(CP10k + 0.5)`.

    Under within-donor tissue exchangeability `H_0`, all `2**N_pairs` sign vectors
    `s in {-1, +1}^N` are equiprobable without parametric or asymptotic assumptions.
    """
    donors = sorted(meta["donor_uid"].unique())
    n_pairs = len(donors)
    if n_pairs < 1:
        raise ValueError("Require >= 1 paired donor for sign-flip test")
    if n_pairs > 16:
        raise ValueError(f"Exhaustive 2**N enumeration capped at N<=16, got N={n_pairs}")

    pt_idx = []
    lm_idx = []
    for d in donors:
        sub = meta.loc[meta["donor_uid"].eq(d)]
        pt_rows = sub.index[sub["tissue"].eq("PT")].tolist()
        lm_rows = sub.index[sub["tissue"].eq("LM")].tolist()
        if len(pt_rows) != 1 or len(lm_rows) != 1:
            raise ValueError(f"Expected exactly one PT and one LM pseudobulk for {d}")
        pt_idx.append(pt_rows[0])
        lm_idx.append(lm_rows[0])

    c_pt = counts.loc[pt_idx].to_numpy(dtype=np.float64)
    c_lm = counts.loc[lm_idx].to_numpy(dtype=np.float64)
    tot_pt = c_pt.sum(axis=1, keepdims=True)
    tot_lm = c_lm.sum(axis=1, keepdims=True)
    if (tot_pt <= 0).any() or (tot_lm <= 0).any():
        raise ValueError("Pseudobulk library has zero total nuclear counts")

    y_pt = np.log2(c_pt * (scale_factor / tot_pt) + pseudocount)
    y_lm = np.log2(c_lm * (scale_factor / tot_lm) + pseudocount)
    diffs = y_lm - y_pt  # shape: (n_pairs, n_genes)

    # Construct exhaustive (2**(n_pairs - 1), n_pairs) half-space signs with first sign fixed +1;
    # since |(-s) . d| == |s . d| identically, half-space enumeration gives the exact same
    # two-sided tail proportion in half the FLOPs and zero +/- symmetry rounding mismatch.
    n_half = 1 << (n_pairs - 1)
    if n_pairs == 1:
        signs = np.ones((1, 1), dtype=np.float64)
    else:
        grid = np.arange(n_half, dtype=np.uint32)[:, None]
        bits = (grid >> np.arange(n_pairs - 2, -1, -1, dtype=np.uint32)) & 1
        rest = np.where(bits == 0, 1.0, -1.0)
        signs = np.concatenate([np.ones((n_half, 1), dtype=np.float64), rest], axis=1)

    n_genes = diffs.shape[1]
    pvals = np.empty(n_genes, dtype=np.float64)
    for start in range(0, n_genes, gene_chunk):
        d_sub = diffs[:, start : start + gene_chunk]
        perm_abs = np.abs(signs @ d_sub)  # row 0 is exact observed (+1, ..., +1)
        obs_abs = perm_abs[0]
        tol = 1e-12 * np.maximum(1.0, obs_abs)
        exceed = (perm_abs >= (obs_abs[None, :] - tol)).sum(axis=0)
        pvals[start : start + d_sub.shape[1]] = exceed.astype(np.float64) / float(n_half)

    return pd.DataFrame(
        {
            "mean_log2_cp10k_diff": diffs.mean(axis=0),
            "median_log2_cp10k_diff": np.median(diffs, axis=0),
            "n_pos_donors_cp10k": (diffs > 0.0).sum(axis=0).astype(np.int32),
            "n_neg_donors_cp10k": (diffs < 0.0).sum(axis=0).astype(np.int32),
            "n_tie_donors_cp10k": (diffs == 0.0).sum(axis=0).astype(np.int32),
            "pvalue_signflip": pvals,
        },
        index=counts.columns,
    )


def collect(records: pd.DataFrame, symbols: Sequence[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Sum verified integer pseudobulk counts across lesions within `(donor_uid, tissue)`."""
    matrix, meta = [], []
    for (donor, tissue), group in records.groupby(["donor_uid", "tissue"], sort=True):
        total = np.zeros(len(symbols), dtype=np.int64)
        for r in group.itertuples():
            counts = sparse.load_npz(ROOT / r.count_matrix)[int(r.matrix_row), :].toarray().ravel()
            if (
                counts.dtype.kind not in "iu"
                or (counts < 0).any()
                or int(counts.sum()) != int(r.common_axis_visible_count_sum)
            ):
                raise ValueError(f"Invalid integer pseudobulk counts for {r.sample_uid} {r.group}")
            total += counts
        key = f"{donor}_{tissue}"
        meta.append(
            {
                "pseudobulk_uid": key,
                "donor_uid": donor,
                "tissue": tissue,
                "sample_uids": " ".join(group.sample_uid),
                "lesion_uids": " ".join(group.lesion_uid),
                "n_lesions": int(group.lesion_uid.nunique()),
                "n_cells": int(group.n_cells.sum()),
                "count_sum": int(total.sum()),
            }
        )
        matrix.append(total)
    meta_df = pd.DataFrame(meta).set_index("pseudobulk_uid")
    meta_df["tissue"] = pd.Categorical(meta_df.tissue, categories=["PT", "LM"])
    return pd.DataFrame(matrix, index=meta_df.index, columns=list(symbols)), meta_df


def fit_pseudobulk_table(
    counts: pd.DataFrame,
    meta: pd.DataFrame,
    symbols: Sequence[str],
    canonical_map: pd.Series,
    *,
    n_cpus: int = 2,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, object]]:
    """Core paired PyDESeq2 + Student-t Wald + exact paired sign-flip estimation."""
    n_pairs = int(meta.donor_uid.nunique())
    if len(meta) != 2 * n_pairs or not meta.groupby("donor_uid", observed=True).tissue.nunique().eq(2).all():
        raise ValueError("One PT and one LM pseudobulk per included patient are required")
    if n_pairs < MIN_PAIRS_FIT:
        raise ValueError(f"Require >= {MIN_PAIRS_FIT} paired donors for NB regression, got {n_pairs}")

    selected = (counts >= 10).sum(axis=0) >= n_pairs
    n_eligible = int(selected.sum())
    if n_eligible == 0:
        raise ValueError("Zero genes passed count eligibility filter")

    with threadpool_limits(limits=1):
        dds = DeseqDataSet(
            counts=counts.loc[:, selected],
            metadata=meta,
            design="~ donor_uid + tissue",
            size_factors_fit_type="poscounts",
            refit_cooks=False,
            n_cpus=n_cpus,
            quiet=True,
        )
        design = dds.obsm["design_matrix"]
        rank = int(np.linalg.matrix_rank(design.to_numpy()))
        residual_df = int(len(design) - rank)
        if list(design.index) != list(meta.index) or rank != design.shape[1] or residual_df != n_pairs - 1:
            raise ValueError("Singular design or unexpected residual degrees of freedom")
        dds.deseq2()
        stats = DeseqStats(
            dds,
            contrast=["tissue", "LM", "PT"],
            alpha=0.05,
            independent_filter=False,
            cooks_filter=False,
            quiet=True,
            n_cpus=n_cpus,
        )
        stats.summary()

    result = stats.results_df.reindex(symbols).copy()
    result.index.name = "gene_symbol"
    result.insert(0, "canonical_symbol", canonical_map.reindex(symbols).to_numpy(dtype=str))
    result["passed_count_filter"] = selected.reindex(symbols).fillna(False).to_numpy(dtype=bool)

    conv_cols = ["_genewise_converged", "_MAP_converged", "_LFC_converged"]
    for col in conv_cols:
        if col in dds.var.columns:
            result[col] = dds.var[col].reindex(symbols)
        else:
            result[col] = np.nan
    for extra_col in ("dispersions", "genewise_dispersions", "fitted_dispersions", "MAP_dispersions"):
        if extra_col in dds.var.columns:
            result[extra_col] = dds.var[extra_col].reindex(symbols)

    gw_ok = result["_genewise_converged"].eq(True)
    map_ok = result["_MAP_converged"].eq(True)
    lfc_ok = result["_LFC_converged"].eq(True)
    stat_vals = result["stat"].to_numpy(dtype=np.float64)
    stat_ok = np.isfinite(stat_vals)
    nb_converged = (result["passed_count_filter"] & gw_ok & map_ok & lfc_ok & stat_ok).to_numpy(dtype=bool)
    result["nb_all_converged"] = nb_converged

    # Preserve unconverged asymptotic output count only as diagnostic scalar; mask gene table p-values!
    n_legacy_unmasked_bh = int((stats.results_df["padj"] < 0.05).sum())

    pval_t = np.full(len(result), np.nan, dtype=np.float64)
    pval_norm = np.full(len(result), np.nan, dtype=np.float64)
    pval_t[nb_converged] = 2.0 * student_t.sf(np.abs(stat_vals[nb_converged]), df=residual_df)
    pval_norm[nb_converged] = 2.0 * norm.sf(np.abs(stat_vals[nb_converged]))

    padj_t = bh_adjust(pval_t)
    padj_norm = bh_adjust(pval_norm)

    # Primary columns use residual-df Student-t Wald on strictly converged genes
    result["pvalue"] = pval_t
    result["padj"] = padj_t
    result["pvalue_t"] = pval_t
    result["padj_t"] = padj_t
    result["pvalue_norm_wald"] = pval_norm
    result["padj_norm_wald"] = padj_norm

    # Exact 2**N paired sign-flip permutation test on nuclear log2(CP10k + 0.5)
    sf_df = exact_paired_signflip_test(counts, meta)
    for col in sf_df.columns:
        result[col] = sf_df[col].reindex(symbols).to_numpy()
    sf_elig_p = np.where(result["passed_count_filter"].to_numpy(dtype=bool), result["pvalue_signflip"].to_numpy(), np.nan)
    result["padj_signflip"] = bh_adjust(sf_elig_p)

    meta_out = meta.copy()
    meta_out["size_factor"] = dds.obs["size_factors"].to_numpy()

    el_var = dds.var
    n_gw_fail = int((el_var["_genewise_converged"] == False).sum()) if "_genewise_converged" in el_var else 0  # noqa: E712
    n_map_fail = int((el_var["_MAP_converged"] == False).sum()) if "_MAP_converged" in el_var else 0  # noqa: E712
    n_lfc_fail = int((el_var["_LFC_converged"] == False).sum()) if "_LFC_converged" in el_var else 0  # noqa: E712
    n_any_fail = int(n_eligible - int(nb_converged.sum()))
    disp_type = str(dds.uns.get("disp_function_type", "unknown"))

    summary = {
        "n_pairs": n_pairs,
        "n_samples": len(meta_out),
        "n_input_genes": len(symbols),
        "n_expression_eligible": n_eligible,
        "disp_function_type": disp_type,
        "n_genewise_unconverged": n_gw_fail,
        "n_map_unconverged": n_map_fail,
        "n_lfc_unconverged": n_lfc_fail,
        "n_any_unconverged_masked": n_any_fail,
        "n_converged_tested": int(nb_converged.sum()),
        "n_gene_p": int(result["pvalue_t"].notna().sum()),
        "n_BH_005": int((result["padj_t"] < 0.05).sum()),
        "n_BH_005_t": int((result["padj_t"] < 0.05).sum()),
        "n_BH_005_signflip": int((result["padj_signflip"] < 0.05).sum()),
        "n_raw_p001_signflip": int((sf_elig_p < 0.01).sum()),
        "n_BH_005_norm_wald_converged": int((result["padj_norm_wald"] < 0.05).sum()),
        "n_BH_005_legacy_unmasked_norm": n_legacy_unmasked_bh,
        "rank": rank,
        "residual_df": residual_df,
        "cooks_filter_mode": "disabled_paired_n2_per_donor",
        "convergence_flags_false": json.dumps(
            {
                "_genewise_converged": n_gw_fail,
                "_MAP_converged": n_map_fail,
                "_LFC_converged": n_lfc_fail,
            }
        ),
    }
    return result, meta_out, design, summary


def fit_and_write(
    records: pd.DataFrame,
    symbols: Sequence[str],
    canonical_map: pd.Series,
    branch: str,
    level: str,
    parent_lineage: str,
    population: str,
    slug: str,
) -> dict[str, object]:
    """Fit one `(population, branch)` contrast and write artifacts under `OUT / slug / branch`."""
    where = OUT / slug / branch
    where.mkdir(parents=True, exist_ok=True)
    n_pairs = int(records["donor_uid"].nunique())
    if n_pairs < MIN_PAIRS_FIT:
        return {
            "level": level,
            "lineage": parent_lineage,
            "population": population,
            "slug": slug,
            "analysis": branch,
            "status": f"skipped_fewer_than_{MIN_PAIRS_FIT}_pairs",
            "n_pairs": n_pairs,
            "n_samples": int(2 * n_pairs),
            "n_input_genes": len(symbols),
            "n_expression_eligible": 0,
            "disp_function_type": "NA",
            "n_genewise_unconverged": 0,
            "n_map_unconverged": 0,
            "n_lfc_unconverged": 0,
            "n_any_unconverged_masked": 0,
            "n_converged_tested": 0,
            "n_gene_p": 0,
            "n_BH_005": 0,
            "n_BH_005_t": 0,
            "n_BH_005_signflip": 0,
            "n_raw_p001_signflip": 0,
            "n_BH_005_norm_wald_converged": 0,
            "n_BH_005_legacy_unmasked_norm": 0,
            "rank": 0,
            "residual_df": max(0, n_pairs - 1),
            "cooks_filter_mode": "disabled_paired_n2_per_donor",
            "convergence_flags_false": "{}",
        }

    counts, meta = collect(records, symbols)
    result, meta_out, design, stats_summary = fit_pseudobulk_table(
        counts, meta, symbols, canonical_map, n_cpus=2
    )
    design.to_csv(where / "design.csv")
    result.to_csv(where / "gene_results.csv.gz")
    meta_out.to_csv(where / "libraries.csv")
    return {
        "level": level,
        "lineage": parent_lineage,
        "population": population,
        "slug": slug,
        "analysis": branch,
        "status": "completed",
        **stats_summary,
    }


def run() -> None:
    start = time.monotonic()
    started = datetime.now(timezone.utc).isoformat()
    symbols = AXIS.read_text(encoding="utf-8").splitlines()
    mapper = get_hgnc_mapper()
    norm_axis = mapper.normalize_gene_axis(symbols)
    canonical_map = pd.Series(
        norm_axis["canonical_symbol"].to_numpy(dtype=str), index=symbols, name="canonical_symbol"
    )

    tables = pd.read_csv(M2)
    tables = tables[(tables.tissue.isin(["PT", "LM"])) & (tables.n_cells >= MIN_CELLS_DEFAULT)].copy()
    manifest = pd.read_csv(ROOT / "data/manifests/m0/sample_manifest_validated.csv")
    if set(manifest.dataset) != {"GSE245552"} or not manifest.patient_mapping_verified.all():
        raise ValueError("Expected source-verified GSE245552 samples")
    table = tables.merge(
        manifest[["sample_uid", "original_specimen_id"]], on="sample_uid", validate="many_to_one"
    )

    OUT.mkdir(parents=True, exist_ok=True)
    all_rows: list[dict[str, object]] = []
    reference: list[pd.DataFrame] = []

    for level, parent_lineage, pop_name, slug in POPULATION_SPECS:
        subset = table[(table.level == level) & (table.lineage == parent_lineage) & (table.group == pop_name)]
        linked = subset.groupby("donor_uid").tissue.nunique()
        paired_donors = linked[linked == 2].index
        subset = subset[subset.donor_uid.isin(paired_donors)].sort_values(
            ["donor_uid", "tissue", "original_specimen_id"]
        )
        for branch_name, rows in [
            ("pooled_lesions", subset),
            ("first_eligible_specimen", subset.drop_duplicates(["donor_uid", "tissue"])),
        ]:
            print(
                f"[{level}] {pop_name} | {branch_name}: {len(rows)} GSMs, "
                f"{rows.donor_uid.nunique()} paired patients",
                flush=True,
            )
            row = fit_and_write(
                rows, symbols, canonical_map, branch_name, level, parent_lineage, pop_name, slug
            )
            all_rows.append(row)
            pd.DataFrame(all_rows).to_csv(OUT / "fit_summary.csv", index=False)
            if row["status"] == "completed":
                table_out = pd.read_csv(OUT / slug / branch_name / "gene_results.csv.gz")
                want_genes = set(REFERENCE_READOUTS[parent_lineage])
                ref = table_out[
                    table_out.canonical_symbol.isin(want_genes) | table_out.gene_symbol.isin(want_genes)
                ].copy()
                ref.insert(0, "analysis", branch_name)
                ref.insert(0, "population", pop_name)
                ref.insert(0, "lineage", parent_lineage)
                ref.insert(0, "level", level)
                reference.append(ref)

    if reference:
        pd.concat(reference, ignore_index=True).to_csv(OUT / "reference_readouts.csv", index=False)

    profile = {
        "started_utc": started,
        "completed_utc": datetime.now(timezone.utc).isoformat(),
        "elapsed_seconds": round(time.monotonic() - start, 3),
        "peak_parent_rss_bytes": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024,
        "peak_child_rss_bytes": resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss * 1024,
        "pydeseq2": version("pydeseq2"),
        "python": sys.version,
        "min_cells_per_library": MIN_CELLS_DEFAULT,
        "populations": [p[2] for p in POPULATION_SPECS],
        "inference_method": "Student_t_df_Npairs_minus_1_plus_exact_2powN_paired_signflip",
        "cooks_filter_mode": "disabled_paired_n2_per_donor",
        "command": " ".join([sys.executable, *sys.argv]),
    }
    (OUT / "profile.json").write_text(json.dumps(profile, indent=2) + "\n", encoding="utf-8")


def self_check() -> None:
    """In-memory verification of convergence masking, t-distribution Wald, and exact 2**N sign-flip."""
    rng = np.random.default_rng(20261121)
    n_pairs = 8
    n_genes = 120
    symbols = [f"GENE{i:03d}" for i in range(n_genes)]
    canonical_map = pd.Series(symbols, index=symbols)
    meta = pd.DataFrame(
        {
            "donor_uid": [f"P{i:02d}" for i in range(n_pairs) for _ in range(2)],
            "tissue": pd.Categorical(["PT", "LM"] * n_pairs, categories=["PT", "LM"]),
        },
        index=[f"P{i:02d}_{t}" for i in range(n_pairs) for t in ("PT", "LM")],
    )
    # Realistic negative-binomial pseudobulk counts (size=5 -> dispersion alpha=0.2) across log-normal means
    gene_means = np.exp(rng.uniform(np.log(25.0), np.log(800.0), size=n_genes))
    donor_mult = np.exp(rng.normal(0.0, 0.25, size=n_pairs))
    raw = np.empty((2 * n_pairs, n_genes), dtype=np.int64)
    size_r = 5.0
    for i in range(n_pairs):
        for t_idx in (0, 1):
            mu = gene_means * donor_mult[i]
            prob = size_r / (size_r + mu)
            raw[2 * i + t_idx] = rng.negative_binomial(size_r, prob)

    # Inject strong deterministic positive effect across all 8 pairs in GENE000
    for i in range(n_pairs):
        raw[2 * i, 0] = 20 + 2 * i
        raw[2 * i + 1, 0] = 220 + 15 * i
    # Inject sparse below-filter gene in last column
    raw[:, -1] = 1
    counts = pd.DataFrame(raw, index=meta.index, columns=symbols)

    res, meta_out, design, summary = fit_pseudobulk_table(
        counts, meta, symbols, canonical_map, n_cpus=1
    )
    assert summary["residual_df"] == n_pairs - 1 == 7
    assert summary["n_converged_tested"] >= 80
    assert not res.loc[symbols[-1], "passed_count_filter"]
    assert np.isnan(res.loc[symbols[-1], "pvalue_t"])
    # Check exact 2**8 sign-flip p-value for unanimous 8/8 positive difference: 2 / 256 == 1 / 128
    assert np.isclose(res.loc["GENE000", "pvalue_signflip"], 2.0 / (2**n_pairs), atol=1e-15)
    assert res.loc["GENE000", "n_pos_donors_cp10k"] == n_pairs
    assert np.isfinite(res.loc["GENE000", "pvalue_t"])
    # Check t-distribution p-value is strictly larger (more conservative) than normal Wald p-value
    conv_idx = res.index[res["nb_all_converged"] & (res["stat"].abs() > 0.5)]
    assert len(conv_idx) >= 20
    assert (res.loc[conv_idx, "pvalue_t"] > res.loc[conv_idx, "pvalue_norm_wald"]).all()
    # Verify any unconverged or ineligible gene is masked to NaN
    unconv = res.index[~res["nb_all_converged"]]
    assert res.loc[unconv, "pvalue"].isna().all() and res.loc[unconv, "padj"].isna().all()
    print(
        f"analyze_paired_counts --self-check PASSED: df={summary['residual_df']}, "
        f"converged={summary['n_converged_tested']}/{summary['n_expression_eligible']} "
        f"(disp={summary['disp_function_type']}), "
        f"GENE000 p_t={res.loc['GENE000', 'pvalue_t']:.3e} vs p_norm={res.loc['GENE000', 'pvalue_norm_wald']:.3e}, "
        f"p_signflip={res.loc['GENE000', 'pvalue_signflip']:.7f} (exact 1/128={1/128:.7f})"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-check", action="store_true", help="Run pure in-memory unit checks without writing files.")
    args = parser.parse_args()
    if args.self_check:
        self_check()
    else:
        run()
