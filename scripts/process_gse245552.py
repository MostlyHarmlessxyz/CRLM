"""Build the discovery-cohort sparse objects and patient-level coverage tables."""
import argparse
import json
import resource
import shlex
import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

import import_gse245552 as pilot
from crlm.annotation_m1 import HELD_OUT_CORROBORATION_GENES, HIERARCHICAL_SUBTYPES, LEVELS, classify
from crlm.hgnc_alias_map import HGNC_PROVENANCE_PATH, HGNC_TSV_PATH, SHARED_NUCLEAR_UNIVERSE

ROOT = pilot.ROOT
MANIFEST = ROOT / "data/manifests/m0/sample_manifest_validated.csv"
BASELINE = ROOT / "artifacts/m1/matrix_audit.csv"
PROCESSED = ROOT / "data/processed/m2"
OUT = ROOT / "artifacts/m2"
COMMON = pilot.AUDIT / "GSE245552_feature_intersection.txt"
NORMAL_EXTRA = pilot.AUDIT / "GSE245552_normal_only_features_353.txt"

MAIN = LEVELS[:3]
RESOLVED_IMMUNE = pilot.RESOLVED_IMMUNE
CORE_THREE_IMMUNE = pilot.CORE_THREE_IMMUNE
MYELOID_FINE_GROUPS = pilot.MYELOID_FINE_GROUPS
LIMITS = (20, 50, 100)

MIN_GENES_PRIMARY = 300
MIN_COUNTS_PRIMARY = 500
MAX_PCT_MT_VISIBLE = 20.0
LOW_MAD_SIGMA = 3.0
MULTIPLET_MAD_SIGMA = 3.5
MULTIPLET_SOFT_GENES = 6000
MULTIPLET_SOFT_COUNTS = 40000
MULTIPLET_HARD_GENES = 8500
MULTIPLET_HARD_COUNTS = 80000
MIN_GROUP_FOR_WITHIN_MAD = 50

# Held-out corroborating markers strictly disjoint from PANELS and AUXILIARY_MARKERS.
CHECK_GENES = HELD_OUT_CORROBORATION_GENES


def cohort():
    frame = pd.read_csv(MANIFEST)
    pairs = [pair for donor in frame.donor_uid.unique() for pair in pilot.select(MANIFEST, donor)]
    if len(pairs) != len(frame) or frame.sample_uid.duplicated().any() or frame.lesion_uid.duplicated().any():
        raise ValueError("GSM/lesion mapping requires a new source review")
    return dict(sorted((row["sample_uid"], (row, record)) for row, record in pairs))


def source_files(record):
    base = ROOT / "data/raw" / pilot.COHORT / record["gsm"]
    files = {k: base / record[k + "_url"].rsplit("/", 1)[-1] for k in ("matrix", "feature", "barcode")}
    for kind, path in files.items():
        if pilot.sha256(path) != record[kind + "_sha256"]:
            raise ValueError(f"Source differs from audited {kind}: {path}")
    return files


def qc_bounds(values, low_sigma=LOW_MAD_SIGMA, high_sigma=MULTIPLET_MAD_SIGMA, cap_high=None):
    logs = np.log1p(np.asarray(values, dtype=np.float64))
    center = float(np.median(logs))
    scale = 1.4826 * float(np.median(np.abs(logs - center)))
    lo = float(np.expm1(center - low_sigma * scale))
    hi = float(np.expm1(center + high_sigma * scale))
    if cap_high is not None:
        hi = min(hi, float(cap_high))
    return lo, hi


def prepare(row, record, files, m1_index):
    gsm = row["sample_uid"]
    if gsm in m1_index.index:
        audit = m1_index.loc[gsm]
        path = ROOT / audit.h5ad_path
        if pilot.sha256(path) != audit.h5ad_sha256:
            raise ValueError(f"M1 source changed: {path}")
        obj = ad.read_h5ad(path)
        obj.var = pilot.enrich_var_hgnc(obj.var)
        if "subtype_hierarchical" not in obj.obs.columns or "myeloid_fine" not in obj.obs.columns:
            labels, details = classify(obj.X, obj.var.gene_symbol.tolist(), pilot.MARKERS)
            for col in details.columns:
                obj.obs[col] = details[col].to_numpy()
            obj.obs["lineage_m1"] = labels
        if "pct_mt_scope" not in obj.obs.columns or "pct_mt_usable_for_qc" not in obj.obs.columns:
            qc_df, _ = pilot.qc_cells(obj.X, obj.var.gene_symbol)
            obj.obs["pct_mt_scope"] = qc_df["pct_mt_scope"].to_numpy()
            obj.obs["pct_mt_usable_for_qc"] = qc_df["pct_mt_usable_for_qc"].to_numpy()
        reference = []  # The M1 hash binds the original 16 source-coordinate checks.
        obj.uns["m1_object_sha256"] = audit.h5ad_sha256
    else:
        counts, genes, barcodes, reference = pilot.read_library(files, record)
        obs, _ = pilot.qc_cells(counts, genes.gene_symbol)
        labels, details = classify(counts, genes.gene_symbol.tolist(), pilot.MARKERS)
        obs = pd.concat([obs, details], axis=1)
        obs["lineage_m1"] = labels
        obs["barcode"] = barcodes
        for key in ("dataset", "donor_uid", "sample_uid", "lesion_uid", "original_specimen_id", "role", "tissue"):
            obs[key] = row[key]
        obs["capture_uid"] = gsm
        obs.index = pd.Index(gsm + ":" + barcodes, name="cell_uid")
        obj = ad.AnnData(X=counts, obs=obs, var=genes)

    obs = obj.obs
    if not obs.sample_uid.eq(gsm).all() or not obs.donor_uid.eq(row["donor_uid"]).all():
        raise ValueError(f"{gsm} object identity differs from manifest")
    obs["lineage_m2"] = obs.lineage_m1.astype(str)
    obs["subtype_m2"] = obs.subtype_hierarchical.astype(str)
    obs["myeloid_fine"] = obs.myeloid_fine.fillna("").astype(str)
    if obs.subtype_m2.str.contains("unresolved", case=False, regex=False).any():
        raise ValueError(f"{gsm} contains unresolved subtype labels after hierarchical disambiguation")

    low_counts, high_counts = qc_bounds(obs.n_counts_visible)
    low_genes, high_genes = qc_bounds(obs.n_genes_visible, cap_high=obj.n_vars)
    obs["qc_joint_low"] = (obs.n_counts_visible < low_counts) & (obs.n_genes_visible < low_genes)

    # Compute lower & upper log1p-MAD bounds within each broad lineage (fallback to library when <50 cells).
    within_low = pd.Series(False, index=obs.index)
    within_high = pd.Series(False, index=obs.index)
    for _, group in obs.groupby("lineage_m2", observed=True):
        if len(group) >= MIN_GROUP_FOR_WITHIN_MAD:
            lo_c, hi_c = qc_bounds(group.n_counts_visible)
            lo_g, hi_g = qc_bounds(group.n_genes_visible, cap_high=obj.n_vars)
        else:
            lo_c, hi_c = low_counts, high_counts
            lo_g, hi_g = low_genes, high_genes
        within_low.loc[group.index] = (group.n_counts_visible < lo_c) & (group.n_genes_visible < lo_g)
        within_high.loc[group.index] = (
            (group.n_counts_visible > hi_c)
            & (group.n_genes_visible > hi_g)
            & ((group.n_genes_visible > MULTIPLET_SOFT_GENES) | (group.n_counts_visible > MULTIPLET_SOFT_COUNTS))
        )

    hard_high = (obs.n_genes_visible > MULTIPLET_HARD_GENES) | (obs.n_counts_visible > MULTIPLET_HARD_COUNTS)
    obs["qc_within_lineage_low"] = within_low
    obs["qc_joint_high"] = within_high | hard_high

    empty = obs.n_counts_visible.eq(0)
    below_floor = (obs.n_genes_visible < MIN_GENES_PRIMARY) | (obs.n_counts_visible < MIN_COUNTS_PRIMARY)
    obs["qc_below_min_floor"] = below_floor

    mt_scope = str(obs["pct_mt_scope"].iloc[0])
    mt_usable = bool(obs["pct_mt_usable_for_qc"].iloc[0])
    high_mt = (obs.pct_mt > MAX_PCT_MT_VISIBLE) if mt_usable else pd.Series(False, index=obs.index)
    obs["qc_high_mt"] = high_mt

    excluded = empty | below_floor | within_low | obs["qc_joint_high"] | high_mt
    obs["qc_retained"] = ~excluded
    obs["qc_exclusion"] = np.select(
        [empty, below_floor, within_low, obs["qc_joint_high"], high_mt],
        ["empty", "below_min_300g_500c", "within_lineage_low", "joint_high_multiplet", "high_pct_mt_gt20"],
        default="",
    )

    # Precomputed sensitivity masks across gene-floor tiers and MT-filter ablation.
    obs["qc_pass_g200"] = ~(
        empty | (obs.n_genes_visible < 200) | (obs.n_counts_visible < 400) | within_low | obs["qc_joint_high"] | high_mt
    )
    obs["qc_pass_g300_primary"] = obs["qc_retained"]
    obs["qc_pass_g500"] = ~(
        empty | (obs.n_genes_visible < 500) | (obs.n_counts_visible < 500) | within_low | obs["qc_joint_high"] | high_mt
    )
    obs["qc_pass_g800"] = ~(
        empty | (obs.n_genes_visible < 800) | (obs.n_counts_visible < 800) | within_low | obs["qc_joint_high"] | high_mt
    )
    obs["qc_pass_no_mt_filter"] = ~(empty | below_floor | within_low | obs["qc_joint_high"])

    obj.uns["qc_limits"] = dict(
        low_counts=low_counts,
        low_genes=low_genes,
        high_counts=high_counts,
        high_genes=high_genes,
        min_genes_primary=MIN_GENES_PRIMARY,
        min_counts_primary=MIN_COUNTS_PRIMARY,
        max_pct_mt_visible=MAX_PCT_MT_VISIBLE,
        multiplet_hard_genes=MULTIPLET_HARD_GENES,
        multiplet_hard_counts=MULTIPLET_HARD_COUNTS,
        pct_mt_scope=mt_scope,
        pct_mt_usable_for_qc=mt_usable,
    )
    obs["published_barcode_prefix"] = obs.barcode.astype(str).str.rsplit("_", n=1).str[0]
    obj.uns["expression"] = {
        "location": "X",
        "scale": "raw_counts",
        "genes_prefiltered": True,
        "gene_universe_id": row["gene_universe_id"],
        "source": record["source_record"],
        "gene_mapping": "exact deposited ID/symbol in var.gene_symbol; deterministic HGNC canonical symbol in var.canonical_symbol",
        "pct_mt_scope": mt_scope,
        "pct_mt_usable_for_qc": mt_usable,
        "capture_protocol": "B07 cohort CD45+ and CD45- magnetic sorting, 2:1 remix",
        "collection_window": "cohort synchronous; exact paired dates unavailable",
    }
    return obj, reference


def build(gsm, row, record, m1_index, code, expected_axes, expected_canon_axes):
    tick = time.perf_counter()
    started = pilot.utc_now()
    files = source_files(record)
    obj, reference = prepare(row, record, files, m1_index)
    if (set(obj.var.gene_symbol) != expected_axes[row["tissue"]]
            or set(obj.var.canonical_symbol) != expected_canon_axes[row["tissue"]]):
        raise ValueError(f"Unexpected measured or canonical gene set for {gsm}")
    obj.uns["provenance_m2"] = dict(code, **{k: record[k] for k in (
        "barcode_url", "feature_url", "matrix_url", "barcode_sha256", "feature_sha256", "matrix_sha256")})
    path = PROCESSED / f"{gsm}.h5ad"
    pilot.write_checked(obj, path, reference)
    metric = {
        "sample_uid": gsm,
        "source_coordinate_checks_including_m1": 16,
        "started_utc": started,
        "ended_utc": pilot.utc_now(),
        "elapsed_seconds": round(time.perf_counter() - tick, 3),
        "process_peak_rss_bytes": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024,
        "raw_bytes": sum(f.stat().st_size for f in files.values()),
        "h5ad_bytes": path.stat().st_size,
        "h5ad_sha256": pilot.sha256(path),
        "command": shlex.join([sys.executable, *sys.argv]),
    }
    (PROCESSED / "runs" / f"{gsm}.json").write_text(json.dumps(metric, indent=2) + "\n")
    print(f"{gsm} {obj.n_obs} input, {obj.obs.qc_retained.sum()} retained, {metric['elapsed_seconds']} s",
          flush=True)


def group_definitions():
    groups = [("lineage", level, level) for level in LEVELS]
    for parent in LEVELS:
        for sub in HIERARCHICAL_SUBTYPES[parent]:
            groups.append(("subtype", parent, sub))
    for mf in MYELOID_FINE_GROUPS:
        groups.append(("myeloid_fine", "Mono/Mac/DC", mf))
    return groups


def sample_report(gsm, row, obj, common, pb_dir=None):
    obs = obj.obs
    sample = {k: row[k] for k in ("donor_uid", "tissue", "lesion_uid", "sample_uid")}
    broad = obs.lineage_m2.astype(str).to_numpy()
    sub = obs.subtype_m2.astype(str).to_numpy()
    mf_arr = obs.myeloid_fine.astype(str).to_numpy()
    keep = obs.qc_retained.to_numpy(dtype=bool)
    nonempty = obs.n_counts_visible.gt(0).to_numpy()
    library_only = nonempty & ~obs.qc_joint_low.to_numpy(dtype=bool)
    pass_g200 = obs.qc_pass_g200.to_numpy(dtype=bool)
    pass_g500 = obs.qc_pass_g500.to_numpy(dtype=bool)
    pass_g800 = obs.qc_pass_g800.to_numpy(dtype=bool)
    pass_no_mt = obs.qc_pass_no_mt_filter.to_numpy(dtype=bool)

    n_total_imm = int((keep & np.isin(broad, RESOLVED_IMMUNE)).sum())
    n_core3_imm = int((keep & np.isin(broad, CORE_THREE_IMMUNE)).sum())

    excl_floor = int(obs.qc_exclusion.eq("below_min_300g_500c").sum())
    excl_wlow = int(obs.qc_exclusion.eq("within_lineage_low").sum())
    excl_jhigh = int(obs.qc_exclusion.eq("joint_high_multiplet").sum())
    excl_hmt = int(obs.qc_exclusion.eq("high_pct_mt_gt20").sum())
    excl_empty = int(obs.qc_exclusion.eq("empty").sum())

    metrics = dict(
        sample,
        input_cells=len(obs),
        retained_cells=int(keep.sum()),
        excluded_total=int((~keep).sum()),
        excluded_empty=excl_empty,
        excluded_joint_low=excl_floor + excl_wlow,
        excluded_below_min_300g_500c=excl_floor,
        excluded_within_lineage_low=excl_wlow,
        excluded_joint_high_multiplet=excl_jhigh,
        excluded_high_pct_mt_gt20=excl_hmt,
        flagged_library_joint_low=int(obs.qc_joint_low.sum()),
        flagged_joint_high=int(obs.qc_joint_high.sum()),
        n_total_resolved_immune=n_total_imm,
        n_core_three_lineage_immune=n_core3_imm,
        nonimmune_retained=int((keep & (broad == "Nonimmune")).sum()),
        mixed_retained=int((keep & (broad == "Mixed")).sum()),
        unknown_retained=int((keep & (broad == "Unknown")).sum()),
        native_genes=obj.n_vars,
        visible_counts=int(obs.n_counts_visible.sum()),
        mt_available=int(obs.pct_mt.notna().sum()),
        **obj.uns["qc_limits"],
    )
    for key in ("n_counts_visible", "n_genes_visible", "pct_mt"):
        for percentile in (.01, .25, .5, .75, .99):
            metrics[f"{key}_p{int(percentile*100):02d}"] = (
                float(obs[key].quantile(percentile)) if obs[key].notna().any() else np.nan
            )

    ids = obj.var_names.get_indexer(common)
    if (ids < 0).any():
        raise ValueError(f"{gsm} missing common-axis genes")
    groups = group_definitions()
    manifest, check_rows = [], []
    weights_i, weights_j = [], []
    for ix, (kind, lineage, name) in enumerate(groups):
        if kind == "lineage":
            labeled = broad == name
        elif kind == "subtype":
            labeled = (broad == lineage) & (sub == name)
        elif kind == "myeloid_fine":
            labeled = (broad == "Mono/Mac/DC") & (sub == "monocyte/macrophage") & (mf_arr == name)
        else:
            raise ValueError(f"Unknown group level: {kind}")

        indexes = np.flatnonzero(labeled & keep)
        weights_i.extend([ix] * len(indexes))
        weights_j.extend(indexes)
        n_grp = len(indexes)
        p_tot = (n_grp / n_total_imm) if (lineage in RESOLVED_IMMUNE and n_total_imm > 0) else np.nan
        p_c3 = (n_grp / n_core3_imm) if (lineage in CORE_THREE_IMMUNE and n_core3_imm > 0) else np.nan
        entry = dict(
            sample,
            level=kind,
            lineage=lineage,
            group=name,
            n_cells=n_grp,
            input_group_cells=int(labeled.sum()),
            n_retained_sample=int(keep.sum()),
            n_total_resolved_immune=n_total_imm,
            n_core_three_lineage_immune=n_core3_imm,
            prop_total_resolved_immune=p_tot,
            prop_core_three_lineage_immune=p_c3,
            n_g200=int((labeled & pass_g200).sum()),
            n_g500=int((labeled & pass_g500).sum()),
            n_g800=int((labeled & pass_g800).sum()),
            n_no_mt_filter=int((labeled & pass_no_mt).sum()),
            n_nonempty=int((labeled & nonempty).sum()),
            n_sample_nonempty=int(nonempty.sum()),
            n_library_qc=int((labeled & library_only).sum()),
            n_sample_library_qc=int(library_only.sum()),
            count_matrix=f"data/processed/m2/pseudobulk/{gsm}.npz",
            matrix_row=ix,
        )
        for limit in LIMITS:
            entry[f"eligible_{limit}"] = n_grp >= limit
        manifest.append(entry)

    weights = sparse.csr_matrix(
        (np.ones(len(weights_j), dtype=np.int64), (weights_i, weights_j)),
        shape=(len(groups), len(obs)),
    )
    counts = (weights @ obj.X)[:, ids].tocsr()
    sums = np.asarray(counts.sum(axis=1)).ravel()

    n_lin = len(LEVELS)
    n_sub = sum(len(HIERARCHICAL_SUBTYPES[l]) for l in LEVELS)
    if (sums[:n_lin].sum() != sums[n_lin:n_lin + n_sub].sum()
            or sum(r["n_cells"] for r in manifest[:n_lin]) != int(keep.sum())
            or sum(r["n_cells"] for r in manifest[n_lin:n_lin + n_sub]) != int(keep.sum())):
        raise ValueError(f"{gsm} nonexhaustive broad/subtype partition")

    for entry, count in zip(manifest, sums):
        entry["common_axis_visible_count_sum"] = int(count)

    # Verify per-parent subtype additivity and myeloid_fine additivity.
    for parent in LEVELS:
        p_row = next(r for r in manifest[:n_lin] if r["group"] == parent)
        s_rows = [r for r in manifest[n_lin:n_lin + n_sub] if r["lineage"] == parent]
        if (sum(r["n_cells"] for r in s_rows) != p_row["n_cells"]
                or sum(r["common_axis_visible_count_sum"] for r in s_rows) != p_row["common_axis_visible_count_sum"]):
            raise ValueError(f"{gsm} subtype partition mismatch inside {parent}")
    mm_row = next(r for r in manifest[n_lin:n_lin + n_sub] if r["group"] == "monocyte/macrophage")
    mf_rows = manifest[n_lin + n_sub:]
    if (sum(r["n_cells"] for r in mf_rows) != mm_row["n_cells"]
            or sum(r["common_axis_visible_count_sum"] for r in mf_rows) != mm_row["common_axis_visible_count_sum"]):
        raise ValueError(f"{gsm} myeloid_fine partition mismatch inside monocyte/macrophage")

    if pb_dir is not None:
        sparse.save_npz(pb_dir / f"{gsm}.npz", counts, compressed=True)

    resolved = keep & ~np.isin(broad, ("Mixed", "Unknown"))
    for lineage, genes in CHECK_GENES.items():
        matches = obj.var_names.get_indexer(genes)
        if (matches < 0).any():
            raise ValueError(f"{gsm} missing prespecified corroboration marker")
        hit = np.asarray(obj.X[:, matches].getnnz(axis=1)) > 0
        target, others = resolved & (broad == lineage), resolved & (broad != lineage)
        check_rows.append(dict(
            sample,
            lineage=lineage,
            heldout_genes=" ".join(genes),
            n_labeled=int(target.sum()),
            n_labeled_detected=int((hit & target).sum()),
            n_other_resolved=int(others.sum()),
            n_other_detected=int((others & hit).sum()),
        ))
    return metrics, manifest, check_rows


def compute_patient_tables(coverage, donors):
    keys = ["donor_uid", "tissue", "level", "lineage", "group"]
    tissue_rows = []
    grouped = coverage.groupby(keys, sort=False)
    for donor in donors:
        for tissue in pilot.TISSUES:
            for kind, lineage, name in group_definitions():
                key = (donor, tissue, kind, lineage, name)
                rows = grouped.get_group(key) if key in grouped.indices else None
                data = dict(zip(keys, key))
                data["n_sampled_lesions"] = 0 if rows is None else len(rows)
                data["cells_all_sampled_lesions"] = np.nan if rows is None else int(rows.n_cells.sum())
                for limit in LIMITS:
                    valid = [] if rows is None else rows.loc[rows.n_cells.ge(limit), "sample_uid"].tolist()
                    data[f"n_eligible_lesions_{limit}"] = len(valid)
                    data[f"eligible_gsms_{limit}"] = "|".join(valid)
                data["status_50"] = (
                    "tissue_unavailable" if rows is None else
                    "no_labeled_cells" if rows.n_cells.sum() == 0 else
                    "below_lesion_threshold" if not data["n_eligible_lesions_50"] else "eligible"
                )
                tissue_rows.append(data)
    patient = pd.DataFrame(tissue_rows)

    paired = []
    for (kind, lineage, name), frame in patient.groupby(["level", "lineage", "group"], sort=False):
        for limit in LIMITS:
            supported = frame.pivot(index="donor_uid", columns="tissue", values=f"n_eligible_lesions_{limit}").gt(0)
            for contrast, tissues in (("PT_LM", ["PT", "LM"]), ("NC_NL", ["NC", "NL"]), ("Q", list(pilot.TISSUES))):
                if not set(tissues).issubset(supported.columns):
                    valid = []
                else:
                    valid = supported.index[supported[tissues].all(axis=1)].tolist()
                paired.append(dict(
                    level=kind,
                    lineage=lineage,
                    group=name,
                    min_cells_per_lesion=limit,
                    contrast=contrast,
                    n_donors=len(valid),
                    donors="|".join(valid),
                ))

    policies = (
        ("primary_min_genes_300", "n_cells"),
        ("min_genes_200", "n_g200"),
        ("min_genes_500", "n_g500"),
        ("min_genes_800_align164", "n_g800"),
        ("no_mt_filter", "n_no_mt_filter"),
        ("nonempty", "n_nonempty"),
        ("library_only", "n_library_qc"),
    )
    comparison = []
    broad = coverage.loc[coverage.level.eq("lineage")]
    for policy, count in policies:
        for level, frame in broad.groupby("lineage", sort=False):
            top = frame.pivot_table(index="donor_uid", columns="tissue", values=count, aggfunc="max")
            for limit in LIMITS:
                for contrast, tissues in (("PT_LM", ["PT", "LM"]), ("NC_NL", ["NC", "NL"]), ("Q", list(pilot.TISSUES))):
                    if not set(tissues).issubset(top.columns):
                        valid = []
                    else:
                        valid = top.index[top[tissues].ge(limit).all(axis=1)].tolist()
                    comparison.append(dict(
                        qc_policy=policy,
                        lineage=level,
                        min_cells_per_lesion=limit,
                        contrast=contrast,
                        n_donors=len(valid),
                        donors="|".join(valid),
                    ))
    return patient, pd.DataFrame(paired), pd.DataFrame(comparison)


def patient_tables(coverage, donors):
    patient, paired, sensitivity = compute_patient_tables(coverage, donors)
    patient.to_csv(OUT / "donor_tissue_lineage_coverage.csv", index=False)
    paired.to_csv(OUT / "paired_support.csv", index=False)
    sensitivity.to_csv(OUT / "qc_sensitivity.csv", index=False)


def qc_plot(qc):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(11, 7), sharey=True, layout="constrained")
    colors = {"PT": "#96402d", "LM": "#305d92", "NC": "#bd8550", "NL": "#47876f"}
    qc = qc.sort_values(["tissue", "sample_uid"]).reset_index(drop=True)
    y = np.arange(len(qc))
    for ax, col, title in zip(axes, ("n_counts_visible", "n_genes_visible"), ("Visible counts per cell", "Detected genes per cell")):
        color = qc.tissue.map(colors)
        ax.hlines(y, qc[col + "_p01"], qc[col + "_p99"], color=color, linewidth=.7)
        ax.hlines(y, qc[col + "_p25"], qc[col + "_p75"], color=color, linewidth=3.1)
        ax.scatter(qc[col + "_p50"], y, color="black", s=6)
        ax.set_xscale("log")
        ax.set_xlabel(title)
        ax.grid(axis="x", alpha=.2)
    axes[0].set_yticks(y, (qc.tissue + " " + qc.sample_uid), fontsize=6.5)
    fig.suptitle("Discovery libraries | P1–P99, interquartile range and median")
    fig.savefig(OUT / "qc_distributions.png", dpi=150)
    plt.close(fig)


def summarize(samples, common, code):
    start = time.perf_counter()
    pb_dir = PROCESSED / "pseudobulk"
    pb_dir.mkdir(exist_ok=True)
    (pb_dir / "gene_symbols.txt").write_text("\n".join(common) + "\n")
    all_qc, all_groups, checks, runs, sources = [], [], [], [], []
    prefix_owner, barcode_owner = {}, set()
    for gsm, (row, record) in samples.items():
        path = PROCESSED / f"{gsm}.h5ad"
        run = json.loads((PROCESSED / "runs" / f"{gsm}.json").read_text())
        if pilot.sha256(path) != run["h5ad_sha256"]:
            raise ValueError(f"Output changed after readback: {gsm}")
        obj = ad.read_h5ad(path)
        if any(obj.uns["provenance_m2"].get(k) != v for k, v in code.items()):
            raise ValueError(f"Different analysis version for {gsm}")
        if obj.n_obs != record["n_cells"] or obj.X.dtype != np.int32 or not sparse.isspmatrix_csr(obj.X):
            raise ValueError(f"Cell axis or raw-count representation changed for {gsm}")
        for prefix in obj.obs.published_barcode_prefix.unique():
            if prefix in prefix_owner:
                raise ValueError(f"Review repeated capture prefix {prefix}: {gsm}, {prefix_owner[prefix]}")
            prefix_owner[prefix] = gsm
        cells = set(obj.obs.barcode)
        if len(cells) != obj.n_obs or not cells.isdisjoint(barcode_owner):
            raise ValueError(f"Repeated full source barcode requires capture review: {gsm}")
        barcode_owner.update(cells)
        metrics, groups, check_rows = sample_report(gsm, row, obj, common, pb_dir=pb_dir)
        metrics["capture_prefixes"] = "|".join(sorted(obj.obs.published_barcode_prefix.unique()))
        all_qc.append(metrics)
        all_groups.extend(groups)
        checks.extend(check_rows)
        runs.append(dict(run, pseudobulk_sha256=pilot.sha256(pb_dir / f"{gsm}.npz")))
        for kind in ("barcode", "feature", "matrix"):
            sources.append(dict(sample_uid=gsm, kind=kind, url=record[kind + "_url"], sha256=record[kind + "_sha256"]))
    qc, groups = pd.DataFrame(all_qc), pd.DataFrame(all_groups)
    qc.to_csv(OUT / "qc_and_exclusions.csv", index=False)
    groups.to_csv(OUT / "pseudobulk_manifest.csv", index=False)
    pd.DataFrame(checks).to_csv(OUT / "annotation_marker_check.csv", index=False)
    pd.DataFrame(runs).to_csv(OUT / "object_manifest.csv", index=False)
    pd.DataFrame(sources).to_csv(OUT / "source_inputs.csv", index=False)
    patient_tables(groups, sorted(qc.donor_uid.unique()))
    qc_plot(qc)
    profile = dict(
        ended_utc=pilot.utc_now(),
        command=shlex.join([sys.executable, *sys.argv]),
        source_versions=code,
        python=sys.version,
        build_seconds_sum=sum(r["elapsed_seconds"] for r in runs),
        summary_seconds=time.perf_counter() - start,
        summary_peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024,
        max_build_peak_rss_bytes=max(r["process_peak_rss_bytes"] for r in runs),
        source_bytes=sum(r["raw_bytes"] for r in runs),
        h5ad_bytes=sum(r["h5ad_bytes"] for r in runs),
        pseudobulk_bytes=sum(p.stat().st_size for p in pb_dir.iterdir()),
    )
    (OUT / "resource_profile.json").write_text(json.dumps(profile, indent=2) + "\n")
    print(f"Accepted {len(qc)} GSMs, {qc.donor_uid.nunique()} patients, {qc.retained_cells.sum()} cells", flush=True)


def run_dry_check(samples, target_gsms, common, expected_axes, expected_canon_axes):
    """Execute full in-memory M2 QC, dual-denominator, and pseudobulk assertions with zero disk writes."""
    t0 = time.perf_counter()
    m1_index = pd.read_csv(BASELINE).set_index("sample_uid")
    all_qc, all_groups, checks = [], [], []
    prefix_owner, barcode_owner = {}, set()
    for gsm in target_gsms:
        if gsm not in samples:
            raise ValueError(f"Unknown GSM in --samples: {gsm}")
        row, record = samples[gsm]
        files = source_files(record)
        obj, _ = prepare(row, record, files, m1_index)
        if (set(obj.var.gene_symbol) != expected_axes[row["tissue"]]
                or set(obj.var.canonical_symbol) != expected_canon_axes[row["tissue"]]):
            raise ValueError(f"Unexpected measured or canonical gene set for {gsm}")
        for prefix in obj.obs.published_barcode_prefix.unique():
            if prefix in prefix_owner:
                raise ValueError(f"Repeated capture prefix {prefix}: {gsm}, {prefix_owner[prefix]}")
            prefix_owner[prefix] = gsm
        cells = set(obj.obs.barcode)
        if len(cells) != obj.n_obs or not cells.isdisjoint(barcode_owner):
            raise ValueError(f"Repeated full source barcode in {gsm}")
        barcode_owner.update(cells)
        metrics, groups, check_rows = sample_report(gsm, row, obj, common, pb_dir=None)
        all_qc.append(metrics)
        all_groups.extend(groups)
        checks.extend(check_rows)

    qc_df = pd.DataFrame(all_qc)
    groups_df = pd.DataFrame(all_groups)
    _, paired_df, sens_df = compute_patient_tables(groups_df, sorted(qc_df.donor_uid.unique()))

    # Verify high_genes cap across every single evaluated library
    if (qc_df["high_genes"] > qc_df["native_genes"]).any():
        raise AssertionError("high_genes exceeded native_genes")

    # Verify dual immune proportion sums == 1.0 on every sample with positive denominator
    broad_df = groups_df.loc[groups_df["level"].eq("lineage")]
    for gsm, grp in broad_df.groupby("sample_uid"):
        if grp["n_total_resolved_immune"].iloc[0] > 0:
            s_tot = grp.loc[grp["lineage"].isin(RESOLVED_IMMUNE), "prop_total_resolved_immune"].sum()
            if not np.isclose(s_tot, 1.0, atol=1e-12):
                raise AssertionError(f"{gsm} prop_total_resolved_immune sum={s_tot} != 1.0")
        if grp["n_core_three_lineage_immune"].iloc[0] > 0:
            s_c3 = grp.loc[grp["lineage"].isin(CORE_THREE_IMMUNE), "prop_core_three_lineage_immune"].sum()
            if not np.isclose(s_c3, 1.0, atol=1e-12):
                raise AssertionError(f"{gsm} prop_core_three_lineage_immune sum={s_c3} != 1.0")
        if grp.loc[~grp["lineage"].isin(RESOLVED_IMMUNE), "prop_total_resolved_immune"].notna().any():
            raise AssertionError(f"{gsm} Nonimmune/Mixed/Unknown received non-NaN immune proportion")

    summary = {
        "mode": "dry_run_in_memory_zero_writes",
        "n_gsms": int(len(qc_df)),
        "n_donors": int(qc_df["donor_uid"].nunique()),
        "input_cells": int(qc_df["input_cells"].sum()),
        "retained_cells": int(qc_df["retained_cells"].sum()),
        "excluded_total": int(qc_df["excluded_total"].sum()),
        "excluded_pct": round(100.0 * float(qc_df["excluded_total"].sum()) / float(qc_df["input_cells"].sum()), 4),
        "exclusions_breakdown": {
            "empty": int(qc_df["excluded_empty"].sum()),
            "below_min_300g_500c": int(qc_df["excluded_below_min_300g_500c"].sum()),
            "within_lineage_low": int(qc_df["excluded_within_lineage_low"].sum()),
            "joint_high_multiplet": int(qc_df["excluded_joint_high_multiplet"].sum()),
            "high_pct_mt_gt20": int(qc_df["excluded_high_pct_mt_gt20"].sum()),
        },
        "denominators_sum": {
            "total_resolved_immune": int(qc_df["n_total_resolved_immune"].sum()),
            "core_three_lineage_immune": int(qc_df["n_core_three_lineage_immune"].sum()),
            "nonimmune_retained": int(qc_df["nonimmune_retained"].sum()),
            "mixed_retained": int(qc_df["mixed_retained"].sum()),
            "unknown_retained": int(qc_df["unknown_retained"].sum()),
        },
        "max_high_genes_threshold": round(float(qc_df["high_genes"].max()), 2),
        "mt_scopes": qc_df.groupby("pct_mt_scope")["sample_uid"].count().to_dict(),
        "paired_PT_LM_50_main": (
            paired_df.loc[
                paired_df["level"].eq("lineage")
                & paired_df["group"].isin(MAIN)
                & paired_df["min_cells_per_lesion"].eq(50)
                & paired_df["contrast"].eq("PT_LM"),
                ["group", "n_donors"],
            ]
            .set_index("group")["n_donors"]
            .to_dict()
        ),
        "paired_Q_50_main": (
            paired_df.loc[
                paired_df["level"].eq("lineage")
                & paired_df["group"].isin(MAIN)
                & paired_df["min_cells_per_lesion"].eq(50)
                & paired_df["contrast"].eq("Q"),
                ["group", "n_donors"],
            ]
            .set_index("group")["n_donors"]
            .to_dict()
        ),
        "qc_sensitivity_PT_LM_50": {
            pol: grp.set_index("lineage")["n_donors"].to_dict()
            for pol, grp in sens_df.loc[
                sens_df["lineage"].isin(MAIN)
                & sens_df["min_cells_per_lesion"].eq(50)
                & sens_df["contrast"].eq("PT_LM")
            ].groupby("qc_policy", sort=False)
        },
        "elapsed_seconds": round(time.perf_counter() - t0, 3),
        "peak_rss_mib": round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0, 2),
    }
    print(json.dumps(summary, indent=2, ensure_ascii=False), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--samples", nargs="+", help="Build specified GSMs; run --reports-only after all are ready.")
    mode.add_argument("--reports-only", action="store_true", help="Rebuild cohort reports from verified M2 objects.")
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Run full QC, annotation, dual-denominator, and pseudobulk checks in memory without writing to M2.",
    )
    args = parser.parse_args()
    if args.dry_run and args.reports_only:
        parser.error("--dry-run cannot be combined with --reports-only")

    samples = cohort()
    common = COMMON.read_text().splitlines()
    normal = set(common) | set(NORMAL_EXTRA.read_text().splitlines())
    expected_axes = {t: set(common) if t in ("PT", "LM") else normal for t in pilot.TISSUES}
    canon_tumor = set(SHARED_NUCLEAR_UNIVERSE["245_tumor_canonical_nuclear"])
    canon_normal = canon_tumor | set(SHARED_NUCLEAR_UNIVERSE["excluded_353_canonical"])
    expected_canon_axes = {t: canon_tumor if t in ("PT", "LM") else canon_normal for t in pilot.TISSUES}

    if args.dry_run:
        target_gsms = args.samples or list(samples.keys())
        run_dry_check(samples, target_gsms, common, expected_axes, expected_canon_axes)
        return

    ad.settings.allow_write_nullable_strings = True
    OUT.mkdir(parents=True, exist_ok=True)
    (PROCESSED / "runs").mkdir(parents=True, exist_ok=True)
    code = {name: pilot.sha256(path) for name, path in {
        "m2_code": Path(__file__),
        "m1_import": Path(pilot.__file__),
        "m1_annotation": ROOT / "src/crlm/annotation_m1.py",
        "hgnc_code": ROOT / "src/crlm/hgnc_alias_map.py",
        "hgnc_table": HGNC_TSV_PATH,
        "hgnc_provenance": HGNC_PROVENANCE_PATH,
        "manifest": MANIFEST,
        "marker_table": pilot.MARKERS,
        "common_genes": COMMON,
        "normal_extra": NORMAL_EXTRA,
    }.items()}
    if not args.reports_only:
        m1_index = pd.read_csv(BASELINE).set_index("sample_uid")
        for gsm in args.samples or samples:
            build(gsm, *samples[gsm], m1_index, code, expected_axes, expected_canon_axes)
    if not args.samples:
        summarize(samples, common, code)


if __name__ == "__main__":
    main()
