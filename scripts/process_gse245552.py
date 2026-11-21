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
from crlm.annotation_m1 import PANELS, classify


ROOT = pilot.ROOT
MANIFEST = ROOT / "data/manifests/m0/sample_manifest_validated.csv"
BASELINE = ROOT / "artifacts/m1/matrix_audit.csv"
PROCESSED = ROOT / "data/processed/m2"
OUT = ROOT / "artifacts/m2"
COMMON = pilot.AUDIT / "GSE245552_feature_intersection.txt"
NORMAL_EXTRA = pilot.AUDIT / "GSE245552_normal_only_features_353.txt"
LEVELS = ("T/NK", "Mono/Mac/DC", "B/Plasma", "Neutrophil", "Mast", "Nonimmune", "Mixed", "Unknown")
MAIN = LEVELS[:3]
LIMITS = (20, 50, 100)
# These corroborating markers never participate in the identity call.
CHECK_GENES = {
    "T/NK": ("CD8A", "CD8B", "CD6", "IL7R", "CCL5"),
    "Mono/Mac/DC": ("LYZ", "MS4A6A", "MAFB", "VCAN", "WDFY4", "SPIB"),
    "B/Plasma": ("IGHA1", "IGHA2", "CD22", "CD72"),
    "Neutrophil": ("CXCL8", "G0S2", "PI3"),
    "Mast": ("TPSD1", "GATA2", "LTC4S"),
    "Nonimmune": ("TFF3", "CEACAM5", "PHGR1", "SPARCL1", "IGFBP7", "MGP"),
}


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


def qc_threshold(values):
    logs = np.log1p(values)
    center = np.median(logs)
    scale = 1.4826 * np.median(np.abs(logs - center))
    return np.expm1(center - 3 * scale), np.expm1(center + 3 * scale)


def prepare(row, record, files, m1_index):
    gsm = row["sample_uid"]
    if gsm in m1_index.index:
        audit = m1_index.loc[gsm]
        path = ROOT / audit.h5ad_path
        if pilot.sha256(path) != audit.h5ad_sha256:
            raise ValueError(f"M1 source changed: {path}")
        obj = ad.read_h5ad(path)
        reference = []  # The hash binds the original 16 source-coordinate checks.
        obj.uns["m1_object_sha256"] = audit.h5ad_sha256
    else:
        counts, genes, barcodes, reference = pilot.read_library(files, record)
        obs, _ = pilot.qc_cells(counts, genes.gene_symbol)
        labels, details = classify(counts, genes.gene_symbol, pilot.MARKERS)
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
    panels = obs.m1_marker_groups.astype(str)
    obs["subtype_m2"] = np.where(panels.str.contains(" | ", regex=False),
                                 obs.lineage_m2 + " unresolved", panels)
    unresolved = obs.lineage_m2.isin(["Unknown", "Mixed"])
    obs.loc[unresolved, "subtype_m2"] = obs.loc[unresolved, "lineage_m2"]
    low_counts, high_counts = qc_threshold(obs.n_counts_visible)
    low_genes, high_genes = qc_threshold(obs.n_genes_visible)
    low = (obs.n_counts_visible < low_counts) & (obs.n_genes_visible < low_genes)
    obs["qc_joint_high"] = (obs.n_counts_visible > high_counts) & (obs.n_genes_visible > high_genes)
    obs["qc_joint_low"] = low
    # Libraries contain lineages with different RNA content. Use the library
    # flag alone for sensitivity and require within-lineage support for exclusion.
    within = pd.Series(False, index=obs.index)
    for _, group in obs.groupby("lineage_m2", observed=True):
        if len(group) >= 50:
            group_counts, _ = qc_threshold(group.n_counts_visible)
            group_genes, _ = qc_threshold(group.n_genes_visible)
            within.loc[group.index] = (group.n_counts_visible < group_counts) & (group.n_genes_visible < group_genes)
    obs["qc_within_lineage_low"] = within
    empty = obs.n_counts_visible.eq(0)
    excluded_low = low & within
    obs["qc_retained"] = ~(excluded_low | empty)
    obs["qc_exclusion"] = np.select([empty, excluded_low], ["empty", "low_counts_and_genes"], default="")
    obj.uns["qc_limits"] = dict(low_counts=low_counts, low_genes=low_genes,
                                high_counts=high_counts, high_genes=high_genes)
    obs["published_barcode_prefix"] = obs.barcode.astype(str).str.rsplit("_", n=1).str[0]
    obj.uns["expression"] = {"location": "X", "scale": "raw_counts", "genes_prefiltered": True,
                             "gene_universe_id": row["gene_universe_id"], "source": record["source_record"],
                             "capture_protocol": "B07 cohort CD45+ and CD45- magnetic sorting, 2:1 remix",
                             "collection_window": "cohort synchronous; exact paired dates unavailable"}
    return obj, reference


def build(gsm, row, record, m1_index, code, expected_axes):
    tick = time.perf_counter()
    started = pilot.utc_now()
    files = source_files(record)
    obj, reference = prepare(row, record, files, m1_index)
    if set(obj.var.gene_symbol) != expected_axes[row["tissue"]]:
        raise ValueError(f"Unexpected measured gene set for {gsm}")
    obj.uns["provenance_m2"] = dict(code, **{k: record[k] for k in (
        "barcode_url", "feature_url", "matrix_url", "barcode_sha256", "feature_sha256", "matrix_sha256")})
    path = PROCESSED / f"{gsm}.h5ad"
    pilot.write_checked(obj, path, reference)
    metric = {"sample_uid": gsm, "source_coordinate_checks_including_m1": 16,
              "started_utc": started, "ended_utc": pilot.utc_now(),
              "elapsed_seconds": round(time.perf_counter() - tick, 3),
              "process_peak_rss_bytes": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024,
              "raw_bytes": sum(f.stat().st_size for f in files.values()),
              "h5ad_bytes": path.stat().st_size, "h5ad_sha256": pilot.sha256(path),
              "command": shlex.join([sys.executable, *sys.argv])}
    (PROCESSED / "runs" / f"{gsm}.json").write_text(json.dumps(metric, indent=2) + "\n")
    print(f"{gsm} {obj.n_obs} input, {obj.obs.qc_retained.sum()} retained, {metric['elapsed_seconds']} s",
          flush=True)


def group_definitions():
    groups = [("lineage", level, level) for level in LEVELS]
    for panel, (lineage, *_) in PANELS.items():
        groups.append(("subtype", lineage, panel))
    for lineage in LEVELS:
        if sum(value[0] == lineage for value in PANELS.values()) > 1:
            groups.append(("subtype", lineage, lineage + " unresolved"))
        elif lineage in ("Mixed", "Unknown"):
            groups.append(("subtype", lineage, lineage))
    return groups


def sample_report(gsm, row, obj, common, pb_dir):
    obs = obj.obs
    sample = {k: row[k] for k in ("donor_uid", "tissue", "lesion_uid", "sample_uid")}
    metrics = dict(sample, input_cells=len(obs), retained_cells=int(obs.qc_retained.sum()),
                   excluded_empty=int(obs.qc_exclusion.eq("empty").sum()),
                   excluded_joint_low=int(obs.qc_exclusion.eq("low_counts_and_genes").sum()),
                   flagged_library_joint_low=int(obs.qc_joint_low.sum()),
                   flagged_joint_high=int(obs.qc_joint_high.sum()),
                   mixed_retained=int((obs.qc_retained & obs.lineage_m2.eq("Mixed")).sum()),
                   unknown_retained=int((obs.qc_retained & obs.lineage_m2.eq("Unknown")).sum()),
                   native_genes=obj.n_vars, visible_counts=int(obs.n_counts_visible.sum()),
                   mt_available=int(obs.pct_mt.notna().sum()), **obj.uns["qc_limits"])
    for key in ("n_counts_visible", "n_genes_visible", "pct_mt"):
        for percentile in (.01, .25, .5, .75, .99):
            metrics[f"{key}_p{int(percentile*100):02d}"] = obs[key].quantile(percentile) if obs[key].notna().any() else np.nan
    broad = obs.lineage_m2.astype(str).to_numpy()
    sub = obs.subtype_m2.astype(str).to_numpy()
    keep = obs.qc_retained.to_numpy()
    nonempty = obs.n_counts_visible.gt(0).to_numpy()
    library_only = nonempty & ~obs.qc_joint_low.to_numpy()
    ids = obj.var_names.get_indexer(common)
    if (ids < 0).any():
        raise ValueError(f"{gsm} missing common-axis genes")
    groups = group_definitions()
    manifest, check_rows = [], []
    weights_i, weights_j = [], []
    for ix, (kind, lineage, name) in enumerate(groups):
        labeled = (broad if kind == "lineage" else sub) == name
        indexes = np.flatnonzero(labeled & keep)
        weights_i.extend([ix] * len(indexes))
        weights_j.extend(indexes)
        entry = dict(sample, level=kind, lineage=lineage, group=name, n_cells=len(indexes),
                     input_group_cells=int(labeled.sum()), n_retained_sample=int(keep.sum()),
                     n_nonempty=int((labeled & nonempty).sum()), n_sample_nonempty=int(nonempty.sum()),
                     n_library_qc=int((labeled & library_only).sum()), n_sample_library_qc=int(library_only.sum()),
                     count_matrix=f"data/processed/m2/pseudobulk/{gsm}.npz", matrix_row=ix)
        for limit in LIMITS:
            entry[f"eligible_{limit}"] = len(indexes) >= limit
        manifest.append(entry)
    weights = sparse.csr_matrix((np.ones(len(weights_j), dtype=np.int64), (weights_i, weights_j)),
                                shape=(len(groups), len(obs)))
    counts = (weights @ obj.X)[:, ids].tocsr()
    sums = np.asarray(counts.sum(axis=1)).ravel()
    if sums[:len(LEVELS)].sum() != sums[len(LEVELS):].sum() or sum(r["n_cells"] for r in manifest[:len(LEVELS)]) != keep.sum():
        raise ValueError(f"{gsm} nonexhaustive broad/subtype partition")
    sparse.save_npz(pb_dir / f"{gsm}.npz", counts, compressed=True)
    for entry, count in zip(manifest, sums):
        entry["common_axis_visible_count_sum"] = int(count)
    resolved = keep & ~np.isin(broad, ("Mixed", "Unknown"))
    for lineage, genes in CHECK_GENES.items():
        matches = obj.var_names.get_indexer(genes)
        if (matches < 0).any():
            raise ValueError(f"{gsm} missing prespecified corroboration marker")
        hit = np.asarray(obj.X[:, matches].getnnz(axis=1)) > 0
        target, others = resolved & (broad == lineage), resolved & (broad != lineage)
        check_rows.append(dict(sample, lineage=lineage, heldout_genes=" ".join(genes),
                               n_labeled=int(target.sum()), n_labeled_detected=int((hit & target).sum()),
                               n_other_resolved=int(others.sum()), n_other_detected=int((others & hit).sum())))
    return metrics, manifest, check_rows


def patient_tables(coverage, donors):
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
                data["status_50"] = ("tissue_unavailable" if rows is None else
                                     "no_labeled_cells" if rows.n_cells.sum() == 0 else
                                     "below_lesion_threshold" if not data["n_eligible_lesions_50"] else "eligible")
                tissue_rows.append(data)
    patient = pd.DataFrame(tissue_rows)
    paired = []
    for (kind, lineage, name), frame in patient.groupby(["level", "lineage", "group"], sort=False):
        for limit in LIMITS:
            supported = frame.pivot(index="donor_uid", columns="tissue", values=f"n_eligible_lesions_{limit}").gt(0)
            for contrast, tissues in (("PT_LM", ["PT", "LM"]), ("NC_NL", ["NC", "NL"]), ("Q", list(pilot.TISSUES))):
                valid = supported.index[supported[tissues].all(axis=1)].tolist()
                paired.append(dict(level=kind, lineage=lineage, group=name, min_cells_per_lesion=limit,
                                   contrast=contrast, n_donors=len(valid), donors="|".join(valid)))
    patient.to_csv(OUT / "donor_tissue_lineage_coverage.csv", index=False)
    pd.DataFrame(paired).to_csv(OUT / "paired_support.csv", index=False)
    comparison = []
    broad = coverage.loc[coverage.level.eq("lineage")]
    for policy, count in (("primary", "n_cells"), ("nonempty", "n_nonempty"), ("library_only", "n_library_qc")):
        for level, frame in broad.groupby("lineage", sort=False):
            top = frame.pivot_table(index="donor_uid", columns="tissue", values=count, aggfunc="max")
            for limit in LIMITS:
                for contrast, tissues in (("PT_LM", ["PT","LM"]), ("NC_NL", ["NC","NL"]), ("Q", list(pilot.TISSUES))):
                    valid = top.index[top[tissues].ge(limit).all(axis=1)].tolist()
                    comparison.append(dict(qc_policy=policy, lineage=level, min_cells_per_lesion=limit,
                                           contrast=contrast, n_donors=len(valid), donors="|".join(valid)))
    pd.DataFrame(comparison).to_csv(OUT / "qc_sensitivity.csv", index=False)


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
        metrics, groups, check_rows = sample_report(gsm, row, obj, common, pb_dir)
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
    profile = dict(ended_utc=pilot.utc_now(), command=shlex.join([sys.executable, *sys.argv]),
                   source_versions=code, python=sys.version, build_seconds_sum=sum(r["elapsed_seconds"] for r in runs),
                   summary_seconds=time.perf_counter()-start, summary_peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,
                   max_build_peak_rss_bytes=max(r["process_peak_rss_bytes"] for r in runs),
                   source_bytes=sum(r["raw_bytes"] for r in runs), h5ad_bytes=sum(r["h5ad_bytes"] for r in runs),
                   pseudobulk_bytes=sum(p.stat().st_size for p in pb_dir.iterdir()))
    (OUT / "resource_profile.json").write_text(json.dumps(profile, indent=2) + "\n")
    print(f"Accepted {len(qc)} GSMs, {qc.donor_uid.nunique()} patients, {qc.retained_cells.sum()} cells", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--samples", nargs="+", help="Build specified GSMs; run --reports-only after all are ready.")
    mode.add_argument("--reports-only", action="store_true", help="Rebuild cohort reports from verified M2 objects.")
    args = parser.parse_args()
    ad.settings.allow_write_nullable_strings = True
    OUT.mkdir(parents=True, exist_ok=True)
    (PROCESSED / "runs").mkdir(parents=True, exist_ok=True)
    samples = cohort()
    common = COMMON.read_text().splitlines()
    normal = set(common) | set(NORMAL_EXTRA.read_text().splitlines())
    expected_axes = {t: set(common) if t in ("PT", "LM") else normal for t in pilot.TISSUES}
    code = {name: pilot.sha256(path) for name, path in {
        "m2_code": Path(__file__), "m1_import": Path(pilot.__file__), "m1_annotation": ROOT/"src/crlm/annotation_m1.py",
        "manifest": MANIFEST, "marker_table": pilot.MARKERS, "common_genes": COMMON, "normal_extra": NORMAL_EXTRA}.items()}
    if not args.reports_only:
        m1_index = pd.read_csv(BASELINE).set_index("sample_uid")
        for gsm in args.samples or samples:
            build(gsm, *samples[gsm], m1_index, code, expected_axes)
    if not args.samples:
        summarize(samples, common, code)


if __name__ == "__main__":
    main()
