"""Score frozen definitions in the original in-vivo UPenn specimens (GSE203608)."""

import csv
import gzip
import hashlib
import json
from pathlib import Path
import resource
import sys
import time

import numpy as np
import pandas as pd
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src"), str(ROOT / "starter")]
from crlm.program_scores import rank_scores
from paired_effects import calculate_patient_contrasts, summarize_patients

DATASET = "GSE203608"
RAW = ROOT / "data/raw" / DATASET
OUT = ROOT / "artifacts/m6" / DATASET
FROZEN = ROOT / "artifacts/m5"
AUDIT = ROOT / "starter/audit_data/late_2026_data"
SITES = {"Colon": "NC", "Liver": "NL", "LiverMet": "LM", "Tumor": "PT"}
MYELOID = {
    "Macrophage",
    "cDC2 (CD1C+)",
    "cDC1 (BATF3+)",
    "LAMP3+ DC (LAMP3+)",
    "pDC (LILRA4+)",
}


def sha(path):
    with path.open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def read_matrix(path, nfeatures, ncells, gene_map, cell_map):
    blocks = []
    checked = 0
    with gzip.open(path, "rb") as f:
        if f.readline().strip() != b"%%MatrixMarket matrix coordinate integer general":
            raise ValueError("Expected author MatrixMarket integer count input")
        header = f.readline()
        while header.startswith(b"%"):
            header = f.readline()
        nf, nc, nnz = map(int, header.split())
        if (nf, nc) != (nfeatures, ncells):
            raise ValueError(
                f"Source annotation axes differ: {(nf, nc)} vs {(nfeatures, ncells)}"
            )
        for chunk in pd.read_csv(
            f,
            sep=r"\s+",
            header=None,
            names=["g", "c", "count"],
            dtype=np.int64,
            chunksize=1_000_000,
        ):
            g = chunk.g.to_numpy() - 1
            c = chunk.c.to_numpy() - 1
            val = chunk["count"].to_numpy()
            if (
                g.min() < 0
                or g.max() >= nf
                or c.min() < 0
                or c.max() >= nc
                or val.min() < 0
                or val.max() > np.iinfo(np.int32).max
            ):
                raise ValueError("Invalid original count or coordinate")
            checked += len(val)
            use = (gene_map[g] >= 0) & (cell_map[c] >= 0) & (val > 0)
            if use.any():
                blocks.append(
                    (cell_map[c[use]], gene_map[g[use]], val[use].astype(np.int32))
                )
        if checked != nnz:
            raise ValueError(f"Read {checked} entries, source declares {nnz}")
    cell, gene, values = (np.concatenate(v) for v in zip(*blocks))
    matrix = sparse.coo_matrix(
        (values, (cell, gene)), shape=(int(cell_map.max()) + 1, int(gene_map.max()) + 1)
    ).tocsr()
    if matrix.nnz != len(values):
        raise ValueError(
            "Duplicate original coordinates among selected cells and features"
        )
    return matrix, nnz


def run():
    start = time.perf_counter()
    registry = pd.read_csv(FROZEN / "frozen_candidate_registry.csv")
    if registry.empty or not registry.status.str.lower().eq("frozen").all():
        raise ValueError("A real frozen registry is required")
    sets = pd.read_csv(FROZEN / "gene_sets/program_genes.csv")
    OUT.mkdir(parents=True, exist_ok=True)
    coverage_dir = OUT / "coverage"
    coverage_dir.mkdir(exist_ok=True)
    cell_file = RAW / "GSE203608_cell_annotation.csv.gz"
    source_meta = json.loads(
        (AUDIT / "GSE203608_download_provenance.json").read_text()
    )[0]
    if sha(cell_file) != source_meta["sha256"]:
        raise ValueError(
            "Original author annotation version differs from identity audit"
        )
    raw_cells = pd.read_csv(cell_file, dtype=str, keep_default_na=False)
    cells = pd.read_csv(
        AUDIT / "GSE203608_original_identity_annotation_subset.csv.gz",
        dtype=str,
        keep_default_na=False,
    )
    pd.testing.assert_frame_equal(
        raw_cells.rename(columns={raw_cells.columns[0]: "barcode"})[cells.columns],
        cells,
    )
    if not cells.barcode.is_unique:
        raise ValueError("Duplicate author barcode")
    in_vivo = cells.SampleType.isin(SITES) & cells.Sample_pass_qc.eq("TRUE")
    cells["lineage"] = "Other/Unknown"
    cells.loc[cells.Cell_type.isin(["Plasma cell", "B cell"]), "lineage"] = "B/Plasma"
    cells.loc[
        cells.Cell_type.eq("Myeloid") & cells.Cell_type_myeloid.isin(MYELOID), "lineage"
    ] = "Mono/Mac/DC"
    cells.loc[cells.Cell_type.eq("T cell"), "lineage"] = "T/NK"
    cells["sample_uid"] = DATASET + ":" + cells.Dataset
    cells["donor_uid"] = DATASET + ":" + cells.Patient
    cells["tissue"] = cells.SampleType.map(SITES)
    samples = cells.loc[
        in_vivo, ["sample_uid", "donor_uid", "tissue"]
    ].drop_duplicates()
    if samples.sample_uid.duplicated().any():
        raise ValueError("A specimen maps to more than one patient or tissue")
    samples = samples.assign(
        lesion_uid=samples.sample_uid,
        dataset=DATASET,
        role="external",
        input_scale="raw_counts",
        patient_mapping_verified=True,
        counts_values_verified=True,
        counts_provenance="Published original gene_by_cell_count_matrix and ordered annotations; checked integer entries",
    )
    coverage = (
        cells.loc[in_vivo]
        .groupby(["sample_uid", "lineage"])
        .size()
        .unstack(fill_value=0)
    )
    coverage.to_csv(coverage_dir / "captured_composition_counts.csv")
    gene_file = RAW / "GSE203608_gene_annotation.csv.gz"
    genes = pd.read_csv(gene_file, dtype=str, keep_default_na=False)
    if genes.gene_id.duplicated().any():
        raise ValueError("Duplicate author gene IDs")
    symbol = genes.gene_short_name
    measurable = symbol.ne("") & ~symbol.duplicated(keep=False)
    gene_map = np.full(len(genes), -1, np.int32)
    gene_map[measurable] = np.arange(measurable.sum())
    gene_axis = pd.Index(symbol[measurable])
    chosen = (in_vivo & cells.lineage.isin(registry.lineage.unique())).to_numpy()
    cell_map = np.full(len(cells), -1, np.int32)
    cell_map[chosen] = np.arange(chosen.sum())
    matrix_file = RAW / "GSE203608_gene_by_cell_count_matrix.txt.gz"
    X, nnz = read_matrix(matrix_file, len(genes), len(cells), gene_map, cell_map)
    selected = cells[chosen].reset_index(drop=True)
    measured = []
    result = []
    for r in registry.itertuples(index=False):
        requested = sets.loc[
            sets.program.eq(r.program) & sets.lineage.eq(r.lineage), "gene_symbol"
        ].tolist()
        keep = [g for g in requested if g in gene_axis]
        ok = (
            len(keep) >= r.min_measured_positive_genes
            and len(keep) >= len(requested) * r.min_measured_fraction_positive
        )
        measured.append(
            dict(
                lineage=r.lineage,
                program=r.program,
                n_frozen=len(requested),
                n_measured=len(keep),
                missing_symbols=" ".join(sorted(set(requested) - set(keep))),
                eligible=ok,
            )
        )
        chosen_rows = np.flatnonzero(selected.lineage.eq(r.lineage))
        scored = selected.iloc[chosen_rows][
            ["sample_uid", "lineage", "Cell_type", "Cell_type_myeloid"]
        ].copy()
        scored["score"] = (
            rank_scores(X[chosen_rows], gene_axis, {r.program: keep})[
                r.program
            ].to_numpy()
            if ok
            else np.nan
        )
        by_sample = scored.groupby("sample_uid").score.agg(["mean", "size"])
        for sample in samples.itertuples(index=False):
            n = (
                int(by_sample.loc[sample.sample_uid, "size"])
                if sample.sample_uid in by_sample.index
                else 0
            )
            value = by_sample.loc[sample.sample_uid, "mean"] if n else np.nan
            result.append(
                dict(
                    sample_uid=sample.sample_uid,
                    lineage=r.lineage,
                    program=r.program,
                    n_cells=n,
                    score=value,
                )
            )
    samples.to_csv(coverage_dir / "sample_manifest.csv", index=False)
    pd.DataFrame(measured).to_csv(coverage_dir / "measurement_coverage.csv", index=False)
    all_scores = pd.DataFrame(result)
    all_scores.to_csv(OUT / "sample_program_measurements.csv", index=False)
    for cutoff in (20, 50, 100):
        scores = all_scores.copy()
        scores.loc[scores.n_cells < cutoff, "score"] = np.nan
        patients, lesions, tissues, _ = calculate_patient_contrasts(samples, scores)
        out = OUT / ("patient_effects" if cutoff == 50 else f"sensitivity_{cutoff}")
        out.mkdir(exist_ok=True)
        for name, fr in [
            ("patient_contrasts", patients),
            ("effects_summary", summarize_patients(patients)),
            ("patient_tissue_means", tissues),
            ("lesion_means", lesions),
        ]:
            fr.to_csv(out / (name + ".csv"), index=False)
        if cutoff == 50:
            summarize_patients(patients[patients.donor_uid.ne(DATASET + ":49")]).to_csv(
                out / "exclude_patient49.csv", index=False
            )
    profile = dict(
        source_metadata_file=cell_file.name,
        metadata_sha256=source_meta["sha256"],
        gene_annotation_sha256=sha(gene_file),
        expression_sha256=sha(matrix_file),
        source_coordinates=nnz,
        source_rows=len(genes),
        source_cells=len(cells),
        unambiguous_symbols=len(gene_axis),
        in_vivo_cells=int(in_vivo.sum()),
        in_vivo_specimens=len(samples),
        donors=int(samples.donor_uid.nunique()),
        elapsed_seconds=time.perf_counter() - start,
        peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024,
    )
    (OUT / "resource_profile.json").write_text(json.dumps(profile, indent=2) + "\n")
    print(json.dumps(profile), flush=True)


if __name__ == "__main__":
    run()
