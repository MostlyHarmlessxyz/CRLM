#!/usr/bin/env python3
"""Offline exact-symbol presence audit; no expression or biological effects are read."""
from pathlib import Path
import csv
import gzip
import hashlib
import itertools
import json

HERE = Path(__file__).resolve().parent
AUDIT = HERE.parent / "audit_data" / "geo"


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path, rows):
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    seed_path = HERE / "literature_program_seeds.csv"
    seeds = read_csv(seed_path)
    if len({(r["seed_set_id"], r["gene_symbol_original"]) for r in seeds}) != len(seeds):
        raise ValueError("Duplicate gene within a seed set; resolve against the original table")
    groups = {}
    for row in seeds:
        groups.setdefault(row["seed_set_id"], set()).add(row["gene_symbol_original"])
    files = {
        "GSE245552_PT_LM_representative_features": AUDIT / "GSM7844812_Pt01_T1_features.tsv.gz",
        "GSE245552_NC_NL_representative_features": AUDIT / "GSM7844838_Pt13_Z1_features.tsv.gz",
    }
    files.update({p.stem: p for p in sorted(AUDIT.glob("GSE335552_*_measured_genes.txt"))})
    if len(files) != 6:
        raise ValueError("Expected two representative feature files and all four spatial panels")
    details, summary, sources = [], [], []
    for name, path in files.items():
        if path.suffix == ".gz":
            with gzip.open(path, "rt", encoding="utf-8") as handle:
                features = [line.rstrip("\n").split("\t")[1] for line in handle if line.strip()]
        else:
            features = [line.strip() for line in path.read_text().splitlines() if line.strip()]
        measured = set(features)
        sources.append({"feature_axis": name, "file": str(path.relative_to(HERE.parent)),
                        "rows": len(features), "unique_original_symbols": len(measured), "sha256": sha256(path)})
        for set_id, genes in groups.items():
            present = genes & measured
            summary.append({"feature_axis": name, "seed_set_id": set_id, "n_seed_genes": len(genes),
                            "n_exact_symbol_present": len(present), "n_exact_symbol_absent": len(genes-present),
                            "fraction_exact_symbol_present": len(present)/len(genes),
                            "absent_original_symbols": ";".join(sorted(genes-present)),
                            "interpretation": "feature_axis_only; no alias mapping, expression, QC or score validation"})
            for gene in sorted(genes):
                details.append({"feature_axis": name, "seed_set_id": set_id,
                                "gene_symbol_original": gene, "exact_symbol_present": gene in measured})
    overlaps = []
    for first, second in itertools.combinations(groups, 2):
        shared, union = groups[first] & groups[second], groups[first] | groups[second]
        overlaps.append({"seed_set_a": first, "seed_set_b": second, "n_shared": len(shared),
                         "jaccard": len(shared)/len(union), "shared_original_symbols": ";".join(sorted(shared))})
    destination = HERE / "coverage"
    destination.mkdir(exist_ok=True)
    write_csv(destination / "seed_coverage_summary.csv", summary)
    write_csv(destination / "seed_coverage_by_gene.csv", details)
    write_csv(destination / "seed_overlap.csv", overlaps)
    audit = {"status": "completed", "scope": "Exact original-symbol matching to provided feature axes only",
             "seed_rows": len(seeds), "seed_unique_symbols": len({r["gene_symbol_original"] for r in seeds}),
             "n_sets": len(groups), "seed_sha256": sha256(seed_path), "source_feature_files": sources,
             "limitations": ["Absent exact symbols may reflect aliases as well as missing features; no alias resolution was performed.",
                             "Presence in a feature axis does not prove detected expression, adequate cells, reliable scoring or biological replication.",
                             "GSE245552 files are representatives of the separately audited identical tumor/adjacent feature sets.",
                             "Unmeasured genes must never be imputed as biological zero."]}
    (destination / "seed_coverage_audit.json").write_text(json.dumps(audit, ensure_ascii=False, indent=2))
    print(json.dumps({"sets": len(groups), "feature_axes": len(files), "summary_rows": len(summary),
                      "gene_presence_rows": len(details), "expression_or_effects_read": False}, ensure_ascii=False))


if __name__ == "__main__":
    main()
