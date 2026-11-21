"""Audit actual donor/large-intestine/liver support in Tabula Sapiens v2.

Default: read the bundled exact metadata subset only, without network or X.
Optional --source-metadata PATH re-extracts it from the official 42 MB gzip CSV.
No small-intestine rows, new cell annotation, program score or expression effect
enters this audit. Fine Mono/Mac/DC and broader author labels stay separate.
"""
from collections import Counter
from pathlib import Path
import argparse
import csv
import gzip
import hashlib
import io
import json

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "normal_tabula"
OUT.mkdir(exist_ok=True)
SOURCE_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE306nnn/GSE306755/suppl/GSE306755_TSP1_30_metadata_min200_2500_decontx_scvi_donorassay_version2d_20241112.csv.gz"
COLS = ["", "donor", "tissue", "anatomical_position", "method", "assay", "sample_id",
        "replicate", "cell_ontology_class", "compartment", "broad_cell_class", "free_annotation",
        "published_2022", "notes", "age", "sex"]
STRICT = {"macrophage", "colon macrophage", "monocyte", "classical monocyte",
          "intermediate monocyte", "non-classical monocyte", "myeloid dendritic cell",
          "plasmacytoid dendritic cell"}


def write_csv(path, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0])
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-metadata", type=Path,
                        help="optional complete official metadata.csv.gz; reads all rows through EOF")
    args = parser.parse_args()
    subset = OUT / "liver_large_intestine_metadata.csv.gz"
    if args.source_metadata:
        selected, donor_tissue = [], Counter()
        with gzip.open(args.source_metadata, "rt") as handle:
            for row in csv.DictReader(handle):
                donor_tissue[(row["donor"], row["tissue"])] += 1
                if row["tissue"] in {"Liver", "Large_Intestine"}:
                    selected.append({k: row[k] for k in COLS})
        text = io.StringIO(newline="")
        writer = csv.DictWriter(text, fieldnames=COLS)
        writer.writeheader()
        writer.writerows(selected)
        body = text.getvalue().encode()
        subset.write_bytes(gzip.compress(body, mtime=0))
        source_bytes = args.source_metadata.read_bytes()
        provenance = {
            "full_metadata_cells": sum(donor_tissue.values()),
            "full_metadata_donors": len({d for d, t in donor_tissue}),
            "full_metadata_tissues": len({t for d, t in donor_tissue}),
            "source_file": args.source_metadata.name, "source_url": SOURCE_URL,
            "source_bytes": len(source_bytes), "source_sha256": hashlib.sha256(source_bytes).hexdigest(),
            "selected_rows": len(selected), "selected_csv_sha256": hashlib.sha256(body).hexdigest(),
            "selected_fields": COLS,
            "scope": "complete original metadata read through gzip EOF; all Liver and Large_Intestine rows retained; no expression read",
        }
        (OUT / "metadata_provenance.json").write_text(json.dumps(provenance, indent=2))
        write_csv(OUT / "donor_tissue_counts.csv", [dict(donor=d, tissue=t, n_cells=n)
                  for (d, t), n in sorted(donor_tissue.items())])
    else:
        provenance = json.loads((OUT / "metadata_provenance.json").read_text())
        body = gzip.decompress(subset.read_bytes())
        assert hashlib.sha256(body).hexdigest() == provenance["selected_csv_sha256"]
        selected = list(csv.DictReader(body.decode().splitlines()))
        assert len(selected) == provenance["selected_rows"]

    assert len({r[""] for r in selected}) == len(selected)
    donor_tissue = Counter((r["donor"], r["tissue"]) for r in selected)
    donor_fine = Counter((r["donor"], r["tissue"], r["cell_ontology_class"]) for r in selected)
    results = []
    for donor, tissue in sorted(donor_tissue):
        values = [r for r in selected if (r["donor"], r["tissue"]) == (donor, tissue)]
        results.append({
            "donor": donor, "tissue": tissue, "n_all": len(values),
            "n_immune": sum(r["compartment"] == "Immune" for r in values),
            "n_strict_monomacdc": sum(r["cell_ontology_class"] in STRICT for r in values),
            "n_broad_myeloid_dc": sum(r["broad_cell_class"] in {"myeloid leukocyte", "dendritic cell"} for r in values),
            "anatomy": ";".join(sorted({r["anatomical_position"] for r in values})),
            "methods": ";".join(sorted({r["method"] for r in values})),
        })
    write_csv(OUT / "liver_colon_donor_support.csv", results)
    write_csv(OUT / "liver_colon_fine_cell_counts.csv", [dict(donor=d, tissue=t, fine_cell_type=c, n_cells=n)
              for (d, t, c), n in sorted(donor_fine.items())])
    stratified = Counter((r["donor"], r["tissue"], r["anatomical_position"], r["method"])
                         for r in selected if r["cell_ontology_class"] in STRICT)
    write_csv(OUT / "strict_myeloid_anatomy_method_counts.csv", [dict(donor=d, tissue=t, anatomy=a, method=m, n_cells=n)
              for (d, t, a, m), n in sorted(stratified.items())])
    liver = {d for d, t in donor_tissue if t == "Liver"}
    colon = {d for d, t in donor_tissue if t == "Large_Intestine"}
    paired = liver & colon
    support = {(r["donor"], r["tissue"]): r for r in results}
    strict_pairs50 = [d for d in paired if all(support[(d, t)]["n_strict_monomacdc"] >= 50 for t in {"Liver", "Large_Intestine"})]
    broad_pairs50 = [d for d in paired if all(support[(d, t)]["n_broad_myeloid_dc"] >= 50 for t in {"Liver", "Large_Intestine"})]
    summary = {
        "source_url": SOURCE_URL, "data_version": "GEO GSE306755 version2d 2024-11-12 snapshot, public 2026-08-26",
        "paper_doi": "10.1016/j.cell.2026.08.010",
        "full_metadata_cells": provenance["full_metadata_cells"],
        "full_metadata_donors": provenance["full_metadata_donors"],
        "full_metadata_tissues": provenance["full_metadata_tissues"],
        "n_selected_liver_large_intestine_cells": len(selected),
        "liver_donors": sorted(liver), "large_intestine_donors": sorted(colon),
        "paired_donors": sorted(paired), "strict_myeloid_definition": sorted(STRICT),
        "paired_strict_myeloid_at50_per_organ": sorted(strict_pairs50),
        "paired_broader_author_myeloid_at50_per_organ": sorted(broad_pairs50),
        "TSP14_liver_strict_myeloid": support[("TSP14", "Liver")]["n_strict_monomacdc"],
        "TSP14_large_intestine_strict_myeloid": support[("TSP14", "Large_Intestine")]["n_strict_monomacdc"],
        "TSP14_large_intestine_broader_author_myeloid": support[("TSP14", "Large_Intestine")]["n_broad_myeloid_dc"],
        "large_intestine_anatomy_labels": sorted({r["anatomical_position"] for r in selected if r["tissue"] == "Large_Intestine"}),
        "rectum_label_observed": any("rect" in r["anatomical_position"].lower() for r in selected if r["tissue"] == "Large_Intestine"),
        "small_intestine_included": False, "expression_read": False,
        "execution_scope": "complete original metadata extraction" if args.source_metadata else "offline exact52krowmetadata subset re-count; no fresh fullmetadata download",
        "interpretation": "No improvement to multiple-donor paired normal liver-colon myeloid inference; original labels only, no reannotation; two methods or colon regions are not independent people",
    }
    assert len(liver) == 4 and len(colon) == 5 and paired == {"TSP14"}
    assert summary["TSP14_liver_strict_myeloid"] == 338
    assert summary["TSP14_large_intestine_strict_myeloid"] == 41
    (OUT / "TabulaSapiens2_support_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
