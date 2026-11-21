"""Recount GSE203608 clinical-tissue metadata, never expression or cell states.

Default: offline, using the two compact, exact metadata-column snapshots.
--make-compact-snapshots: read already-downloaded official CSV.gz files locally
and rebuild the snapshots. Neither mode downloads data or reads expression.
"""
from collections import Counter, defaultdict
from pathlib import Path
import argparse
import csv
import gzip
import hashlib
import io
import json

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "late_2026_data"
SOURCES = [
    "GSE203608_cell_annotation.csv.gz",
    "GSE203608_GSM7294726-GSM7294729_cell_annotation.csv.gz",
]
FIELDS = ["barcode", "Dataset", "Patient", "SampleType", "Cell_type",
          "Cell_type_myeloid", "Sample_pass_qc"]
INVIVO = {"Colon": "NC", "Tumor": "PT", "Liver": "NL", "LiverMet": "LM"}
MYELOID = {"Macrophage", "cDC2 (CD1C+)", "cDC1 (BATF3+)",
           "LAMP3+ DC (LAMP3+)", "pDC (LILRA4+)"}


def sha(body):
    return hashlib.sha256(body).hexdigest()


def snapshot_name(i):
    return "GSE203608_" + ("original" if i == 0 else "expanded") + "_identity_annotation_subset.csv.gz"


def make_snapshots():
    provenance = []
    for i, source in enumerate(SOURCES):
        body = (DATA / source).read_bytes()
        records = list(csv.DictReader(io.StringIO(gzip.decompress(body).decode())))
        buffer = io.StringIO(newline="")
        writer = csv.DictWriter(buffer, fieldnames=FIELDS)
        writer.writeheader()
        for row in records:
            writer.writerow({key: row.get("" if key == "barcode" else key, "") for key in FIELDS})
        selected = buffer.getvalue().encode()
        target = DATA / snapshot_name(i)
        target.write_bytes(gzip.compress(selected, mtime=0))
        provenance.append({
            "source": source,
            "source_url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE203nnn/GSE203608/suppl/"+source,
            "source_bytes": len(body), "source_sha256": sha(body),
            "snapshot": target.name, "snapshot_bytes": target.stat().st_size,
            "snapshot_sha256": sha(target.read_bytes()),
            "snapshot_uncompressed_sha256": sha(selected), "rows": len(records),
            "scope": "All original rows; selected metadata columns only; no filtering or patient renaming. Missing columns in expanded file remain blank. No expression or Macrophage_state column retained.",
        })
    (DATA / "GSE203608_identity_subset_provenance.json").write_text(json.dumps(provenance, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--make-compact-snapshots", action="store_true")
    args = parser.parse_args()
    if args.make_compact_snapshots:
        make_snapshots()
    provenance = json.loads((DATA / "GSE203608_identity_subset_provenance.json").read_text())
    tables = []
    for proof in provenance:
        body = (DATA / proof["snapshot"]).read_bytes()
        assert sha(body) == proof["snapshot_sha256"]
        decoded = gzip.decompress(body)
        assert sha(decoded) == proof["snapshot_uncompressed_sha256"]
        records = list(csv.DictReader(io.StringIO(decoded.decode())))
        assert len(records) == proof["rows"]
        assert len({r["barcode"] for r in records}) == len(records)
        tables.append(records)
    original, expanded = tables
    old = {r["barcode"]: r for r in original}
    new = {r["barcode"]: r for r in expanded}
    shared = set(old) & set(new)
    assert set(old) <= set(new)
    common_fields = ["Dataset", "Patient", "SampleType", "Cell_type"]
    assert all(all(old[b][k] == new[b][k] for k in common_fields) for b in shared)

    clinical = [r for r in original if r["SampleType"] in INVIVO]
    counts = Counter((r["Patient"], INVIVO[r["SampleType"]], r["Dataset"]) for r in clinical)
    myeloid = Counter((r["Patient"], INVIVO[r["SampleType"]], r["Dataset"])
                      for r in clinical if r["Cell_type_myeloid"] in MYELOID)
    support = [{"patient_code": p, "tissue": t, "author_dataset": s,
                "metadata_cells": n, "author_annotated_myeloid_cells": myeloid[(p,t,s)]}
               for (p,t,s), n in sorted(counts.items())]
    with (DATA / "GSE203608_patient_tissue_support.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=support[0]); writer.writeheader(); writer.writerows(support)
    sites = defaultdict(set)
    supports = defaultdict(dict)
    for row in support:
        sites[row["patient_code"]].add(row["tissue"])
        supports[row["patient_code"]][row["tissue"]] = row["author_annotated_myeloid_cells"]
    pairs = sorted(p for p,s in sites.items() if {"PT","LM"} <= s)
    quads = sorted(p for p,s in sites.items() if {"PT","LM","NC","NL"} <= s)
    result = {
        "scope": "Official deposited metadata recounted offline; no expression rows read, no reannotation, no biological effects or state-fraction comparisons.",
        "original_metadata_rows": len(original), "expanded_metadata_rows": len(expanded),
        "shared_barcodes": len(shared), "original_is_subset_of_expanded": set(old)<=set(new),
        "shared_identity_annotation_fields_identical": True,
        "new_expanded_rows": len(set(new)-set(old)),
        "new_expanded_SampleType_counts": dict(Counter(r["SampleType"] for r in expanded if r["barcode"] not in old)),
        "in_vivo_cells": len(clinical), "in_vivo_patients": len(sites),
        "in_vivo_specimens": len(counts),
        "in_vivo_cells_by_tissue": dict(Counter(INVIVO[r["SampleType"]] for r in clinical)),
        "in_vivo_specimens_by_tissue": dict(Counter(row["tissue"] for row in support)),
        "PT_LM_paired_patient_codes": pairs, "four_tissue_patient_codes": quads,
        "myeloid_definition": sorted(MYELOID),
        "myeloid_excluded_labels": ["Unannotated","Macrophage-T cell doublet","Mast cell"],
        "myeloid_cells_in_clinical_tissue": sum(myeloid.values()),
        "paired_myeloid_support": {p: {t:supports[p][t] for t in ["PT","LM"]} for p in pairs},
        "myeloid_pairs_at50_per_tissue": [p for p in pairs if min(supports[p]["PT"],supports[p]["LM"])>=50],
        "myeloid_pairs_at100_per_tissue": [p for p in pairs if min(supports[p]["PT"],supports[p]["LM"])>=100],
        "myeloid_quads_at50_per_tissue": [p for p in quads if min(supports[p].values())>=50],
        "full_expression_content_verified": False,
        "clinical_note": "Paper reports16patients,14MSS2MSI overall; metadata clinical-tissue subset has15. Do not replace actual subset denominator by paper total. Cultures/co-cultures and HTO libraries are not extra clinical subjects.",
        "status": "passed",
    }
    assert len(original)==126014 and len(expanded)==146988 and len(clinical)==59661
    assert len(sites)==15 and pairs==["09","27","44","49"] and quads==["09"]
    assert all(r["SampleType"] not in INVIVO for r in expanded if r["barcode"] not in old)
    (DATA / "GSE203608_metadata_audit_summary.json").write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
