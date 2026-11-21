"""Parse original GSE245552 S2/S3 DOCX tables and retain their exact scope.

S2 units are millions, not sequenced-cell/QC counts. S3 gene symbols come from
the final `gene` column; its first column contains R-generated row labels.
Default is offline using bundled original supplements and feature universe.
"""
from collections import Counter
from pathlib import Path
import csv
import hashlib
import json
import re
import xml.etree.ElementTree as ET
import zipfile

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "papers/GSE245552_supplement"
NS = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main"}


def read_docx(name):
    path = SOURCE / name
    with zipfile.ZipFile(path) as archive:
        xml = archive.read("word/document.xml")
    tree = ET.fromstring(xml)
    tables = []
    for table in tree.findall(".//w:tbl", NS):
        tables.append([[" ".join(x.text or "" for x in cell.findall(".//w:t", NS))
                        for cell in row.findall("w:tc", NS)] for row in table.findall("w:tr", NS)])
    return tree, tables, {"file": str(path.relative_to(ROOT)), "bytes": path.stat().st_size,
                          "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def write_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0])
        writer.writeheader()
        writer.writerows(rows)


def main():
    tree, tables, s2_source = read_docx("can-23-2123_table_s2_suppst2.docx")
    s2_source.update(doi="10.1158/0008-5472.25607325.v1",
                     api_url="https://api.figshare.com/v2/articles/25607325",
                     download_url="https://ndownloader.figshare.com/files/45666639")
    header_row = tree.find(".//w:tbl/w:tr", NS)
    superscript_sixes = []
    for run in header_row.findall(".//w:r", NS):
        va = run.find("w:rPr/w:vertAlign", NS)
        value = "".join(t.text or "" for t in run.findall("w:t", NS))
        if value == "6" and va is not None and va.get("{"+NS["w"]+"}val") == "superscript":
            superscript_sixes.append(value)
    assert len(superscript_sixes) == 2, "Original CD45+/- headers must show x10 superscript6"
    s2 = []
    for row in tables[0][1:]:
        assert len(row) == 5
        s2.append({
            "patient_id": re.sub(r"\s+", "", row[0]),
            "author_specimen_id": re.sub(r"\s+", "", row[1]), "author_sample_type": row[2],
            "CD45_positive_original_value": re.sub(r"\s+", "", row[3]),
            "CD45_negative_original_value": re.sub(r"\s+", "", row[4]),
            "original_unit": "x10^6 cells;not scRNA barcode/QC cells", "source_doi": s2_source["doi"],
        })
    assert len(s2) == 39 and len({r["patient_id"] for r in s2}) == 16
    write_csv(ROOT / "GSE245552_original_S2_cell_yields.csv", s2)

    _, tables, s3_source = read_docx("can-23-2123_table_s3_suppst3.docx")
    s3_source.update(doi="10.1158/0008-5472.25607322.v1",
                     api_url="https://api.figshare.com/v2/articles/25607322",
                     download_url="https://ndownloader.figshare.com/files/45666636")
    measured = set((ROOT / "GSE245552_feature_intersection.txt").read_text().splitlines())
    markers = []
    for index, table in enumerate(tables, 1):
        assert table[0][-1] == "gene"
        for row in table[1:]:
            assert len(row) == len(table[0])
            data = dict(zip(table[0], row))
            markers.append({
                "table_index": index, "author_row_label": row[0],
                "author_group": data.get("cluster", data.get("subcluster", "")), "gene": data["gene"],
                "p_val_author_cell_level": data["p_val"], "avg_log2FC_author": data["avg_log2FC"],
                "pct_1": data["pct.1"], "pct_2": data["pct.2"],
                "p_val_adj_author_cell_level": data["p_val_adj"], "source_doi": s3_source["doi"],
                "present_in_public_tumor_features": data["gene"] in measured,
            })
    write_csv(ROOT / "GSE245552_author_marker_tables.csv", markers)
    summary = {
        "S2": {"source": s2_source, "specimen_rows": len(s2),
               "patients": len({r["patient_id"] for r in s2}),
               "sample_types": dict(Counter(r["author_sample_type"] for r in s2)),
               "original_unit": "x10^6 cells; both superscript6 units confirmed in original DOCX XML",
               "limitation": "not scRNA barcode/QC counts; cell-counting stage and physical enrichment protocol not established from this table; exact ordinal correspondence of multiplePT/LM lesions to GEO_Tnumbers not inferred"},
        "S3": {"source": s3_source, "tables": len(tables), "marker_rows": len(markers),
               "unique_genes": len({r["gene"] for r in markers}),
               "groups": len({r["author_group"] for r in markers}),
               "unmeasured_genes": sorted({r["gene"] for r in markers if not r["present_in_public_tumor_features"]}),
               "gene_field": "last gene column; first column is not a safe gene-symbol field",
               "limitation": "author annotation markers, not new patient-level effects or verified interchangeable cross-cohort cell states"},
        "execution_scope": "offline original DOCX parsing and gene-availability check; no expression analysis",
    }
    assert (summary["S3"]["marker_rows"], summary["S3"]["unique_genes"], summary["S3"]["groups"]) == (600, 390, 30)
    assert not summary["S3"]["unmeasured_genes"]
    (ROOT / "GSE245552_S2_S3_provenance.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2))
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
