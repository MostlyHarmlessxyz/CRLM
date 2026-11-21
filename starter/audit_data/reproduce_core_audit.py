"""Reaggregate the bundled CRLM input evidence without network access.

This script does NOT rescan the original 1 GB of expression matrices. It checks
the retained full-scan records, original small feature files, source/clinical
mapping, and saved metadata group tables. The output names each evidence scope.
For a new complete matrix scan use audit_GSE245552_all_matrices.py --force.
All paths are relative to this script, so the bundle can be moved elsewhere.
"""
from collections import Counter, defaultdict
from pathlib import Path
import csv
import gzip
import hashlib
import json
import zipfile
import xml.etree.ElementTree as ET
from audit_GSE245552_supplement_tables import read_docx, NS

ROOT = Path(__file__).resolve().parent
GEO = ROOT / "geo"


def rows(name):
    with (ROOT / name).open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def read_json(name):
    return json.loads((ROOT / name).read_text())


def read_identity_subset(name):
    provenance = read_json("geo/identity_subset_provenance.json")[name]
    body = gzip.decompress((GEO / name).read_bytes())
    assert hashlib.sha256(body).hexdigest() == provenance["derived_uncompressed_sha256"]
    values = list(csv.DictReader(body.decode().splitlines()))
    assert len(values) == provenance["n_rows"]
    assert len({r["barcode"] for r in values}) == len(values)
    return values


def groups_summary(group_counts, tumour_sites, four_sites):
    by_patient = defaultdict(set)
    by_tissue = Counter()
    for key, n in group_counts.items():
        patient, tissue = key.split("|", 1)
        by_patient[patient].add(tissue)
        by_tissue[tissue] += int(n)
    return {
        "n_patient_codes": len(by_patient),
        "n_patient_tissue_groups": len(group_counts),
        "n_cells": sum(by_tissue.values()),
        "n_paired_codes": sum(tumour_sites <= s for s in by_patient.values()),
        "n_four_tissue_codes": sum(four_sites <= s for s in by_patient.values()),
        "cells_by_tissue": dict(by_tissue),
    }


def main():
    result = {
        "execution_scope": "offline reaggregation of bundled audit evidence; no expression matrix downloaded or newly scanned",
        "historical_scan_scope": "39 original GSE245552 MTX streams were fully read to EOF in the recorded audit; their SHA256 and per-file results are retained",
    }

    manifest = rows("GSE245552_sample_manifest_audited.csv")
    starter = rows("GSE245552_sample_manifest.csv")
    audits = read_json("GSE245552_all_matrix_audit.json")
    expected = read_json("GSE245552_matrix_audit_summary.json")
    by_gsm = {r["gsm"]: r for r in audits}
    assert len(by_gsm) == len(audits) == len(manifest) == len(starter) == 39
    sites = defaultdict(set)
    for row in manifest:
        check = by_gsm[row["gsm"]]
        assert check["status"] == "verified" and check["full_stream_to_eof"]
        assert check["integer_nonnegative"] and check["dimensions_match"]
        assert check["coordinates_unique"]
        assert check["entries_checked"] == check["declared_nnz"]
        assert check["n_cells"] == int(row["input_barcodes"])
        assert (check["patient_id"], check["tissue"]) == (row["patient_id"], row["site"])
        assert len(check["matrix_sha256"]) == 64
        sites[row["patient_id"]].add(row["site"])
    for row in starter:
        assert row["counts_values_verified"].lower() == "true"
        assert row["role"] == "development"
        if row["tissue"] in {"PT", "LM"}:
            assert row["pct_mt_recomputable"] == "false"
            assert row["mitochondrial_gene_rows_available"] == "0"
    aggregate = {
        "n_libraries": len(audits), "n_donors": len(sites),
        "n_paired_PT_LM": sum({"PT", "LM"} <= x for x in sites.values()),
        "n_complete_four": sum({"PT", "LM", "NC", "NL"} <= x for x in sites.values()),
        "libraries_by_tissue": dict(Counter(r["site"] for r in manifest)),
        "total_input_barcodes": sum(r["n_cells"] for r in audits),
        "total_matrix_entries_checked": sum(r["entries_checked"] for r in audits),
        "total_compressed_matrix_bytes": sum(r["matrix_bytes"] for r in audits),
        "total_visible_molecules": sum(r["total_molecules"] for r in audits),
        "verification_scope": "reaggregated retained complete-scan records, not a new matrix-content scan",
    }
    for key in ["n_libraries", "n_donors", "total_input_barcodes", "total_matrix_entries_checked", "total_compressed_matrix_bytes"]:
        assert aggregate[key] == expected[key]
    assert (aggregate["n_donors"], aggregate["n_paired_PT_LM"], aggregate["n_complete_four"]) == (16, 13, 1)

    # These two representative files are original GEO feature files, not
    # gene lists reconstructed from a plot. All 39 file hashes are retained.
    feature_sets = {}
    for gsm in ["GSM7844812", "GSM7844838"]:
        check = by_gsm[gsm]
        path = GEO / check["feature_url"].rsplit("/", 1)[-1]
        compressed = path.read_bytes()
        assert hashlib.sha256(compressed).hexdigest() == check["feature_sha256"]
        body = gzip.decompress(compressed)
        assert hashlib.sha256(body).hexdigest() == expected["feature_file_sha256"][gsm]
        feature_sets[gsm] = {line.split("\t")[0] for line in body.decode().splitlines()}
    tumour = feature_sets["GSM7844812"]
    adjacent = feature_sets["GSM7844838"]
    assert len(tumour) == 29569 and len(adjacent) == 29922 and tumour < adjacent
    assert len(adjacent - tumour) == 353
    assert not any(g.startswith("MT-") for g in tumour)
    assert sum(g.startswith("MT-") for g in adjacent) == 13
    aggregate["gene_universe"] = {
        "tumour_genes": len(tumour), "adjacent_genes": len(adjacent),
        "genes_absent_from_tumour": len(adjacent-tumour),
        "tumour_mt_rows": 0, "adjacent_mt_rows": 13,
        "scope": "two original feature files read and SHA256 checked; per-file feature hashes for all39 retained",
    }

    ns = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main"}
    with zipfile.ZipFile(ROOT / "papers/can-23-2123_table_s1_suppst1.docx") as archive:
        tree = ET.fromstring(archive.read("word/document.xml"))
    table = tree.find(".//w:tbl", ns)
    cells = [[" ".join(x.text or "" for x in cell.findall(".//w:t", ns))
              for cell in row.findall("w:tc", ns)] for row in table.findall("w:tr", ns)]
    clinical = {r[0]: dict(zip(cells[0], r)) for r in cells[1:] if r and r[0].startswith("Pt")}
    assert len(clinical) == 18 and all(r["MSI status"] == "MSS" for r in clinical.values())
    assert set(sites) <= set(clinical)
    assert set(clinical)-set(sites) == {"Pt19", "Pt20"}
    aggregate["original_S1_patients"] = len(clinical)
    aggregate["all_S1_MSS"] = True
    aggregate["treatment"] = "unknown; no treatment column in original S1 or GEO sample metadata"
    supplement_proof = read_json("GSE245552_S2_S3_provenance.json")
    s2_tree, s2_tables, s2_source = read_docx("can-23-2123_table_s2_suppst2.docx")
    assert s2_source["sha256"] == supplement_proof["S2"]["source"]["sha256"]
    s2_people = {"".join(row[0].split()) for row in s2_tables[0][1:]}
    assert len(s2_tables[0])-1 == 39 and s2_people == set(sites)
    units6 = 0
    for run in s2_tree.find(".//w:tbl/w:tr", NS).findall(".//w:r", NS):
        align = run.find("w:rPr/w:vertAlign", NS)
        value = "".join(x.text or "" for x in run.findall("w:t", NS))
        units6 += int(value == "6" and align is not None and align.get("{"+NS["w"]+"}val") == "superscript")
    assert units6 == 2
    _, s3_tables, s3_source = read_docx("can-23-2123_table_s3_suppst3.docx")
    assert s3_source["sha256"] == supplement_proof["S3"]["source"]["sha256"]
    assert all(t[0][-1] == "gene" for t in s3_tables)
    s3_genes = {r[-1] for table in s3_tables for r in table[1:]}
    s3_rows = sum(len(table)-1 for table in s3_tables)
    assert s3_rows == 600 and len(s3_genes) == 390 and s3_genes <= tumour
    aggregate["original_S2_S3"] = {
        "S2_specimens": 39, "S2_donors": 16, "S2_units": "millions (x10^6); not scRNA QC cells",
        "S3_marker_rows": s3_rows, "S3_unique_genes": len(s3_genes),
        "S3_all_genes_measured_in_tumour": True, "S3_gene_field": "last gene column, not first row-label column",
        "scope": "original DOCX files SHA-checked and read offline; no expression or biological effects analysed",
    }
    result["GSE245552"] = aggregate

    # Exact three-column subsets preserve every original metadata row. Source
    # hashes are recorded, and the full metadata can be fetched separately.
    meta = read_json("data_audit_reproduced.json")
    ids164 = read_identity_subset("GSE164522_metadata_identity_columns.csv.gz")
    groups164 = Counter(r["patient_code"]+"|"+r["tissue"] for r in ids164)
    assert dict(groups164) == meta["GSE164522_metadata"]["counts_by_patient_tissue"]
    summary164 = groups_summary(groups164,
                                {"primary tumor", "metastasis tumor"},
                                {"primary tumor", "metastasis tumor", "primary normal", "metastasis normal"})
    summary164["four_tissue_cells"] = sum(v for k, v in summary164["cells_by_tissue"].items() if k not in {"PBMC", "lymph node"})
    summary164["scope"] = "complete original metadata identity columns recounted; all6originalexpressionheaders reread and matched; no expression gene rows read"
    assert (summary164["n_patient_codes"], summary164["four_tissue_cells"]) == (10, 163347)
    header_audit = read_json("geo/GSE164522_all_expression_headers_audit.json")
    header_counts = {}
    for report in header_audit["files"]:
        body = gzip.decompress((GEO / ("GSE164522_"+report["tissue_code"]+"_expression_header.csv.gz")).read_bytes())
        assert hashlib.sha256(body).hexdigest() == report["header_sha256"]
        observed = next(csv.reader([body.decode().rstrip("\r\n")]))[1:]
        wanted = {r["barcode"].replace("-", ".") for r in ids164 if r["tissue"] == report["metadata_tissue"]}
        assert len(observed) == len(set(observed)) and set(observed) == wanted
        header_counts[report["tissue_code"]] = len(observed)
    assert len(header_counts) == 6 and sum(header_counts.values()) == 257843
    summary164["all_expression_header_counts"] = header_counts
    result["GSE164522"] = summary164
    ids225 = read_identity_subset("GSE225857_immune_identity_columns.csv.gz")
    groups225 = Counter(r["patient_code"]+"|"+r["tissue"] for r in ids225)
    assert dict(groups225) == meta["GSE225857_immune_metadata"]["counts_by_patient_tissue"]
    summary225 = groups_summary(groups225,
                                {"CCL", "LCL"}, {"CCL", "LCL", "CNL", "LNL"})
    assert (summary225["n_patient_codes"], summary225["n_patient_tissue_groups"], summary225["n_cells"]) == (7, 30, 196473)
    assert (summary225["n_paired_codes"], summary225["n_four_tissue_codes"]) == (6, 5)
    summary225["scope"] = "complete original metadata identity columns recounted; patient codes are NOT verified independent individuals"
    summary225["mapping_status"] = "unresolved: clinical supplement has6codedpatients, expression includesextra0816"
    summary225["cells_after_unjustified_removal_of_s0816"] = sum(1 for r in ids225 if r["patient_code"] != "s0816")
    assert summary225["cells_after_unjustified_removal_of_s0816"] == 188057
    ids225_nonimmune = read_identity_subset("GSE225857_nonimmune_identity_columns.csv.gz")
    summary225["nonimmune_cells"] = len(ids225_nonimmune)
    assert summary225["nonimmune_cells"] == 41892
    result["GSE225857"] = summary225

    small = rows("small_file_cell_counts.csv")
    spatial = [r for r in small if r["dataset"] == "GSE335552"]
    panels = []
    for row in spatial:
        genes = set((GEO / ("GSE335552_"+row["title"]+"_measured_genes.txt")).read_text().splitlines())
        assert len(genes) == int(row["genes"])
        panels.append(genes)
    result["GSE335552"] = {
        "samples": len(spatial), "input_segmented_cells": sum(int(r["cells"]) for r in spatial),
        "genes_per_panel": [len(s) for s in panels],
        "gene_intersection": len(set.intersection(*panels)), "gene_union": len(set.union(*panels)),
        "scope": "retained exact gene lists extracted from4originalH5 and saved per-H5 barcode counts; original H5 not bundled",
    }
    assert result["GSE335552"]["input_segmented_cells"] == 646711
    assert result["GSE335552"]["gene_intersection"] == 309
    assert all(not ({"APOE", "APOC1", "MCAM", "TNFSF12", "TNFRSF12A"} & s) for s in panels)

    organ_rows = rows("E-MTAB-11536_donor_organ_cell_counts.csv")
    lineage_rows = rows("E-MTAB-11536_liver_colon_lineage_support.csv")
    obs_body = gzip.decompress((GEO / "E-MTAB-11536_obs_selected.json.gz").read_bytes())
    assert hashlib.sha256(obs_body).hexdigest() == read_json("geo/E-MTAB-11536_access_probe.json")["obs_json_sha256"]
    obs = json.loads(obs_body)
    recalculated_organs = Counter(zip(obs["Donor"], obs["Organ"]))
    assert recalculated_organs == {(r["donor"], r["organ"]): int(r["n_immune_cells"]) for r in organ_rows}
    per_site = {(r["donor"], r["organ"]): int(r["n_cells"]) for r in lineage_rows if r["lineage"] == "Myeloid"}
    recalculated_myeloid = Counter((d, o) for d, o, c in zip(obs["Donor"], obs["Organ"], obs["Majority_voting_CellTypist_high"])
                                  if o in {"LIV", "CAE", "TCL", "SCL"} and c in {"Monocytes", "Macrophages", "DC", "pDC"})
    assert recalculated_myeloid == per_site
    liver_donors = {r["donor"] for r in organ_rows if r["organ"] == "LIV"}
    large_bowel = {"CAE", "TCL", "SCL"}
    paired50 = []
    for donor in liver_donors:
        bowel_total = sum(per_site.get((donor, organ), 0) for organ in large_bowel)
        if bowel_total >= 50 and per_site.get((donor, "LIV"), 0) >= 50:
            paired50.append(donor)
    result["E-MTAB-11536"] = {
        "immune_cells": sum(int(r["n_immune_cells"]) for r in organ_rows),
        "total_donors": len({r["donor"] for r in organ_rows}),
        "liver_donors": sorted(liver_donors),
        "transverse_sigmoid_donors": sorted({r["donor"] for r in organ_rows if r["organ"] in {"TCL", "SCL"}}),
        "myeloid_pairs_at50_even_if_all_large_bowel_sites_pooled": paired50,
        "scope": "complete329762authorobs rows recounted from exact retained columns; expression not read, no new cell annotation",
    }
    assert result["E-MTAB-11536"]["immune_cells"] == 329762
    assert len(liver_donors) == 7 and not paired50

    tabula_provenance = read_json("normal_tabula/metadata_provenance.json")
    tabula_saved = read_json("normal_tabula/TabulaSapiens2_support_summary.json")
    tabula_body = gzip.decompress((ROOT / "normal_tabula/liver_large_intestine_metadata.csv.gz").read_bytes())
    assert hashlib.sha256(tabula_body).hexdigest() == tabula_provenance["selected_csv_sha256"]
    tabula_rows = list(csv.DictReader(tabula_body.decode().splitlines()))
    tabula_liver = {r["donor"] for r in tabula_rows if r["tissue"] == "Liver"}
    tabula_colon = {r["donor"] for r in tabula_rows if r["tissue"] == "Large_Intestine"}
    strict_types = set(tabula_saved["strict_myeloid_definition"])
    tabula_strict = Counter((r["donor"], r["tissue"]) for r in tabula_rows if r["cell_ontology_class"] in strict_types)
    assert len(tabula_rows) == 52298 and len(tabula_liver) == 4 and len(tabula_colon) == 5
    assert tabula_liver & tabula_colon == {"TSP14"}
    assert tabula_strict[("TSP14", "Liver")] == 338 and tabula_strict[("TSP14", "Large_Intestine")] == 41
    result["TabulaSapiens2_GSE306755"] = {
        "actual_recounted_liver_large_intestine_metadata_rows": len(tabula_rows),
        "liver_donors": sorted(tabula_liver), "large_intestine_donors": sorted(tabula_colon),
        "paired_donors": sorted(tabula_liver & tabula_colon),
        "TSP14_strict_myeloid_liver": 338, "TSP14_strict_myeloid_large_intestine": 41,
        "strict_paired_donors_at50_per_organ": [],
        "scope": "exact selected metadata recounted offline; no full42MBmetadata or expression downloaded in this command",
    }

    result["status"] = "passed"
    (ROOT / "core_audit_reproduced.json").write_text(json.dumps(result, ensure_ascii=False, indent=2))
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
