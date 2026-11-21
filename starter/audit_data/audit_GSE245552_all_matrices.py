"""Audit every deposited GSE245552 count matrix, at two downloads maximum.

This is an input-capability audit, not biological reanalysis. Matrix files are
fully streamed through EOF and removed after successful verification by default.
Feature/barcode files, SHA256s, per-sample JSON and failure logs are retained.
No FASTQ, cell typing, program scoring or candidate effect estimation is done.
Network downloads/full scans now require --force; by default only a complete
local audit cache is reaggregated. The compact bundle's offline entry point is
reproduce_core_audit.py, which does not need every original feature/barcode file.
"""
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from itertools import islice
from pathlib import Path
import argparse
import csv
import gzip
import hashlib
import json
import traceback
import urllib.error
import urllib.request
import numpy as np

ROOT = Path(__file__).resolve().parent
GEO = ROOT / "geo"
OUT = ROOT / "geo/GSE245552_matrix_audit"
OUT.mkdir(parents=True, exist_ok=True)
Q = [0, .01, .05, .25, .5, .75, .95, .99, 1]
QL = ["min", "p01", "p05", "p25", "median", "p75", "p95", "p99", "max"]

def sha(path):
    with path.open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()

def fetch(url):
    destination = GEO / url.rsplit("/", 1)[-1]
    if destination.exists():
        return destination
    temporary = destination.with_suffix(destination.suffix + ".partial")
    errors = []
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=60) as response, temporary.open("wb") as out:
                expected = response.headers.get("Content-Length")
                copied = 0
                while chunk := response.read(1048576):
                    out.write(chunk)
                    copied += len(chunk)
            if expected is not None and copied != int(expected):
                raise ValueError(f"Truncated download {copied} != {expected}")
            temporary.replace(destination)
            return destination
        except (OSError, urllib.error.URLError, ValueError) as error:
            errors.append(f"attempt {attempt+1}: {type(error).__name__}: {error}")
    raise RuntimeError("; ".join(errors))

def quantiles(values):
    return dict(zip(QL, map(float, np.quantile(values, Q))))

def audit(row, keep=False, force=False):
    gsm = row["gsm"]
    report_path = OUT / (gsm + ".json")
    if report_path.exists() and not force:
        previous = json.loads(report_path.read_text())
        if previous.get("status") == "verified":
            return previous
    report = {"gsm": gsm, "specimen": row["sample_id"], "patient_id": row["patient_id"],
              "tissue": row["site"], "started_utc": datetime.now(timezone.utc).isoformat(),
              "status": "failed", "source_record": row["GEO_source_url"],
              "source_provenance": "Original GEO GSM: Cell Ranger v5.0.1 demultiplexing, barcoded processing, gene counting and aggregation; author-deposited MTX"}
    try:
        urls = row["files"].split(";")
        barcode_url = next(u for u in urls if u.endswith("_barcodes.tsv.gz"))
        feature_url = next(u for u in urls if u.endswith("_features.tsv.gz"))
        matrix_url = next(u for u in urls if u.endswith("_matrix.mtx.gz"))
        barcode_path, feature_path = fetch(barcode_url), fetch(feature_url)
        with gzip.open(barcode_path, "rt") as handle:
            barcodes = [line.rstrip("\r\n") for line in handle]
        with gzip.open(feature_path, "rt") as handle:
            features = [line.rstrip("\r\n").split("\t") for line in handle]
        feature_ids = [f[0] for f in features]
        feature_symbols = [f[1] if len(f) > 1 else f[0] for f in features]
        report.update(n_barcodes=len(barcodes), unique_barcodes=len(set(barcodes)),
                      n_features=len(features), unique_feature_ids=len(set(feature_ids)),
                      unique_feature_symbols=len(set(feature_symbols)),
                      feature_types=sorted({f[2] if len(f)>2 else "unspecified" for f in features}),
                      barcode_sha256=sha(barcode_path), feature_sha256=sha(feature_path),
                      matrix_url=matrix_url, feature_url=feature_url, barcode_url=barcode_url)
        if len(barcodes) != len(set(barcodes)) or len(features) != len(set(feature_ids)):
            raise ValueError("Duplicate barcode or feature ID")
        matrix_path = fetch(matrix_url)
        report.update(matrix_bytes=matrix_path.stat().st_size, matrix_sha256=sha(matrix_path))
        with gzip.open(matrix_path, "rt") as stream:
            header = stream.readline().strip()
            if header != "%%MatrixMarket matrix coordinate integer general":
                raise ValueError(f"Unexpected original matrix format: {header}")
            line = stream.readline()
            while line.startswith("%"):
                line = stream.readline()
            ng, nc, expected_nnz = map(int, line.split())
            report.update(matrix_header=header, n_genes=ng, n_cells=nc, declared_nnz=expected_nnz)
            if (ng, nc) != (len(features), len(barcodes)):
                raise ValueError("Matrix dimensions do not match features/barcodes")
            if nc != int(row["input_barcodes"]):
                raise ValueError("Barcode number differs from previous audited manifest")
            total_counts = np.zeros(nc, dtype=np.int64)
            genes_detected = np.zeros(nc, dtype=np.int64)
            entries = zero_entries = 0
            smallest = None
            largest = 0
            last_coordinate = -1
            for chunk in iter(lambda: list(islice(stream, 200000)), []):
                # loadtxt enforces three columns per physical line. Parse float64
                # then check integrality explicitly rather than silently truncate.
                values = np.loadtxt(chunk, dtype=np.float64, ndmin=2)
                if values.shape != (len(chunk), 3):
                    raise ValueError("Malformed MTX record")
                if not np.isfinite(values).all() or not np.equal(values, np.floor(values)).all():
                    raise ValueError("Nonfinite or noninteger matrix entries/indices")
                if (values > 2**53).any():
                    raise ValueError("Values exceed exact integer range for verification")
                values = values.astype(np.int64)
                ri, ci, v = values[:, 0], values[:, 1], values[:, 2]
                if (ri < 1).any() or (ri > ng).any() or (ci < 1).any() or (ci > nc).any() or (v < 0).any():
                    raise ValueError("Negative values or out-of-range matrix indices")
                coordinate = (ci-1) * ng + (ri-1)
                # Deposited matrices are in strict column/row order. This also
                # verifies no duplicate coordinates, so detected-gene counts are valid.
                if coordinate[0] <= last_coordinate or (np.diff(coordinate) <= 0).any():
                    raise ValueError("Duplicate or unsorted coordinate: requires separate review")
                last_coordinate = int(coordinate[-1])
                entries += len(values)
                zero_entries += int((v == 0).sum())
                smallest = int(v.min()) if smallest is None else min(smallest, int(v.min()))
                largest = max(largest, int(v.max()))
                total_counts += np.bincount(ci-1, weights=v, minlength=nc).astype(np.int64)
                genes_detected += np.bincount(ci[v>0]-1, minlength=nc)
        if entries != expected_nnz:
            raise ValueError(f"Entries {entries} != declared {expected_nnz}")
        report.update(status="verified", entries_checked=entries, zero_entries=zero_entries,
                      minimum_value=smallest, maximum_value=largest, full_stream_to_eof=True,
                      integer_nonnegative=True, dimensions_match=True, coordinates_unique=True,
                      column_row_sorted=True, total_molecules=int(total_counts.sum()),
                      zero_count_cells=int((total_counts==0).sum()),
                      library_size_quantiles=quantiles(total_counts),
                      genes_detected_quantiles=quantiles(genes_detected),
                      verification_scope="entire author-deposited matrix and barcode/feature files; not upstream FASTQ or biological QC",
                      completed_utc=datetime.now(timezone.utc).isoformat())
        report_path.write_text(json.dumps(report, indent=2))
        if not keep:
            matrix_path.unlink()
        print(f"VERIFIED {gsm} {row['sample_id']} {nc} cells {entries} entries", flush=True)
        return report
    except Exception as error:
        report.update(error=f"{type(error).__name__}: {error}", traceback=traceback.format_exc(),
                      completed_utc=datetime.now(timezone.utc).isoformat())
        report_path.write_text(json.dumps(report, indent=2))
        print(f"FAILED {gsm} {report['error']}", flush=True)
        return report

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, choices=[1,2], default=2)
    parser.add_argument("--keep-matrices", action="store_true")
    parser.add_argument("--force", action="store_true", help="explicitly permit GEO downloads and fully rescan every matrix")
    args = parser.parse_args()
    with (ROOT/"GSE245552_sample_manifest_audited.csv").open(encoding="utf-8-sig", newline="") as f:
        rows = list(csv.DictReader(f))
    if not args.force:
        complete_cache = True
        for row in rows:
            report_path = OUT / (row["gsm"]+".json")
            if not report_path.exists():
                complete_cache = False
                break
            record = json.loads(report_path.read_text())
            if record.get("status") != "verified" or any(not (GEO/record[k].rsplit("/",1)[-1]).exists()
                                                        for k in ["feature_url", "barcode_url"]):
                complete_cache = False
                break
        if not complete_cache:
            raise SystemExit("No complete local scan cache. Offline bundle review: python reproduce_core_audit.py. "
                             "A new full GEO download/scan requires explicit --force (about 1 GB). No files downloaded.")
    results = []
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(audit, row, args.keep_matrices, args.force) for row in rows]
        for future in as_completed(futures):
            results.append(future.result())
    results.sort(key=lambda r:r["gsm"])
    (ROOT/"GSE245552_all_matrix_audit.json").write_text(json.dumps(results, indent=2))
    flattened = []
    keys = ["gsm","specimen","patient_id","tissue","status","n_genes","n_cells","declared_nnz",
            "entries_checked","matrix_bytes","matrix_sha256","integer_nonnegative","dimensions_match",
            "coordinates_unique","total_molecules","zero_count_cells","unique_feature_symbols","error"]
    for result in results:
        item = {k: result.get(k, "") for k in keys}
        for metric in ["library_size", "genes_detected"]:
            for q in QL:
                item[metric+"_"+q] = result.get(metric+"_quantiles",{}).get(q,"")
        flattened.append(item)
    with (ROOT/"GSE245552_all_matrix_audit.csv").open("w",newline="",encoding="utf-8-sig") as f:
        writer=csv.DictWriter(f,fieldnames=flattened[0]);writer.writeheader();writer.writerows(flattened)
    verified = [r for r in results if r["status"] == "verified"]
    if verified:
        feature_hashes, feature_sets, mt_counts = {}, {}, {}
        raw_barcodes, repeated_barcodes = set(), 0
        for result in verified:
            feature_body = gzip.decompress((GEO/result["feature_url"].rsplit("/",1)[-1]).read_bytes())
            feature_hashes[result["gsm"]] = hashlib.sha256(feature_body).hexdigest()
            feature_sets[result["gsm"]] = {line.split("\t")[0] for line in feature_body.decode().splitlines()}
            mt_counts[result["gsm"]] = sum(g.startswith("MT-") for g in feature_sets[result["gsm"]])
            ids = set(gzip.decompress((GEO/result["barcode_url"].rsplit("/",1)[-1]).read_bytes()).decode().splitlines())
            repeated_barcodes += len(ids & raw_barcodes)
            raw_barcodes.update(ids)
        common = set.intersection(*feature_sets.values())
        union = set.union(*feature_sets.values())
        tumour_sets = [feature_sets[r["gsm"]] for r in verified if r["tissue"] in {"PT","LM"}]
        normal_sets = [feature_sets[r["gsm"]] for r in verified if r["tissue"] in {"NC","NL"}]
        tumour_union = set.union(*tumour_sets) if tumour_sets else set()
        normal_union = set.union(*normal_sets) if normal_sets else set()
        normal_only = normal_union - tumour_union
        (ROOT/"GSE245552_feature_intersection.txt").write_text("\n".join(sorted(common))+"\n")
        (ROOT/"GSE245552_normal_only_features.txt").write_text("\n".join(sorted(normal_only))+"\n")
        if len(normal_only) == 353:
            (ROOT/"GSE245552_normal_only_features_353.txt").write_text("\n".join(sorted(normal_only))+"\n")
        median_range = lambda key: [min((r[key]["median"],r["specimen"]) for r in verified), max((r[key]["median"],r["specimen"]) for r in verified)]
        summary = {
            "n_libraries": len(results), "n_donors": len({r["patient_id"] for r in verified}),
            "verified_libraries": len(verified), "failed_libraries": len(results)-len(verified),
            "total_input_barcodes": sum(r["n_cells"] for r in verified),
            "total_matrix_entries_checked": sum(r["entries_checked"] for r in verified),
            "total_compressed_matrix_bytes": sum(r["matrix_bytes"] for r in verified),
            "total_molecules": sum(r["total_molecules"] for r in verified),
            "zero_count_cells": sum(r["zero_count_cells"] for r in verified),
            "n_genes_per_file": sorted({r["n_genes"] for r in verified}),
            "unique_decompressed_feature_file_hashes": len(set(feature_hashes.values())),
            "feature_file_sha256": feature_hashes, "mt_gene_rows_per_sample": mt_counts,
            "feature_intersection_count": len(common), "feature_union_count": len(union),
            "tumor_only_feature_count": len(tumour_union-normal_union),
            "normal_only_feature_count": len(normal_only),
            "sample_median_library_size_range": median_range("library_size_quantiles"),
            "sample_median_detected_genes_range": median_range("genes_detected_quantiles"),
            "cross_sample_barcodes_seen_again": repeated_barcodes,
            "globally_unique_raw_barcode_strings": len(raw_barcodes),
            "matrix_files_removed_after_audit": not args.keep_matrices,
            "source": "Original GEO Cell Ranger v5.0.1 declarations plus complete successful MTX streams and matching feature/barcode files",
            "feature_universe_confound": "Tumour matrices lack genes present in adjacent-normal matrices, including MT- and many ribosomal genes; absent rows are not biological zeros",
            "limitation": "Input verification only; no cell filtering/annotation/program effect analysis; deposited counts are not complete unfiltered Cell Ranger outputs; tumour percent.mt cannot be recomputed",
        }
        (ROOT/"GSE245552_matrix_audit_summary.json").write_text(json.dumps(summary,indent=2))
    print(json.dumps({"execution_scope":"new complete matrix scan" if args.force else "retained complete-scan cache reaggregation; no fresh matrix scan",
                      "samples":len(results),"verified":sum(x["status"]=="verified" for x in results),
                      "failed":sum(x["status"]!="verified" for x in results),
                      "verified_input_cells":sum(x.get("n_cells",0) for x in results if x["status"]=="verified")}))

if __name__ == "__main__":
    main()
