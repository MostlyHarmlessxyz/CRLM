"""Match every public GSE164522 expression column header to original metadata.

Reads only the first uncompressed line of each of six GEO expression files.
No gene rows, program scores or biological effects are analysed. At most 2
requests run together, and each response is closed immediately after its header.
Header SHA256 is explicitly NOT the checksum of the entire expression file.
"""
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
from urllib.request import urlopen
import csv
import gzip
import hashlib
import json

ROOT = Path(__file__).resolve().parent
GEO = ROOT / "geo"
URL_BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE164nnn/GSE164522/suppl/"
TISSUES = {"PT": "primary tumor", "MT": "metastasis tumor", "PN": "primary normal",
           "MN": "metastasis normal", "LN": "lymph node", "PBMC": "PBMC"}


class CountedReader:
    def __init__(self, source):
        self.source, self.consumed = source, 0

    def read(self, n=-1):
        if n < 0 or self.consumed + n > 8 * 1048576:
            raise ValueError("Bounded header-read limit exceeded")
        data = self.source.read(n)
        self.consumed += len(data)
        return data


def main():
    metadata_path = GEO / "GSE164522_CRLM_metadata.csv.gz"
    with gzip.open(metadata_path, "rt") as handle:
        metadata = list(csv.DictReader(handle))
    metadata_sha = hashlib.sha256(metadata_path.read_bytes()).hexdigest()
    expected = {code: {r[""].replace("-", ".") for r in metadata if r["tissue"] == tissue}
                for code, tissue in TISSUES.items()}
    patient_ids = {code: sorted({r["patient"] for r in metadata if r["tissue"] == tissue})
                   for code, tissue in TISSUES.items()}

    def audit(code):
        url = URL_BASE + f"GSE164522_CRLM_{code}_expression.csv.gz"
        result = {"tissue_code": code, "metadata_tissue": TISSUES[code], "url": url,
                  "metadata_file_sha256": metadata_sha,
                  "scope": "complete column header only; no gene row and no whole-expression-file checksum",
                  "patient_mapping_source": "original complete metadata, not inferred from expression values",
                  "patient_ids": patient_ids[code], "status": "failed"}
        try:
            with urlopen(url, timeout=60) as response:
                counted = CountedReader(response)
                with gzip.GzipFile(fileobj=counted) as stream:
                    header = stream.readline(32 * 1048576)
                result["compressed_bytes_consumed"] = counted.consumed
                result["http_content_length"] = response.headers.get("Content-Length")
            if not header.endswith(b"\n"):
                raise ValueError("Complete header newline not reached")
            columns = next(csv.reader([header.decode().rstrip("\r\n")]))[1:]
            observed = set(columns)
            wanted = expected[code]
            result.update(
                header_columns=len(columns), unique_header_columns=len(observed),
                expected_metadata_columns=len(wanted), intersection=len(wanted & observed),
                expression_only=len(observed-wanted), metadata_only=len(wanted-observed),
                expression_only_examples=sorted(observed-wanted)[:5],
                metadata_only_examples=sorted(wanted-observed)[:5],
                header_uncompressed_bytes=len(header),
                header_sha256=hashlib.sha256(header).hexdigest(),
                id_normalization="metadata barcode '-' replaced by '.' to match R-written matrix names",
                checked_utc=datetime.now(timezone.utc).isoformat(),
            )
            if len(columns) != len(observed) or observed != wanted:
                raise ValueError("Header mismatch or duplicate column")
            (GEO / f"GSE164522_{code}_expression_header.csv.gz").write_bytes(gzip.compress(header, mtime=0))
            result["status"] = "matched"
        except Exception as error:
            result.update(error=f"{type(error).__name__}: {error}", checked_utc=datetime.now(timezone.utc).isoformat())
        print(json.dumps(result), flush=True)
        return result

    with ThreadPoolExecutor(max_workers=2) as pool:
        reports = list(pool.map(audit, TISSUES))
    result = {
        "metadata_cells": len(metadata),
        "metadata_patient_ids": sorted({r["patient"] for r in metadata}),
        "metadata_sha256": metadata_sha,
        "all_headers_match": all(r["status"] == "matched" for r in reports),
        "total_matched_header_cells": sum(r.get("header_columns", 0) for r in reports if r["status"] == "matched"),
        "total_compressed_bytes_consumed": sum(r.get("compressed_bytes_consumed", 0) for r in reports),
        "files": reports,
        "limitation": "Column identity audit only; neither full matrix values nor expression processing provenance for each tissue are verified here",
    }
    (GEO / "GSE164522_all_expression_headers_audit.json").write_text(json.dumps(result, ensure_ascii=False, indent=2))
    print(json.dumps({k:v for k,v in result.items() if k!="files"}), flush=True)


if __name__ == "__main__":
    main()
