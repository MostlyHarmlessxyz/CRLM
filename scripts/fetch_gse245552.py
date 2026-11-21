"""Materialize audited GSM files from the GEO series archive or the individual URLs."""
import argparse
import csv
import hashlib
import json
import os
import subprocess
import tarfile
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
AUDIT = ROOT / "starter/audit_data/geo"
RECORDS = AUDIT / "GSE245552_matrix_audit"
RAW = ROOT / "data/raw/GSE245552"
KINDS = ("barcode", "feature", "matrix")


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def from_archive(path, records):
    entries = {r[f"{k}_url"].rsplit("/", 1)[1]: (r["gsm"], r[f"{k}_sha256"])
               for r in records for k in KINDS}
    available = path.stat().st_size
    with tarfile.open(path, "r:") as archive:
        member = archive.next()
        while member is not None and member.offset_data + member.size <= available:
            name = Path(member.name).name
            if name in entries and member.isfile():
                gsm, expected = entries[name]
                dest = RAW / gsm / name
                dest.parent.mkdir(exist_ok=True)
                if not dest.exists():
                    with archive.extractfile(member) as source, tempfile.NamedTemporaryFile(dir=dest.parent) as tmp:
                        while block := source.read(1 << 20):
                            tmp.write(block)
                        tmp.flush()
                        if digest(Path(tmp.name)) != expected:
                            raise ValueError(f"Archive member differs from audited input: {name}")
                        os.link(tmp.name, dest)
            if archive.offset + 512 > available:
                break
            member = archive.next()


def fetch(record):
    folder = RAW / record["gsm"]
    folder.mkdir(exist_ok=True)
    result = []
    for kind in KINDS:
        url = record[f"{kind}_url"]
        path = folder / url.rsplit("/", 1)[1]
        expected = record[f"{kind}_sha256"]
        origin = "local"
        if not path.exists() and kind == "feature":
            for ref in AUDIT.glob("GSM*_features.tsv.gz"):
                if digest(ref) == expected:
                    os.link(ref, path)
                    origin = "identical archived feature file"
                    break
        if not path.exists():
            with tempfile.NamedTemporaryFile(dir=folder) as tmp:
                run = subprocess.run(["curl", "-fsSL", "--connect-timeout", "10", "--max-time", "120",
                                      "--retry", "2", "--retry-all-errors", "-o", tmp.name, url])
                if run.returncode:
                    print("pending", record["gsm"], kind, run.returncode, flush=True)
                    return result
                if digest(Path(tmp.name)) != expected:
                    raise ValueError(f"Downloaded file differs from audited input: {path}")
                os.link(tmp.name, path)
            origin = "GEO GSM"
        if digest(path) != expected:
            raise ValueError(f"Local file differs from audited input: {path}")
        result.append(dict(sample_uid=record["gsm"], kind=kind, path=str(path.relative_to(ROOT)),
                           url=url, bytes=path.stat().st_size, sha256=expected, acquisition=origin))
    print("ready",record["gsm"],flush=True)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, help="Use available complete members of a GEO RAW.tar first")
    args = parser.parse_args()
    RAW.mkdir(parents=True, exist_ok=True)
    records = [json.loads(p.read_text()) for p in sorted(RECORDS.glob("GSM*.json"))]
    if args.archive:
        from_archive(args.archive, records)
    with ThreadPoolExecutor(max_workers=2) as pool:
        rows = [row for batch in pool.map(fetch, records) for row in batch]
    output = ROOT / "artifacts/m2/source_fetch.csv"
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["sample_uid","kind","path","url","bytes","sha256","acquisition"])
        writer.writeheader()
        writer.writerows(rows)
    if len(rows) != 3 * len(records):
        raise SystemExit("Some GEO files remain pending. Rerun with the series archive or accessible GSM files.")
