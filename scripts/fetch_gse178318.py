#!/usr/bin/env python3
"""Download GSE178318 supplementary single-cell UMI files in verified resumable HTTP ranges.

Prospective zero-exposure audit guarantee:
- Verifies that `artifacts/m5/frozen_inputs_sha256.txt` and `artifacts/m5/output_sha256.txt`
  already match locked Git commit `5b01467` BEFORE downloading any expression bytes.
- Downloads `GSE178318_barcodes.tsv.gz`, `GSE178318_genes.tsv.gz`, and `GSE178318_matrix.mtx.gz`
  into `data/raw/GSE178318/` with byte-exact Content-Length verification, gzip stream check,
  MatrixMarket dimension cross-check, and SHA256 recording in `artifacts/m6/GSE178318/source_inputs.csv`.
"""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import gzip
import hashlib
from pathlib import Path
import shutil
import subprocess
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "scripts") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts"))

from validate_gse164522 import range_file  # noqa: E402

DATASET = "GSE178318"
BASE_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE178nnn/GSE178318/suppl"
RAW_DIR = ROOT / "data/raw" / DATASET
OUT_DIR = ROOT / "artifacts/m6" / DATASET
M5_DIR = ROOT / "artifacts/m5"
LOCKED_M5_COMMIT = "5b01467"

EXPECTED_FILES = (
    {
        "file_name": "GSE178318_barcodes.tsv.gz",
        "bytes": 592450,
        "last_modified_http": "Wed, 16 Jun 2021 07:37:47 GMT",
        "kind": "barcodes",
    },
    {
        "file_name": "GSE178318_genes.tsv.gz",
        "bytes": 264505,
        "last_modified_http": "Wed, 16 Jun 2021 07:37:50 GMT",
        "kind": "genes",
    },
    {
        "file_name": "GSE178318_matrix.mtx.gz",
        "bytes": 545942023,
        "last_modified_http": "Wed, 16 Jun 2021 08:29:56 GMT",
        "kind": "matrix",
    },
)


def sha256_file(path: Path) -> str:
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _verify_sha_entries(manifest_path: Path) -> int:
    count = 0
    for line in manifest_path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line:
            continue
        expected, rel = line.split(None, 1)
        rel = rel.strip().lstrip("*")
        target = ROOT / rel
        if not target.is_file():
            raise FileNotFoundError(f"Missing frozen M5 input/output file: {rel}")
        actual = sha256_file(target)
        if actual != expected:
            raise ValueError(f"SHA256 mismatch for {rel}: got {actual}, expected {expected}")
        count += 1
    return count


def verify_pre_download_m5_lock() -> tuple[str, str]:
    """Assert M5 freeze manifests match `5b01467` AND all 72+19 underlying frozen files match SHA256."""
    frozen_in = M5_DIR / "frozen_inputs_sha256.txt"
    frozen_out = M5_DIR / "output_sha256.txt"
    if not frozen_in.is_file() or not frozen_out.is_file():
        raise FileNotFoundError("Missing M5 freeze SHA256 manifests; refusing external download")

    committed_in = subprocess.check_output(
        ["git", "show", f"{LOCKED_M5_COMMIT}:artifacts/m5/frozen_inputs_sha256.txt"],
        cwd=ROOT,
    )
    committed_out = subprocess.check_output(
        ["git", "show", f"{LOCKED_M5_COMMIT}:artifacts/m5/output_sha256.txt"],
        cwd=ROOT,
    )
    if frozen_in.read_bytes() != committed_in:
        raise ValueError("artifacts/m5/frozen_inputs_sha256.txt differs from locked commit 5b01467")
    if frozen_out.read_bytes() != committed_out:
        raise ValueError("artifacts/m5/output_sha256.txt differs from locked commit 5b01467")
    n_in = _verify_sha_entries(frozen_in)
    n_out = _verify_sha_entries(frozen_out)
    if n_in != 72 or n_out != 19:
        raise ValueError(f"Unexpected M5 manifest entry counts: frozen_inputs={n_in}, output={n_out}")
    return sha256_file(frozen_in), sha256_file(frozen_out)


def _file_mtime_utc(path: Path) -> str:
    return datetime.fromtimestamp(path.stat().st_mtime, tz=timezone.utc).isoformat(timespec="seconds")


def load_recorded_timestamps(out_csv: Path) -> dict[str, tuple[str, str]]:
    """Return `{file_name: (sha256, downloaded_utc)}` from existing `source_inputs.csv` if present."""
    if not out_csv.is_file():
        return {}
    prev = pd.read_csv(out_csv, dtype=str, keep_default_na=False)
    if not {"file_name", "sha256", "downloaded_utc"}.issubset(prev.columns):
        return {}
    return {
        str(r.file_name): (str(r.sha256), str(r.downloaded_utc))
        for r in prev.itertuples(index=False)
        if str(r.downloaded_utc).strip()
    }


def download_small_file(
    url: str,
    target: Path,
    expected_size: int,
    recorded: dict[str, tuple[str, str]],
) -> str:
    if target.is_file() and target.stat().st_size == expected_size:
        rec = recorded.get(target.name)
        if rec is not None:
            actual_sha = sha256_file(target)
            if actual_sha != rec[0]:
                raise ValueError(
                    f"Cached {target.name} SHA256 {actual_sha} differs from recorded {rec[0]}; refusing silent overwrite"
                )
            return rec[1]
        return _file_mtime_utc(target)
    tmp = target.with_name(target.name + ".part")
    subprocess.run(
        [
            "curl",
            "-fsSL",
            "--retry",
            "5",
            "--retry-all-errors",
            "--connect-timeout",
            "20",
            "-o",
            str(tmp),
            url,
        ],
        check=True,
    )
    if tmp.stat().st_size != expected_size:
        raise ValueError(
            f"Size mismatch for {target.name}: got {tmp.stat().st_size}, expected {expected_size}"
        )
    ts = utc_now()
    tmp.replace(target)
    return ts


def download_range_matrix(
    url: str,
    target: Path,
    expected_size: int,
    recorded: dict[str, tuple[str, str]],
) -> str:
    if target.is_file() and target.stat().st_size == expected_size:
        rec = recorded.get(target.name)
        if rec is not None:
            actual_sha = sha256_file(target)
            if actual_sha != rec[0]:
                raise ValueError(
                    f"Cached {target.name} SHA256 {actual_sha} differs from recorded {rec[0]}; refusing silent overwrite"
                )
            return rec[1]
        return _file_mtime_utc(target)
    prefix = target.with_name(target.name + ".part")
    start = prefix.stat().st_size if prefix.is_file() else 0
    parts = target.parent / "ranges"
    parts.mkdir(exist_ok=True, parents=True)
    chunk_bytes = 8 << 20  # 8 MiB
    jobs = []
    for first in range(start, expected_size, chunk_bytes):
        last = min(first + chunk_bytes, expected_size) - 1
        jobs.append((url, first, last, expected_size, parts / f"{first}-{last}.bytes"))

    with ThreadPoolExecutor(max_workers=6) as pool:
        for idx, _ in enumerate(pool.map(range_file, jobs), 1):
            if idx % 10 == 0 or idx == len(jobs):
                print(f"[fetch_gse178318] {idx} / {len(jobs)} ranges downloaded", flush=True)

    joined = target.with_name(target.name + ".assembling")
    with joined.open("wb") as output:
        for part_path in ([prefix] if start else []) + [job[-1] for job in jobs]:
            with part_path.open("rb") as stream:
                shutil.copyfileobj(stream, output, length=1 << 20)
    if joined.stat().st_size != expected_size:
        raise ValueError(
            f"Assembled matrix size {joined.stat().st_size} differs from expected {expected_size}"
        )
    ts = utc_now()
    joined.replace(target)
    if start and prefix.exists():
        prefix.unlink()
    for job in jobs:
        if job[-1].exists():
            job[-1].unlink()
    if parts.exists() and not any(parts.iterdir()):
        parts.rmdir()
    return ts


def verify_and_summarize(raw_dir: Path, timestamps: dict[str, str], m5_in_sha: str, m5_out_sha: str) -> pd.DataFrame:
    barcodes_path = raw_dir / "GSE178318_barcodes.tsv.gz"
    genes_path = raw_dir / "GSE178318_genes.tsv.gz"
    matrix_path = raw_dir / "GSE178318_matrix.mtx.gz"

    with gzip.open(barcodes_path, "rt", encoding="utf-8") as handle:
        barcodes = [line.rstrip("\r\n") for line in handle if line.strip()]
    with gzip.open(genes_path, "rt", encoding="utf-8") as handle:
        genes = [line.rstrip("\r\n") for line in handle if line.strip()]

    with gzip.open(matrix_path, "rb") as handle:
        banner = handle.readline().strip()
        if banner != b"%%MatrixMarket matrix coordinate integer general":
            raise ValueError(f"Unexpected MatrixMarket banner in {matrix_path.name}: {banner!r}")
        line = handle.readline()
        while line.startswith(b"%"):
            line = handle.readline()
        n_genes_hdr, n_cells_hdr, nnz_hdr = map(int, line.split())

    if n_genes_hdr != len(genes) or n_cells_hdr != len(barcodes):
        raise ValueError(
            f"MatrixMarket header ({n_genes_hdr}, {n_cells_hdr}) disagrees with "
            f"genes ({len(genes)}) and barcodes ({len(barcodes)})"
        )

    summary_map = {
        "GSE178318_barcodes.tsv.gz": f"n_cells={len(barcodes)}",
        "GSE178318_genes.tsv.gz": f"n_genes={len(genes)}",
        "GSE178318_matrix.mtx.gz": f"shape={n_genes_hdr}x{n_cells_hdr};nnz={nnz_hdr}",
    }

    rows = []
    for spec in EXPECTED_FILES:
        fname = spec["file_name"]
        fpath = raw_dir / fname
        actual_bytes = fpath.stat().st_size
        if actual_bytes != spec["bytes"]:
            raise ValueError(f"Unexpected byte size for {fname}: {actual_bytes} != {spec['bytes']}")
        rows.append(
            {
                "file_name": fname,
                "path": str(fpath.relative_to(ROOT)),
                "url": f"{BASE_URL}/{fname}",
                "bytes": actual_bytes,
                "last_modified_http": spec["last_modified_http"],
                "lines_or_nnz": summary_map[fname],
                "sha256": sha256_file(fpath),
                "downloaded_utc": timestamps[fname],
                "pre_download_m5_commit": LOCKED_M5_COMMIT,
                "pre_download_m5_frozen_inputs_sha256": m5_in_sha,
                "pre_download_m5_output_sha256": m5_out_sha,
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    m5_in_sha, m5_out_sha = verify_pre_download_m5_lock()
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out_csv = OUT_DIR / "source_inputs.csv"
    recorded = load_recorded_timestamps(out_csv)

    timestamps: dict[str, str] = {}
    for spec in EXPECTED_FILES:
        fname = spec["file_name"]
        url = f"{BASE_URL}/{fname}"
        target = RAW_DIR / fname
        if spec["kind"] == "matrix":
            timestamps[fname] = download_range_matrix(url, target, spec["bytes"], recorded)
        else:
            timestamps[fname] = download_small_file(url, target, spec["bytes"], recorded)
        print(f"[fetch_gse178318] ready: {fname} ({target.stat().st_size:,} bytes) at {timestamps[fname]}", flush=True)

    df = verify_and_summarize(RAW_DIR, timestamps, m5_in_sha, m5_out_sha)
    df.to_csv(out_csv, index=False)
    print(df[["file_name", "bytes", "lines_or_nnz", "sha256", "downloaded_utc"]].to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
