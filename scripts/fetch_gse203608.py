"""Download the original UPenn matrix in verified resumable HTTP ranges."""

import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from validate_gse164522 import range_file


def main():
    raw = Path(__file__).resolve().parents[1] / "data/raw/GSE203608"
    name = "GSE203608_gene_by_cell_count_matrix.txt.gz"
    url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE203nnn/GSE203608/suppl/" + name
    size = 995467292
    raw.mkdir(exist_ok=True, parents=True)
    for annotation in (
        "GSE203608_cell_annotation.csv.gz",
        "GSE203608_gene_annotation.csv.gz",
    ):
        path = raw / annotation
        if not path.is_file():
            tmp = raw / (annotation + ".part")
            subprocess.run(
                [
                    "curl",
                    "-fsSL",
                    "--retry",
                    "4",
                    "-o",
                    str(tmp),
                    url.rsplit("/", 1)[0] + "/" + annotation,
                ],
                check=True,
            )
            tmp.replace(path)
    final = raw / name
    if final.exists() and final.stat().st_size == size:
        return
    prefix = raw / (name + ".part")
    start = prefix.stat().st_size if prefix.exists() else 0
    parts = raw / "ranges"
    parts.mkdir(exist_ok=True)
    jobs = []
    for first in range(start, size, 8 << 20):
        last = min(first + (8 << 20), size) - 1
        jobs.append((url, first, last, size, parts / f"{first}-{last}.bytes"))
    with ThreadPoolExecutor(max_workers=4) as pool:
        for i, _ in enumerate(pool.map(range_file, jobs), 1):
            if i % 10 == 0:
                print(f"{i} / {len(jobs)} ranges", flush=True)
    joined = raw / (name + ".assembling")
    with joined.open("wb") as output:
        for path in ([prefix] if start else []) + [j[-1] for j in jobs]:
            with path.open("rb") as stream:
                shutil.copyfileobj(stream, output, length=1 << 20)
    if joined.stat().st_size != size:
        raise ValueError("Original matrix size differs from the GEO response")
    joined.replace(final)
    if start:
        prefix.unlink()
    for job in jobs:
        job[-1].unlink()
    print(f"{final.name}: {size} bytes", flush=True)


if __name__ == "__main__":
    main()
