#!/usr/bin/env python3
"""Read-only host inventory unless --output is explicitly supplied. No installs."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from importlib import metadata
import json
import os
from pathlib import Path
import platform
import shutil
import sys

ROOT = Path(__file__).resolve().parents[1]
GIB = 1024 ** 3


def host_memory():
    """Best-effort physical memory; cgroup/VM limits can be lower."""
    meminfo = Path("/proc/meminfo")
    if meminfo.exists():
        values = {}
        for line in meminfo.read_text().splitlines():
            key, value = line.split(":", 1)
            if key in ("MemTotal", "MemAvailable"):
                values[key] = int(value.split()[0]) * 1024
        return {"source": "/proc/meminfo", "total_bytes": values.get("MemTotal"), "available_bytes": values.get("MemAvailable"), "note": "OS-visible physical memory, not a guarantee of process/cgroup allocation."}
    if platform.system() == "Windows":
        import ctypes
        from ctypes import wintypes

        class MemoryStatus(ctypes.Structure):
            _fields_ = [("dwLength", wintypes.DWORD), ("dwMemoryLoad", wintypes.DWORD)] + [(key, ctypes.c_ulonglong) for key in ("ullTotalPhys", "ullAvailPhys", "ullTotalPageFile", "ullAvailPageFile", "ullTotalVirtual", "ullAvailVirtual", "ullAvailExtendedVirtual")]
        status = MemoryStatus()
        status.dwLength = ctypes.sizeof(status)
        if ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(status)):
            return {"source": "GlobalMemoryStatusEx", "total_bytes": status.ullTotalPhys, "available_bytes": status.ullAvailPhys}
    if hasattr(os, "sysconf"):
        try:
            return {"source": "os.sysconf", "total_bytes": os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES"), "available_bytes": None}
        except (ValueError, OSError):
            pass
    return {"source": "unavailable", "total_bytes": None, "available_bytes": None}


def cpu_name():
    value = platform.processor()
    if not value and Path("/proc/cpuinfo").exists():
        for line in Path("/proc/cpuinfo").read_text().splitlines():
            if line.startswith("model name"):
                value = line.partition(":")[2].strip()
                break
    return value or None


def installed_versions():
    result = {}
    for name in ("numpy", "pandas", "scipy", "anndata", "scanpy", "torch", "scgpt"):
        try:
            result[name] = metadata.version(name)
        except metadata.PackageNotFoundError:
            result[name] = None
    return result


def probe_torch():
    try:
        import torch
        result = {"torch_version": torch.__version__, "cuda_available": bool(torch.cuda.is_available()), "xpu_available": bool(hasattr(torch, "xpu") and torch.xpu.is_available())}
        if result["cuda_available"]:
            result["cuda_devices"] = [torch.cuda.get_device_name(i) for i in range(torch.cuda.device_count())]
        if result["xpu_available"]:
            result["xpu_devices"] = [torch.xpu.get_device_name(i) for i in range(torch.xpu.device_count())]
        result["scope"] = "Backend availability only; not a scGPT compatibility, memory or performance test."
        return result
    except Exception as error:
        return {"probe_error": f"{type(error).__name__}: {error}", "scope": "No packages were installed or modified."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--torch-probe", action="store_true", help="Explicitly import installed PyTorch and query CUDA/XPU; does not run a model.")
    parser.add_argument("--output", help="Explicit output JSON path; relative paths are resolved against the project root. Writes only when supplied.")
    args = parser.parse_args()
    config = json.loads((ROOT / "configs/project.json").read_text(encoding="utf-8"))
    disk = shutil.disk_usage(ROOT)
    cells, genes, dim = 190190, 29569, 512
    report = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "This runtime only. It is the user's workstation only if run on that workstation.",
        "project_root": str(ROOT), "python": sys.version.split()[0], "python_executable": sys.executable,
        "platform": platform.platform(), "cpu_name": cpu_name(), "logical_cpu_count": os.cpu_count(),
        "memory": host_memory(), "project_filesystem_disk": {"total_bytes": disk.total, "free_bytes": disk.free, "free_gib": round(disk.free / GIB, 2)},
        "user_reported_target_hardware_unverified": config["user_reported_hardware"],
        "installed_packages": installed_versions(), "R_executable": shutil.which("R"),
        "torch_probe": probe_torch() if args.torch_probe else {"status": "not_requested", "use": "--torch-probe to query an already installed backend"},
        "arithmetic_examples_not_benchmarks": {
            "shared_gene_dense_float32_gib": round(cells * genes * 4 / GIB, 3),
            "shared_gene_dense_float64_gib": round(cells * genes * 8 / GIB, 3),
            "embedding_example_512d_float32_gib": round(cells * dim * 4 / GIB, 3),
            "basis": {"cells_before_QC": cells, "shared_genes": genes, "example_embedding_dim": dim},
            "excludes": "Copies, sparse index arrays, model weights, activations, metadata, temporary files and other cohorts. Actual scGPT embedding width follows the selected checkpoint."
        },
        "next_step": "Use sparse per-sample import on CPU; profile actual model inputs before choosing GPU/cloud resources. This report does not complete M0."
    }
    text = json.dumps(report, ensure_ascii=False, indent=2)
    if args.output:
        target = Path(args.output)
        if not target.is_absolute():
            target = ROOT / target
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(text + "\n", encoding="utf-8")
    print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
