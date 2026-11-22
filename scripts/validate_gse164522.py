"""Prepare the four GSE164522 tissues and apply frozen V1.1 discovery programs.

Supports:
1. Official HGNC symbol canonicalization onto `245_164_shared_nuclear` (18,815 nuclear genes),
   preserving symbols already on the shared nuclear axis without context-free re-mapping.
2. Frozen detection-matched background controls (`matched_controls.csv` inside `freeze_dir`)
   with strict `requires_debiased` enforcement (zero silent fallback or external refitting):
   - Primary metric: `debiased_ucell_scores` (`U'(P) - U'(P_ctrl)`).
   - Parallel sensitivity metrics: `log1p_cp10k_debiased` (re-calibrated to CP10k on the shared
     nuclear axis), `ucell_raw`, and `legacy_positive` (`positive_midrank_1500`).
3. Both broad lineages (`T/NK`, `Mono/Mac/DC`, `B/Plasma`) and fine subtypes (`CD4 T cell`,
   `CD8 T cell`, `NK cell`, `Macrophage`, `Monocyte`, `monocyte/macrophage`, `B cell`, `plasma cell`)
   evaluated at both `>=50` (primary) and `>=20` cell thresholds.
4. Cross-lineage `CD3E+` doublet removal in `Mono/Mac/DC` (`hC36/41/47/56`, 3,662 cells) and
   `B/Plasma` (`hC75_B-CD3E`, 3,174 cells) in the primary clean-author branch, PLUS a unified
   marker-rule re-annotation branch (`annotation_m1.classify`) and dual immune denominators
   (`total_resolved_immune` including `ILC3`/`Mast`, and `core_three_lineage_immune`).
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import resource
import shlex
import shutil
import subprocess
import sys
import time
from typing import Iterable, Mapping, Sequence

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))
if str(ROOT / "starter") not in sys.path:
    sys.path.insert(1, str(ROOT / "starter"))

from adjust_registered_tests import adjust_registered  # noqa: E402
from composition_decomposition import decompose_pair  # noqa: E402
from crlm.annotation_m1 import classify  # noqa: E402
from crlm.hgnc_alias_map import (  # noqa: E402
    align_axis_to_canonical_targets,
    get_hgnc_mapper,
    load_shared_nuclear_universe,
    normalize_gene_axis,
)
from crlm.model_bags import cell_seed  # noqa: E402
from crlm.program_scores import (  # noqa: E402
    controls_from_dataframe,
    debiased_ucell_scores,
    legacy_positive_rank_scores,
    log1p_cp10k_scores,
    ucell_scores,
)
from paired_effects import calculate_patient_contrasts, summarize_patients  # noqa: E402

GSE = "GSE164522"
RAW = ROOT / "data/raw" / GSE
OUT = ROOT / "artifacts/m6" / GSE
OBJECTS = ROOT / "data/processed/m6" / GSE
GEO = ROOT / "starter/audit_data/geo"
METADATA = RAW / "GSE164522_CRLM_metadata.csv.gz"
FREEZE = ROOT / "artifacts/m5"

ORGANS = {"PT": "PT", "MT": "LM", "PN": "NC", "MN": "NL"}
MAJOR = {
    "Bcell": "B/Plasma",
    "Monocyte": "Mono/Mac/DC",
    "Macrophage": "Mono/Mac/DC",
    "DC": "Mono/Mac/DC",
    "CD4 T": "T/NK",
    "CD8 T": "T/NK",
    "NK cell": "T/NK",
    "gdT": "T/NK",
    "NKT": "T/NK",
    "MAST": "Mast",
    "ILC3": "ILC3",
    "CD45-": "CD45-",
}
EXPECTED_GLOBAL = {
    "T/NK": "Lymphoid-T-NK",
    "B/Plasma": "Lymphoid-B",
    "Mono/Mac/DC": "Myeloid",
}
PARENT_LINEAGE_MAP = {
    "T/NK": "T/NK",
    "CD4 T cell": "T/NK",
    "CD8 T cell": "T/NK",
    "Double_neg_or_unsplit_T": "T/NK",
    "NK cell": "T/NK",
    "Mono/Mac/DC": "Mono/Mac/DC",
    "monocyte/macrophage": "Mono/Mac/DC",
    "Macrophage": "Mono/Mac/DC",
    "Monocyte": "Mono/Mac/DC",
    "conventional DC": "Mono/Mac/DC",
    "plasmacytoid DC": "Mono/Mac/DC",
    "mature DC": "Mono/Mac/DC",
    "B/Plasma": "B/Plasma",
    "B cell": "B/Plasma",
    "plasma cell": "B/Plasma",
    "Mast": "Mast",
    "mast cell": "Mast",
    "Neutrophil": "Neutrophil",
    "neutrophil": "Neutrophil",
}
AUTHOR_RESOLVED_IMMUNE = ("T/NK", "Mono/Mac/DC", "B/Plasma", "Mast", "ILC3")
RULE_RESOLVED_IMMUNE = ("T/NK", "Mono/Mac/DC", "B/Plasma", "Neutrophil", "Mast")
CORE_THREE_IMMUNE = ("T/NK", "Mono/Mac/DC", "B/Plasma")
FOCAL_SUBTYPES_BY_PARENT = {
    "T/NK": ("CD4 T cell", "CD8 T cell", "NK cell"),
    "Mono/Mac/DC": ("Macrophage", "Monocyte", "monocyte/macrophage"),
    "B/Plasma": ("plasma cell", "B cell"),
}

# Exact M3 LIVER_RBC ambient hepatocyte/erythrocyte set (intentionally excludes myeloid core APOE/APOC1).
LIVER_RBC = {
    "ALB", "APOA1", "APOA2", "APOC3", "APOH", "TTR",
    "FGA", "FGB", "FGG", "HP", "HPR", "HBA1", "HBA2", "HBB",
}
SCORE_COLUMNS = ["sample_uid", "lineage", "program", "score", "n_cells"]
SEED = 20261121


def sha(path: Path) -> str:
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def canonicalize_preserving_axis(symbols: Iterable[object], axis_set: set[str]) -> list[str]:
    """Canonicalize gene symbols while preserving any symbol already on `axis_set` verbatim."""
    mapper = get_hgnc_mapper()
    out: list[str] = []
    seen: set[str] = set()
    for raw in symbols:
        if raw is None:
            continue
        s = str(raw).strip()
        if not s:
            continue
        canon = s if s in axis_set else mapper.resolve_symbol(s).canonical_symbol
        if canon and canon not in seen:
            seen.add(canon)
            out.append(canon)
    return out


def registry_requires_debiased(registry: pd.DataFrame) -> bool:
    if "scoring_method" not in registry.columns or registry.empty:
        return True
    methods = set(registry["scoring_method"].astype(str).str.strip())
    return any(m != "positive_midrank_1500" for m in methods)


def map_author_annotations(frame: pd.DataFrame) -> pd.DataFrame:
    """Attach broad author lineage, cross-lineage CD3E+ doublet flag, and harmonized fine subtype."""
    out = frame.copy()
    out["lineage"] = out["celltype_major"].map(MAJOR)
    if out["lineage"].isna().any():
        bad = sorted(out.loc[out["lineage"].isna(), "celltype_major"].unique())
        raise ValueError(f"Unsupported GSE164522 celltype_major labels: {bad}")

    sub_str = out["celltype_sub"].astype(str)
    major_str = out["celltype_major"].astype(str)
    is_cd3e_doublet = out["lineage"].isin(["Mono/Mac/DC", "B/Plasma"]) & sub_str.str.contains("CD3E", regex=False)
    out["is_cd3e_doublet"] = is_cd3e_doublet
    out["clean_author_qc"] = ~is_cd3e_doublet

    fine = np.full(len(out), "Unknown", dtype=object)
    fine[major_str.eq("CD4 T")] = "CD4 T cell"
    fine[major_str.eq("CD8 T")] = "CD8 T cell"
    fine[major_str.eq("NK cell")] = "NK cell"
    fine[major_str.isin(["gdT", "NKT"])] = "Double_neg_or_unsplit_T"
    fine[major_str.eq("Macrophage")] = "Macrophage"
    fine[major_str.eq("Monocyte")] = "Monocyte"
    is_dc = major_str.eq("DC")
    fine[is_dc] = "conventional DC"
    fine[is_dc & sub_str.str.contains("LAMP3", regex=False)] = "mature DC"
    fine[is_dc & (sub_str.str.contains("pDC", regex=False) | sub_str.str.contains("LILRA4", regex=False))] = "plasmacytoid DC"
    is_b = major_str.eq("Bcell")
    fine[is_b] = "B cell"
    fine[is_b & sub_str.str.contains("plasmaB", regex=False)] = "plasma cell"
    fine[major_str.eq("MAST")] = "mast cell"
    fine[major_str.eq("ILC3")] = "ILC3"
    fine[major_str.eq("CD45-")] = "Nonimmune"
    fine[is_cd3e_doublet.to_numpy()] = "Cross_lineage_doublet"
    out["author_subtype"] = pd.Series(fine, index=out.index, dtype=str)
    return out


def match_population_mask(
    target: str,
    broad: np.ndarray,
    subtype: np.ndarray,
    myeloid_fine: np.ndarray | None = None,
) -> np.ndarray:
    """Match either a broad lineage (`T/NK`, `Mono/Mac/DC`, `B/Plasma`) or fine subtype."""
    if target in ("T/NK", "Mono/Mac/DC", "B/Plasma", "Mast", "Neutrophil", "ILC3", "Nonimmune"):
        return broad == target
    if target == "monocyte/macrophage":
        return (broad == "Mono/Mac/DC") & np.isin(subtype, ("monocyte/macrophage", "Macrophage", "Monocyte"))
    if target in ("Macrophage", "Monocyte"):
        if myeloid_fine is not None:
            return (broad == "Mono/Mac/DC") & ((subtype == target) | (myeloid_fine == target))
        return (broad == "Mono/Mac/DC") & (subtype == target)
    parent = PARENT_LINEAGE_MAP.get(target)
    if parent is not None:
        return (broad == parent) & (subtype == target)
    return (broad == target) | (subtype == target)


def leave_one_out(patients: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (lineage, program), group in patients.groupby(["lineage", "program"]):
        available = group.dropna(subset=["D_T"])
        for patient in available.donor_uid:
            others = available[available.donor_uid != patient].D_T
            rows.append(
                (
                    lineage,
                    program,
                    patient,
                    len(others),
                    float(others.mean()) if len(others) else np.nan,
                )
            )
    return pd.DataFrame(
        rows, columns=["lineage", "program", "omitted_patient", "n_remaining", "mean_D_T"]
    )


def decompose(manifest: pd.DataFrame, subtypes: pd.DataFrame) -> pd.DataFrame:
    """Lesion-equal-weighted symmetric composition/state decomposition via starter.decompose_pair."""
    usable = subtypes[subtypes.n_cells >= 50].merge(
        manifest[["sample_uid", "donor_uid", "lesion_uid", "tissue"]],
        on="sample_uid",
        validate="many_to_one",
    )
    usable = usable[usable.tissue.isin(["PT", "LM"])].copy()
    if usable.empty:
        return pd.DataFrame(
            columns=[
                "donor_uid", "lineage", "program", "min_subtype_cells", "status",
                "observed_mean_PT", "observed_mean_LM", "observed_total_difference",
                "full_composition", "full_within_type_state",
                "common_retained_mass_PT", "common_retained_mass_LM",
                "common_composition", "common_within_type_state", "unidentified_subtypes",
            ]
        )
    usable["proportion"] = usable.n_subtype / usable.n_cells
    usable["score_mass"] = usable.proportion * usable.mean_score
    rows = []
    for (donor_uid, lineage, program), group in usable.groupby(["donor_uid", "lineage", "program"]):
        frame = []
        observed: dict[str, float] = {}
        names = sorted(group.subtype.unique())
        for tissue, tissue_rows in group.groupby("tissue"):
            n_lesions = tissue_rows.lesion_uid.nunique()
            agg = (
                tissue_rows.groupby("subtype")[["n_subtype", "proportion", "score_mass"]]
                .sum()
                .reindex(names, fill_value=0)
            )
            p = agg.proportion / n_lesions
            mu = (agg.score_mass / agg.proportion).where(p > 0)
            observed[tissue] = float(agg.score_mass.sum() / n_lesions)
            frame.extend(
                dict(
                    tissue=tissue,
                    subtype=s,
                    proportion=float(p.at[s]),
                    mean_score=float(mu.at[s]) if pd.notna(mu.at[s]) else np.nan,
                    n_subtype=int(agg.n_subtype.at[s]),
                )
                for s in names
            )
        values = pd.DataFrame(frame)
        unmasked_pt = observed.get("PT")
        unmasked_lm = observed.get("LM")
        unmasked_diff = (
            (unmasked_lm - unmasked_pt)
            if (unmasked_pt is not None and unmasked_lm is not None)
            else None
        )
        for support_min in (1, 20):
            checked = values.copy()
            checked.loc[checked.n_subtype < support_min, "mean_score"] = np.nan
            pair = decompose_pair(checked)
            common = pair.get("common_support") or {}
            full = pair.get("full_decomposition") or {}
            cond = common.get("renormalized_components") or {}
            rows.append(
                {
                    "donor_uid": donor_uid,
                    "lineage": lineage,
                    "program": program,
                    "min_subtype_cells": support_min,
                    "status": pair["status"],
                    "observed_mean_PT": unmasked_pt if pair["status"] != "missing_tissue" else None,
                    "observed_mean_LM": unmasked_lm if pair["status"] != "missing_tissue" else None,
                    "observed_total_difference": unmasked_diff if pair["status"] != "missing_tissue" else None,
                    "full_composition": full.get("composition"),
                    "full_within_type_state": full.get("within_type_state"),
                    "common_retained_mass_PT": common.get("retained_mass_PT"),
                    "common_retained_mass_LM": common.get("retained_mass_LM"),
                    "common_composition": cond.get("composition"),
                    "common_within_type_state": cond.get("within_type_state"),
                    "unidentified_subtypes": "|".join(pair.get("unidentified_subtypes") or []),
                }
            )
    return pd.DataFrame(rows)


def load_frozen_matched_controls(
    freeze_dir: Path,
    keys: Iterable[tuple[str, str]],
    registry: pd.DataFrame,
    shared_axis: Sequence[str] | None = None,
    *,
    strict_keys: Iterable[tuple[str, str]] | None = None,
) -> dict[tuple[str, str], dict[str, list[str]]]:
    """Load frozen development matched_controls.csv from `freeze_dir` without external refitting."""
    freeze_resolved = Path(freeze_dir).resolve()
    root_resolved = ROOT.resolve()

    # Priority 1: directly inside freeze_dir so custom/archived --freeze-dir never cross-reads live M5.
    candidates_paths: list[Path] = [
        freeze_resolved / "gene_sets/matched_controls.csv",
        freeze_resolved / "matched_controls.csv",
    ]
    if "matched_controls_file" in registry.columns:
        for rel in registry["matched_controls_file"].dropna().unique():
            rel_s = str(rel).strip()
            if not rel_s:
                continue
            if rel_s.startswith("artifacts/m5/"):
                candidates_paths.append(freeze_resolved / Path(rel_s).relative_to("artifacts/m5"))
            else:
                candidates_paths.append((root_resolved / rel_s).resolve())

    mc_path = next((p for p in candidates_paths if p.is_file()), None)
    requires_debiased = registry_requires_debiased(registry)
    if mc_path is None:
        if requires_debiased:
            raise FileNotFoundError(
                f"V1.1 debiased external scoring requires frozen matched_controls.csv under {freeze_resolved}; "
                "refitting control bins on external cohorts is prohibited"
            )
        return {}

    axis_set = (
        set(str(g) for g in shared_axis)
        if shared_axis is not None
        else set(load_shared_nuclear_universe("245_164_shared_nuclear"))
    )
    mc_df = pd.read_csv(mc_path, dtype=str, keep_default_na=False)
    mc_df["replicate"] = mc_df["replicate"].astype(int)
    mapper = get_hgnc_mapper()
    uniq_genes = set(mc_df["target_gene"]).union(mc_df["control_gene"])
    canon_lookup = {
        g: (g if g in axis_set else mapper.resolve_symbol(g).canonical_symbol)
        for g in uniq_genes
    }
    mc_df["target_gene"] = mc_df["target_gene"].map(canon_lookup)
    mc_df["control_gene"] = mc_df["control_gene"].map(canon_lookup)

    key_list = list(keys)
    # Also attempt to load dedicated M3 `__remove_liver_rbc_HLA_DRB345` control sets when present.
    extra_variant_keys = [(lin, f"{prog}__remove_liver_rbc_HLA_DRB345") for lin, prog in key_list]

    resolved: dict[tuple[str, str], dict[str, list[str]]] = {}
    has_lin = "lineage" in mc_df.columns
    for lin, prog in key_list + extra_variant_keys:
        if has_lin:
            sub = mc_df.loc[mc_df["lineage"].eq(lin) & mc_df["program"].eq(prog)]
            if sub.empty and lin in PARENT_LINEAGE_MAP:
                sub = mc_df.loc[mc_df["lineage"].eq(PARENT_LINEAGE_MAP[lin]) & mc_df["program"].eq(prog)]
        else:
            sub = mc_df.loc[mc_df["program"].eq(prog)]
        if not sub.empty:
            parsed = controls_from_dataframe(sub[["program", "target_gene", "replicate", "control_gene"]])
            if prog in parsed and parsed[prog]:
                resolved[(lin, prog)] = parsed[prog]

    if requires_debiased:
        must_have = list(strict_keys) if strict_keys is not None else key_list
        missing_ctrls = [f"{lin}::{prog}" for lin, prog in must_have if (lin, prog) not in resolved]
        if missing_ctrls:
            raise ValueError(
                f"Missing frozen matched_controls in {mc_path} for required programs {missing_ctrls}; "
                "silent fallback to non-debiased scores is prohibited"
            )
    return resolved


def sources():
    audit = json.loads((GEO / "GSE164522_all_expression_headers_audit.json").read_text())
    return {r["tissue_code"]: r for r in audit["files"] if r["tissue_code"] in ORGANS}


def range_file(job):
    url, first, last, size, dest = job
    length = last - first + 1
    if dest.exists() and dest.stat().st_size == length:
        return
    headers = dest.with_suffix(".headers")
    tmp = dest.with_suffix(".partial")
    subprocess.run(
        [
            "curl", "-fsSL", "--retry", "5", "--retry-all-errors", "--connect-timeout", "12",
            "--speed-limit", "1024", "--speed-time", "60", "--max-time", "360",
            "--range", f"{first}-{last}", "-D", str(headers), "-o", str(tmp), url,
        ],
        check=True,
    )
    expected = f"content-range: bytes {first}-{last}/{size}"
    if expected not in headers.read_text().lower().splitlines() or tmp.stat().st_size != length:
        raise ValueError(f"GEO range reply differs for {dest.name}")
    tmp.replace(dest)
    headers.unlink()


def fetch(records, out_dir: Path = OUT):
    RAW.mkdir(parents=True, exist_ok=True)
    if not METADATA.exists():
        tmp = METADATA.with_suffix(".gz.part")
        base = next(iter(records.values()))["url"].rsplit("/", 1)[0]
        subprocess.run(
            ["curl", "-fsSL", "--retry", "4", "--retry-all-errors", "-o", str(tmp), base + "/" + METADATA.name],
            check=True,
        )
        if sha(tmp) != next(iter(records.values()))["metadata_file_sha256"]:
            raise ValueError("Downloaded metadata differs from archived GEO metadata")
        tmp.replace(METADATA)
    ranges = RAW / "ranges"
    ranges.mkdir(exist_ok=True)
    jobs, files = [], {}
    for code, r in records.items():
        size = int(r["http_content_length"])
        path = RAW / r["url"].rsplit("/", 1)[1]
        if path.exists() and path.stat().st_size == size:
            continue
        prefix = Path(str(path) + ".part")
        start = prefix.stat().st_size if prefix.exists() else 0
        if start > size:
            raise ValueError(f"Invalid prefix for {code}")
        block = 8 << 20
        parts = []
        for left in range(start, size, block):
            right = min(left + block, size) - 1
            part = ranges / f"{code}_{left}_{right}.bytes"
            parts.append(part)
            jobs.append((code, (r["url"], left, right, size, part)))
        files[code] = (prefix, path, parts, size)

    def join(code):
        prefix, path, parts, size = files[code]
        if all(p.is_file() for p in parts) and not path.is_file():
            with prefix.open("ab") as out:
                for part in parts:
                    with part.open("rb") as stream:
                        shutil.copyfileobj(stream, out, length=1 << 20)
            if prefix.stat().st_size != size:
                raise ValueError(f"Invalid assembled size for {code}")
            prefix.replace(path)
            for p in parts:
                p.unlink()
            print(f"{now()} {code} complete {size} bytes", flush=True)

    with ThreadPoolExecutor(max_workers=8) as executor:
        futures = {executor.submit(range_file, job): code for code, job in jobs}
        for fut in as_completed(futures):
            fut.result()
            join(futures[fut])
    for code in files:
        join(code)
    summary = []
    for code, r in records.items():
        file = RAW / r["url"].rsplit("/", 1)[1]
        summary.append(
            dict(
                tissue=ORGANS[code],
                path=str(file.relative_to(ROOT)),
                url=r["url"],
                bytes=file.stat().st_size,
                sha256=sha(file),
            )
        )
    summary.append(
        dict(
            tissue="metadata",
            path=str(METADATA.relative_to(ROOT)),
            url=next(iter(records.values()))["url"].rsplit("/", 1)[0] + "/" + METADATA.name,
            bytes=METADATA.stat().st_size,
            sha256=sha(METADATA),
        )
    )
    pd.DataFrame(summary).to_csv(out_dir / "source_inputs.csv", index=False)


def metadata(records):
    expected = next(iter(records.values()))["metadata_file_sha256"]
    if sha(METADATA) != expected:
        raise ValueError("Complete author metadata differs from the GEO audit")
    frame = pd.read_csv(METADATA, index_col=0)
    frame.index.name = "source_barcode"
    identity = pd.read_csv(GEO / "GSE164522_metadata_identity_columns.csv.gz", index_col=0)
    if (
        not frame.index.is_unique
        or not frame.index.equals(identity.index)
        or not np.array_equal(
            frame[["patient", "tissue"]].to_numpy(),
            identity[["patient_code", "tissue"]].to_numpy(),
        )
    ):
        raise ValueError("Cell/donor/tissue mapping differs from the archived original identity table")
    tissues = {r["metadata_tissue"]: ORGANS[code] for code, r in records.items()}
    frame = frame.loc[frame.tissue.isin(tissues)].copy()
    frame["author_tissue"] = frame.tissue
    frame["tissue"] = frame.author_tissue.map(tissues)
    frame["source_barcode"] = frame.index
    frame["expression_column"] = frame.index.str.replace("-", ".", regex=False)
    frame = map_author_annotations(frame)
    if not frame.expression_column.is_unique:
        raise ValueError("Normalized barcode collision in GSE164522")
    frame["donor_uid"] = GSE + ":" + frame.patient
    frame["sample_uid"] = frame.donor_uid + ":" + frame["sample"].str.replace(" ", "_", regex=False)
    frame["lesion_uid"] = frame.sample_uid
    return frame.set_index("expression_column", drop=False)


def header(stream):
    return next(csv.reader([stream.readline().decode("utf-8-sig").rstrip("\r\n")]))[1:]


def write_identity(frame: pd.DataFrame, records, out_dir: Path = OUT) -> pd.DataFrame:
    cov = out_dir / "coverage"
    cov.mkdir(parents=True, exist_ok=True)
    checks = []
    for code, record in records.items():
        raw = RAW / record["url"].rsplit("/", 1)[1]
        path = raw if raw.is_file() else GEO / f"GSE164522_{code}_expression_header.csv.gz"
        with gzip.open(path, "rb") as stream:
            line = stream.readline()
        if hashlib.sha256(line).hexdigest() != record["header_sha256"]:
            raise ValueError(f"Archived and selected expression header differ for {code}")
        cols = pd.Index(next(csv.reader([line.decode("utf-8-sig").rstrip("\r\n")]))[1:])
        meta = frame.index[frame.tissue.eq(ORGANS[code])]
        if not cols.is_unique or len(meta) != len(cols) or not cols.isin(meta).all():
            raise ValueError(f"Header mismatch for {code}")
        checks.append(
            dict(
                code=code,
                tissue=ORGANS[code],
                n_header_cells=len(cols),
                n_matched=len(meta),
                header_source=str(path.relative_to(ROOT)),
                header_sha256=record["header_sha256"],
            )
        )
    pd.DataFrame(checks).to_csv(cov / "header_audit.csv", index=False)
    keys = ["donor_uid", "sample_uid", "lesion_uid", "patient", "tissue", "author_tissue", "sample", "ID"]
    manifest = frame.groupby(keys, observed=True).size().rename("input_cells").reset_index()
    clean_counts = (
        frame.loc[frame["clean_author_qc"]]
        .groupby("sample_uid", observed=True)
        .size()
        .rename("clean_cells")
    )
    doublet_counts = (
        frame.loc[frame["is_cd3e_doublet"]]
        .groupby("sample_uid", observed=True)
        .size()
        .rename("excluded_cd3e_doublets")
    )
    manifest = (
        manifest.merge(clean_counts, on="sample_uid", how="left")
        .merge(doublet_counts, on="sample_uid", how="left")
    )
    manifest["clean_cells"] = manifest["clean_cells"].fillna(0).astype(int)
    manifest["excluded_cd3e_doublets"] = manifest["excluded_cd3e_doublets"].fillna(0).astype(int)
    if manifest.sample_uid.duplicated().any() or manifest.ID.duplicated().any():
        raise ValueError("Original author ID and patient/specimen mapping disagree")
    manifest["dataset"] = GSE
    manifest["role"] = "external"
    manifest["input_scale"] = "normalized"
    manifest["patient_mapping_verified"] = True
    manifest["treatment"] = "treatment-naive at cohort level in B02"
    manifest["episode"] = "source patient and specimen; exact collection dates unavailable"
    manifest.to_csv(cov / "sample_manifest.csv", index=False)
    frame[
        [
            "source_barcode", "expression_column", *keys,
            "lineage", "author_subtype", "is_cd3e_doublet", "clean_author_qc",
            "celltype_global", "celltype_major", "celltype_sub",
        ]
    ].to_csv(cov / "cell_identity.csv.gz", index=False, compression={"method": "gzip", "mtime": 0})
    frame[
        ["celltype_global", "celltype_major", "lineage", "author_subtype", "is_cd3e_doublet", "celltype_sub"]
    ].drop_duplicates().sort_values(["lineage", "celltype_major", "celltype_sub"]).to_csv(
        cov / "author_lineage_mapping.csv", index=False
    )
    return manifest


def parse_published(path: Path, cells: pd.DataFrame):
    """Read one numeric row at a time; only nonzero values enter CSR storage."""
    genes, blocks, positions, values, offsets = [], [], [], [], [0]

    def flush(n_cells: int):
        if values:
            block = sparse.csr_matrix(
                (np.concatenate(values), np.concatenate(positions), np.asarray(offsets)),
                shape=(len(values), n_cells),
            )
            blocks.append(block)
            values.clear()
            positions.clear()
            offsets[:] = [0]

    with gzip.open(path, "rb") as stream:
        actual = pd.Index(header(stream))
        if not actual.is_unique or len(cells) != len(actual) or not actual.isin(cells.index).all():
            raise ValueError(f"{path.name}: author cells and expression columns disagree")
        for line in stream:
            symbol, sep, numeric = line.partition(b",")
            if not sep:
                raise ValueError(f"Missing gene row in {path.name}")
            genes.append(symbol.strip(b'"').decode())
            vector = np.fromstring(numeric, dtype=np.float64, sep=",")
            if len(vector) != len(actual) or not np.isfinite(vector).all() or vector.min() < 0:
                raise ValueError(f"{path.name}: invalid published numeric row {len(genes)}")
            keep = np.flatnonzero(vector).astype(np.int32)
            positions.append(keep)
            values.append(vector[keep])
            offsets.append(offsets[-1] + len(keep))
            if len(values) == 256:
                flush(len(actual))
        flush(len(actual))
    if not genes or "" in genes:
        raise ValueError(f"{path.name}: missing published gene axis")
    return sparse.vstack(blocks, format="csr").T.tocsr(), genes, actual


def convert(code: str, frame: pd.DataFrame, record, objects_dir: Path = OBJECTS, out_dir: Path = OUT):
    path = RAW / record["url"].rsplit("/", 1)[1]
    if path.stat().st_size != int(record["http_content_length"]):
        raise ValueError(f"{path.name}: incomplete GEO file")
    compressed_sha = sha(path)
    dest = objects_dir / f"{code}.h5ad"
    if dest.exists() and (objects_dir / f"{code}.json").exists():
        saved = ad.read_h5ad(dest, backed="r")
        same = (
            saved.uns.get("expression_sha256") == compressed_sha
            and saved.uns.get("metadata_sha256") == record["metadata_file_sha256"]
            and "clean_author_qc" in saved.obs.columns
            and "canonical_symbol" in saved.var.columns
        )
        saved.file.close()
        if same:
            return
    start = time.perf_counter()
    subset = frame.loc[frame.tissue.eq(ORGANS[code])]
    matrix, genes, columns = parse_published(path, subset)
    obs = subset.loc[columns].copy()
    obs.index = pd.Index(GSE + ":" + obs.source_barcode.to_numpy(str), name="cell_uid")
    gene_index = pd.Index(genes, name="source_gene")
    norm_df = normalize_gene_axis(genes, excluded_canonical=set(load_shared_nuclear_universe("excluded_353_canonical")))
    var = pd.DataFrame(
        {
            "gene_symbol": genes,
            "unambiguous_symbol": ~gene_index.duplicated(keep=False),
            "canonical_symbol": norm_df["canonical_symbol"].to_numpy(),
            "hgnc_id": norm_df["hgnc_id"].to_numpy(),
            "hgnc_match_type": norm_df["match_type"].to_numpy(),
            "is_unique_measurable": norm_df["is_unique_measurable"].to_numpy(),
            "in_excluded_353": norm_df["in_excluded_353"].to_numpy(),
        },
        index=gene_index,
    )
    obj = ad.AnnData(
        X=matrix,
        obs=obs,
        var=var,
        uns={
            "expression_scale": "author processed continuous; published value in X",
            "expression_sha256": compressed_sha,
            "metadata_sha256": record["metadata_file_sha256"],
            "source_url": record["url"],
            "tissue_code": code,
        },
    )
    obj.var_names_make_unique(join="::duplicate-row-")
    scales = np.empty(len(obs))
    for first in range(0, len(obs), 512):
        block = matrix[first : first + 512].copy()
        block.data = np.expm1(block.data)
        scales[first : first + len(block.indptr) - 1] = np.asarray(block.sum(axis=1)).ravel()
    check = pd.DataFrame(
        {
            "sample_uid": obs.sample_uid.to_numpy(),
            "inverse_sum_diagnostic": scales,
            "within_one_of_10000": np.abs(scales - 10000) <= 1,
            "within_9500_10500": (scales >= 9500) & (scales <= 10500),
        }
    )
    check = check.loc[np.diff(matrix.indptr) > 0]
    scale_summary = (
        check.groupby("sample_uid")
        .agg(
            n_nonempty=("inverse_sum_diagnostic", "size"),
            min_nonempty_inverse_sum=("inverse_sum_diagnostic", "min"),
            median_nonempty_inverse_sum=("inverse_sum_diagnostic", "median"),
            max_nonempty_inverse_sum=("inverse_sum_diagnostic", "max"),
            fraction_within_one_of_10000=("within_one_of_10000", "mean"),
            fraction_within_five_percent_of_10000=("within_9500_10500", "mean"),
        )
        .reindex(sorted(obs.sample_uid.unique()))
    )
    scale_summary["n_nonempty"] = scale_summary.n_nonempty.fillna(0).astype(int)
    (out_dir / "coverage").mkdir(parents=True, exist_ok=True)
    scale_summary.to_csv(out_dir / "coverage" / f"{code}_scale_check.csv")
    out = dest.with_suffix(".partial.h5ad")
    obj.write_h5ad(out, compression="lzf")
    reread = ad.read_h5ad(out, backed="r")
    checked = 0
    try:
        if (
            not reread.obs_names.equals(obj.obs_names)
            or not reread.var_names.equals(obj.var_names)
            or reread.X.dtype != np.float64
        ):
            raise ValueError(f"{code}: output axes or dtype differ")
        for row in (0, len(obs) // 2, len(obs) - 1):
            if (matrix[row] != reread.X[row]).nnz:
                raise ValueError(f"{code}: readback scores differ")
        with gzip.open(path, "rt", newline="") as stream:
            table = csv.reader(stream)
            source_columns = next(table)[1:]
            for gene_no in range(min(8, len(genes))):
                gene, *v = next(table)
                if reread.var.gene_symbol.iloc[gene_no] != gene:
                    raise ValueError(f"{code}: source feature mismatch")
                positive = next((j for j, value in enumerate(v) if float(value) > 0), 0)
                for j in (0, len(obs) // 2, len(obs) - 1, positive):
                    if source_columns[j] != reread.obs.expression_column.iloc[j] or float(v[j]) != reread.X[j, gene_no]:
                        raise ValueError(f"{code}: published value/column mismatch")
                    checked += 1
    finally:
        reread.file.close()
    out.replace(dest)
    entry = dict(
        code=code,
        tissue=ORGANS[code],
        input_bytes=path.stat().st_size,
        source_sha256=compressed_sha,
        genes=len(genes),
        cells=len(obs),
        nnz=matrix.nnz,
        source_entries_checked=checked,
        ambiguous_symbols=sorted(set(gene_index[gene_index.duplicated(keep=False)])),
        conversion_seconds=round(time.perf_counter() - start, 2),
        h5ad_bytes=dest.stat().st_size,
        h5ad_sha256=sha(dest),
        completed_utc=now(),
    )
    (objects_dir / f"{code}.json").write_text(json.dumps(entry, indent=2) + "\n")
    print(f"Converted {code}: {len(genes)} x {len(obs)}", flush=True)


def frozen_definitions(freeze_dir: Path = FREEZE, objects_dir: Path = OBJECTS, out_dir: Path = OUT):
    freeze_dir = Path(freeze_dir).resolve()
    objects_dir = Path(objects_dir).resolve()
    out_dir = Path(out_dir).resolve()

    reg = pd.read_csv(freeze_dir / "frozen_candidate_registry.csv")
    if reg.empty or not reg.status.astype(str).str.lower().eq("frozen").all():
        raise ValueError("A frozen candidate registry is required before external validation")
    genes = pd.read_csv(freeze_dir / "gene_sets/program_genes.csv")
    model_path = freeze_dir / "gene_sets/model_program_genes.csv"
    model_genes = pd.read_csv(model_path) if model_path.is_file() else pd.DataFrame(columns=["lineage", "program", "gene_symbol"])

    axis = list(load_shared_nuclear_universe("245_164_shared_nuclear"))
    for code in ORGANS:
        obj = ad.read_h5ad(objects_dir / f"{code}.h5ad", backed="r")
        try:
            al = align_axis_to_canonical_targets(obj.var["gene_symbol"].tolist(), axis)
            if not al["present"].all():
                missing = al.loc[~al["present"], "canonical_symbol"].tolist()[:10]
                raise ValueError(f"GSE164522 {code}: missing expected shared nuclear genes {missing}")
        finally:
            obj.file.close()

    cov = out_dir / "coverage"
    cov.mkdir(parents=True, exist_ok=True)
    (cov / "external_rank_background.txt").write_text("\n".join(axis) + "\n")
    background = set(axis)

    candidate_keys = set(zip(reg.lineage.astype(str), reg.program.astype(str)))
    model_keys = set(zip(model_genes.lineage.astype(str), model_genes.program.astype(str)))
    sec_path = freeze_dir / "secondary_sensitivity_candidates.csv"
    sec_keys: set[tuple[str, str]] = set()
    if sec_path.is_file():
        sec_df = pd.read_csv(sec_path)
        sec_keys = set(zip(sec_df.lineage.astype(str), sec_df.program.astype(str)))
    if ("B/Plasma", "HALLMARK_UNFOLDED_PROTEIN_RESPONSE") in model_keys:
        sec_keys.add(("B/Plasma", "HALLMARK_UNFOLDED_PROTEIN_RESPONSE"))
    eval_keys = candidate_keys | (sec_keys & model_keys)
    all_keys = sorted(candidate_keys | model_keys)
    requires_debiased = registry_requires_debiased(reg)

    reg_lookup = {(str(r.lineage), str(r.program)): r for r in reg.itertuples(index=False)}
    programs: dict[tuple[str, str], list[str]] = {}
    coverage = []
    for lineage, program in all_keys:
        source = genes if (lineage, program) in candidate_keys else model_genes
        raw_g = source.loc[source.program.eq(program) & source.lineage.eq(lineage), "gene_symbol"].tolist()
        if not raw_g:
            raise ValueError(f"Invalid frozen members: {lineage}::{program}")
        canon_g = canonicalize_preserving_axis(raw_g, background)
        if (lineage, program) in candidate_keys and (lineage, program) in model_keys:
            model_raw = model_genes.loc[model_genes.program.eq(program) & model_genes.lineage.eq(lineage), "gene_symbol"].tolist()
            if set(canon_g) != set(canonicalize_preserving_axis(model_raw, background)):
                raise ValueError(f"Frozen program differs between biology and model: {program}")
        r_meta = reg_lookup.get((lineage, program))
        min_pos = int(r_meta.min_measured_positive_genes) if r_meta is not None else 5
        min_frac = float(r_meta.min_measured_fraction_positive) if r_meta is not None else 0.8
        measured = [s for s in canon_g if s in background]
        sufficient = len(measured) >= min_pos and len(measured) / len(canon_g) >= min_frac
        programs[(lineage, program)] = measured if sufficient else []
        coverage.append(
            dict(
                lineage=lineage,
                program=program,
                primary_candidate=(lineage, program) in candidate_keys,
                secondary_or_falsification=(lineage, program) in (eval_keys - candidate_keys),
                model_input=(lineage, program) in model_keys,
                n_axis_genes=len(axis),
                frozen_members_raw=len(raw_g),
                frozen_members_canonical=len(canon_g),
                measured_members=len(measured),
                measured_fraction=len(measured) / len(canon_g),
                eligible=sufficient,
                absent_members="|".join(sorted(set(canon_g) - set(measured))),
            )
        )

    eligible_scored_keys = [k for k in all_keys if programs[k]]
    controls_map = load_frozen_matched_controls(
        freeze_dir,
        all_keys,
        reg,
        shared_axis=axis,
        strict_keys=eligible_scored_keys,
    )
    cov_df = pd.DataFrame(coverage)
    cov_df["has_frozen_matched_controls"] = [
        (r["lineage"], r["program"]) in controls_map for _, r in cov_df.iterrows()
    ]
    cov_df.to_csv(cov / "candidate_measurement.csv", index=False)
    return axis, programs, controls_map, candidate_keys, eval_keys, model_keys, requires_debiased


def score_anndata_dict(
    h5ads_by_code: Mapping[str, ad.AnnData],
    meta_frame: pd.DataFrame,
    axis: Sequence[str],
    programs: Mapping[tuple[str, str], list[str]],
    controls_map: Mapping[tuple[str, str], Mapping[str, Sequence[str]]],
    candidate_keys: set[tuple[str, str]],
    *,
    requires_debiased: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Score GSE164522 across clean-author, raw-author, and unified annotation_m1 branches."""
    if requires_debiased and not controls_map and any(programs.values()):
        raise ValueError(
            "requires_debiased=True but controls_map is empty; silent fallback to non-debiased scores is prohibited"
        )

    values = []
    concordance_rows = []
    comp_rows = []
    meta_by_bc = meta_frame.set_index("source_barcode", drop=False)

    for code, obj in h5ads_by_code.items():
        al = align_axis_to_canonical_targets(obj.var["gene_symbol"].tolist(), axis)
        if not al["present"].all():
            raise ValueError(f"GSE164522 {code}: missing shared nuclear feature")
        positions = al["source_index"].to_numpy(dtype=np.int32)

        # Canonicalize full axis for unified two-stage rule annotation (`annotation_m1.classify`).
        full_norm = normalize_gene_axis(obj.var["gene_symbol"].tolist())
        rule_labels, rule_details = classify(obj.X, full_norm["canonical_symbol"].tolist())
        rule_broad = np.asarray(rule_labels, dtype=str)
        rule_sub = rule_details["subtype_hierarchical"].astype(str).to_numpy()
        rule_mf = rule_details["myeloid_fine"].fillna("").astype(str).to_numpy()

        bcs = obj.obs["source_barcode"].astype(str).to_numpy()
        aligned_meta = meta_by_bc.loc[bcs]
        auth_broad = aligned_meta["lineage"].astype(str).to_numpy()
        auth_sub = aligned_meta["author_subtype"].astype(str).to_numpy()
        clean_qc = aligned_meta["clean_author_qc"].to_numpy(dtype=bool)
        is_doublet = aligned_meta["is_cd3e_doublet"].to_numpy(dtype=bool)

        ctab = (
            pd.DataFrame(
                {
                    "tissue": ORGANS[code],
                    "author_lineage": auth_broad,
                    "author_subtype": auth_sub,
                    "clean_author_qc": clean_qc,
                    "rule_lineage_m1": rule_broad,
                    "rule_subtype_m1": rule_sub,
                }
            )
            .value_counts()
            .rename("n_cells")
            .reset_index()
        )
        concordance_rows.append(ctab)

        sample_arr = aligned_meta["sample_uid"].astype(str).to_numpy()
        for samp in np.unique(sample_arr):
            sm = sample_arr == samp
            n_cap = int(sm.sum())
            n_clean = int((sm & clean_qc).sum())
            n_auth_tot_imm = int((sm & clean_qc & np.isin(auth_broad, AUTHOR_RESOLVED_IMMUNE)).sum())
            n_auth_c3_imm = int((sm & clean_qc & np.isin(auth_broad, CORE_THREE_IMMUNE)).sum())
            n_rule_tot_imm = int((sm & clean_qc & np.isin(rule_broad, RULE_RESOLVED_IMMUNE)).sum())
            n_rule_c3_imm = int((sm & clean_qc & np.isin(rule_broad, CORE_THREE_IMMUNE)).sum())
            for lin_name in ("T/NK", "Mono/Mac/DC", "B/Plasma", "Mast", "ILC3", "Neutrophil"):
                n_a = int((sm & clean_qc & (auth_broad == lin_name)).sum())
                n_r = int((sm & clean_qc & (rule_broad == lin_name)).sum())
                comp_rows.append(
                    {
                        "sample_uid": samp,
                        "tissue": ORGANS[code],
                        "lineage": lin_name,
                        "n_captured": n_cap,
                        "n_clean_author_sample": n_clean,
                        "n_author_clean_cells": n_a,
                        "n_author_total_resolved_immune": n_auth_tot_imm,
                        "n_author_core_three_immune": n_auth_c3_imm,
                        "prop_author_total_resolved_immune": (
                            n_a / n_auth_tot_imm if n_auth_tot_imm > 0 and lin_name in AUTHOR_RESOLVED_IMMUNE else np.nan
                        ),
                        "prop_author_core_three_immune": (
                            n_a / n_auth_c3_imm if n_auth_c3_imm > 0 and lin_name in CORE_THREE_IMMUNE else np.nan
                        ),
                        "prop_author_captured": (n_a / n_clean if n_clean > 0 else np.nan),
                        "n_rule_m1_cells": n_r,
                        "n_rule_total_resolved_immune": n_rule_tot_imm,
                        "n_rule_core_three_immune": n_rule_c3_imm,
                        "prop_rule_total_resolved_immune": (
                            n_r / n_rule_tot_imm if n_rule_tot_imm > 0 and lin_name in RULE_RESOLVED_IMMUNE else np.nan
                        ),
                        "prop_rule_core_three_immune": (
                            n_r / n_rule_c3_imm if n_rule_c3_imm > 0 and lin_name in CORE_THREE_IMMUNE else np.nan
                        ),
                    }
                )

        targets = sorted({key[0] for key in programs if programs[key]})
        for target_lin in targets:
            active_keys = [
                k
                for k, members in programs.items()
                if k[0] == target_lin and members and ((k in candidate_keys) or (code in ("PT", "MT")))
            ]
            if not active_keys:
                continue

            auth_raw_match = match_population_mask(target_lin, auth_broad, auth_sub)
            auth_clean_match = auth_raw_match & clean_qc
            rule_match = match_population_mask(target_lin, rule_broad, rule_sub, rule_mf) & clean_qc
            union_mask = auth_raw_match | rule_match
            rows = np.flatnonzero(union_mask)
            if len(rows) == 0:
                continue

            measured = obj.X[rows][:, positions]
            prog_dict: dict[str, list[str]] = {}
            ctrl_dict: dict[str, Mapping[str, Sequence[str]]] = {}
            drop_map: dict[str, str] = {}

            for k in active_keys:
                p = k[1]
                members = programs[k]
                prog_dict[p] = members
                if k in controls_map:
                    ctrl_dict[p] = controls_map[k]
                elif requires_debiased:
                    raise ValueError(
                        f"Missing frozen matched_controls for {k[0]}::{p}; silent fallback is prohibited"
                    )
                if k in candidate_keys:
                    filt = [g for g in members if g not in LIVER_RBC and g not in ("HLA-DRB3", "HLA-DRB4", "HLA-DRB5")]
                    if len(filt) >= 5:
                        dkey = f"{p}::drop_liver_Hb_DRB345"
                        prog_dict[dkey] = filt
                        drop_map[p] = dkey
                        m3_var_key = (k[0], f"{p}__remove_liver_rbc_HLA_DRB345")
                        if m3_var_key in controls_map:
                            ctrl_dict[dkey] = controls_map[m3_var_key]
                        elif k in controls_map:
                            ctrl_dict[dkey] = controls_map[k]

            leg_df = legacy_positive_rank_scores(measured, axis, prog_dict, max_rank=1500)
            if requires_debiased:
                if len(ctrl_dict) != len(prog_dict):
                    missing_p = sorted(set(prog_dict) - set(ctrl_dict))
                    raise ValueError(f"Missing frozen matched_controls for {missing_p}")
                deb_df, raw_u_df, _ = debiased_ucell_scores(
                    measured, axis, prog_dict, matched_controls=ctrl_dict, max_rank=1500, return_components=True
                )
                log_deb_df = log1p_cp10k_scores(
                    measured,
                    axis,
                    prog_dict,
                    matched_controls=ctrl_dict,
                    input_scale="log1p",
                    renormalize_subaxis=True,
                )
            else:
                raw_u_df = ucell_scores(measured, axis, prog_dict, max_rank=1500)
                log_deb_df = log1p_cp10k_scores(
                    measured, axis, prog_dict, debiased=False, input_scale="log1p", renormalize_subaxis=True
                )
                deb_df = leg_df.copy()

            base_obs = pd.DataFrame(
                {
                    "cell_uid": GSE + ":" + bcs[rows],
                    "source_barcode": bcs[rows],
                    "donor_uid": aligned_meta["donor_uid"].to_numpy(str)[rows],
                    "sample_uid": sample_arr[rows],
                    "tissue": aligned_meta["tissue"].to_numpy(str)[rows],
                    "lineage": target_lin,
                    "author_lineage": auth_broad[rows],
                    "author_subtype": auth_sub[rows],
                    "is_cd3e_doublet": is_doublet[rows],
                    "in_clean_author": auth_clean_match[rows],
                    "in_raw_author": auth_raw_match[rows],
                    "in_rule_m1": rule_match[rows],
                    "rule_lineage_m1": rule_broad[rows],
                    "rule_subtype_m1": rule_sub[rows],
                    "rule_myeloid_fine_m1": rule_mf[rows],
                    "celltype_sub": aligned_meta["celltype_sub"].to_numpy(str)[rows],
                    "celltype_major": aligned_meta["celltype_major"].to_numpy(str)[rows],
                    "celltype_global": aligned_meta["celltype_global"].to_numpy(str)[rows],
                    "n_genes_common": np.diff(measured.indptr),
                }
            )
            for k in active_keys:
                p = k[1]
                dkey = drop_map.get(p)
                primary_vec = deb_df[p].to_numpy() if requires_debiased else leg_df[p].to_numpy()
                drop_vec = (
                    (deb_df[dkey].to_numpy() if requires_debiased else leg_df[dkey].to_numpy())
                    if dkey is not None
                    else np.full(len(rows), np.nan)
                )
                values.append(
                    base_obs.assign(
                        program=p,
                        score=primary_vec,
                        ucell_debiased=deb_df[p].to_numpy(),
                        ucell_raw=raw_u_df[p].to_numpy(),
                        log1p_cp10k_debiased=log_deb_df[p].to_numpy(),
                        legacy_positive=leg_df[p].to_numpy(),
                        drop_liver_Hb_DRB345=drop_vec,
                    )
                )

    empty_cols = [
        "cell_uid", "source_barcode", "donor_uid", "sample_uid", "tissue", "lineage",
        "author_lineage", "author_subtype", "is_cd3e_doublet", "in_clean_author",
        "in_raw_author", "in_rule_m1", "rule_lineage_m1", "rule_subtype_m1",
        "rule_myeloid_fine_m1", "celltype_sub", "celltype_major", "celltype_global",
        "n_genes_common", "program", "score", "ucell_debiased", "ucell_raw",
        "log1p_cp10k_debiased", "legacy_positive", "drop_liver_Hb_DRB345",
    ]
    cells_df = pd.concat(values, ignore_index=True) if values else pd.DataFrame(columns=empty_cols)
    conc_df = pd.concat(concordance_rows, ignore_index=True) if concordance_rows else pd.DataFrame()
    comp_df = pd.DataFrame(comp_rows)
    return cells_df, conc_df, comp_df


def measure(
    axis,
    programs,
    controls_map,
    candidates,
    meta_frame,
    *,
    requires_debiased: bool = True,
    objects_dir: Path = OBJECTS,
    out_dir: Path = OUT,
):
    h5ads = {code: ad.read_h5ad(objects_dir / f"{code}.h5ad") for code in ORGANS}
    cells_df, conc_df, comp_df = score_anndata_dict(
        h5ads,
        meta_frame,
        axis,
        programs,
        controls_map,
        candidates,
        requires_debiased=requires_debiased,
    )
    cov = out_dir / "coverage"
    cov.mkdir(parents=True, exist_ok=True)
    conc_df.to_csv(cov / "author_vs_rule_annotation_concordance.csv", index=False)
    comp_df.to_csv(cov / "sample_immune_composition.csv", index=False)
    return cells_df, comp_df, h5ads


def model_inputs(
    cells: pd.DataFrame,
    keys: set[tuple[str, str]],
    support: pd.Series,
    *,
    h5ads_by_code: Mapping[str, ad.AnnData] | None = None,
    meta_frame: pd.DataFrame | None = None,
    axis: Sequence[str] | None = None,
    freeze_dir: Path = FREEZE,
    out_dir: Path = OUT,
):
    """Export clean PT/LM program inputs aligned to `models/m4/full/<lineage>/matched_controls.csv`."""
    rows = []
    freeze_resolved = Path(freeze_dir).resolve()
    model_dirs: dict[str, Path] = {}
    frozen_models_csv = freeze_resolved / "frozen_models.csv"
    if frozen_models_csv.is_file():
        fm = pd.read_csv(frozen_models_csv)
        for r in fm.loc[fm["accepted"] & fm["method"].eq("program_elastic")].itertuples(index=False):
            pdir = (ROOT / str(r.path)).resolve().parent
            if (pdir / "program_defs.json").is_file() and (pdir / "matched_controls.csv").is_file():
                model_dirs[str(r.lineage)] = pdir

    # When M4 model directories with dedicated program_defs.json + matched_controls.csv are present
    # alongside PT/MT AnnData objects, score bags identically to `predict_gse164522.score_frozen_program_bag`.
    if model_dirs and h5ads_by_code is not None and meta_frame is not None and axis is not None:
        axis_list = list(axis)
        meta_by_bc = meta_frame.set_index("source_barcode", drop=False)
        for code in ("PT", "MT"):
            if code not in h5ads_by_code:
                continue
            obj = h5ads_by_code[code]
            al = align_axis_to_canonical_targets(obj.var["gene_symbol"].tolist(), axis_list)
            shared_cols = al["source_index"].to_numpy(dtype=np.int64)
            bcs = obj.obs["source_barcode"].astype(str).to_numpy()
            aligned_meta = meta_by_bc.loc[bcs]
            clean_qc = aligned_meta["clean_author_qc"].to_numpy(dtype=bool)
            lin_arr = aligned_meta["lineage"].astype(str).to_numpy()
            samp_arr = aligned_meta["sample_uid"].astype(str).to_numpy()

            for samp in sorted(np.unique(samp_arr)):
                sm = samp_arr == samp
                donor_uid = str(aligned_meta.loc[bcs[np.flatnonzero(sm)[0]], "donor_uid"])
                tissue = ORGANS[code]
                for lineage, pdir in sorted(model_dirs.items()):
                    clean_idx = np.flatnonzero(sm & clean_qc & (lin_arr == lineage))
                    if len(clean_idx) < 50:
                        continue
                    uids = np.asarray([f"{GSE}:{b}" for b in bcs[clean_idx]], dtype=str)
                    clean_idx = clean_idx[np.argsort(uids)]
                    # Recalibrate expm1 -> CP10k -> log1p on shared nuclear axis exactly as predict_gse164522
                    sub_mat = obj.X[clean_idx][:, shared_cols].tocsr().copy()
                    sub_mat.sum_duplicates()
                    sub_mat.eliminate_zeros()
                    sub_mat.data = np.expm1(sub_mat.data.astype(np.float64))
                    nuc_sums = np.asarray(sub_mat.sum(axis=1)).ravel()
                    row_nnz = np.diff(sub_mat.indptr)
                    sub_mat.data = np.log1p(sub_mat.data * np.repeat(10000.0 / nuc_sums, row_nnz))

                    spec = json.loads((pdir / "program_defs.json").read_text(encoding="utf-8"))
                    m_ctrls = controls_from_dataframe(pd.read_csv(pdir / "matched_controls.csv"))
                    c_scores = debiased_ucell_scores(
                        sub_mat,
                        axis_list,
                        spec["programs"],
                        matched_controls=m_ctrls,
                        max_rank=int(spec.get("max_rank", 1500)),
                    )
                    if len(clean_idx) > 256:
                        rng = np.random.default_rng(cell_seed(f"{samp}|{lineage}"))
                        pick256 = np.sort(rng.choice(len(clean_idx), 256, replace=False))
                    else:
                        pick256 = np.arange(len(clean_idx))
                    for program in sorted(spec["programs"]):
                        rows.append(
                            dict(
                                sample_uid=samp,
                                lineage=lineage,
                                program=program,
                                score=float(c_scores[program].mean()),
                                score_sub256=float(c_scores[program].iloc[pick256].mean()),
                                n_cells=len(clean_idx),
                                n_population=len(clean_idx),
                                donor_uid=donor_uid,
                                tissue=tissue,
                                n_captured=int(support[samp]),
                            )
                        )
    else:
        clean_cells = cells.loc[cells["in_clean_author"]]
        for lineage, program in sorted(keys):
            values = clean_cells.loc[clean_cells.lineage.eq(lineage) & clean_cells.program.eq(program)]
            for sample_uid, group in values.groupby("sample_uid"):
                if group.tissue.iloc[0] not in ("PT", "LM") or len(group) < 50:
                    continue
                ordered = group.sort_values("cell_uid")
                if len(ordered) > 256:
                    rng = np.random.default_rng(cell_seed(f"{sample_uid}|{lineage}"))
                    sub256 = ordered.iloc[rng.choice(len(ordered), 256, replace=False)]
                else:
                    sub256 = ordered
                rows.append(
                    dict(
                        sample_uid=sample_uid,
                        lineage=lineage,
                        program=program,
                        score=float(ordered.score.mean()),
                        score_sub256=float(sub256.score.mean()),
                        n_cells=len(ordered),
                        n_population=len(ordered),
                        donor_uid=group.donor_uid.iloc[0],
                        tissue=group.tissue.iloc[0],
                        n_captured=int(support[sample_uid]),
                    )
                )

    target = out_dir / "model_predictions"
    target.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(
        rows,
        columns=[
            "sample_uid", "lineage", "program", "score", "score_sub256",
            "n_cells", "n_population", "donor_uid", "tissue", "n_captured",
        ],
    ).to_csv(target / "external_model_inputs.csv", index=False)


def summarize(
    manifest: pd.DataFrame,
    meta: pd.DataFrame,
    cells: pd.DataFrame,
    programs: Mapping[tuple[str, str], list[str]],
    comp_df: pd.DataFrame,
    out_dir: Path = OUT,
):
    effects = out_dir / "patient_effects"
    effects.mkdir(parents=True, exist_ok=True)
    sizes = meta.groupby("sample_uid").size()
    clean_sizes = meta.loc[meta["clean_author_qc"]].groupby("sample_uid").size()

    clean_cells = cells.loc[cells["in_clean_author"]]
    clean_counts = clean_cells.groupby(["sample_uid", "lineage", "program"]).size()
    rule_counts = cells.loc[cells["in_rule_m1"]].groupby(["sample_uid", "lineage", "program"]).size()

    coverage = []
    for sample in manifest.sample_uid:
        n_cap = int(sizes.get(sample, 0))
        n_cln = int(clean_sizes.get(sample, 0))
        for lineage, p in programs:
            n = int(clean_counts.get((sample, lineage, p), 0))
            n_rule = int(rule_counts.get((sample, lineage, p), 0))
            entry = dict(
                sample_uid=sample,
                lineage=lineage,
                program=p,
                n_cells=n,
                n_cells_rule_m1=n_rule,
                n_clean_sample=n_cln,
                n_captured=n_cap,
                captured_fraction=(n / n_cap if n_cap else np.nan),
            )
            for minimum in (20, 50, 100):
                entry[f"status_at_{minimum}"] = (
                    "insufficient_gene_coverage" if not programs[(lineage, p)] else
                    "no_lineage_cells" if not n else
                    "below_specimen_cell_threshold" if n < minimum else "eligible"
                )
            coverage.append(entry)
    support = pd.DataFrame(coverage)
    (out_dir / "coverage").mkdir(parents=True, exist_ok=True)
    support.to_csv(out_dir / "coverage/sample_lineage_support.csv", index=False)

    policies = [
        ("main", "in_clean_author", "score", 50),
        ("n20", "in_clean_author", "score", 20),
        ("n100", "in_clean_author", "score", 100),
        ("log1p_cp10k_debiased", "in_clean_author", "log1p_cp10k_debiased", 50),
        ("ucell_raw", "in_clean_author", "ucell_raw", 50),
        ("legacy_positive", "in_clean_author", "legacy_positive", 50),
        ("rule_annotation_m1", "in_rule_m1", "score", 50),
        ("rule_annotation_m1_n20", "in_rule_m1", "score", 20),
        ("rule_annotation_m1_log1p_debiased", "in_rule_m1", "log1p_cp10k_debiased", 50),
        ("include_CD3E_doublets_raw_author", "in_raw_author", "score", 50),
        ("major_global_concordance", "in_clean_author", "score", 50),
        ("drop_liver_Hb_DRB345", "in_clean_author", "drop_liver_Hb_DRB345", 50),
    ]
    sensitive = []
    base_grid = support[["sample_uid", "lineage", "program"]].copy()
    for label, mask_col, score_col, threshold in policies:
        subset = cells.loc[cells[mask_col]]
        if label == "major_global_concordance":
            exp_glob = subset["lineage"].map(lambda l: EXPECTED_GLOBAL.get(PARENT_LINEAGE_MAP.get(l, l), ""))
            subset = subset.loc[exp_glob.eq("") | subset["celltype_global"].eq(exp_glob)]
        agg = subset.groupby(["sample_uid", "lineage", "program"])[score_col].agg(["mean", "size"]).reset_index()
        table = base_grid.merge(agg, how="left", on=["sample_uid", "lineage", "program"], validate="one_to_one")
        table["n_cells"] = table["size"].fillna(0).astype(int)
        table["score"] = table["mean"].where(table["n_cells"] >= threshold)
        patients, lesions, means, _ = calculate_patient_contrasts(manifest, table[SCORE_COLUMNS])
        summaries = summarize_patients(patients, seed=SEED)
        summaries.insert(0, "analysis", label)
        sensitive.append(summaries)
        if label == "main":
            patients.to_csv(effects / "patient_contrasts.csv", index=False)
            means.to_csv(effects / "patient_tissue_means.csv", index=False)
            lesions.to_csv(effects / "lesion_means.csv", index=False)
            table[SCORE_COLUMNS].to_csv(effects / "sample_program_scores.csv", index=False)
            summaries.drop(columns="analysis").to_csv(effects / "effects_summary.csv", index=False)
            leave_one_out(patients).to_csv(effects / "leave_one_patient_out.csv", index=False)
        else:
            patients.to_csv(effects / f"patients_{label}.csv", index=False)
    pd.concat(sensitive, ignore_index=True).to_csv(effects / "sensitivity_summary.csv", index=False)

    # Focal fine-subtype contrasts at both >=50 (primary) and >=20 (recovering smaller fine populations e.g. Monocyte)
    focal_summaries, focal_patients = [], []
    for min_cells in (50, 20):
        for branch_name, mask_col, sub_col, mf_col in (
            ("clean_author", "in_clean_author", "author_subtype", None),
            ("rule_annotation_m1", "in_rule_m1", "rule_subtype_m1", "rule_myeloid_fine_m1"),
        ):
            b_cells = cells.loc[cells[mask_col]]
            for (lin, prog), _ in programs.items():
                parent_lin = PARENT_LINEAGE_MAP.get(lin, lin)
                sub_targets = FOCAL_SUBTYPES_BY_PARENT.get(parent_lin, ())
                p_cells = b_cells.loc[b_cells["lineage"].eq(lin) & b_cells["program"].eq(prog)]
                if p_cells.empty or not sub_targets:
                    continue
                for st in sub_targets:
                    if st == "monocyte/macrophage":
                        smask = p_cells[sub_col].isin(["monocyte/macrophage", "Macrophage", "Monocyte"])
                    elif st in ("Macrophage", "Monocyte") and mf_col is not None:
                        smask = p_cells[sub_col].eq(st) | p_cells[mf_col].eq(st)
                    else:
                        smask = p_cells[sub_col].eq(st)
                    s_agg = p_cells.loc[smask].groupby("sample_uid")["score"].agg(["mean", "size"]).reset_index()
                    s_tab = pd.DataFrame({"sample_uid": manifest["sample_uid"], "lineage": st, "program": prog}).merge(
                        s_agg, on="sample_uid", how="left"
                    )
                    s_tab["n_cells"] = s_tab["size"].fillna(0).astype(int)
                    s_tab["score"] = s_tab["mean"].where(s_tab["n_cells"] >= min_cells)
                    p_df, _, _, _ = calculate_patient_contrasts(manifest, s_tab[SCORE_COLUMNS])
                    s_df = summarize_patients(p_df, seed=SEED)
                    p_df.insert(0, "min_cells_threshold", min_cells)
                    p_df.insert(0, "parent_candidate_lineage", lin)
                    p_df.insert(0, "annotation_branch", branch_name)
                    s_df.insert(0, "min_cells_threshold", min_cells)
                    s_df.insert(0, "parent_candidate_lineage", lin)
                    s_df.insert(0, "annotation_branch", branch_name)
                    focal_patients.append(p_df)
                    focal_summaries.append(s_df)
    if focal_summaries:
        pd.concat(focal_summaries, ignore_index=True).to_csv(effects / "focal_subtype_effects_summary.csv", index=False)
        pd.concat(focal_patients, ignore_index=True).to_csv(effects / "focal_subtype_patient_contrasts.csv", index=False)

    # Dual immune denominator composition contrasts (total_resolved_immune incl. ILC3/Mast + core_three_lineage_immune)
    if not comp_df.empty:
        comp_frames = []
        for denom_name, prop_col, n_col in (
            ("author_total_resolved_immune", "prop_author_total_resolved_immune", "n_author_clean_cells"),
            ("author_core_three_immune", "prop_author_core_three_immune", "n_author_clean_cells"),
            ("rule_m1_total_resolved_immune", "prop_rule_total_resolved_immune", "n_rule_m1_cells"),
            ("rule_m1_core_three_immune", "prop_rule_core_three_immune", "n_rule_m1_cells"),
        ):
            sub_c = comp_df.loc[comp_df[prop_col].notna(), ["sample_uid", "lineage", prop_col, n_col]].copy()
            if sub_c.empty:
                continue
            sub_c["program"] = denom_name
            sub_c = sub_c.rename(columns={prop_col: "score", n_col: "n_cells"})
            c_pat, _, _, _ = calculate_patient_contrasts(manifest, sub_c[SCORE_COLUMNS])
            comp_frames.append(c_pat)
        if comp_frames:
            pd.concat(comp_frames, ignore_index=True).to_csv(effects / "patient_lineage_composition.csv", index=False)

    def strata(frame: pd.DataFrame, labels, filename: str):
        if frame.empty:
            return
        work = frame.assign(subtype=labels)
        sub = work.groupby(["sample_uid", "lineage", "program", "subtype"]).score.agg(["size", "mean"]).reset_index()
        tot = work.groupby(["sample_uid", "lineage", "program"]).size().rename("n_cells").reset_index()
        sub = sub.rename(columns={"size": "n_subtype", "mean": "mean_score"}).merge(
            tot, on=["sample_uid", "lineage", "program"]
        )
        sub.to_csv(effects / (filename + "_sample.csv"), index=False)
        decompose(manifest, sub).to_csv(effects / (filename + ".csv"), index=False)

    strata(clean_cells, clean_cells.author_subtype, "patient_subtype_decomposition")
    b_clean = clean_cells.loc[clean_cells.author_lineage.eq("B/Plasma")]
    if not b_clean.empty:
        strata(
            b_clean,
            np.where(b_clean.author_subtype.eq("plasma cell"), "Plasma", "B"),
            "patient_B_plasma_decomposition",
        )
    if not clean_cells.empty:
        bins = np.digitize(clean_cells.n_genes_common.to_numpy(), [500, 1000, 2000])
        strata(clean_cells, bins.astype(str), "patient_detection_depth_decomposition")


def verify_freeze_integrity(freeze_dir: Path, freeze_commit: str | None = None) -> int:
    """Verify SHA256 integrity of `freeze_dir` files (resolving relative paths strictly)."""
    freeze_resolved = Path(freeze_dir).resolve()
    root_resolved = ROOT.resolve()
    prog_scores_copy = freeze_resolved / "program_scores.py"
    if not prog_scores_copy.is_file():
        raise FileNotFoundError(f"Missing frozen program_scores.py in {freeze_resolved}")
    if sha(prog_scores_copy) != sha(root_resolved / "src/crlm/program_scores.py"):
        raise ValueError("Current program_scores.py differs from the frozen M5 copy")

    checked = 0
    sha_manifest = freeze_resolved / "frozen_inputs_sha256.txt"
    if sha_manifest.is_file():
        for line in sha_manifest.read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if not line:
                continue
            digest_hex, rel_path = line.split(None, 1)
            rel_clean = rel_path.strip()
            if rel_clean == "artifacts/m5" or rel_clean.startswith("artifacts/m5/"):
                target = (freeze_resolved / Path(rel_clean).relative_to("artifacts/m5")).resolve()
                in_freeze = True
            else:
                target = (root_resolved / rel_clean).resolve()
                in_freeze = target.is_relative_to(freeze_resolved)
            if in_freeze:
                if not target.is_file():
                    raise FileNotFoundError(f"Frozen M5 file declared in manifest is missing: {rel_clean} -> {target}")
                if sha(target) != digest_hex:
                    raise ValueError(f"Frozen M5 checksum mismatch for {rel_clean}")
                checked += 1
        if checked == 0:
            raise ValueError(f"Zero freeze_dir files matched in {sha_manifest}")

    if freeze_commit and freeze_commit != "SHA_LOCKED":
        frozen_files = [
            "artifacts/m5/frozen_candidate_registry.csv",
            "artifacts/m5/gene_sets/program_genes.csv",
            "artifacts/m5/registered_tests.csv",
            "artifacts/m5/program_scores.py",
            "artifacts/m5/validation_rules.md",
        ]
        for opt in ("artifacts/m5/gene_sets/model_program_genes.csv", "artifacts/m5/gene_sets/matched_controls.csv"):
            if (root_resolved / opt).is_file():
                frozen_files.append(opt)
        for name in frozen_files:
            frozen_bytes = subprocess.check_output(["git", "show", f"{freeze_commit}:{name}"], cwd=root_resolved)
            if frozen_bytes != (root_resolved / name).read_bytes():
                raise ValueError(f"Working and frozen git commit {freeze_commit} disagree for {name}")
    return checked


def run_selfcheck() -> None:
    """Pure in-memory synthetic & contract verification of V1.1 rules (zero project disk writes)."""
    import os
    import tempfile
    from crlm.annotation_m1 import AUXILIARY_MARKERS, PANELS

    rng = np.random.default_rng(SEED)
    full_shared = list(load_shared_nuclear_universe("245_164_shared_nuclear"))
    full_shared_set = set(full_shared)
    # 1. Verify axis-preserving canonicalization keeps PRR34-AS1 / LINC01358 intact
    for tricky in ("PRR34-AS1", "LINC01358", "APOE", "APOC1"):
        if tricky in full_shared_set:
            assert canonicalize_preserving_axis([tricky], full_shared_set) == [tricky]
    assert "APOE" not in LIVER_RBC and "APOC1" not in LIVER_RBC

    # 2. Verify relative freeze_dir resolution, missing-file guard, and missing-control ValueError guard
    with tempfile.TemporaryDirectory() as tmp:
        tdir = Path(tmp)
        (tdir / "program_scores.py").write_bytes((ROOT / "src/crlm/program_scores.py").read_bytes())
        (tdir / "gene_sets").mkdir()
        mc_sample = pd.DataFrame(
            [{"lineage": "Mono/Mac/DC", "program": "SYN_SAM", "target_gene": "APOE", "replicate": 0, "control_gene": "A1BG"}]
        )
        mc_file = tdir / "gene_sets/matched_controls.csv"
        mc_sample.to_csv(mc_file, index=False)
        (tdir / "frozen_inputs_sha256.txt").write_text(
            f"{sha(mc_file)}  artifacts/m5/gene_sets/matched_controls.csv\n", encoding="utf-8"
        )
        rel_or_abs = Path(os.path.relpath(tdir, Path.cwd()))
        assert verify_freeze_integrity(rel_or_abs, "SHA_LOCKED") == 1
        # Missing control for a second candidate MUST raise ValueError (no silent fallback)
        reg_two = pd.DataFrame(
            [
                {"lineage": "Mono/Mac/DC", "program": "SYN_SAM", "scoring_method": "ucell_debiased"},
                {"lineage": "T/NK", "program": "MISSING_PROG", "scoring_method": "ucell_debiased"},
            ]
        )
        raised_val = False
        try:
            load_frozen_matched_controls(
                tdir, [("Mono/Mac/DC", "SYN_SAM"), ("T/NK", "MISSING_PROG")], reg_two, shared_axis=full_shared
            )
        except ValueError:
            raised_val = True
        assert raised_val, "Expected ValueError when a candidate program is missing from matched_controls.csv"

    axis = full_shared[:200]
    marker_inject = [
        "CD3D", "CD3E", "CD4", "CD8A", "TRAC", "IL7R",
        "MS4A1", "CD79A", "MZB1", "JCHAIN", "SDC1", "TNFRSF17",
        "C1QA", "C1QB", "CD14", "CD68", "CSF1R", "SPP1", "MARCO", "PTPRC", "KIT", "APOE", "APOC1",
    ]
    for i, mg in enumerate(marker_inject):
        if mg not in axis:
            axis[i] = mg
    all_panel_genes = {g for _, a, s in PANELS.values() for g in (a + " " + s).split()} | set(AUXILIARY_MARKERS)
    non_panel = [g for g in axis if g not in all_panel_genes]
    prog_genes = list(dict.fromkeys(["APOE", "APOC1", *non_panel]))[:10]
    ctrl_pool = [g for g in non_panel if g not in prog_genes][:100]
    controls_map = {
        ("Mono/Mac/DC", "SYN_SAM"): {
            g: [ctrl_pool[(idx * 10 + r) % len(ctrl_pool)] for r in range(10)]
            for idx, g in enumerate(prog_genes)
        },
        ("CD4 T cell", "SYN_TH17"): {
            g: [ctrl_pool[(idx * 10 + r + 5) % len(ctrl_pool)] for r in range(10)]
            for idx, g in enumerate(prog_genes)
        },
    }
    programs = {
        ("Mono/Mac/DC", "SYN_SAM"): prog_genes,
        ("CD4 T cell", "SYN_TH17"): prog_genes,
    }

    h5ads = {}
    meta_rows = []
    col_idx = {g: i for i, g in enumerate(axis)}
    for code, tissue in ORGANS.items():
        n_c = 170  # 60 Macrophage + 60 CD4 T + 10 ILC3 + 40 B-CD3E doublet
        mat = sparse.random(n_c, len(axis), density=0.15, format="csr", random_state=rng)
        mat.data = np.log1p(rng.integers(1, 20, size=mat.nnz) * 10.0)
        mat_lil = mat.tolil()
        for g in all_panel_genes:
            if g in col_idx:
                mat_lil[:, col_idx[g]] = 0.0
        for r in range(60):
            for g in ("C1QA", "C1QB", "CD68", "CSF1R", "APOE", "APOC1"):
                mat_lil[r, col_idx[g]] = 3.5
            if tissue == "LM":
                for g in prog_genes[:5]:
                    mat_lil[r, col_idx[g]] = 4.5
        for r in range(60, 120):
            for g in ("CD3D", "CD3E", "CD4", "TRAC"):
                mat_lil[r, col_idx[g]] = 3.5
        for r in range(120, 130):
            mat_lil[r, col_idx["KIT"]] = 3.5
        for r in range(130, 170):
            for g in ("CD3D", "CD3E", "MS4A1", "CD79A"):
                mat_lil[r, col_idx[g]] = 3.5
        mat_csr = mat_lil.tocsr()

        bcs = [f"{code}_cell_{i}" for i in range(n_c)]
        majors = ["Macrophage"] * 60 + ["CD4 T"] * 60 + ["ILC3"] * 10 + ["Bcell"] * 40
        subs = (
            ["hC54_Mph-SPP1"] * 60
            + ["hC23_CD4_Th17-IL23R"] * 60
            + ["hC58_ILC3-KIT"] * 10
            + ["hC75_B-CD3E"] * 40
        )
        globs = ["Myeloid"] * 60 + ["Lymphoid-T-NK"] * 60 + ["Myeloid"] * 10 + ["Lymphoid-B"] * 40
        obs_df = pd.DataFrame(
            {
                "source_barcode": bcs,
                "donor_uid": f"{GSE}:p01",
                "sample_uid": f"{GSE}:p01:{code}",
                "lesion_uid": f"{GSE}:p01:{code}",
                "patient": "p01",
                "tissue": tissue,
                "author_tissue": code,
                "sample": code,
                "ID": code,
                "celltype_major": majors,
                "celltype_sub": subs,
                "celltype_global": globs,
            },
            index=pd.Index([f"{GSE}:{b}" for b in bcs], name="cell_uid"),
        )
        obs_df = map_author_annotations(obs_df)
        meta_rows.append(obs_df)
        var_df = pd.DataFrame({"gene_symbol": axis}, index=pd.Index(axis, name="source_gene"))
        h5ads[code] = ad.AnnData(X=mat_csr, obs=obs_df, var=var_df)

    meta_all = pd.concat(meta_rows, ignore_index=True)
    assert int(meta_all["is_cd3e_doublet"].sum()) == 160
    cells_df, conc_df, comp_df = score_anndata_dict(
        h5ads, meta_all, axis, programs, controls_map, set(programs), requires_debiased=True
    )
    assert not cells_df.empty and not conc_df.empty and not comp_df.empty
    assert not cells_df.loc[cells_df["in_clean_author"], "is_cd3e_doublet"].any()
    assert int(cells_df.loc[cells_df["lineage"].eq("CD4 T cell") & cells_df["in_clean_author"]].shape[0]) == 240
    assert int(cells_df.loc[cells_df["lineage"].eq("CD4 T cell") & cells_df["in_rule_m1"]].shape[0]) == 240
    # Verify dual immune proportions sum to 1.0 EVEN WHEN ILC3 cells are present!
    for _, grp in comp_df.groupby("sample_uid"):
        assert np.isclose(grp["prop_author_total_resolved_immune"].dropna().sum(), 1.0)
        assert np.isclose(grp["prop_author_core_three_immune"].dropna().sum(), 1.0)
        assert np.isclose(grp["prop_rule_total_resolved_immune"].dropna().sum(), 1.0)
        assert np.isclose(grp["prop_rule_core_three_immune"].dropna().sum(), 1.0)
    print(
        json.dumps(
            {
                "selfcheck": "passed",
                "script": "validate_gse164522.py",
                "scored_rows": int(len(cells_df)),
                "cd3e_doublets_stripped": int(meta_all["is_cd3e_doublet"].sum()),
                "ilc3_composition_sum_verified": True,
            }
        ),
        flush=True,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("step", choices=("fetch", "identity", "convert", "score", "selfcheck"))
    parser.add_argument("--tissue-code", choices=tuple(ORGANS), nargs="+")
    parser.add_argument("--freeze-commit", help="Git revision (or SHA_LOCKED) containing the frozen M5 registry and tests")
    parser.add_argument("--freeze-dir", type=Path, default=FREEZE)
    parser.add_argument("--out-dir", type=Path, default=OUT)
    parser.add_argument("--objects-dir", type=Path, default=OBJECTS)
    args = parser.parse_args()

    if args.step == "selfcheck":
        run_selfcheck()
        return

    out_dir: Path = args.out_dir.resolve()
    freeze_dir: Path = args.freeze_dir.resolve()
    objects_dir: Path = args.objects_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    ad.settings.allow_write_nullable_strings = True
    records = sources()
    if args.step == "fetch":
        fetch(records, out_dir=out_dir)
        return
    frame = metadata(records)
    manifest = write_identity(frame, records, out_dir=out_dir)
    if args.step == "identity":
        print(
            f"Author identity: {len(frame)} cells ({int(frame.clean_author_qc.sum())} clean, "
            f"{int(frame.is_cd3e_doublet.sum())} CD3E+ doublets), "
            f"{len(manifest)} specimens, {manifest.donor_uid.nunique()} donors"
        )
        return
    if not args.freeze_commit:
        parser.error("--freeze-commit is required for expression processing (pass git commit or SHA_LOCKED)")
    verify_freeze_integrity(freeze_dir, args.freeze_commit)
    objects_dir.mkdir(parents=True, exist_ok=True)
    if args.step == "convert":
        for code in args.tissue_code or ORGANS:
            convert(code, frame, records[code], objects_dir=objects_dir, out_dir=out_dir)
        return

    start = time.perf_counter()
    axis, programs, controls_map, candidate_keys, eval_keys, model_keys, req_deb = frozen_definitions(
        freeze_dir=freeze_dir, objects_dir=objects_dir, out_dir=out_dir
    )
    cells, comp_df, h5ads = measure(
        axis,
        programs,
        controls_map,
        eval_keys,
        frame,
        requires_debiased=req_deb,
        objects_dir=objects_dir,
        out_dir=out_dir,
    )
    cells.to_csv(objects_dir / "frozen_cell_scores.csv.gz", index=False, compression={"method": "gzip", "mtime": 0})
    model_inputs(
        cells,
        model_keys,
        frame.groupby("sample_uid").size(),
        h5ads_by_code=h5ads,
        meta_frame=frame,
        axis=axis,
        freeze_dir=freeze_dir,
        out_dir=out_dir,
    )
    is_eval = pd.MultiIndex.from_frame(cells[["lineage", "program"]]).isin(list(eval_keys))
    summarize(
        manifest,
        frame,
        cells.loc[is_eval],
        {k: v for k, v in programs.items() if k in eval_keys},
        comp_df,
        out_dir=out_dir,
    )
    tests = pd.read_csv(freeze_dir / "registered_tests.csv", keep_default_na=False)
    plan = tests.loc[tests.cohort.eq(GSE)]
    summary = pd.read_csv(out_dir / "patient_effects/effects_summary.csv")
    adjusted, extra, families, audit = adjust_registered(plan, summary)
    target = out_dir / "registered_tests"
    target.mkdir(exist_ok=True)
    for name, table in (
        ("registered_results", adjusted),
        ("family_summary", families),
        ("unplanned_results", extra),
    ):
        table.to_csv(target / f"{name}.csv", index=False)
    (target / "audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    profile = dict(
        completed_utc=now(),
        command=shlex.join([sys.executable, *sys.argv]),
        freeze_commit=args.freeze_commit,
        seconds=time.perf_counter() - start,
        peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024,
        code_sha256=sha(Path(__file__)),
        rank_code_sha256=sha(ROOT / "src/crlm/program_scores.py"),
        hgnc_code_sha256=sha(ROOT / "src/crlm/hgnc_alias_map.py"),
        annotation_code_sha256=sha(ROOT / "src/crlm/annotation_m1.py"),
        n_rank_genes=len(axis),
        datasets=[json.loads((objects_dir / f"{c}.json").read_text()) for c in ORGANS if (objects_dir / f"{c}.json").is_file()],
    )
    (out_dir / "run_profile.json").write_text(json.dumps(profile, indent=2) + "\n")


if __name__ == "__main__":
    main()
