#!/usr/bin/env python3
"""Validate sample identities, roles and declared expression capabilities.

This checks declarations for consistency, not the truth of an author's identity
mapping or the contents/normalization of an expression matrix.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from _common import bool_column, cli_error, read_csv, require_columns, require_nonempty, sha256_file, write_json

REQUIRED = ["dataset", "donor_uid", "sample_uid", "tissue", "lesion_uid", "role", "input_scale", "patient_mapping_verified"]
ROLES = {"development", "external", "exploratory", "excluded"}
SCALES = {"raw_counts", "normalized", "log_normalized", "corrected", "unknown"}
TISSUES = {"PT", "LM", "NC", "NL"}


def validate_frame(frame):
    require_columns(frame, REQUIRED, "sample manifest")
    require_nonempty(frame, REQUIRED, "sample manifest")
    result = frame.copy()
    for field, allowed in [("role", ROLES), ("input_scale", SCALES), ("tissue", TISSUES)]:
        invalid = sorted(set(result[field]) - allowed)
        if invalid:
            raise ValueError(f"{field}: unsupported values {invalid}; allowed={sorted(allowed)}")
    result["patient_mapping_verified"] = bool_column(result["patient_mapping_verified"], "patient_mapping_verified")
    if result["sample_uid"].duplicated().any():
        duplicated = result.loc[result["sample_uid"].duplicated(False), "sample_uid"].unique().tolist()
        raise ValueError(f"sample_uid must identify exactly one manifest row; duplicate/conflicting samples: {duplicated[:10]}")
    for donor, group in result.groupby("donor_uid", sort=True):
        if group["role"].nunique() != 1:
            raise ValueError(f"donor {donor} crosses roles {sorted(group['role'].unique())}; all tissues/lesions must share one role")
        if group["patient_mapping_verified"].nunique() != 1:
            raise ValueError(f"donor {donor}: inconsistent identity verification")
    # A lesion cannot change donor or organ/tissue across technical samples.
    for lesion, group in result.groupby("lesion_uid", sort=True):
        if group["donor_uid"].nunique() != 1 or group["tissue"].nunique() != 1:
            raise ValueError(f"lesion_uid {lesion} maps to different donors/tissues; use globally unique lesion IDs")
    if "source_audit_scale" in result:
        declared = result["source_audit_scale"]
        invalid = set(declared) - SCALES - {""}
        if invalid:
            raise ValueError(f"source_audit_scale: unsupported values {sorted(invalid)}")
        contradiction = (result.input_scale == "raw_counts") & declared.isin(SCALES - {"raw_counts", "unknown"})
        if contradiction.any():
            raise ValueError("Processed source_audit_scale cannot be declared raw_counts")
    if "can_nb_pseudobulk" in result:
        claimed = bool_column(result["can_nb_pseudobulk"], "can_nb_pseudobulk")
        if (claimed & result.input_scale.ne("raw_counts")).any():
            raise ValueError("Processed/unknown expression cannot claim can_nb_pseudobulk=true")
    # Never infer raw counts from float/int dtype or from a filename.
    provenance = result.get("counts_provenance", pd.Series("", index=result.index)).astype(str)
    verified_values = bool_column(result["counts_values_verified"], "counts_values_verified") if "counts_values_verified" in result else pd.Series(False, index=result.index)
    can_nb = result.input_scale.eq("raw_counts") & provenance.ne("") & verified_values
    if "can_nb_pseudobulk" in result and (claimed & ~can_nb).any():
        raise ValueError("A raw-count NB claim requires counts_provenance and counts_values_verified=true; no counts are reconstructed here")
    result["can_nb_pseudobulk"] = can_nb
    result["eligible_for_paired_analysis"] = result.patient_mapping_verified & result.role.ne("excluded")
    warnings = []
    unverified = sorted(result.loc[~result.patient_mapping_verified, "donor_uid"].unique())
    if unverified:
        warnings.append(f"Unverified donor codes are ineligible for paired conclusions: {unverified}")
    unsupported_raw = sorted(result.loc[result.input_scale.eq("raw_counts") & ~can_nb, "sample_uid"].unique())
    if unsupported_raw:
        warnings.append(f"Raw-count labels without recorded verification are NOT NB-eligible: {unsupported_raw}")
    multisource = result.groupby("donor_uid").dataset.nunique()
    repeated = multisource[multisource > 1].index.tolist()
    if repeated:
        warnings.append(f"Donors appear in multiple datasets; deduplicate and do not claim dataset independence: {repeated}")
    donors = []
    for donor, group in result.groupby("donor_uid", sort=True):
        donors.append({
            "donor_uid": donor,
            "datasets": json.dumps(sorted(group.dataset.unique()), ensure_ascii=False),
            "role": group.role.iloc[0],
            "patient_mapping_verified": bool(group.patient_mapping_verified.iloc[0]),
            "sample_uids": json.dumps(sorted(group.sample_uid), ensure_ascii=False),
            "tissues": json.dumps(sorted(group.tissue.unique())),
            "n_samples": len(group),
        })
    return result, pd.DataFrame(donors), warnings


def run(manifest, outdir):
    samples, donors, warnings = validate_frame(read_csv(manifest))
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    samples.to_csv(outdir / "sample_manifest_validated.csv", index=False)
    donors.to_csv(outdir / "donor_manifest_validated.csv", index=False)
    report = {
        "status": "passed_with_warnings" if warnings else "passed",
        "source_sha256": sha256_file(manifest),
        "n_samples": len(samples), "n_donor_codes": len(donors),
        "n_verified_donors": int(donors.patient_mapping_verified.sum()),
        "warnings": warnings,
        "limitation": "Checks declared provenance and internal consistency; cannot establish true patient identity or raw-count provenance without source audit.",
    }
    write_json(outdir / "manifest_validation.json", report)
    print(json.dumps(report, ensure_ascii=False))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()
    cli_error(parser, lambda: run(args.manifest, args.outdir))


if __name__ == "__main__":
    main()
