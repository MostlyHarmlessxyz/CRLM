#!/usr/bin/env python3
"""Match a frozen test registry to paired summaries and apply family-wise Holm.

Missing planned tests retain raw_p=NA and use 1 only in the adjustment vector.
Unplanned results are exported separately (or rejected), never added to a family.
Statistical passage is not a declaration of complete biological replication.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binomtest

from _common import bool_column, cli_error, read_csv, require_columns, require_nonempty, sha256_file, write_json

PLAN_REQUIRED = ["test_id", "family", "cohort", "role", "lineage", "program", "contrast",
                 "expected_direction", "direction_rule", "alpha", "test_method", "adjustment",
                 "analysis_tier", "freeze_status", "freeze_timestamp_utc", "freeze_version", "synthetic_example"]
KEY = ["cohort", "role", "lineage", "program", "contrast"]
RESULT_KEY = ["dataset", "role", "lineage", "program", "contrast"]
RESULT_REQUIRED = RESULT_KEY + ["n_patients", "status", "mean_difference", "n_positive", "n_negative", "n_ties", "n_sign_test", "sign_p_two_sided"]
CONTRASTS = {"D_T", "D_N", "D_T_Q", "D_N_Q", "I"}


def number(value, field, missing=False, integer=False):
    text = str(value).strip()
    if text.lower() in {"", "na", "nan", "none", "null"}:
        if missing:
            return None
        raise ValueError(f"{field}: missing numeric value")
    try:
        result = float(text)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{field}: invalid numeric value {value!r}") from error
    if not np.isfinite(result) or (integer and (result < 0 or result != int(result))):
        raise ValueError(f"{field}: require finite {'nonnegative integer' if integer else 'number'}")
    return int(result) if integer else result


def validate_plan(frame):
    require_columns(frame, PLAN_REQUIRED, "frozen test registry")
    require_nonempty(frame, PLAN_REQUIRED, "frozen test registry")
    if frame[PLAN_REQUIRED].isna().any().any():
        raise ValueError("Frozen registry required fields may not contain NA")
    plan = frame.copy()
    if plan.test_id.duplicated().any() or plan.duplicated(KEY).any():
        raise ValueError("test_id and cohort/role/lineage/program/contrast keys must each be unique; no duplicate hypothesis across families")
    enums = {"contrast": CONTRASTS, "role": {"development", "external", "exploratory"},
             "expected_direction": {"positive", "negative"}, "direction_rule": {"mean_and_sign_majority"},
             "test_method": {"exact_sign_two_sided"}, "adjustment": {"Holm"},
             "analysis_tier": {"primary", "secondary", "background", "diagnostic"}, "freeze_status": {"frozen"}}
    for field, values in enums.items():
        invalid = set(plan[field]) - values
        if invalid:
            raise ValueError(f"{field}: unsupported/unfrozen values {sorted(invalid)}; supported={sorted(values)}")
    if ((plan.analysis_tier != "diagnostic") & (plan.role != "external")).any():
        raise ValueError("Primary, secondary and background external-evidence families require role=external; development results may only be diagnostic here")
    plan["alpha"] = [number(x, "alpha") for x in plan.alpha]
    if not plan.alpha.between(0, 1, inclusive="neither").all():
        raise ValueError("Require 0 < alpha < 1")
    plan["synthetic_example"] = bool_column(plan.synthetic_example, "synthetic_example")
    if plan.synthetic_example.nunique() != 1:
        raise ValueError("Do not mix synthetic and real planned tests")
    for stamp in plan.freeze_timestamp_utc:
        try:
            parsed = datetime.fromisoformat(stamp.replace("Z", "+00:00"))
        except ValueError as error:
            raise ValueError("freeze_timestamp_utc must be an ISO-8601 UTC timestamp") from error
        if parsed.tzinfo is None or parsed.utcoffset() != timedelta(0):
            raise ValueError("freeze_timestamp_utc must explicitly use UTC (Z or +00:00)")
    for family, group in plan.groupby("family", sort=True):
        for field in ["alpha", "adjustment", "test_method", "analysis_tier", "freeze_version", "freeze_timestamp_utc"]:
            if group[field].nunique() != 1:
                raise ValueError(f"Family {family}: inconsistent {field}")
    return plan


def validate_results(frame):
    missing = set(RESULT_REQUIRED) - set(frame.columns)
    if missing:
        raise ValueError(f"paired effect summaries: missing columns {sorted(missing)}")
    # An empty, correctly headed summary is allowed: every planned test stays missing.
    if frame.empty:
        return frame.copy()
    require_nonempty(frame, RESULT_KEY + ["status"], "paired effect summaries")
    if frame[RESULT_KEY + ["status"]].isna().any().any():
        raise ValueError("Summary keys/status may not contain NA")
    if frame.duplicated(RESULT_KEY).any():
        raise ValueError("Duplicate summary keys: matching must be strictly one-to-one, including unplanned rows")
    if not set(frame.contrast).issubset(CONTRASTS):
        raise ValueError("Unknown contrast in paired effect summaries")
    result = frame.copy()
    counts = ["n_patients", "n_positive", "n_negative", "n_ties", "n_sign_test"]
    for field in counts:
        result[field] = [number(x, field, integer=True) for x in result[field]]
    result["mean_difference"] = [number(x, "mean_difference", missing=True) for x in result.mean_difference]
    result["sign_p_two_sided"] = [number(x, "sign_p_two_sided", missing=True) for x in result.sign_p_two_sided]
    for row in result.itertuples(index=False):
        if row.n_positive + row.n_negative != row.n_sign_test or row.n_sign_test + row.n_ties != row.n_patients:
            raise ValueError("Direction counts, ties, n_sign_test and n_patients disagree")
        p_missing = pd.isna(row.sign_p_two_sided)
        mean_missing = pd.isna(row.mean_difference)
        if row.n_patients == 0:
            if not p_missing or not mean_missing:
                raise ValueError("N=0 requires missing raw p and mean; never manufacture an observed zero")
            continue
        if p_missing or mean_missing or not 0 <= row.sign_p_two_sided <= 1:
            raise ValueError("Observed contrasts require a finite mean and a raw p in [0,1]")
        expected = binomtest(row.n_positive, row.n_sign_test, 0.5, alternative="two-sided").pvalue if row.n_sign_test else 1.0
        if not np.isclose(row.sign_p_two_sided, expected, rtol=1e-9, atol=1e-12):
            raise ValueError("Input raw p does not match the registered exact two-sided sign test and direction counts; retain full precision")
    return result


def holm_adjust(p_values):
    p = np.asarray(p_values, dtype=float)
    if p.ndim != 1 or len(p) == 0 or not np.isfinite(p).all() or ((p < 0) | (p > 1)).any():
        raise ValueError("Holm requires a nonempty finite p vector in [0,1]")
    order = np.argsort(p, kind="stable")
    adjusted = np.minimum(1.0, np.maximum.accumulate((len(p) - np.arange(len(p))) * p[order]))
    output = np.empty_like(p)
    output[order] = adjusted
    return output


def direction_fields(row, expected_direction):
    if row is None or row["n_patients"] == 0:
        return {"direction_evaluable": False, "direction_status": "no_observed_contrast",
                "mean_direction": None, "nonzero_majority_direction": None,
                "mean_direction_matches": None, "majority_direction_matches": None, "direction_matches": None}
    if row["n_sign_test"] == 0:
        return {"direction_evaluable": False, "direction_status": "all_ties_no_direction_information",
                "mean_direction": None, "nonzero_majority_direction": None,
                "mean_direction_matches": None, "majority_direction_matches": None, "direction_matches": None}
    mean_direction = "positive" if row["mean_difference"] > 0 else "negative" if row["mean_difference"] < 0 else "zero"
    majority = "positive" if row["n_positive"] > row["n_negative"] else "negative" if row["n_positive"] < row["n_negative"] else "tied_majority"
    mean_match, majority_match = mean_direction == expected_direction, majority == expected_direction
    return {"direction_evaluable": True,
            "direction_status": "matches_registered_direction" if mean_match and majority_match else "opposite_zero_or_mean_majority_disagreement",
            "mean_direction": mean_direction, "nonzero_majority_direction": majority,
            "mean_direction_matches": mean_match, "majority_direction_matches": majority_match,
            "direction_matches": mean_match and majority_match}


def adjust_registered(plan_frame, summary_frame, unplanned_policy="separate"):
    if unplanned_policy not in {"separate", "error"}:
        raise ValueError("unplanned_policy must be separate or error")
    plan, summaries = validate_plan(plan_frame), validate_results(summary_frame)
    lookup = {tuple(row[field] for field in RESULT_KEY): row.to_dict() for _, row in summaries.iterrows()}
    planned_keys = {tuple(row[field] for field in KEY) for _, row in plan.iterrows()}
    exploratory_mask = [tuple(row[field] for field in RESULT_KEY) not in planned_keys for _, row in summaries.iterrows()]
    exploratory = summaries.loc[exploratory_mask].copy() if len(summaries) else summaries.copy()
    if len(exploratory) and unplanned_policy == "error":
        raise ValueError(f"{len(exploratory)} unplanned summary rows; they cannot enter the frozen families")
    exploratory["registry_exclusion_reason"] = "not_in_frozen_plan"
    exploratory["registry_analysis_tier"] = "unplanned_exploratory"
    rows = []
    for _, planned in plan.iterrows():
        key = tuple(planned[field] for field in KEY)
        observed = lookup.get(key)
        raw_p = None if observed is None or pd.isna(observed["sign_p_two_sided"]) else float(observed["sign_p_two_sided"])
        output = planned.to_dict()
        output.update(matched_summary_row=observed is not None,
                      result_status="missing_summary_row" if observed is None else observed["status"],
                      raw_p=raw_p, adjustment_only_p=1.0 if raw_p is None else raw_p,
                      planned_test_missing=raw_p is None)
        for field in ["n_patients", "n_positive", "n_negative", "n_ties", "n_sign_test", "mean_difference"]:
            output[field] = None if observed is None else observed[field]
        output.update(direction_fields(observed, planned.expected_direction))
        rows.append(output)
    adjusted = pd.DataFrame(rows)
    family_rows = []
    for family, group in adjusted.groupby("family", sort=True):
        idx = group.index
        adjusted.loc[idx, "holm_adjusted_p"] = holm_adjust(group.adjustment_only_p.to_numpy())
        adjusted.loc[idx, "family_size_planned"] = len(group)
        family_rows.append({"family": family, "analysis_tier": group.analysis_tier.iloc[0], "alpha": float(group.alpha.iloc[0]),
                            "n_planned_tests": len(group), "n_raw_p_available": int(group.raw_p.notna().sum()),
                            "n_missing_planned_tests": int(group.planned_test_missing.sum()),
                            "n_matched_summary_rows": int(group.matched_summary_row.sum()),
                            "n_all_ties": int(group.direction_status.eq("all_ties_no_direction_information").sum())})
    adjusted["family_size_planned"] = adjusted.family_size_planned.astype(int)
    adjusted["holm_reject"] = adjusted.raw_p.notna() & (adjusted.holm_adjusted_p <= adjusted.alpha)
    adjusted["statistical_direction_criterion_met"] = adjusted.holm_reject & adjusted.direction_matches.eq(True)
    adjusted["biological_replication_status"] = "not_assessed_by_this_tool"
    audit = {"n_planned_tests": len(plan), "n_families": plan.family.nunique(), "n_unplanned_results": len(exploratory),
             "n_missing_planned_tests": int(adjusted.planned_test_missing.sum()), "unplanned_policy": unplanned_policy,
             "synthetic_example": bool(plan.synthetic_example.iloc[0]),
             "notes": ["Missing planned tests keep raw_p missing and use p=1 only for Holm adjustment; family size never shrinks.",
                       "Direction requires both mean difference and nonzero-patient sign majority to match the frozen expected direction; all ties are not evaluable.",
                       "Primary/secondary/background labels come from the frozen plan; separate families do not imply study-wide error control.",
                       "A Holm rejection plus direction match is a statistical criterion, not proof of full biological replication, independence or mechanism.",
                       "A frozen timestamp declaration is checked for format, not historical truth; use the external-access log and archived registry hash."]}
    return adjusted, exploratory, pd.DataFrame(family_rows), audit


def run(plan_path, summary_path, outdir, unplanned_policy="separate"):
    adjusted, exploratory, families, audit = adjust_registered(read_csv(plan_path), read_csv(summary_path), unplanned_policy)
    destination = Path(outdir)
    destination.mkdir(parents=True, exist_ok=True)
    adjusted.to_csv(destination / "registered_test_results.csv", index=False)
    exploratory.to_csv(destination / "unplanned_exploratory_results.csv", index=False)
    families.to_csv(destination / "registered_families.csv", index=False)
    audit.update(plan_sha256=sha256_file(plan_path), summaries_sha256=sha256_file(summary_path))
    write_json(destination / "registered_test_audit.json", audit)
    print(f"Adjusted {len(adjusted)} planned tests in {len(families)} frozen families; {len(exploratory)} unplanned results exported separately.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", required=True)
    parser.add_argument("--summaries", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--unplanned-policy", choices=["separate", "error"], default="separate")
    args = parser.parse_args()
    cli_error(parser, lambda: run(args.plan, args.summaries, args.outdir, args.unplanned_policy))


if __name__ == "__main__":
    main()
