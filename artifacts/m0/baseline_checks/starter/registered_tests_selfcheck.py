#!/usr/bin/env python3
"""Synthetic, independent checks of registry matching and multiplicity control."""
from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd

from _common import read_csv, write_json
from adjust_registered_tests import PLAN_REQUIRED, RESULT_REQUIRED, adjust_registered, holm_adjust
from paired_effects import summarize_difference


def synthetic_fixture():
    """Create artificial outcomes, never accessed from a biological cohort."""
    candidates = [("Myeloid", "fixed_demo"), ("T_NK", "toy_T_program"), ("B_Plasma", "toy_B_program")]
    definitions = [("PRIMARY_DT_164_V1", "primary", "SYN_GSE164522", "D_T"),
                   ("BACKGROUND_I_164_V1", "background", "SYN_GSE164522", "I"),
                   ("SECONDARY_DT_178_V1", "secondary", "SYN_GSE178318", "D_T")]
    plan_rows = []
    for family, tier, cohort, contrast in definitions:
        for index, (lineage, program) in enumerate(candidates, 1):
            plan_rows.append(dict(test_id=f"{family}_C{index}", family=family, cohort=cohort,
                                  role="external", lineage=lineage, program=program, contrast=contrast,
                                  expected_direction="positive", direction_rule="mean_and_sign_majority",
                                  alpha=0.05, test_method="exact_sign_two_sided", adjustment="Holm",
                                  analysis_tier=tier, freeze_status="frozen", freeze_timestamp_utc="2026-11-21T00:00:00Z",
                                  freeze_version="SYNTHETIC_TOY_V1", synthetic_example=True))
    plan = pd.DataFrame(plan_rows, columns=PLAN_REQUIRED)
    # None: absent summary row; []: submitted but no complete patient contrast.
    artificial_differences = [[1.0]*10, [-1.0]*10, None,
                              [1.0]*6, [0.0]*6, [],
                              [1.0]*6, [1.0]*6, [1.0, -1.0, 1.0, -1.0, 0.0, 0.0]]
    observed = []
    for row, values in zip(plan.to_dict("records"), artificial_differences):
        if values is None:
            continue
        summary = {"mean_difference": None, "sign_p_two_sided": None,
                   **summarize_difference(values, seed=17, n_bootstrap=200)}
        observed.append(dict(dataset=row["cohort"], role=row["role"], lineage=row["lineage"],
                             program=row["program"], contrast=row["contrast"], **summary))
    observed.append(dict(dataset="SYN_GSE164522", role="external", lineage="Other",
                         program="unplanned_program", contrast="D_T",
                         **summarize_difference([1.0]*20, seed=23, n_bootstrap=200)))
    return plan, pd.DataFrame(observed)


def main():
    checks = []

    def check(name, condition):
        if not bool(condition):
            raise AssertionError(name)
        checks.append({"name": name, "passed": True})

    def rejects(name, operation):
        try:
            operation()
        except ValueError:
            check(name, True)
        else:
            raise AssertionError(name + ": unsafe input accepted")

    plan, summary = synthetic_fixture()
    output, exploratory, families, audit = adjust_registered(plan, summary)
    lookup = output.set_index("test_id")
    check("hand_computed_Holm_with_missing_p_one", np.allclose(holm_adjust([0.03, 1.0, 0.01]), [0.06, 1.0, 0.03]))
    check("nine_planned_tests_three_families", len(output) == 9 and len(families) == 3 and output.family_size_planned.eq(3).all())
    first = lookup.loc["PRIMARY_DT_164_V1_C1"]
    absent = lookup.loc["PRIMARY_DT_164_V1_C3"]
    check("missing_planned_test_retains_NA_raw_p_and_adjustment_one", pd.isna(absent.raw_p) and absent.adjustment_only_p == 1 and absent.planned_test_missing)
    check("missing_planned_test_keeps_Holm_family_size", np.isclose(first.holm_adjusted_p, 3*0.001953125) and not np.isclose(first.holm_adjusted_p, 2*0.001953125))
    background = lookup.loc["BACKGROUND_I_164_V1_C1"]
    check("different_families_are_not_accidentally_pooled", np.isclose(background.holm_adjusted_p, 3*0.03125))
    opposite = lookup.loc["PRIMARY_DT_164_V1_C2"]
    check("significant_opposite_direction_fails_direction_criterion", opposite.holm_reject and not opposite.direction_matches and not opposite.statistical_direction_criterion_met)
    check("matching_direction_and_Holm_rejection_are_recorded", first.holm_reject and first.direction_matches and first.statistical_direction_criterion_met)
    all_ties = lookup.loc["BACKGROUND_I_164_V1_C2"]
    check("all_ties_have_p_one_and_no_evaluable_direction", all_ties.raw_p == 1 and not all_ties.direction_evaluable and pd.isna(all_ties.direction_matches) and not all_ties.statistical_direction_criterion_met)
    no_patients = lookup.loc["BACKGROUND_I_164_V1_C3"]
    check("N_zero_is_distinct_from_missing_summary_row", no_patients.matched_summary_row and not absent.matched_summary_row and pd.isna(no_patients.raw_p) and no_patients.n_patients == 0)
    check("unplanned_results_are_exported_only_to_exploratory", len(exploratory) == 1 and exploratory.registry_analysis_tier.eq("unplanned_exploratory").all() and "unplanned_program" not in set(output.program))
    without_extra = summary.loc[summary.program.ne("unplanned_program")]
    without_output = adjust_registered(plan, without_extra)[0]
    check("extreme_unplanned_p_cannot_change_planned_adjustment", np.array_equal(output.holm_adjusted_p, without_output.holm_adjusted_p))
    shuffled = adjust_registered(plan.sample(frac=1, random_state=21), summary.sample(frac=1, random_state=8))[0]
    check("input_row_order_does_not_change_adjustment", output.set_index("test_id").holm_adjusted_p.sort_index().equals(shuffled.set_index("test_id").holm_adjusted_p.sort_index()))
    check("statistical_passage_never_declares_complete_biological_replication", output.biological_replication_status.eq("not_assessed_by_this_tool").all())
    empty = pd.DataFrame(columns=RESULT_REQUIRED)
    no_results = adjust_registered(plan, empty)[0]
    check("all_planned_tests_missing_keep_each_family_and_p_one", len(no_results) == 9 and no_results.raw_p.isna().all() and no_results.holm_adjusted_p.eq(1).all() and no_results.family_size_planned.eq(3).all() and not no_results.holm_reject.any())

    bad = plan.copy(); bad.loc[1, "test_id"] = bad.loc[0, "test_id"]
    rejects("duplicate_test_id_rejected", lambda: adjust_registered(bad, summary))
    duplicate = plan.iloc[[0]].copy(); duplicate["test_id"] = "DIFFERENT_ID"; duplicate["family"] = "ANOTHER_FAMILY"
    rejects("duplicate_hypothesis_across_families_rejected", lambda: adjust_registered(pd.concat([plan, duplicate], ignore_index=True), summary))
    rejects("duplicate_summary_key_rejected", lambda: adjust_registered(plan, pd.concat([summary, summary.iloc[[0]]], ignore_index=True)))
    rejects("unplanned_result_strict_policy_rejected", lambda: adjust_registered(plan, summary, "error"))
    bad = plan.copy(); bad.loc[1, "alpha"] = 0.01
    rejects("inconsistent_family_alpha_rejected", lambda: adjust_registered(bad, summary))
    bad = plan.copy(); bad.loc[0, "freeze_status"] = "DRAFT_NOT_FROZEN"
    rejects("unfrozen_plan_rejected", lambda: adjust_registered(bad, summary))
    bad = plan.copy(); bad.loc[0, "freeze_timestamp_utc"] = "2026-11-21T00:00:00"
    rejects("ambiguous_timestamp_without_timezone_rejected", lambda: adjust_registered(bad, summary))
    bad = summary.copy(); bad.loc[0, "sign_p_two_sided"] *= 0.5
    rejects("one_sided_p_masquerading_as_two_sided_rejected", lambda: adjust_registered(plan, bad))
    bad = summary.copy(); bad.loc[0, "n_patients"] += 1
    rejects("inconsistent_patient_sign_counts_rejected", lambda: adjust_registered(plan, bad))
    bad = plan.copy(); bad.loc[0, "role"] = "development"
    rejects("development_result_cannot_be_primary_external_evidence", lambda: adjust_registered(bad, summary))
    bad = plan.copy(); bad.loc[0, "synthetic_example"] = False
    rejects("synthetic_and_real_registry_mixture_rejected", lambda: adjust_registered(bad, summary))

    one = plan.iloc[[0]].copy()
    discordant = pd.DataFrame([dict(dataset=one.cohort.iloc[0], role="external", lineage=one.lineage.iloc[0],
                                   program=one.program.iloc[0], contrast="D_T",
                                   **summarize_difference([1.0]*9+[-20.0], seed=1, n_bootstrap=200))])
    discordant_result = adjust_registered(one, discordant)[0].iloc[0]
    check("mean_and_sign_majority_disagreement_is_explicit", discordant_result.majority_direction_matches and not discordant_result.mean_direction_matches and not discordant_result.direction_matches)
    negative_plan = plan.iloc[[1]].copy(); negative_plan["expected_direction"] = "negative"
    negative = adjust_registered(negative_plan, summary.loc[summary.program.eq("toy_T_program") & summary.contrast.eq("D_T") & summary.dataset.eq("SYN_GSE164522")])[0].iloc[0]
    check("prespecified_negative_direction_is_supported", negative.direction_matches and negative.statistical_direction_criterion_met)

    here = Path(__file__).resolve().parent
    with tempfile.TemporaryDirectory() as directory:
        tmp = Path(directory)
        plan.to_csv(tmp / "plan.csv", index=False)
        summary.to_csv(tmp / "summaries.csv", index=False)
        command = [sys.executable, str(here / "adjust_registered_tests.py"), "--plan", str(tmp / "plan.csv"),
                   "--summaries", str(tmp / "summaries.csv"), "--outdir", str(tmp / "outputs")]
        subprocess.run(command, capture_output=True, text=True, check=True)
        check("CLI_writes_all_four_auditable_outputs", all((tmp / "outputs" / filename).exists() for filename in ["registered_test_results.csv", "unplanned_exploratory_results.csv", "registered_families.csv", "registered_test_audit.json"]))
        saved = read_csv(tmp / "outputs/registered_test_results.csv")
        check("CLI_csv_preserves_missing_raw_p_as_blank", saved.set_index("test_id").loc["PRIMARY_DT_164_V1_C3", "raw_p"] == "")
        persisted_audit = json.loads((tmp / "outputs/registered_test_audit.json").read_text())
        check("CLI_records_input_hashes_and_synthetic_status", len(persisted_audit["plan_sha256"]) == 64 and len(persisted_audit["summaries_sha256"]) == 64 and persisted_audit["synthetic_example"] is True)
        rejected = subprocess.run(command + ["--unplanned-policy", "error"], capture_output=True, text=True)
        check("CLI_strict_unplanned_policy_exits_two", rejected.returncode == 2 and "unplanned" in rejected.stderr)
        # Real interface regression: retain only the EXISTING wholly missing
        # synthetic program, then run both CLIs. No observed program is present
        # to create the mean/p-value columns incidentally.
        existing_scores = read_csv(here / "synthetic_demo/program_scores.csv")
        existing_scores.loc[existing_scores.program.eq("missing_demo")].to_csv(tmp / "only_missing_scores.csv", index=False)
        missing_plan_rows = []
        for contrast in ["D_T", "D_N", "D_T_Q", "D_N_Q", "I"]:
            row = plan.iloc[0].to_dict()
            row.update(test_id="ALL_MISSING_" + contrast, family="SYN_ALL_UNMEASURABLE_DIAGNOSTIC",
                       cohort="DEV", role="development", lineage="Myeloid", program="missing_demo",
                       contrast=contrast, analysis_tier="diagnostic")
            missing_plan_rows.append(row)
        pd.DataFrame(missing_plan_rows).to_csv(tmp / "only_missing_plan.csv", index=False)
        subprocess.run([sys.executable, str(here / "paired_effects.py"), "--manifest", str(here / "synthetic_demo/sample_manifest.csv"),
                        "--scores", str(tmp / "only_missing_scores.csv"), "--n-bootstrap", "100", "--outdir", str(tmp / "only_missing_effects")],
                       capture_output=True, text=True, check=True)
        subprocess.run([sys.executable, str(here / "adjust_registered_tests.py"), "--plan", str(tmp / "only_missing_plan.csv"),
                        "--summaries", str(tmp / "only_missing_effects/effects_summary.csv"), "--outdir", str(tmp / "only_missing_registered")],
                       capture_output=True, text=True, check=True)
        missing_results = pd.read_csv(tmp / "only_missing_registered/registered_test_results.csv")
        check("entirely_unmeasurable_existing_program_survives_paired_to_registry_CLIs",
              len(missing_results) == 5 and missing_results.matched_summary_row.eq(True).all()
              and missing_results.n_patients.eq(0).all() and missing_results.raw_p.isna().all()
              and missing_results.adjustment_only_p.eq(1).all() and missing_results.holm_adjusted_p.eq(1).all()
              and missing_results.family_size_planned.eq(5).all() and missing_results.planned_test_missing.eq(True).all()
              and missing_results.holm_reject.eq(False).all())
    result = {"status": "passed", "n_checks": len(checks), "checks": checks,
              "scope": "Synthetic registry, Holm, direction and CLI checks only; no biological data accessed.",
              "demo_expected_values": {"primary_C1_raw_p": float(first.raw_p), "primary_C1_holm_p": float(first.holm_adjusted_p),
                                       "background_C1_holm_p": float(background.holm_adjusted_p), "n_missing_planned_tests": audit["n_missing_planned_tests"]}}
    write_json(here / "registered_tests_selfcheck_results.json", result)
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
