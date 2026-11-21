#!/usr/bin/env python3
"""Small adversarial checks using synthetic data only; no biological results."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd
import scipy

from _common import read_csv, write_json
from composition_decomposition import analyze_frame, symmetric_components
from loss_counterexample import review_loss
from make_splits import build_splits
from paired_effects import calculate_patient_contrasts, summarize_difference, summarize_patients
from validate_manifest import validate_frame

ROOT = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default=str(ROOT / "selfcheck_results.json"))
    args = parser.parse_args()
    passed = []

    def check(name, condition):
        if not condition:
            raise AssertionError(name)
        passed.append(name)

    def rejected(name, callback):
        try:
            callback()
        except ValueError:
            passed.append(name)
        else:
            raise AssertionError(f"Did not reject: {name}")

    manifest = read_csv(ROOT / "synthetic_demo/sample_manifest.csv")
    scores = read_csv(ROOT / "synthetic_demo/program_scores.csv")
    validated, donors, warnings = validate_frame(manifest)
    check("Processed input is accepted but cannot enter NB models", not validated.can_nb_pseudobulk.any())
    check("Unresolved identities are recorded as pair-ineligible", not validated.loc[validated.donor_uid.eq("UNRESOLVED_U1"), "eligible_for_paired_analysis"].any())
    mixed_role = manifest.copy()
    mixed_role.loc[0, "role"] = "external"
    rejected("Patient role leakage is rejected", lambda: validate_frame(mixed_role))
    conflicting = pd.concat([manifest, manifest.iloc[[0]].assign(tissue="LM")], ignore_index=True)
    rejected("Duplicate or conflicting sample labels are rejected", lambda: validate_frame(conflicting))
    false_raw = manifest.copy()
    false_raw.loc[0, "input_scale"] = "raw_counts"
    rejected("Processed source cannot masquerade as raw counts", lambda: validate_frame(false_raw))
    false_nb = manifest.copy()
    false_nb["can_nb_pseudobulk"] = "false"
    false_nb.loc[0, "can_nb_pseudobulk"] = "true"
    rejected("NB capability cannot be asserted for processed expression", lambda: validate_frame(false_nb))
    split = build_splits(donors, n_folds=2, seed=7)
    check("Splits are deterministic for a fixed donor manifest/seed", split == build_splits(donors, 2, 7))
    holdout = set(sum(split["holdout_donors"].values(), []))
    test_occurrences = []
    for fold in split["folds"]:
        check(f"Fold {fold['fold']} keeps all donor groups disjoint", not (set(fold["train_donors"]) & set(fold["test_donors"])) and not (holdout & set(fold["train_donors"] + fold["test_donors"])))
        for donor in fold["test_donors"]:
            expected = set(validated.loc[validated.donor_uid.eq(donor), "sample_uid"])
            check(f"Fold {fold['fold']} retains every lesion/tissue of {donor}", expected.issubset(fold["test_sample_uids"]) and expected.isdisjoint(fold["train_sample_uids"]))
        test_occurrences.extend(fold["test_donors"])
    check("Every development donor is held out exactly once", sorted(test_occurrences) == sorted(split["development_donors"]))
    false_identity = donors.copy()
    false_identity.loc[false_identity.donor_uid.eq("DEV_D1"), "patient_mapping_verified"] = False
    rejected("Unverified development identities cannot generate valid folds", lambda: build_splits(false_identity, 2, 7))
    bad_scores = pd.concat([scores, pd.DataFrame([{"sample_uid": "U1_PT", "lineage": "Myeloid", "program": "fixed_demo", "score": "1", "n_cells": "100"}])], ignore_index=True)
    rejected("Unverified patient codes cannot enter paired effect estimates", lambda: calculate_patient_contrasts(manifest, bad_scores))
    rejected("Duplicated score rows cannot overweight a sample", lambda: calculate_patient_contrasts(manifest, pd.concat([scores, scores.iloc[[0]]], ignore_index=True)))
    patients, _, _, _ = calculate_patient_contrasts(manifest, scores)
    fixed = patients[patients.program.eq("fixed_demo")].set_index("donor_uid")
    check("Technical samples and unequal lesion counts preserve lesion-equal D1 effect", np.isclose(fixed.loc["DEV_D1", "A_LM"], 5) and np.isclose(fixed.loc["DEV_D1", "D_T"], 4) and fixed.loc["DEV_D1", "n_lesions_LM"] == 2)
    check("Patients have equal influence despite huge cell-count differences", np.isclose(fixed.loc[["DEV_D1", "DEV_D2"], "D_T"].mean(), 1))
    check("Four-tissue contrast uses same patient and correct sign", np.isclose(fixed.loc["DEV_D1", "I"], 2) and np.isclose(fixed.loc["DEV_D2", "I"], -2))
    check("Missing normal tissue never becomes a zero interaction", pd.isna(fixed.loc["DEV_D3", "I"]))
    check("Four-tissue background components use exactly the interaction subset", fixed.D_T_Q.notna().equals(fixed.I.notna()) and fixed.D_N_Q.notna().equals(fixed.I.notna()) and pd.isna(fixed.loc["DEV_D3", "D_T_Q"]))
    check("Same-subset background means preserve difference identity", np.isclose(fixed.D_T_Q.mean() - fixed.D_N_Q.mean(), fixed.I.mean()))
    check("Missing LM score never becomes a zero tumor difference", pd.isna(fixed.loc["EXT_E3", "D_T"]))
    summary = summarize_patients(patients, n_bootstrap=200)
    empty = summary[summary.program.eq("missing_demo")]
    check("Completely missing submitted candidate remains visible with N=0", len(empty) == 5 and empty.n_patients.eq(0).all())
    ties = summarize_difference([0, 0, 0, 0], seed=1, n_bootstrap=200)
    check("All ties give effective N=0, p=1 and no directional information", ties["n_sign_test"] == 0 and ties["n_ties"] == 4 and ties["sign_p_two_sided"] == 1 and ties["positive_fraction_exact_ci_low"] == 0 and ties["positive_fraction_exact_ci_high"] == 1)
    exact = summarize_difference([1] * 6, seed=1, n_bootstrap=200)
    check("Six concordant patients have exact p=0.03125 and wide CI", np.isclose(exact["sign_p_two_sided"], 0.03125) and np.isclose(exact["positive_fraction_exact_ci_low"], 0.025 ** (1 / 6)))
    comp = analyze_frame(read_csv(ROOT / "synthetic_demo/composition.csv"))
    full = comp[0]["full_decomposition"]
    check("Known composition/state example matches hand-derived values", np.isclose(full["total_difference"], 3.2) and np.isclose(full["composition"], 1.8) and np.isclose(full["within_type_state"], 1.4))
    missing = comp[1]
    check("Missing support withholds full decomposition without zero imputation", missing["full_decomposition"] is None and np.isclose(missing["observed_total_difference"], 1) and np.isclose(missing["common_support"]["renormalized_components"]["total_difference"], 0) and np.isclose(missing["common_support"]["outside_mass_LM"], 0.5))
    rejected("Low-level decomposition rejects unknown subtype means", lambda: symmetric_components([1, 0], [.5, .5], [1, None], [1, 3]))
    loss = review_loss()["rows"]
    check("Normal-gap reversal improves original objective without changing tumor BCE", loss[2]["original_plus_symmetric_penalty"] < loss[1]["original_plus_symmetric_penalty"] and loss[0]["tumor_BCE_unchanged"] == loss[2]["tumor_BCE_unchanged"])
    check("Symmetric-only normal penalty has no reward for reversing its sign", loss[0]["symmetric_penalty_without_difference_reward"] == loss[2]["symmetric_penalty_without_difference_reward"] > loss[1]["symmetric_penalty_without_difference_reward"])
    with tempfile.TemporaryDirectory(prefix="crlm_selfcheck_") as temp:
        temp = Path(temp)

        def run(script, *arguments):
            subprocess.run([sys.executable, str(ROOT / script), *map(str, arguments)], check=True, capture_output=True, text=True, timeout=30)

        run("validate_manifest.py", "--manifest", ROOT / "synthetic_demo/sample_manifest.csv", "--outdir", temp / "manifest")
        run("make_splits.py", "--donor-manifest", temp / "manifest/donor_manifest_validated.csv", "--n-folds", 2, "--output", temp / "folds.json")
        run("paired_effects.py", "--manifest", temp / "manifest/sample_manifest_validated.csv", "--scores", ROOT / "synthetic_demo/program_scores.csv", "--n-bootstrap", 200, "--outdir", temp / "effects")
        run("composition_decomposition.py", "--input", ROOT / "synthetic_demo/composition.csv", "--output", temp / "composition.json")
        run("loss_counterexample.py", "--output", temp / "loss.json")
        check("All five CLI entry points run end to end", (temp / "effects/effects_summary.csv").exists() and json.loads((temp / "folds.json").read_text())["n_folds"] == 2)
    report = {"status": "passed", "n_checks": len(passed), "checks": passed,
              "data": "Synthetic examples only; no biological performance was tested.",
              "environment": {"python": sys.version.split()[0], "numpy": np.__version__, "pandas": pd.__version__, "scipy": scipy.__version__}}
    write_json(args.output, report)
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
