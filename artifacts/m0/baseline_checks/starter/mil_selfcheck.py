#!/usr/bin/env python3
"""Meaningful numerical and leakage checks for the small frozen-embedding MIL."""
from __future__ import annotations

import copy
import json
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np

from _common import read_csv, write_json
from linear_attention_mil import (build_design, fit_model, fit_scaler, load_bags, load_model,
                                  loss_and_gradient, predict_bags, save_model,
                                  select_training_rows, validate_provenance)
from make_splits import build_splits
from make_synthetic_embeddings import generate
from validate_manifest import validate_frame


def main():
    here = Path(__file__).resolve().parent
    manifest = here / "synthetic_demo/sample_manifest.csv"
    frame = read_csv(manifest)
    checks, gradient_results = [], []

    def check(name, condition):
        if not condition:
            raise AssertionError(name)
        checks.append(name)

    def rejects(name, function):
        try:
            function()
        except (ValueError, KeyError):
            checks.append(name)
            return
        raise AssertionError("Expected rejection: " + name)

    with tempfile.TemporaryDirectory(prefix="crlm_mil_") as temporary:
        tmp = Path(temporary)
        arrays, provenance = generate(manifest, tmp)
        rows, planned, _ = select_training_rows(frame)
        design = build_design(rows)
        bags = load_bags(tmp / "cell_embeddings.npz", design.samples)
        validate_provenance(provenance, "SyntheticImmune", planned, 5)
        check("only_verified_development_tumor_bags_selected", set(rows.role) == {"development"} and set(rows.tissue) == {"PT", "LM"} and set(rows.donor_uid) == {"DEV_D1", "DEV_D2", "DEV_D3", "DEV_D4"})
        weights = design.weight_table.set_index("sample_uid").bce_weight
        check("unequal_technical_replicates_and_lesions_get_correct_BCE_weights", np.isclose(weights["D1_PT"], 1/8) and np.isclose(weights["D1_LM_a"], 1/48) and np.isclose(weights["D1_LM_2"], 1/16))
        check("each_patient_has_equal_total_BCE_weight", np.allclose(design.weight_table.groupby("donor_uid").bce_weight.sum(), 0.25))
        check("each_patient_tumor_site_has_equal_BCE_weight", np.allclose(design.weight_table.groupby(["donor_uid", "tissue"]).bce_weight.sum(), 0.125))
        fake_logits = np.array([9.0 if s == "D1_LM_2" else 1.0 for s in design.samples])
        check("paired_logits_aggregate_lesions_before_patient_difference", np.isclose((design.pair_matrix @ fake_logits)[design.pair_donors.index("DEV_D1")], 4.0))
        for pooling in ["mean", "attention"]:
            d = bags[0].shape[1]
            size = d + 1 if pooling == "mean" else 2*d + 1
            theta = np.random.default_rng(444).normal(0, 0.25, size)
            for pair_weight in [0.0, 0.7]:
                objective, analytical = loss_and_gradient(theta, bags, design, pooling, pair_weight, 0.13)
                h = 1e-6
                numerical = np.empty_like(theta)
                for index in range(len(theta)):
                    step = np.zeros_like(theta)
                    step[index] = h
                    numerical[index] = (loss_and_gradient(theta + step, bags, design, pooling, pair_weight, 0.13)[0] - loss_and_gradient(theta - step, bags, design, pooling, pair_weight, 0.13)[0]) / (2*h)
                error = float(np.max(np.abs(numerical - analytical)))
                gradient_results.append({"pooling": pooling, "lambda_pair": pair_weight, "maximum_absolute_error": error, "step": h})
                check(f"finite_difference_gradient_{pooling}_pair{pair_weight}", error < 1e-7)
            zero = np.zeros(size)
            loss0, _ = loss_and_gradient(zero, bags, design, pooling, 0.7, 0.13)
            check(f"zero_logit_loss_matches_BCE_plus_pair_{pooling}", np.isclose(loss0, 1.7*np.log(2)))
        leaked = frame.copy()
        leaked.loc[leaked.sample_uid.eq("E1_PT"), "donor_uid"] = "DEV_D1"
        rejects("patient_role_leakage_rejected_before_filtering", lambda: select_training_rows(leaked))
        unverified = frame.copy()
        unverified.loc[unverified.sample_uid.eq("U1_PT"), "role"] = "development"
        rejects("unverified_development_identity_rejected", lambda: select_training_rows(unverified))
        _, donors, _ = validate_frame(frame)
        splits = build_splits(donors, n_folds=2)
        split_rows, split_donors, _ = select_training_rows(frame, splits, fold=0)
        check("fold_uses_all_eligible_tumor_bags_of_only_train_donors", set(split_rows.donor_uid) == set(splits["folds"][0]["train_donors"]) and not set(split_rows.donor_uid) & set(splits["holdout_donors"]["external"]))
        corrupt_split = copy.deepcopy(splits)
        corrupt_split["folds"][0]["train_donors"].append("EXT_E1")
        rejects("external_donor_in_supplied_split_rejected", lambda: select_training_rows(frame, corrupt_split, fold=0))
        corrupt_split = copy.deepcopy(splits)
        corrupt_split["folds"][0]["train_sample_uids"].pop()
        rejects("lesion_or_normal_sample_split_inconsistency_rejected", lambda: select_training_rows(frame, corrupt_split, fold=0))
        bad_provenance = dict(provenance, preprocessing_fit_scope="train_fold_only", preprocessing_training_donors=planned)
        rejects("full_development_PCA_not_accepted_as_fold_fitted", lambda: validate_provenance(bad_provenance, "SyntheticImmune", split_donors, 5))
        bad_provenance = dict(provenance, gene_coverage_audited=False)
        rejects("missing_gene_coverage_audit_rejected", lambda: validate_provenance(bad_provenance, "SyntheticImmune", planned, 5))
        no_pairs = rows.loc[rows.sample_uid.isin(["D1_PT", "D2_LM"])].reset_index(drop=True)
        no_pair_design = build_design(no_pairs)
        rejects("positive_pair_loss_without_pairs_rejected", lambda: fit_model([arrays[s] for s in no_pair_design.samples], no_pair_design))
        rejects("missing_embedding_bag_rejected", lambda: load_bags(tmp / "cell_embeddings.npz", ["MISSING"]))
        rejects("normal_tissues_cannot_receive_binary_tumor_site_BCE", lambda: build_design(frame.loc[frame.sample_uid.isin(["D1_PT", "D1_NL"])].reset_index(drop=True)))
        corrupted = dict(arrays)
        for sample in frame.loc[frame.role.ne("development"), "sample_uid"]:
            corrupted[sample] = np.array([[np.nan]])
        np.savez_compressed(tmp / "corrupt_external.npz", **corrupted)
        isolated = load_bags(tmp / "corrupt_external.npz", design.samples)
        check("training_loader_does_not_read_external_arrays", all(np.array_equal(a, b) for a, b in zip(bags, isolated)))
        center, scale = fit_scaler(bags, design.bce_weights, True)
        manual_center = sum(q*h.mean(axis=0) for q,h in zip(design.bce_weights, bags))
        check("scaler_is_training_patient_weighted", np.allclose(center, manual_center) and not np.allclose(center, np.concatenate(bags).mean(axis=0)))
        doubled = [np.repeat(h, 3, axis=0) for h in bags]
        c2, s2 = fit_scaler(doubled, design.bce_weights, True)
        check("duplicating_cells_does_not_change_scaler_or_patient_weights", np.allclose(center, c2) and np.allclose(scale, s2))
        for pooling in ["mean", "attention"]:
            model, fit = fit_model(bags, design, pooling, standardize=True, seed=123)
            check(f"synthetic_fit_converges_{pooling}", fit["converged"])
            metadata = {"lineage": "SyntheticImmune", "fit": fit, "training_samples": design.samples}
            pred, attention = predict_bags(model, metadata, bags, design.samples)
            one, _ = predict_bags(model, metadata, [bags[0]], [design.samples[0]])
            check(f"single_bag_inference_equals_batch_{pooling}", one.logit.iloc[0] == pred.logit.iloc[0])
            shuffled = [h[::-1] for h in bags]
            reverse, _ = predict_bags(model, metadata, shuffled, design.samples)
            check(f"cell_order_invariant_{pooling}", np.allclose(pred.logit, reverse.logit, atol=1e-12, rtol=1e-12))
            more, _ = predict_bags(model, metadata, [bags[0], np.full((3, 5), 1e6)], [design.samples[0], "UNRELATED"])
            check(f"unrelated_inference_bag_cannot_change_first_prediction_{pooling}", more.logit.iloc[0] == one.logit.iloc[0])
            check(f"attention_is_normalized_nonnegative_{pooling}", all(np.isclose(a.sum(), 1.0) and (a >= 0).all() for a in attention.values()))
            save_model(tmp / f"{pooling}.npz", model, metadata)
            restored, restored_meta = load_model(tmp / f"{pooling}.npz")
            loaded, _ = predict_bags(restored, restored_meta, bags, design.samples)
            check(f"saved_model_reproduces_logits_{pooling}", np.array_equal(pred.logit.to_numpy(), loaded.logit.to_numpy()))
            repeat_model, _ = fit_model(isolated, design, pooling, standardize=True, seed=123)
            check(f"same_seed_and_changed_external_data_leave_fit_identical_{pooling}", np.array_equal(model["theta"], repeat_model["theta"]) and np.array_equal(model["center"], repeat_model["center"]))
        # Exercise user-facing CLIs, including a file with deliberately invalid external arrays.
        for pooling in ["mean", "attention"]:
            command = [sys.executable, str(here / "linear_attention_mil.py"), "train", "--manifest", str(manifest),
                       "--embeddings", str(tmp / "corrupt_external.npz"), "--provenance", str(tmp / "embedding_provenance.json"),
                       "--lineage", "SyntheticImmune", "--pooling", pooling, "--standardize", "--output", str(tmp / f"cli_{pooling}.npz")]
            subprocess.run(command, capture_output=True, text=True, check=True)
            check(f"train_CLI_works_without_external_loading_{pooling}", (tmp / f"cli_{pooling}.npz").exists())
        subprocess.run([sys.executable, str(here / "linear_attention_mil.py"), "predict", "--model", str(tmp / "cli_attention.npz"),
                        "--embeddings", str(tmp / "cell_embeddings.npz"), "--sample", "E1_PT", "--output", str(tmp / "one_external.csv")], capture_output=True, text=True, check=True)
        check("independent_external_predict_CLI", read_csv(tmp / "one_external.csv").sample_uid.tolist() == ["E1_PT"])
    result = {"status": "passed", "n_checks": len(checks), "checks": checks, "gradient_checks": gradient_results,
              "scope": "Synthetic numerical, weighting, identity, role and inference isolation checks; no real biological performance claim."}
    write_json(here / "mil_selfcheck_results.json", result)
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
