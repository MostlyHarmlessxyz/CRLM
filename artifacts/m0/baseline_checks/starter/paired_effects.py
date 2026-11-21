#!/usr/bin/env python3
"""Patient-level effects for already fixed sample x lineage x program scores.

Technical samples within lesion -> equal lesion means within donor/site ->
patient contrasts. Cells and lesions never count as independent patients.
Mean-effect percentile bootstrap is approximate; the sign test/interval use
the exact binomial model after excluding ties.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binomtest

from _common import cli_error, read_csv, require_columns, require_nonempty, sha256_file, stable_seed, write_json
from validate_manifest import validate_frame

GROUP = ["dataset", "donor_uid", "role", "lineage", "program"]


def calculate_patient_contrasts(manifest, scores, skip_ineligible=False):
    manifest, _, _ = validate_frame(manifest)
    require_columns(scores, ["sample_uid", "lineage", "program", "score"], "scores")
    require_nonempty(scores, ["sample_uid", "lineage", "program"], "scores")
    keys = ["sample_uid", "lineage", "program"]
    if scores.duplicated(keys).any():
        raise ValueError("Scores must have one row per sample_uid x lineage x program; duplicate rows could overweight a sample")
    unknown = sorted(set(scores.sample_uid) - set(manifest.sample_uid))
    if unknown:
        raise ValueError(f"Scores have samples absent from manifest: {unknown[:10]}")
    allowed_columns = keys + ["score", "n_cells"]
    unsupported = set(scores.columns) - set(allowed_columns)
    if unsupported:
        raise ValueError(f"Keep sample metadata in manifest, not scores; unsupported score columns: {sorted(unsupported)}")
    score_values = pd.to_numeric(scores["score"].replace("", np.nan), errors="raise")
    if np.isinf(score_values.to_numpy(dtype=float)).any():
        raise ValueError("Infinite scores are invalid")
    scores = scores.copy()
    scores["score"] = score_values
    joined = scores.merge(manifest, on="sample_uid", how="left", validate="many_to_one")
    invalid = ~joined.eligible_for_paired_analysis
    if invalid.any() and not skip_ineligible:
        bad = sorted(joined.loc[invalid, "donor_uid"].unique())
        raise ValueError(f"Scores include excluded or identity-unverified donors {bad}; curate input or explicitly --skip-ineligible")
    audit = {
        "n_input_score_rows": len(joined),
        "n_ineligible_score_rows": int(invalid.sum()),
        "ineligible_donors": sorted(joined.loc[invalid, "donor_uid"].unique()),
        "n_missing_score_rows": int(joined.score.isna().sum()),
        "weighting": "Equal technical-sample means within lesion, then equal lesion means within donor/tissue, then equal patients. n_cells is never an inference weight.",
    }
    eligible = joined.loc[~invalid].copy()
    if eligible.empty:
        raise ValueError("No identity-verified, non-excluded donor scores")
    usable = eligible.loc[eligible.score.notna()].copy()
    lesion_keys = GROUP + ["tissue", "lesion_uid"]
    lesions = usable.groupby(lesion_keys, as_index=False, sort=True).agg(score=("score", "mean"), n_technical_samples=("sample_uid", "nunique"))
    tissue = lesions.groupby(GROUP + ["tissue"], as_index=False, sort=True).agg(score=("score", "mean"), n_lesions=("lesion_uid", "nunique"))
    expected_index = pd.MultiIndex.from_frame(eligible[GROUP].drop_duplicates())
    values = tissue.pivot(index=GROUP, columns="tissue", values="score").reindex(index=expected_index, columns=["PT", "LM", "NC", "NL"])
    counts = tissue.pivot(index=GROUP, columns="tissue", values="n_lesions").reindex(index=expected_index, columns=["PT", "LM", "NC", "NL"]).fillna(0).astype(int)
    patients = values.rename(columns=lambda x: "A_" + x).join(counts.rename(columns=lambda x: "n_lesions_" + x)).reset_index()
    patients["D_T"] = patients.A_LM - patients.A_PT
    patients["D_N"] = patients.A_NL - patients.A_NC
    patients["I"] = patients.D_T - patients.D_N
    # Background interpretation must compare the SAME complete four-site
    # patient set. All-pair D_T and all-normal-pair D_N may have different N.
    complete_four = patients["I"].notna()
    patients["D_T_Q"] = patients["D_T"].where(complete_four)
    patients["D_N_Q"] = patients["D_N"].where(complete_four)
    patients["complete_four_tissue"] = complete_four
    audit["n_unique_eligible_donors"] = int(patients.donor_uid.nunique())
    audit["n_patient_program_rows"] = len(patients)
    return patients, lesions, tissue, audit


def summarize_difference(values, seed, n_bootstrap=10000, tie_tolerance=0.0):
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]
    if n_bootstrap < 100:
        raise ValueError("Use at least 100 bootstrap replicates (default 10000)")
    if tie_tolerance < 0 or not np.isfinite(tie_tolerance):
        raise ValueError("tie_tolerance must be finite and >= 0, fixed before viewing results")
    n = len(x)
    if not n:
        # Preserve the complete output schema even when EVERY submitted
        # program/contrast is unmeasurable. Downstream registries must be able
        # to retain these planned tests without inventing observed zero effects.
        return {"n_patients": 0, "status": "no_complete_contrast",
                "mean_difference": None, "median_difference": None,
                "bootstrap_mean_ci_low": None, "bootstrap_mean_ci_high": None,
                "n_positive": 0, "n_negative": 0, "n_ties": 0, "n_sign_test": 0,
                "positive_fraction_nonzero": None, "positive_fraction_ties_half_descriptive": None,
                "positive_fraction_exact_ci_low": None, "positive_fraction_exact_ci_high": None,
                "sign_p_two_sided": None}
    positive = int(np.sum(x > tie_tolerance))
    negative = int(np.sum(x < -tie_tolerance))
    ties = n - positive - negative
    effective_n = positive + negative
    if effective_n:
        result = binomtest(positive, effective_n, p=0.5, alternative="two-sided")
        interval = result.proportion_ci(confidence_level=0.95, method="exact")
        pvalue, lower, upper = float(result.pvalue), float(interval.low), float(interval.high)
    else:
        pvalue, lower, upper = 1.0, 0.0, 1.0  # No nonzero directions: no evidence, not a successful replication.
    if n >= 2:
        rng = np.random.default_rng(seed)
        means = np.mean(x[rng.integers(0, n, size=(n_bootstrap, n))], axis=1)
        ci_low, ci_high = map(float, np.quantile(means, [0.025, 0.975]))
    else:
        ci_low, ci_high = None, None
    return {"n_patients": n, "status": "all_ties_no_direction_information" if effective_n == 0 else ("one_patient_no_bootstrap_ci" if n == 1 else "estimated"),
            "mean_difference": float(x.mean()), "median_difference": float(np.median(x)),
            "bootstrap_mean_ci_low": ci_low, "bootstrap_mean_ci_high": ci_high,
            "n_positive": positive, "n_negative": negative, "n_ties": ties, "n_sign_test": effective_n,
            "positive_fraction_nonzero": positive / effective_n if effective_n else None,
            "positive_fraction_ties_half_descriptive": (positive + 0.5 * ties) / n,
            "positive_fraction_exact_ci_low": lower, "positive_fraction_exact_ci_high": upper,
            "sign_p_two_sided": pvalue}


def summarize_patients(patients, seed=20261121, n_bootstrap=10000, tie_tolerance=0.0):
    rows = []
    group_keys = ["dataset", "role", "lineage", "program"]
    for key, group in patients.groupby(group_keys, sort=True):
        if group.donor_uid.duplicated().any():
            raise ValueError("A donor occurs twice in a program summary")
        for contrast in ["D_T", "D_N", "D_T_Q", "D_N_Q", "I"]:
            summary = summarize_difference(group[contrast], stable_seed(seed, *key, contrast), n_bootstrap, tie_tolerance)
            rows.append({**dict(zip(group_keys, key)), "contrast": contrast, **summary})
    return pd.DataFrame(rows)


def run(manifest_path, scores_path, outdir, seed=20261121, n_bootstrap=10000, tie_tolerance=0.0, skip_ineligible=False):
    patients, lesions, tissue, audit = calculate_patient_contrasts(read_csv(manifest_path), read_csv(scores_path), skip_ineligible)
    summary = summarize_patients(patients, seed, n_bootstrap, tie_tolerance)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    for name, frame in [("lesion_means", lesions), ("patient_tissue_means", tissue), ("patient_contrasts", patients), ("effects_summary", summary)]:
        frame.to_csv(outdir / (name + ".csv"), index=False)
    audit.update({"manifest_sha256": sha256_file(manifest_path), "scores_sha256": sha256_file(scores_path), "seed": seed,
                  "n_bootstrap": n_bootstrap, "tie_tolerance": tie_tolerance,
                  "ci_caution": "Percentile bootstrap is approximate, particularly for tiny N. The exact binomial CI concerns nonzero direction probability, not mean effect.",
                  "multiplicity": "Raw two-sided sign p-values only. Freeze candidate families before choosing any multiple-testing correction.",
                  "common_four_tissue_subset": "D_T_Q, D_N_Q and I use exactly the same complete four-tissue patients; do not subtract overall D_T and D_N means with different coverage.",
                  "causal_caution": "I is the observed organ-by-tumor contrast, not an identified causal metastasis effect."})
    write_json(outdir / "analysis_audit.json", audit)
    print(f"Wrote {len(patients)} patient-program rows and {len(summary)} cohort-specific effect summaries.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--scores", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--seed", type=int, default=20261121)
    parser.add_argument("--n-bootstrap", type=int, default=10000)
    parser.add_argument("--tie-tolerance", type=float, default=0.0)
    parser.add_argument("--skip-ineligible", action="store_true")
    args = parser.parse_args()
    cli_error(parser, lambda: run(args.manifest, args.scores, args.outdir, args.seed, args.n_bootstrap, args.tie_tolerance, args.skip_ineligible))


if __name__ == "__main__":
    main()
