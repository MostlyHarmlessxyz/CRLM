#!/usr/bin/env python3
"""Objective-function counterexample, not a biological simulation or model fit."""
from __future__ import annotations

import argparse
import math

from _common import cli_error, write_json


def softplus(x):
    return max(x, 0) + math.log1p(math.exp(-abs(x)))


def review_loss():
    rows = []
    tumor_gap, margin, penalty_weight = 2.0, 1.0, 0.05
    for normal_gap in [2.0, 0.0, -2.0]:
        original = softplus(margin - tumor_gap + normal_gap)
        symmetric = penalty_weight * normal_gap**2
        rows.append({
            "s_PT": -1.0, "s_LM": 1.0, "s_NC": 0.0, "s_NL": normal_gap,
            "tumor_BCE_unchanged": softplus(-1.0),
            "D_T": tumor_gap, "D_N": normal_gap, "I": tumor_gap - normal_gap,
            "original_background_loss": original,
            "original_plus_symmetric_penalty": original + symmetric,
            "symmetric_penalty_without_difference_reward": symmetric,
        })
    return {"type": "mathematical_counterexample_no_patient_data", "margin": margin, "lambda_normal": penalty_weight,
            "rows": rows,
            "conclusion": "Adding a small symmetric square penalty to the original difference reward still favors a negative normal gap. Removing the difference reward and using the square alone removes this directional incentive; it does not establish causal deconfounding."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    cli_error(parser, lambda: write_json(args.output, review_loss()))
    print("Wrote the loss counterexample; no training or patient data were used.")


if __name__ == "__main__":
    main()
