#!/usr/bin/env python3
"""Symmetric composition/state decomposition of patient-level linear means.

Requires PT and LM proportions and subtype means for the SAME fixed program.
Missing conditional means are never filled with zero. For incomplete support,
the full decomposition is withheld and common-support results are identified
as a different, conditional estimand. No classifier outputs are decomposed.
"""
from __future__ import annotations

import argparse

import numpy as np
import pandas as pd

from _common import cli_error, read_csv, require_columns, require_nonempty, sha256_file, write_json


def symmetric_components(pi_pt, pi_lm, mu_pt, mu_lm):
    arrays = [np.asarray(a, dtype=float) for a in [pi_pt, pi_lm, mu_pt, mu_lm]]
    if not all(a.ndim == 1 and a.size == arrays[0].size for a in arrays) or not arrays[0].size:
        raise ValueError("Equal, nonempty one-dimensional vectors are required")
    if not all(np.isfinite(a).all() for a in arrays):
        raise ValueError("Missing/non-finite subtype means: full components are not identifiable")
    pp, pm, up, um = arrays
    if (pp < 0).any() or (pm < 0).any() or not np.isclose(pp.sum(), 1) or not np.isclose(pm.sum(), 1):
        raise ValueError("Proportions must be nonnegative and sum to one")
    total = float(np.dot(pm, um) - np.dot(pp, up))
    composition = float(np.dot(pm - pp, (um + up) / 2))
    state = float(np.dot((pm + pp) / 2, um - up))
    return {"total_difference": total, "composition": composition, "within_type_state": state,
            "identity_residual": total - composition - state}


def decompose_pair(group):
    if set(group.tissue) != {"PT", "LM"}:
        return {"status": "missing_tissue", "full_decomposition": None, "common_support": None}
    pt = group[group.tissue == "PT"].set_index("subtype")
    lm = group[group.tissue == "LM"].set_index("subtype")
    for tissue, data in [("PT", pt), ("LM", lm)]:
        if not np.isclose(data.proportion.sum(), 1.0, rtol=0.0, atol=1e-6):
            raise ValueError(f"{tissue}: supplied proportions must sum to 1; received {data.proportion.sum()}")
    subtypes = sorted(set(pt.index) | set(lm.index))
    pt = pt.reindex(subtypes)
    lm = lm.reindex(subtypes)
    pp, pm = pt.proportion.to_numpy(dtype=float), lm.proportion.to_numpy(dtype=float)
    up, um = pt.mean_score.to_numpy(dtype=float), lm.mean_score.to_numpy(dtype=float)
    for tissue, proportion, mean in [("PT", pp, up), ("LM", pm, um)]:
        if ((proportion == 0) & np.isfinite(mean)).any():
            raise ValueError(f"{tissue}: zero-proportion subtype must have a missing conditional mean, not an invented/imputed value")
    # Zero mass contributes zero to an observed marginal mean; its unknown
    # conditional mean remains unknown and is NEVER changed to zero.
    def observed_marginal(proportion, mean):
        if not np.isfinite(proportion).all():
            return None
        positive = proportion > 0
        if not np.isfinite(mean[positive]).all():
            return None
        return float(np.dot(proportion[positive], mean[positive]))
    mean_pt, mean_lm = observed_marginal(pp, up), observed_marginal(pm, um)
    observed_total = mean_lm - mean_pt if mean_pt is not None and mean_lm is not None else None
    both_zero = (pp == 0) & (pm == 0)
    active = ~both_zero
    complete = np.isfinite(pp[active]).all() and np.isfinite(pm[active]).all() and np.isfinite(up[active]).all() and np.isfinite(um[active]).all()
    full = symmetric_components(pp[active], pm[active], up[active], um[active]) if complete else None
    common = (pp > 0) & (pm > 0) & np.isfinite(up) & np.isfinite(um)
    mass_pt, mass_lm = float(np.nansum(pp[common])), float(np.nansum(pm[common]))
    common_result = None
    if common.any() and mass_pt > 0 and mass_lm > 0:
        common_result = {
            "subtypes": [s for s, keep in zip(subtypes, common) if keep],
            "retained_mass_PT": mass_pt, "retained_mass_LM": mass_lm,
            "outside_mass_PT": 1 - mass_pt, "outside_mass_LM": 1 - mass_lm,
            "renormalized_components": symmetric_components(pp[common] / mass_pt, pm[common] / mass_lm, up[common], um[common]),
            "estimand": "Conditional on common measured support; cannot replace the total population difference.",
        }
    reasons = []
    for index, subtype in enumerate(subtypes):
        if both_zero[index]:
            continue
        if not all(np.isfinite(a[index]) for a in [pp, pm, up, um]):
            reasons.append(subtype)
    return {"status": "complete_support" if complete else "full_components_not_identifiable",
            "observed_mean_PT": mean_pt, "observed_mean_LM": mean_lm, "observed_total_difference": observed_total,
            "unidentified_subtypes": reasons, "full_decomposition": full, "common_support": common_result}


def analyze_frame(frame):
    required = ["dataset", "donor_uid", "lineage", "program", "tissue", "subtype", "proportion", "mean_score"]
    require_columns(frame, required, "composition input")
    require_nonempty(frame, required[:-1], "composition input")
    if not set(frame.tissue).issubset({"PT", "LM"}):
        raise ValueError("This entry point takes PT/LM only; submit any other contrast separately with an explicitly documented mapping")
    keys = ["dataset", "donor_uid", "lineage", "program", "tissue", "subtype"]
    if frame.duplicated(keys).any():
        raise ValueError("Duplicate donor/tissue/subtype/program rows; aggregate to the declared patient-level estimand first")
    data = frame.copy()
    data["proportion"] = pd.to_numeric(data.proportion, errors="raise")
    data["mean_score"] = pd.to_numeric(data.mean_score.replace("", np.nan), errors="raise")
    if not np.isfinite(data.proportion).all() or (data.proportion < 0).any() or (data.proportion > 1).any():
        raise ValueError("Each proportion must be finite and between 0 and 1")
    if np.isinf(data.mean_score).any():
        raise ValueError("Infinite mean_score is invalid; missing means must be blank")
    groups = ["dataset", "donor_uid", "lineage", "program"]
    results = []
    for key, group in data.groupby(groups, sort=True):
        try:
            result = decompose_pair(group)
        except ValueError as error:
            raise ValueError(f"{dict(zip(groups, key))}: {error}") from error
        results.append({**dict(zip(groups, key)), **result})
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    def work():
        results = analyze_frame(read_csv(args.input))
        write_json(args.output, {"source_sha256": sha256_file(args.input), "results": results,
                                "scope": "Algebraic decomposition of linear patient/site means, not causal effects or nonlinear prediction scores."})
        print(f"Wrote {len(results)} patient-level decompositions; unsupported full components were withheld.")

    cli_error(parser, work)


if __name__ == "__main__":
    main()
