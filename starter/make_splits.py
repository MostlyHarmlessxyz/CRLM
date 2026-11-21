#!/usr/bin/env python3
"""Create deterministic donor folds from validate_manifest.py's donor table.

External, exploratory and excluded roles never enter cross-validation.
This is a single outer split plan; model selection needs donor-grouped inner
validation or a configuration fixed before these outer results are inspected.
"""
from __future__ import annotations

import argparse
import json

import numpy as np

from _common import bool_column, cli_error, read_csv, require_columns, require_nonempty, sha256_file, write_json
from validate_manifest import ROLES


def build_splits(frame, n_folds=5, seed=20261121):
    required = ["donor_uid", "role", "patient_mapping_verified", "sample_uids"]
    require_columns(frame, required, "donor manifest")
    require_nonempty(frame, required, "donor manifest")
    if frame.donor_uid.duplicated().any():
        raise ValueError("donor_manifest must have exactly one row per global donor_uid")
    if not set(frame.role).issubset(ROLES):
        raise ValueError("Unknown role in donor manifest")
    verified = bool_column(frame.patient_mapping_verified, "patient_mapping_verified")
    if ((frame.role == "development") & ~verified).any():
        raise ValueError("Unverified development donor identities cannot be used to make independent folds")
    sample_by_donor, owner = {}, {}
    for row in frame.itertuples(index=False):
        samples = json.loads(row.sample_uids)
        if not isinstance(samples, list) or not samples or not all(isinstance(s, str) and s for s in samples):
            raise ValueError(f"{row.donor_uid}: sample_uids must be a nonempty JSON string array")
        if len(samples) != len(set(samples)):
            raise ValueError(f"{row.donor_uid}: duplicate sample_uids")
        for sample in samples:
            if sample in owner:
                raise ValueError(f"sample_uid {sample} belongs to more than one donor")
            owner[sample] = row.donor_uid
        sample_by_donor[row.donor_uid] = sorted(samples)
    donors = sorted(frame.loc[frame.role == "development", "donor_uid"])
    if not 2 <= n_folds <= len(donors):
        raise ValueError(f"Require 2 <= n_folds <= number of verified development donors ({len(donors)})")
    rng = np.random.default_rng(seed)
    shuffled = np.asarray(donors, dtype=object)[rng.permutation(len(donors))]
    assignment = {str(donor): idx % n_folds for idx, donor in enumerate(shuffled)}
    sample_list = lambda ids: sorted(s for d in ids for s in sample_by_donor[d])
    folds = []
    for fold in range(n_folds):
        test = sorted(d for d in donors if assignment[d] == fold)
        train = sorted(set(donors) - set(test))
        folds.append({"fold": fold, "train_donors": train, "test_donors": test,
                      "train_sample_uids": sample_list(train), "test_sample_uids": sample_list(test)})
    holdouts = {role: sorted(frame.loc[frame.role == role, "donor_uid"]) for role in sorted(ROLES - {"development"})}
    return {"seed": int(seed), "n_folds": n_folds, "group_unit": "donor_uid",
            "development_donors": donors, "holdout_donors": holdouts,
            "donor_test_fold": assignment, "folds": folds,
            "note": "All samples and lesions of a donor move together. No external donor enters CV. Donor identities must already be source-verified."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--donor-manifest", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--n-folds", type=int, default=5)
    parser.add_argument("--seed", type=int, default=20261121)
    args = parser.parse_args()

    def work():
        result = build_splits(read_csv(args.donor_manifest), args.n_folds, args.seed)
        result["source_sha256"] = sha256_file(args.donor_manifest)
        write_json(args.output, result)
        print(f"Wrote {args.n_folds} donor folds; {len(result['development_donors'])} development donors; external donors excluded.")

    cli_error(parser, work)


if __name__ == "__main__":
    main()
