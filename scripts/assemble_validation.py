"""Apply the archived test plan to source-specific patient estimates."""

from pathlib import Path
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "starter"))
from adjust_registered_tests import run


def main():
    out = ROOT / "artifacts/m6"
    plan = ROOT / "artifacts/m5/registered_tests.csv"
    tables = []
    patient_tables = []
    for accession in ("GSE164522", "GSE203608", "GSE178318"):
        folder = out / accession / "patient_effects"
        if accession == "GSE178318" and not (folder / "effects_summary.csv").exists():
            continue
        summary = pd.read_csv(folder / "effects_summary.csv")
        patients = pd.read_csv(folder / "patient_contrasts.csv")
        keys = ["dataset", "donor_uid", "lineage", "program"]
        if patients.duplicated(keys).any() or set(patients.dataset) != {accession}:
            raise ValueError(f"Invalid external patient identity for {accession}")
        q = patients[patients.I.notna()]
        if not (
            q[["A_PT", "A_LM", "A_NC", "A_NL"]].notna().all().all()
            and np.allclose(q.I, q.D_T_Q - q.D_N_Q)
            and np.allclose(q.D_T_Q, q.D_T)
            and np.allclose(q.D_N_Q, q.D_N)
        ):
            raise ValueError(
                f"Four-tissue contrasts differ from their common patients in {accession}"
            )
        tables.append(summary)
        patient_tables.append(patients)
    summary_path = out / "effects_summary.csv"
    pd.concat(tables).to_csv(summary_path, index=False)
    pd.concat(patient_tables).to_csv(out / "patient_contrasts.csv", index=False)
    run(plan, summary_path, out / "registered_tests")
    adjusted = pd.read_csv(out / "registered_tests/registered_test_results.csv")
    assert len(adjusted) == len(pd.read_csv(plan)) and adjusted.test_id.is_unique
    assert adjusted.groupby("family").size().eq(2).all()
    for accession in ("GSE203608", "GSE178318"):
        folder = out / accession / "registered_tests"
        folder.mkdir(parents=True, exist_ok=True)
        adjusted.loc[adjusted.cohort.eq(accession)].to_csv(folder / "registered_results.csv", index=False)
    print(
        adjusted[
            [
                "cohort",
                "contrast",
                "program",
                "n_patients",
                "mean_difference",
                "direction_status",
                "raw_p",
                "holm_adjusted_p",
            ]
        ].to_string(index=False)
    )


if __name__ == "__main__":
    main()
