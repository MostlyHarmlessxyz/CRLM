"""Patient/site summaries and patient-clustered paired comparisons of OOF scores."""

from itertools import combinations

import numpy as np
import pandas as pd
from sklearn.metrics import (
    average_precision_score,
    balanced_accuracy_score,
    roc_auc_score,
)


def summarize(oof, seed):
    keys = ["method", "lineage", "seed", "fold", "donor_uid", "tissue"]
    lesions = oof.groupby(keys + ["lesion_uid"], as_index=False).logit.mean()
    patients = lesions.groupby(keys, as_index=False).logit.mean()
    metrics = []
    for (method, lineage, model_seed), group in patients.groupby(
        ["method", "lineage", "seed"]
    ):
        y = group.tissue.eq("LM")
        pairs = (
            group.pivot(index="donor_uid", columns="tissue", values="logit")
            .reindex(columns=["PT", "LM"])
            .dropna()
        )
        diff = pairs.LM - pairs.PT
        metrics.append(
            dict(
                method=method,
                lineage=lineage,
                seed=model_seed,
                patients=group.donor_uid.nunique(),
                n_patient_sites=len(group),
                paired_n=len(diff),
                mean_paired_logit=diff.mean(),
                fraction_positive_pairs=(diff > 0).mean(),
                auroc=roc_auc_score(y, group.logit) if y.nunique() == 2 else np.nan,
                average_precision=average_precision_score(y, group.logit)
                if y.any()
                else np.nan,
                balanced_accuracy=balanced_accuracy_score(y, group.logit >= 0)
                if y.nunique() == 2
                else np.nan,
            )
        )
    contrasts = []
    for lineage, group in patients[patients.seed == seed].groupby("lineage"):
        wide = group.pivot(
            index=["donor_uid", "tissue"], columns="method", values="logit"
        )
        for first, second in combinations(sorted(wide.columns), 2):
            if not (
                "PCA_mean" in (first, second)
                or {first, second} == {"scGPT_mean", "scGPT_attention"}
            ):
                continue
            common = wide[[first, second]].dropna().reset_index()
            pos, neg = common.tissue.eq("LM"), common.tissue.eq("PT")
            if not pos.any() or not neg.any():
                continue
            donor_codes, donors = pd.factorize(common.donor_uid, sort=True)
            n = len(donors)
            rng = np.random.default_rng(seed)
            weights = rng.multinomial(n, np.full(n, 1 / n), size=2000)
            wp, wn = weights[:, donor_codes[pos]], weights[:, donor_codes[neg]]
            den = wp.sum(1) * wn.sum(1)
            keep = den > 0
            comparisons = []
            for m in (first, second):
                diffs = (
                    common.loc[pos, m].to_numpy()[:, None]
                    - common.loc[neg, m].to_numpy()
                )
                comparisons.append((diffs > 0).astype(float) + 0.5 * (diffs == 0))
            delta = comparisons[1] - comparisons[0]
            samples = np.einsum("bi,ij,bj->b", wp[keep], delta, wn[keep]) / den[keep]
            lo, hi = np.quantile(samples, [0.025, 0.975])
            contrasts.append(
                dict(
                    lineage=lineage,
                    first=first,
                    second=second,
                    n_patients=n,
                    n_sites=len(common),
                    auc_second_minus_first=delta.mean(),
                    lower95=lo,
                    upper95=hi,
                    valid_resamples=int(keep.sum()),
                )
            )
    return patients, pd.DataFrame(metrics), pd.DataFrame(contrasts)
