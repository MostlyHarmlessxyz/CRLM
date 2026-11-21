"""Patient-level comparison of direct MIL adapters to their matched mean baselines."""
from pathlib import Path

import numpy as np
import pandas as pd
from crlm.benchmark_statistics import summarize


def paired_auc(patient, baseline, direct, seed):
    matched = patient.pivot(index=["donor_uid","tissue"],columns="method",values="logit")[[baseline,direct]].dropna().reset_index()
    donor, donors = pd.factorize(matched.donor_uid,sort=True)
    y = matched.tissue.eq("LM")
    kernel = []
    for method in (baseline,direct):
        diff=matched.loc[y,method].to_numpy()[:,None]-matched.loc[~y,method].to_numpy()
        kernel.append((diff>0).astype(float)+.5*(diff==0))
    delta=kernel[1]-kernel[0]
    n=len(donors)
    weights=np.random.default_rng(seed).multinomial(n,np.full(n,1/n),size=2000)
    pos,neg=weights[:,donor[y]],weights[:,donor[~y]]
    den=pos.sum(1)*neg.sum(1)
    valid=den>0
    boot=np.einsum("bi,ij,bj->b",pos[valid],delta,neg[valid])/den[valid]
    lo,hi=np.quantile(boot,[.025,.975])
    return {"reference":baseline,"method":direct,"auc_delta":delta.mean(),"lower95":lo,"upper95":hi,
            "n_patients":n,"n_sites":len(matched),"n_boot":len(boot)}


def run():
    root=Path(__file__).resolve().parents[2]
    stage=root/"artifacts/m4"
    dest=stage/"direct_prior"
    seed=20261121
    orig=pd.concat([pd.read_csv(stage/rep/"oof_predictions.csv") for rep in ("PCA","scGPT")])
    orig=orig[orig.method.isin(["PCA_mean","scGPT_mean"]) & (orig.seed==seed)]
    direct=pd.concat([pd.read_csv(dest/rep/"oof_predictions.csv") for rep in ("PCA32","scGPT32","scGPT512")])
    samples=pd.concat([orig,direct],ignore_index=True)
    if samples.duplicated(["method","lineage","sample_uid","seed"]).any():
        raise ValueError("Repeated OOF sample in comparison")
    pat,met,_=summarize(samples,seed)
    pat.to_csv(dest/"patient_site_logits.csv",index=False)
    met.to_csv(dest/"patient_metrics.csv",index=False)
    compared=[]
    for lineage,group in pat.groupby("lineage"):
        for rep,baseline in (("PCA32","PCA_mean"),("scGPT32","scGPT32_mean"),("scGPT512","scGPT_mean")):
            for head in (("HA","CTA") if rep=="scGPT512" else ("HA","CTA","MixMIL")):
                compared.append({"lineage":lineage,**paired_auc(group,baseline,f"{rep}_{head}",seed)})
    pd.DataFrame(compared).to_csv(dest/"matched_mean_comparisons.csv",index=False)


if __name__=="__main__":
    run()
