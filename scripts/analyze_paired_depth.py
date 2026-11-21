"""Development sensitivity of primary program effects to library-depth distributions."""
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
from crlm.program_scores import rank_scores


def run():
    out=ROOT/"artifacts/m3"
    members=pd.read_csv(out/"program_definitions_draft/program_genes.csv")
    members=members[(members.role=="complete_program") & members.on_discovery_axis]
    axis=(out/"program_definitions_draft/ranking_background.txt").read_text().splitlines()
    manifest=pd.read_csv(ROOT/"data/manifests/m0/sample_manifest_validated.csv")
    summary=[]
    for r in manifest[manifest.tissue.isin(["PT","LM"])].itertuples(index=False):
        d=ad.read_h5ad(ROOT/f"data/processed/m2/{r.sample_uid}.h5ad")
        cols=pd.Index(d.var.gene_symbol).get_indexer(axis)
        if (cols<0).any():
            raise ValueError("Missing discovery background")
        for lineage, defs in members.groupby("lineage"):
            ix=np.flatnonzero(d.obs.qc_retained & (d.obs.lineage_m2==lineage))
            if len(ix)<50:
                continue
            programs=defs.groupby("program").gene_symbol.apply(list).to_dict()
            score=rank_scores(d.X[ix,:][:,cols],axis,programs)
            bucket=np.digitize(d.obs.n_counts_visible.to_numpy()[ix],[1000,3000,10000])
            for name in programs:
                for group in range(4):
                    subset=bucket==group
                    summary.append((r.donor_uid,r.tissue,r.lesion_uid,r.sample_uid,lineage,name,group,
                                    int(subset.sum()),float(subset.mean()),float(score.loc[subset,name].sum()/len(ix))))
    cells=pd.DataFrame(summary,columns=["donor_uid","tissue","lesion_uid","sample_uid","lineage","program",
                                         "depth_bin","n_cells","proportion","joint_score"])
    if manifest.lesion_uid.duplicated().any():
        raise ValueError("Aggregate technical specimens to lesions before standardizing")
    patient=cells.groupby(["donor_uid","tissue","lineage","program","depth_bin"]).agg(
        n_cells=("n_cells","sum"),p=("proportion","mean"),joint=("joint_score","mean")).reset_index()
    rows=[]
    for (donor,lineage,program),g in patient.groupby(["donor_uid","lineage","program"]):
        if set(g.tissue)!={"PT","LM"}:
            continue
        pt=g[g.tissue=="PT"].set_index("depth_bin")
        lm=g[g.tissue=="LM"].set_index("depth_bin")
        supported=(pt.n_cells>=10)&(lm.n_cells>=10)
        overlap=np.minimum(pt.p,lm.p).where(supported,0)
        mass=overlap.sum()
        differences=(lm.joint/lm.p)-(pt.joint/pt.p)
        contrast=(differences[supported]*overlap[supported]).sum()/mass if mass>0 else np.nan
        rows.append((donor,lineage,program,lm.joint.sum()-pt.joint.sum(),contrast,mass,int(supported.sum())))
    result=pd.DataFrame(rows,columns=["donor_uid","lineage","program","observed_D_T","depth_standardized_D_T",
                                       "overlap_mass","shared_bins"])
    primary=pd.read_csv(out/"patient_effects/patient_contrasts.csv")
    required=primary[primary.program.isin(members.program) & primary.D_T.notna()]
    proof=result.merge(required,on=["donor_uid","lineage","program"],validate="one_to_one")
    if len(proof)!=len(required) or len(proof)!=len(result) or not np.allclose(proof.observed_D_T,proof.D_T,rtol=0,atol=1e-12):
        raise ValueError("Paired support or observed effect differs from primary analysis")
    result.to_csv(out/"patient_depth_standardized.csv",index=False)
    print(result.groupby("program")["depth_standardized_D_T"].agg(["count","mean","median"]).to_string())


if __name__=="__main__":
    run()
