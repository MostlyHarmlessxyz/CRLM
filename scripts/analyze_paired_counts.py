"""Paired pseudobulk negative-binomial contrasts in the GSE245552 development set."""
import json
import resource
import sys
import time
from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path

import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from scipy import sparse
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[1]
M2 = ROOT/"artifacts/m2/pseudobulk_manifest.csv"
AXIS = ROOT/"data/processed/m2/pseudobulk/gene_symbols.txt"
OUT = ROOT/"artifacts/m3/paired_counts"
PARENTS = {"T/NK":"T_NK","Mono/Mac/DC":"Mono_Mac_DC","B/Plasma":"B_Plasma"}
REFERENCE_READOUTS = {"T/NK":["STAT1","IFIT1","GZMB","NKG7"],
                      "Mono/Mac/DC":["SPP1","CXCL9","CXCL10","HLA-DRA"],
                      "B/Plasma":["XBP1","JCHAIN","MZB1","HSPA5"]}


def collect(records, symbols):
    matrix, meta = [], []
    for (donor, tissue), group in records.groupby(["donor_uid","tissue"],sort=True):
        total=np.zeros(len(symbols),dtype=np.int64)
        for r in group.itertuples():
            counts=sparse.load_npz(ROOT/r.count_matrix)[int(r.matrix_row),:].toarray().ravel()
            if counts.dtype.kind not in "iu" or (counts<0).any() or counts.sum()!=r.common_axis_visible_count_sum:
                raise ValueError(f"Invalid counts for {r.sample_uid} {r.group}")
            total+=counts
        key=f"{donor}_{tissue}"
        meta.append({"pseudobulk_uid":key,"donor_uid":donor,"tissue":tissue,
                     "sample_uids":" ".join(group.sample_uid),"lesion_uids":" ".join(group.lesion_uid),
                     "n_cells":int(group.n_cells.sum()),"count_sum":int(total.sum())})
        matrix.append(total)
    meta=pd.DataFrame(meta).set_index("pseudobulk_uid")
    meta["tissue"]=pd.Categorical(meta.tissue,categories=["PT","LM"])
    return pd.DataFrame(matrix,index=meta.index,columns=symbols), meta


def fit(records, symbols, branch, lineage):
    where=OUT/PARENTS[lineage]/branch
    where.mkdir(parents=True,exist_ok=True)
    counts,meta=collect(records,symbols)
    n_pairs=meta.donor_uid.nunique()
    if len(meta)!=2*n_pairs or not meta.groupby("donor_uid",observed=True).tissue.nunique().eq(2).all():
        raise ValueError("One PT and one LM pseudobulk per included patient are required")
    # Preselected expression filter; preserve flags and gene-level NA in the output.
    selected=(counts>=10).sum(axis=0)>=n_pairs
    with threadpool_limits(limits=1):
        dds=DeseqDataSet(counts=counts.loc[:,selected],metadata=meta,design="~ donor_uid + tissue",
                         size_factors_fit_type="poscounts",refit_cooks=False,n_cpus=2,quiet=True)
        design=dds.obsm["design_matrix"]
        rank=np.linalg.matrix_rank(design.to_numpy())
        if list(design.index)!=list(meta.index) or rank!=design.shape[1] or len(design)-rank<=0:
            raise ValueError("Singular design or misordered pseudobulk libraries")
        design.to_csv(where/"design.csv")
        dds.deseq2()
        stats=DeseqStats(dds,contrast=["tissue","LM","PT"],alpha=.05,independent_filter=False,
                         cooks_filter=True,quiet=True,n_cpus=2)
        stats.summary()
    result=stats.results_df.reindex(symbols)
    result.index.name="gene_symbol"
    result["passed_count_filter"]=selected.reindex(symbols).to_numpy()
    converged_cols=[c for c in dds.var.columns if "converged" in c]
    for col in converged_cols+["dispersions","_pvalue_cooks_outlier"]:
        result[col]=dds.var[col].reindex(symbols)
    result.to_csv(where/"gene_results.csv.gz")
    meta["size_factor"]=dds.obs["size_factors"].to_numpy()
    meta.to_csv(where/"libraries.csv")
    convergence={c:int((dds.var[c]==False).sum()) for c in converged_cols}
    return {"lineage":lineage,"analysis":branch,"n_pairs":n_pairs,"n_samples":len(meta),
            "n_input_genes":len(symbols),"n_expression_eligible":int(selected.sum()),
            "n_gene_p":int(result.pvalue.notna().sum()),"n_BH_005":int((result.padj<.05).sum()),
            "rank":int(rank),"residual_df":int(len(design)-rank),"convergence_flags_false":json.dumps(convergence)}


def run():
    start=time.monotonic()
    started=datetime.now(timezone.utc).isoformat()
    symbols=AXIS.read_text().splitlines()
    tables=pd.read_csv(M2)
    tables=tables[(tables.level=="lineage") & (tables.tissue.isin(["PT","LM"])) & (tables.n_cells>=50)]
    manifest=pd.read_csv(ROOT/"data/manifests/m0/sample_manifest_validated.csv")
    if set(manifest.dataset)!={"GSE245552"} or not manifest.patient_mapping_verified.all():
        raise ValueError("Expected source-verified GSE245552 samples")
    table=tables.merge(manifest[["sample_uid","original_specimen_id"]],on="sample_uid",validate="many_to_one")
    all_rows=[]
    reference=[]
    for parent in PARENTS:
        subset=table[table.group==parent]
        linked=subset.groupby("donor_uid").tissue.nunique()
        subset=subset[subset.donor_uid.isin(linked[linked==2].index)].sort_values(["donor_uid","tissue","original_specimen_id"])
        for name, rows in [("pooled_lesions",subset),("first_eligible_specimen",subset.drop_duplicates(["donor_uid","tissue"]))]:
            print(parent,name,len(rows),"GSM",rows.donor_uid.nunique(),"patients",flush=True)
            row=fit(rows,symbols,name,parent)
            all_rows.append(row)
            pd.DataFrame(all_rows).to_csv(OUT/"fit_summary.csv",index=False)
            table_out=pd.read_csv(OUT/PARENTS[parent]/name/"gene_results.csv.gz")
            ref=table_out[table_out.gene_symbol.isin(REFERENCE_READOUTS[parent])].copy()
            ref.insert(0,"analysis",name)
            ref.insert(0,"lineage",parent)
            reference.append(ref)
    pd.concat(reference,ignore_index=True).to_csv(OUT/"reference_readouts.csv",index=False)
    profile={"started_utc":started,"elapsed_seconds":time.monotonic()-start,
             "peak_parent_rss_bytes":resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,
             "peak_child_rss_bytes":resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss*1024,
             "pydeseq2":version("pydeseq2"),"python":sys.version,"command":" ".join([sys.executable,*sys.argv])}
    (OUT/"profile.json").write_text(json.dumps(profile,indent=2)+"\n")


if __name__=="__main__":
    run()
