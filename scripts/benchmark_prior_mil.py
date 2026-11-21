"""Direct HA, CTA and full Gaussian variational MixMIL comparators on M4 bags."""
import argparse
import json
import resource
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT/"src"),str(ROOT/"starter")]
from crlm.model_bags import LINEAGES, SEED, sample_bags
from crlm.prior_mil import BagBatch, HierarchicalMIL, MixedMIL, official_posterior, train_fixed
from linear_attention_mil import build_design, fit_model, fit_scaler, predict_bags, save_model

OUT = ROOT/"artifacts/m4/direct_prior"
MODELS = ROOT/"models/m4_direct"
CACHE = Path.home()/".cache/crlm-m4-priors"
MIXMIL_REV = "16d92349ccf19b73c1dccd167934fcc3202b6be4"


def load_inputs(rep, rows):
    published=pd.read_csv(ROOT/"artifacts/m4/selected_cells.csv.gz",index_col="cell_uid")
    groups={k:{} for k in LINEAGES}
    for gsm in rows.sample_uid:
        bags=sample_bags(ROOT/f"data/processed/m2/{gsm}.h5ad")
        cols=None
        if rep.startswith("scGPT"):
            with np.load(ROOT/f"embeddings/m4/scgpt/{gsm}.npz") as f:
                cols=pd.Index(f["cell_uid"])
                enc=f["embedding"]
        for lineage,bag in bags.items():
            expected=published[(published.sample_uid==gsm)&(published.lineage_m2==lineage)]
            if not bag.obs.index.equals(expected.index):
                raise ValueError(f"{gsm} {lineage}: cell selection differs from M4")
            if rep=="PCA32":
                x=bag.counts.astype(np.float32)
                x.data*=np.repeat(10000/np.asarray(x.sum(axis=1)).ravel(),np.diff(x.indptr))
                np.log1p(x.data,out=x.data)
            else:
                ix=cols.get_indexer(bag.obs.index)
                if (ix<0).any() or not cols.is_unique:
                    raise ValueError(f"{gsm}: missing or ambiguous scGPT cell identities")
                x=enc[ix]
            groups[lineage][gsm]=(x,bag.genes,bag.obs)
    return groups


def project(group, design, rep, lineage, fold, destination):
    if rep=="PCA32":
        src=ROOT/f"models/m4/{lineage.replace('/','_')}/fold{fold}/pca.npz"
        with np.load(src) as f:
            if set(f["training_samples"])!=set(design.samples):
                raise ValueError("Saved PCA has a different training fold")
            genes,center,components=f["genes"],f["center"],f["components"]
        matrix={}
        for gsm,(x,symbols,_) in group.items():
            indices=symbols.get_indexer(genes)
            if (indices<0).any():
                raise ValueError("Missing a fitted PCA feature")
            matrix[gsm]=np.asarray((x[:,indices].toarray()-center)@components.T,dtype=np.float32)
    elif rep=="scGPT32":
        xs=np.concatenate([group[s][0] for s in design.samples]).astype(np.float64)
        center=xs.mean(axis=0)
        xs-=center
        _,vectors=np.linalg.eigh(xs.T@xs)
        components=vectors[:,-32:][:,::-1].T
        # Fix arbitrary eigenvector signs for a portable projection.
        signs=np.sign(components[np.arange(len(components)),np.abs(components).argmax(axis=1)])
        components*=signs[:,None]
        matrix={k:np.asarray((x-center)@components.T,dtype=np.float32) for k,(x,_,_) in group.items()}
        np.savez_compressed(destination/"projection.npz", center=center, components=components,
                            training_samples=np.asarray(design.samples,dtype=str))
    else:
        matrix={s:x for s,(x,_,_) in group.items()}
    loc,scale=fit_scaler([matrix[s].astype(float) for s in design.samples], design.bce_weights, True)
    np.savez_compressed(destination/"scaler.npz",center=loc,scale=scale,
                        training_samples=np.asarray(design.samples,dtype=str))
    return {s:np.asarray((a-loc)/scale,dtype=np.float32) for s,a in matrix.items()}


def run(rep, repo):
    start=time.perf_counter()
    source_rev=subprocess.check_output(["git","-C",str(repo),"rev-parse","HEAD"],text=True).strip()
    if source_rev!=MIXMIL_REV:
        raise ValueError(f"Expected MixMIL checkout {MIXMIL_REV}, found {source_rev}")
    posterior=official_posterior(repo)
    out=OUT/rep
    out.mkdir(parents=True,exist_ok=True)
    rows=pd.read_csv(ROOT/"data/manifests/m0/sample_manifest_validated.csv").query("role=='development' and tissue in ['PT','LM']")
    folds=json.loads((ROOT/"data/manifests/m0/donor_folds_pre_qc.json").read_text())["folds"]
    groups=load_inputs(rep,rows)
    torch.set_num_threads(2)
    runs,predictions,cell_parts=[],[],[]
    for lineage,group in groups.items():
        ready=rows[rows.sample_uid.isin(group)]
        for split in folds:
            fold=split["fold"]
            train=ready[ready.donor_uid.isin(split["train_donors"])].reset_index(drop=True)
            test=ready[ready.donor_uid.isin(split["test_donors"])].reset_index(drop=True)
            if set(train.donor_uid)&set(test.donor_uid) or set(train.sample_uid)&set(test.sample_uid):
                raise ValueError("Patient or sample leakage")
            design=build_design(train)
            if not design.pair_donors or test.empty:
                raise ValueError("A study fold lacks training pairs or test tissue")
            dest=MODELS/rep/lineage.replace("/","_")/f"fold{fold}"
            dest.mkdir(parents=True,exist_ok=True)
            h=project(group,design,rep,lineage,fold,dest)
            ids=design.samples
            test_ids=test.sample_uid.tolist()
            trmat=[h[s] for s in ids]
            train_batch=BagBatch(trmat,[group[s][2].subtype_m2.astype(str).to_numpy() for s in ids])
            test_batch=BagBatch([h[s] for s in test_ids],[group[s][2].subtype_m2.astype(str).to_numpy() for s in test_ids])
            kinds=["CTA","HA"]
            if rep!="scGPT512":
                kinds.append("MixMIL")
            if rep=="scGPT32":
                kinds+=["mean","attention"]
            for kind in kinds:
                began=time.perf_counter()
                method=f"{rep}_{kind}"
                common={"method":method,"lineage":lineage,"fold":fold,"seed":SEED,"n_train_bags":len(train),
                        "n_train_donors":train.donor_uid.nunique(),"n_test_bags":len(test),"configs_tried":1}
                if kind in ("mean","attention"):
                    model,fit=fit_model([x.astype(float) for x in trmat],design,pooling=kind,lambda_pair=.1,
                                        lambda_l2=.01,seed=SEED,standardize=False,maxiter=500)
                    meta={"fit":fit,"lineage":lineage,"training_samples":ids,"training_donors":sorted(train.donor_uid.unique())}
                    save_model(dest/f"{kind}.npz",model,meta)
                    pred,_=predict_bags(model,meta,[h[s].astype(float) for s in test_ids],test_ids)
                    logits=pred.logit.to_numpy()
                    n_params=fit["n_parameters"]
                    steps=fit["iterations"]
                    status="converged" if fit["converged"] else fit["message"]
                    obj=fit["objective"] if "objective" in fit else None
                else:
                    torch.manual_seed(SEED)
                    width=trmat[0].shape[1]
                    model=(MixedMIL(trmat,design.labels,posterior) if kind=="MixMIL" else HierarchicalMIL(width,kind=="HA"))
                    steps=500 if kind=="MixMIL" else 100
                    obj=train_fixed(model,train_batch,design,train.donor_uid.nunique(),steps)
                    n_params=sum(p.numel() for p in model.parameters() if p.requires_grad)
                    torch.save({"state_dict":model.state_dict(),"width":width,"method":method,"seed":SEED,
                                "training_samples":ids,"training_donors":sorted(train.donor_uid.unique())},
                               dest/f"{kind}.pt")
                    with torch.inference_mode():
                        lin,att,signed,bias=model(test_batch)
                        logits=lin.numpy()
                        # Every specimen can be scored alone using the frozen training parameters.
                        first=test_ids[0]
                        single=BagBatch([h[first]],[group[first][2].subtype_m2.astype(str).to_numpy()])
                        np.testing.assert_allclose(model(single)[0].numpy(),logits[:1],rtol=2e-5,atol=2e-5)
                        at=att.numpy()
                        si=signed.numpy()
                        owner=test_batch.owner.numpy()
                        cells=pd.concat([group[s][2][["sample_uid","subtype_m2"]] for s in test_ids]).copy()
                        np.testing.assert_allclose(np.bincount(owner,weights=at),np.ones(len(test_ids)),atol=2e-5)
                        np.testing.assert_allclose(np.bincount(owner,weights=si)+bias.item(),logits,rtol=2e-5,atol=2e-5)
                        cell_parts.append(cells.assign(method=method,lineage=lineage,fold=fold,
                                                       attention=at,signed_logit_contribution=si,bias=bias.item()))
                    status="completed_fixed_steps"
                if not np.isfinite(logits).all():
                    raise FloatingPointError(f"{method}: non-finite OOF logit")
                runs.append(dict(common,steps=steps,n_parameters=n_params,seconds=time.perf_counter()-began,status=status,objective=obj))
                if status in ("converged","completed_fixed_steps"):
                    predictions.append(test[["sample_uid","donor_uid","lesion_uid","tissue"]].assign(
                        method=method,lineage=lineage,fold=fold,seed=SEED,logit=logits))
                print(f"{rep} {lineage} fold={fold} {kind}: {status} ({time.perf_counter()-began:.2f}s)",flush=True)
            pd.DataFrame(runs).to_csv(out/"benchmark_runs.csv",index=False)
            pd.concat(predictions,ignore_index=True).to_csv(out/"oof_predictions.csv",index=False)
    pd.concat(cell_parts).to_csv(out/"cell_contributions.csv.gz",compression={"method":"gzip","mtime":0})
    (out/"profile.json").write_text(json.dumps({"mixmil_commit":source_rev,"torch":torch.__version__,
        "torch_threads":torch.get_num_threads(),"seconds":time.perf_counter()-start,
        "peak_rss_bytes":resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,
        "command":" ".join([sys.executable,*sys.argv])},indent=2)+"\n")


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--representation",choices=["PCA32","scGPT32","scGPT512"],required=True)
    parser.add_argument("--mixmil-repo",type=Path,default=CACHE/"mixmil")
    args=parser.parse_args()
    with threadpool_limits(limits=2):
        run(args.representation,args.mixmil_repo)
