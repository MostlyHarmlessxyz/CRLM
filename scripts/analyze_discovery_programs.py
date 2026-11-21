"""Fixed cell programs, lesion-weighted patient effects and composition in GSE245552."""
import argparse
import hashlib
import itertools
import json
import os
import resource
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src"), str(ROOT / "starter")]
from crlm.program_scores import rank_scores
from paired_effects import calculate_patient_contrasts, summarize_patients
from composition_decomposition import decompose_pair
from _common import stable_seed

OUT = ROOT / "artifacts/m3"
DEFS = OUT / "program_definitions_draft"
MANIFEST = ROOT / "data/manifests/m0/sample_manifest_validated.csv"
M2 = ROOT / "data/processed/m2"
COMMON = M2 / "pseudobulk/gene_symbols.txt"
M1_AXIS = ROOT / "starter/audit_data/GSE245552_feature_intersection.txt"
SEEDS = ROOT / "starter/literature_seeds/literature_program_seeds.csv"
BACK = SEEDS.parent / "literature_program_sources/B06_S6_background_sets.csv"
HALLMARK = DEFS / "h.all.v2025.1.Hs.symbols.gmt"
SEED = 20261121
PARENTS = ("T/NK", "Mono/Mac/DC", "B/Plasma")
LABELS = PARENTS + ("Mast", "Neutrophil", "Mixed", "Unknown", "Nonimmune")
HALLMARK_LINEAGE = {"HALLMARK_INTERFERON_GAMMA_RESPONSE": "T/NK",
                    "HALLMARK_UNFOLDED_PROTEIN_RESPONSE": "B/Plasma"}
LIVER_RBC = {"ALB", "APOA1", "APOA2", "TTR", "FGA", "FGB", "FGG", "HP", "HBA1", "HBA2", "HBB"}
SCORE_COLUMNS = ["sample_uid", "lineage", "program", "score", "n_cells"]


def sha(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def definitions():
    """Return all source membership rows, including missing original symbols."""
    axis = M1_AXIS.read_text().splitlines()
    if COMMON.exists() and set(COMMON.read_text().splitlines()) != set(axis):
        raise ValueError("M2 common gene axis disagrees with M1")
    src = pd.read_csv(SEEDS)
    rows = []
    for row in src.itertuples(index=False):
        if row.seed_set_id.startswith("B19_S4_"):
            role = "complete_program"
        elif row.seed_set_id.startswith("B06_S5_"):
            role = "pathway_overlap"
        else:
            continue
        rows.append(("Mono/Mac/DC", row.seed_set_id, row.gene_symbol_original, role, row.source_doi, row.source_location))
    for row in pd.read_csv(BACK).itertuples(index=False):
        rows.append(("Mono/Mac/DC", row.set_id, row.gene_symbol_original, "noncancer_background",
                     row.source_doi, row.source_location))
    for line in HALLMARK.read_text().splitlines():
        prog, url, *members = line.rstrip().split("\t")
        if prog in HALLMARK_LINEAGE:
            rows.extend((HALLMARK_LINEAGE[prog], prog, member, "complete_program",
                         "10.1016/j.cels.2015.12.004", f"MSigDB 2025.1.Hs {url}") for member in members)
    genes = pd.DataFrame(rows, columns=["lineage","program","gene_symbol","role","source_doi","source_location"])
    if set(HALLMARK_LINEAGE) - set(genes.program) or genes.duplicated(["program","gene_symbol"]).any():
        raise ValueError("Incomplete or duplicated original program definitions")
    genes["on_discovery_axis"] = genes.gene_symbol.isin(axis)
    genes.to_csv(DEFS / "program_genes.csv", index=False)
    (DEFS / "ranking_background.txt").write_text("\n".join(axis) + "\n")
    return genes, axis


def program_variants(genes):
    info, sets, overlaps = [], {}, []
    base = {k: set(v.gene_symbol) for k,v in genes[genes.on_discovery_axis].groupby("program")}
    role = genes.groupby("program").role.first().to_dict()
    for a,b in itertools.combinations(base, 2):
        shared = base[a] & base[b]
        overlaps.append({"program_a": a, "program_b": b, "n_shared": len(shared),
                         "jaccard": len(shared)/len(base[a] | base[b]), "shared_genes": " ".join(sorted(shared))})
    for name, group in genes.groupby("program",sort=True):
        members = base.get(name,set())
        lineage = group.lineage.iloc[0]
        available = len(members) >= (5 if role[name]=="complete_program" else 3) and len(members) >= .8*len(group)
        def add(variant, selected):
            key = name if variant == "original" else f"{name}__{variant}"
            info.append({"score_id":key,"program":name,"lineage":lineage,"role":role[name],"variant":variant,
                         "n_scored_genes":len(selected),"eligible":bool(selected) and available})
            if selected and available:
                sets[key] = sorted(selected)
        add("original", members)
        if role[name] == "noncancer_background":
            continue
        shared = members & set().union(*(g for p,g in base.items() if p!=name and role[p]!="noncancer_background"))
        add("unique",members-shared)
        add("shared",shared)
        add("remove_liver_rbc_HLA_DRB345",members-LIVER_RBC-{"HLA-DRB3","HLA-DRB4","HLA-DRB5"})
    pd.DataFrame(overlaps).to_csv(DEFS/"gene_overlap.csv", index=False)
    return pd.DataFrame(info), sets


def matched_controls(manifest, axis, genes, info, sets):
    """Estimate detection on an equal-sized (at most 200) lineage cell sample per GSM."""
    rng = np.random.default_rng(SEED)
    detection = {k:np.zeros(len(axis),dtype=float) for k in PARENTS}
    libraries = {k:0 for k in PARENTS}
    for gsm in manifest.sample_uid:
        d = ad.read_h5ad(M2/f"{gsm}.h5ad")
        cols = pd.Index(d.var.gene_symbol).get_indexer(axis)
        if (cols < 0).any():
            raise ValueError(f"{gsm}: absent ranking background")
        for lineage in PARENTS:
            ix = np.flatnonzero(d.obs.qc_retained & (d.obs.lineage_m2==lineage))
            if not len(ix):
                continue
            sel = np.sort(rng.choice(ix,min(200,len(ix)),replace=False))
            prob = np.asarray((d.X[sel,:][:,cols]>0).sum(axis=0)).ravel()/len(sel)
            detection[lineage] += prob
            libraries[lineage] += 1
    axis_arr = np.asarray(axis)
    controls, rows = [], []
    for row in info[(info.role=="complete_program") & (info.variant=="original") & info.eligible].itertuples():
        candidates = ~np.isin(axis_arr, genes.gene_symbol.unique())
        frequency = detection[row.lineage] / libraries[row.lineage]
        # Use nearest-rank bins to match sparse genes, breaking frequency ties by symbol order.
        order = np.lexsort((axis_arr, frequency))
        bins = np.empty(len(axis),int)
        bins[order] = np.arange(len(axis))*30//len(axis)
        members = sets[row.score_id]
        members_at = pd.Index(axis).get_indexer(members)
        for draw in range(1,11):
            selected = []
            for bin_id, n_genes in zip(*np.unique(bins[members_at],return_counts=True)):
                pool = np.flatnonzero(candidates & (bins==bin_id))
                selected.extend(axis_arr[rng.choice(pool,int(n_genes),replace=False)].tolist())
            name = f"{row.program}__matched{draw:02d}"
            sets[name] = sorted(selected)
            controls.append(dict(score_id=name,program=row.program,lineage=row.lineage,role=row.role,
                                 variant=f"matched{draw:02d}",n_scored_genes=len(selected),eligible=True))
            rows.extend((name,gene) for gene in sorted(selected))
    pd.DataFrame(rows,columns=["score_id","gene_symbol"]).to_csv(DEFS/"matched_control_genes.csv",index=False)
    return pd.concat([info,pd.DataFrame(controls)],ignore_index=True)


def score_libraries(manifest, axis, info, sets):
    sample, comp, subtype, tech, precision = [], [], [], [], []
    for gsm in manifest.sample_uid:
        d = ad.read_h5ad(M2/f"{gsm}.h5ad")
        if set(d.obs.sample_uid)!= {gsm} or not set(d.obs.lineage_m2).issubset(LABELS):
            raise ValueError(f"Unexpected M2 identity or annotation in {gsm}")
        keep = d.obs.qc_retained.to_numpy(dtype=bool)
        labels = d.obs.lineage_m2.to_numpy(dtype=str)
        cols = pd.Index(d.var.gene_symbol).get_indexer(axis)
        if (cols < 0).any():
            raise ValueError(f"{gsm} missing a common feature")
        for lineage in LABELS:
            n = int((keep & (labels==lineage)).sum())
            comp.append((gsm,lineage,n,int(keep.sum()),n/keep.sum()))
        for lineage in PARENTS:
            idx = np.flatnonzero((labels==lineage) & (d.obs.n_counts_visible.to_numpy()>0))
            specs = info[info.lineage==lineage]
            selected = {p:sets[p] for p in specs.loc[specs.eligible,"score_id"]}
            scores = rank_scores(d.X[idx,:][:,cols],axis,selected)
            kept = keep[idx]
            strict = kept & ~d.obs.qc_joint_low.to_numpy(dtype=bool)[idx]
            surviving = scores.loc[kept]
            subtype_labels = d.obs.subtype_m2.to_numpy(dtype=str)[idx][kept]
            cells = d.obs.iloc[idx].loc[kept]
            for item in specs.itertuples(index=False):
                if item.score_id in scores and len(surviving):
                    mean = float(surviving[item.score_id].mean())
                    nonempty_mean = float(scores[item.score_id].mean())
                else:
                    mean = np.nan
                    nonempty_mean = float(scores[item.score_id].mean()) if item.score_id in scores else np.nan
                strict_mean=scores.loc[strict,item.score_id].mean() if item.score_id in scores else np.nan
                sample.append((gsm,lineage,item.score_id,mean,len(surviving),nonempty_mean,len(scores),strict_mean,int(strict.sum())))
                if item.variant!="original" or item.role!="complete_program" or not len(surviving):
                    continue
                obs_scores = surviving[item.score_id].to_numpy()
                if len(obs_scores)>=2:
                    order=np.random.default_rng(stable_seed(SEED,gsm,lineage)).permutation(len(obs_scores))
                    half=len(order)//2
                    first,second=obs_scores[order[:half]].mean(),obs_scores[order[half:]].mean()
                    precision.append((gsm,lineage,item.program,len(obs_scores),first,second,abs(first-second)))
                for label in sorted(set(subtype_labels)):
                    mask = subtype_labels==label
                    subtype.append((gsm,lineage,item.program,label,int(mask.sum()),len(surviving),float(obs_scores[mask].mean())))
                for var in ["n_counts_visible","n_genes_visible"]:
                    vals = cells[var].to_numpy()
                    corr = spearmanr(vals,obs_scores).statistic if np.unique(vals).size>1 and np.unique(obs_scores).size>1 else np.nan
                    tech.append((gsm,lineage,item.program,var,float(np.median(vals)),float(corr),len(vals)))
        print(f"{gsm} {keep.sum()} QC-retained cells; completed program scoring",flush=True)
    return (pd.DataFrame(sample,columns=["sample_uid","lineage","program","unfiltered_score","n_cells","nonempty_score","n_nonempty","strict_score","n_strict"]),
            pd.DataFrame(comp,columns=["sample_uid","lineage","n_cells","n_retained","proportion"]),
            pd.DataFrame(subtype,columns=["sample_uid","lineage","program","subtype","n_subtype","n_cells","mean_score"]),
            pd.DataFrame(tech,columns=["sample_uid","lineage","program","qc_metric","median_qc","cell_spearman","n_cells"]),
            pd.DataFrame(precision,columns=["sample_uid","lineage","program","n_cells","half_a","half_b","abs_gap"]))


def patient_outputs(manifest, scores, where, *, summary=True):
    patients, lesions, tissue, audit = calculate_patient_contrasts(manifest, scores[SCORE_COLUMNS])
    where.mkdir(parents=True,exist_ok=True)
    patients.to_csv(where/"patient_contrasts.csv",index=False)
    tissue.to_csv(where/"patient_tissue_means.csv",index=False)
    if summary:
        lesions.to_csv(where/"lesion_means.csv",index=False)
        result = summarize_patients(patients,seed=SEED)
        result.to_csv(where/"effects_summary.csv",index=False)
    else:
        result = None
    return patients,result


def leave_one_out(patients):
    rows = []
    for (lineage,program), group in patients.groupby(["lineage","program"]):
        available = group.dropna(subset="D_T")
        for patient in available.donor_uid:
            others = available[available.donor_uid!=patient].D_T
            rows.append((lineage,program,patient,len(others),others.mean()))
    return pd.DataFrame(rows,columns=["lineage","program","omitted_patient","n_remaining","mean_D_T"])


def decompose(manifest, subtypes):
    usable = subtypes[subtypes.n_cells>=50].merge(
        manifest[["sample_uid","donor_uid","lesion_uid","tissue"]],on="sample_uid",validate="many_to_one")
    usable = usable[usable.tissue.isin(["PT","LM"])].copy()
    # Each source row covers a unique GSM/lesion; a missing subtype has zero mass.
    if manifest.duplicated("lesion_uid").any():
        raise ValueError("Separate technical captures must first be averaged within each lesion")
    usable["proportion"] = usable.n_subtype / usable.n_cells
    usable["score_mass"] = usable.proportion * usable.mean_score
    rows = []
    keys = ["donor_uid","lineage","program"]
    for grouping, group in usable.groupby(keys):
        frame = []
        observed = {}
        names = sorted(group.subtype.unique())
        for tissue, tissue_rows in group.groupby("tissue"):
            n_lesions = tissue_rows.lesion_uid.nunique()
            aggregated = tissue_rows.groupby("subtype")[["n_subtype","proportion","score_mass"]].sum().reindex(names,fill_value=0)
            p = aggregated.proportion/n_lesions
            mu = (aggregated.score_mass/aggregated.proportion).where(p>0)
            observed[tissue] = float(aggregated.score_mass.sum()/n_lesions)
            frame.extend(dict(tissue=tissue,subtype=s,proportion=p.at[s],mean_score=mu.at[s],n_subtype=aggregated.n_subtype.at[s]) for s in names)
        values = pd.DataFrame(frame)
        for support in (1,20):
            checked = values.copy()
            checked.loc[checked.n_subtype<support,"mean_score"] = np.nan
            pair = decompose_pair(checked)
            common = pair.get("common_support") or {}
            full = pair.get("full_decomposition") or {}
            cond = common.get("renormalized_components",{})
            actual = observed.get("LM",np.nan)-observed.get("PT",np.nan)
            if full and (abs(full["identity_residual"])>1e-10 or
                         not np.isclose(actual,full["total_difference"],rtol=0,atol=1e-10)):
                raise ValueError(f"Marginal mean decomposition failed for {grouping}")
            rows.append(dict(zip(keys,grouping)) | {"min_subtype_cells_per_tissue":support,
                "status":pair["status"],"observed_total_D_T":actual,"full_composition":full.get("composition"),
                "full_within_subtype_state":full.get("within_type_state"),
                "common_mass_PT":common.get("retained_mass_PT"),"common_mass_LM":common.get("retained_mass_LM"),
                "common_total_D_T":cond.get("total_difference"),"common_composition":cond.get("composition"),
                "common_within_subtype_state":cond.get("within_type_state"),
                "common_subtypes":" | ".join(common.get("subtypes",[]))})
    return pd.DataFrame(rows)


def state_summaries(manifest, raw, info, where):
    rows = []
    base_ids = set(info.loc[info.variant=="original","score_id"])
    for name, threshold, value, ncol, mask in (
        ("main",50,"unfiltered_score","n_cells",raw.program.isin(base_ids)),
        ("n20",20,"unfiltered_score","n_cells",raw.program.isin(base_ids)),
        ("n100",100,"unfiltered_score","n_cells",raw.program.isin(base_ids)),
        ("before_lower_tail_qc",50,"nonempty_score","n_nonempty",raw.program.isin(base_ids)),
        ("strict_library_qc",50,"strict_score","n_strict",raw.program.isin(base_ids)),
        ("gene_sensitivity",50,"unfiltered_score","n_cells",~raw.program.isin(base_ids)),
    ):
        current = raw[mask].copy()
        current["score"] = current[value].where(current[ncol]>=threshold)
        current["n_cells"] = current[ncol]
        if name=="main":
            current[SCORE_COLUMNS].to_csv(OUT/"sample_program_scores.csv",index=False)
        patients,summary = patient_outputs(manifest,current,where if name=="main" else where/name,summary=name=="main")
        for (lineage,program), group in patients.groupby(["lineage","program"]):
            val=group.D_T.dropna()
            rows.append((name,lineage,program,len(val),val.mean(),val.median(),int((val>0).sum()),int((val<0).sum())))
        if name=="main":
            loo=leave_one_out(patients)
            loo.to_csv(where/"leave_one_patient_out.csv",index=False)
    return pd.DataFrame(rows,columns=["analysis","lineage","program","n_pairs","mean_D_T","median_D_T","positive","negative"]),loo


def shortlist(info, sensitivity, loo):
    outcomes=[]
    for r in info[(info.role=="complete_program") & (info.variant=="original")].itertuples():
        metrics=sensitivity[(sensitivity.program==r.program) & (sensitivity.lineage==r.lineage)].set_index("analysis")
        main=metrics.loc["main"]
        d=np.sign(main.mean_D_T)
        concordant=main.positive if d==1 else main.negative if d==-1 else 0
        n_nonzero=main.positive+main.negative
        share=concordant/n_nonzero if n_nonzero else 0.
        influence=loo[(loo.program==r.program)&(loo.lineage==r.lineage)].mean_D_T
        reasons=[]
        if main.n_pairs<8: reasons.append("fewer_than_8_pairs")
        if not n_nonzero or share<2/3: reasons.append("direction_support_below_two_thirds")
        if d==0 or len(influence)<2 or not (np.sign(influence)==d).all(): reasons.append("leave_one_out_direction")
        for subset in ("n20","n100"):
            if metrics.at[subset,"n_pairs"]>=4 and np.sign(metrics.at[subset,"mean_D_T"])!=d:
                reasons.append(f"{subset}_direction")
        outcomes.append({"lineage":r.lineage,"program":r.program,"n_pairs":int(main.n_pairs),
                         "mean_D_T":main.mean_D_T,"median_D_T":main.median_D_T,"direction":d,
                         "n_nonzero":int(n_nonzero),"direction_fraction":share,"loo_min":influence.min(),
                         "loo_max":influence.max(),"screen_pass":not reasons,"screen_reason":" | ".join(reasons)})
    result=pd.DataFrame(outcomes).sort_values(["screen_pass","direction_fraction","n_pairs","program"],ascending=[False,False,False,True])
    result["shortlisted"]=result.screen_pass & (np.arange(len(result))<3)
    return result


def run():
    start=time.monotonic()
    started=datetime.now(timezone.utc).isoformat()
    manifest = pd.read_csv(MANIFEST)
    if len(manifest)!=39 or set(manifest.dataset)!={"GSE245552"} or set(manifest.role)!={"development"}:
        raise ValueError("Expected the validated GSE245552 discovery manifest")
    if not (ROOT/"artifacts/m2/pseudobulk_manifest.csv").exists():
        raise FileNotFoundError("M2 has not published the complete pseudobulk manifest")
    genes,axis=definitions()
    info,sets=program_variants(genes)
    info=matched_controls(manifest,axis,genes,info,sets)
    info.to_csv(DEFS/"scored_definitions.csv",index=False)
    raw,comp,subtype,tech,precision=score_libraries(manifest,axis,info,sets)
    composition=OUT/"composition_results"
    composition.mkdir(exist_ok=True)
    raw.to_csv(OUT/"sample_program_measurements.csv.gz",index=False)
    comp.to_csv(composition/"sample_proportions.csv",index=False)
    tech.to_csv(OUT/"depth_sensitivity.csv",index=False)
    precision.to_csv(OUT/"split_half_precision.csv",index=False)
    subtype.to_csv(composition/"subtype_contributions.csv",index=False)
    decompose(manifest,subtype).to_csv(composition/"patient_subtype_decomposition.csv",index=False)
    comp_inputs=comp.assign(program="captured_cell_fraction",score=comp.proportion)
    patient_outputs(manifest,comp_inputs,composition)
    where=OUT/"patient_effects"
    sensitivity,loo=state_summaries(manifest,raw,info,where)
    sensitivity.to_csv(where/"sensitivity_summary.csv",index=False)
    result=shortlist(info,sensitivity,loo)
    patient_rows=pd.read_csv(where/"patient_contrasts.csv")
    depth_rows=[]
    for metric, statistical in tech[tech.n_cells>=50].groupby("qc_metric"):
        qc=statistical.merge(manifest[["sample_uid","donor_uid","tissue","lesion_uid"]],on="sample_uid",validate="many_to_one")
        qc=qc.groupby(["donor_uid","lineage","program","tissue"])["median_qc"].mean().unstack("tissue")
        delta=(qc.LM-qc.PT).rename("depth_D_T").reset_index()
        joined=patient_rows.merge(delta,on=["donor_uid","lineage","program"]).dropna(subset=["depth_D_T","D_T"])
        for (lineage,program),g in joined.groupby(["lineage","program"]):
            rho=spearmanr(g.depth_D_T,g.D_T).statistic if min(g.depth_D_T.nunique(),g.D_T.nunique())>1 else np.nan
            depth_rows.append({"qc_metric":metric,"lineage":lineage,"program":program,"n_pairs":len(g),"spearman_depth_delta_program_delta":rho})
    pd.DataFrame(depth_rows).to_csv(OUT/"paired_depth_sensitivity.csv",index=False)
    result.to_csv(OUT/"candidate_screen.csv",index=False)
    report=["# 发现集候选草案", "",
        "本表由 M3 登记的四个完整文献程序产生。患者内差为 LM 减 PT，单位为 0–1 秩分数。", "",
        "| 谱系 | 程序 | 完整患者对 | 平均差 | 同向患者 | 留一范围 | 短名单 |",
        "|---|---|---:|---:|---:|---|---|"]
    reason_labels={"fewer_than_8_pairs":"配对人数低于 8",
                   "direction_support_below_two_thirds":"同向患者少于三分之二",
                   "leave_one_out_direction":"留一患者后方向变化",
                   "n20_direction":"20 细胞门槛方向变化","n100_direction":"100 细胞门槛方向变化"}
    for r in result.itertuples(index=False):
        conclusion="入选" if r.shortlisted else "、".join(reason_labels[k] for k in r.screen_reason.split(" | ") if k) or "排序靠后"
        report.append(f"| {r.lineage} | {r.program} | {r.n_pairs} | {r.mean_D_T:.5f} | "
                      f"{round(r.n_nonzero*r.direction_fraction)}/{r.n_nonzero} | "
                      f"{r.loo_min:.5f} 至 {r.loo_max:.5f} | {conclusion} |")
    report.extend(["","每个候选的门槛 20/100、宽松 QC、去共享成员、去肝/血红蛋白及 DRB3/4/5、"
        "检测率匹配随机集合的结果见[敏感性结果](patient_effects/sensitivity_summary.csv)。"
        "同一父谱系内的[患者分解](composition_results/patient_subtype_decomposition.csv)列出共同支持质量。"
        "冻前方向按患者均差确定，模型增量归 M4，确认家族归 M5。",
        "", "测序深度、细胞组成和伪批量的独立分析见[阶段报告](stage_report.md)。"])
    (OUT/"candidate_evidence_draft.md").write_text("\n".join(report)+"\n")
    inputs=[Path(__file__),ROOT/"src/crlm/program_scores.py",MANIFEST,DEFS/"development_rules.md",
            SEEDS,BACK,HALLMARK,M1_AXIS,COMMON,ROOT/"starter/paired_effects.py",ROOT/"starter/composition_decomposition.py"]
    for gsm in manifest.sample_uid:
        inputs.append(M2/f"{gsm}.h5ad")
    profile={"started_utc":started,"ended_utc":datetime.now(timezone.utc).isoformat(),"seed":SEED,
             "command":" ".join([sys.executable,*sys.argv]),"inputs_sha256":{str(p.relative_to(ROOT)):sha(p) for p in inputs},
             "peak_rss_bytes":resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,
             "elapsed_seconds":time.monotonic()-start,"discovery_samples":len(manifest),
             "retained_cells":int(comp.n_cells.sum()),"development_patients":manifest.donor_uid.nunique(),
             "shortlisted":result.loc[result.shortlisted,["lineage","program","direction"]].to_dict("records")}
    (OUT/"resource_profile.json").write_text(json.dumps(profile,indent=2)+"\n")
    print(result[["lineage","program","n_pairs","mean_D_T","shortlisted"]].to_string(index=False),flush=True)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--definitions-only",action="store_true")
    args=parser.parse_args()
    if args.definitions_only:
        program_variants(definitions()[0])
    else:
        run()
