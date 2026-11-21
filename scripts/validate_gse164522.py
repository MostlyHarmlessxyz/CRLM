"""Prepare the four GSE164522 tissues and apply frozen discovery programs."""
import argparse
import csv
import gzip
import hashlib
import json
import resource
import shlex
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src"), str(ROOT / "starter")]
from crlm.program_scores import rank_scores
from crlm.model_bags import cell_seed
from adjust_registered_tests import adjust_registered
from paired_effects import calculate_patient_contrasts, summarize_patients
from analyze_discovery_programs import LIVER_RBC, decompose, leave_one_out

GSE = "GSE164522"
RAW = ROOT / "data/raw" / GSE
OUT = ROOT / "artifacts/m6" / GSE
OBJECTS = ROOT / "data/processed/m6" / GSE
GEO = ROOT / "starter/audit_data/geo"
METADATA = RAW / "GSE164522_CRLM_metadata.csv.gz"
FREEZE = ROOT / "artifacts/m5"
ORGANS = {"PT": "PT", "MT": "LM", "PN": "NC", "MN": "NL"}
MAJOR = {"Bcell": "B/Plasma", "Monocyte": "Mono/Mac/DC", "Macrophage": "Mono/Mac/DC", "DC": "Mono/Mac/DC",
         "CD4 T": "T/NK", "CD8 T": "T/NK", "NK cell": "T/NK", "gdT": "T/NK", "NKT": "T/NK",
         "MAST": "Mast", "ILC3": "ILC3", "CD45-": "CD45-"}
EXPECTED_GLOBAL = {"T/NK": "Lymphoid-T-NK", "B/Plasma": "Lymphoid-B", "Mono/Mac/DC": "Myeloid"}
SCORE_COLUMNS = ["sample_uid", "lineage", "program", "score", "n_cells"]
SEED = 20261121


def sha(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def now():
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def sources():
    audit = json.loads((GEO / "GSE164522_all_expression_headers_audit.json").read_text())
    return {r["tissue_code"]: r for r in audit["files"] if r["tissue_code"] in ORGANS}


def range_file(job):
    url, first, last, size, dest = job
    length = last-first+1
    if dest.exists() and dest.stat().st_size == length:
        return
    headers = dest.with_suffix(".headers")
    tmp = dest.with_suffix(".partial")
    subprocess.run(["curl","-fsSL","--retry","5","--retry-all-errors","--connect-timeout","12",
                    "--speed-limit","1024","--speed-time","60","--max-time","360",
                    "--range",f"{first}-{last}","-D",str(headers),"-o",str(tmp),url],check=True)
    expected = f"content-range: bytes {first}-{last}/{size}"
    if expected not in headers.read_text().lower().splitlines() or tmp.stat().st_size != length:
        raise ValueError(f"GEO range reply differs for {dest.name}")
    tmp.replace(dest)
    headers.unlink()


def fetch(records):
    RAW.mkdir(parents=True, exist_ok=True)
    if not METADATA.exists():
        tmp = METADATA.with_suffix(".gz.part")
        base = next(iter(records.values()))["url"].rsplit("/",1)[0]
        subprocess.run(["curl","-fsSL","--retry","4","--retry-all-errors","-o",str(tmp),base+"/"+METADATA.name],check=True)
        if sha(tmp)!=next(iter(records.values()))["metadata_file_sha256"]:
            raise ValueError("Downloaded metadata differs from archived GEO metadata")
        tmp.replace(METADATA)
    ranges=RAW/"ranges"
    ranges.mkdir(exist_ok=True)
    jobs, files = [], {}
    for code, r in records.items():
        size=int(r["http_content_length"])
        path=RAW/r["url"].rsplit("/",1)[1]
        if path.exists() and path.stat().st_size==size:
            continue
        prefix=Path(str(path)+".part")
        start=prefix.stat().st_size if prefix.exists() else 0
        if start>size:
            raise ValueError(f"Invalid prefix for {code}")
        block=8<<20
        parts=[]
        for left in range(start,size,block):
            right=min(left+block,size)-1
            part=ranges/f"{code}_{left}_{right}.bytes"
            parts.append(part)
            jobs.append((code,(r["url"],left,right,size,part)))
        files[code]=(prefix,path,parts,size)
    def join(code):
        prefix,path,parts,size=files[code]
        if all(p.is_file() for p in parts) and not path.is_file():
            with prefix.open("ab") as out:
                for part in parts:
                    with part.open("rb") as stream:
                        shutil.copyfileobj(stream,out,length=1<<20)
            if prefix.stat().st_size!=size:
                raise ValueError(f"Invalid assembled size for {code}")
            prefix.replace(path)
            for p in parts:
                p.unlink()
            print(f"{now()} {code} complete {size} bytes",flush=True)
    with ThreadPoolExecutor(max_workers=8) as executor:
        futures={executor.submit(range_file,job):code for code,job in jobs}
        for fut in as_completed(futures):
            fut.result()
            join(futures[fut])
    for code in files:
        join(code)
    summary=[]
    for code,r in records.items():
        file=RAW/r["url"].rsplit("/",1)[1]
        summary.append(dict(tissue=ORGANS[code],path=str(file.relative_to(ROOT)),url=r["url"],
                            bytes=file.stat().st_size,sha256=sha(file)))
    summary.append(dict(tissue="metadata",path=str(METADATA.relative_to(ROOT)),
                        url=next(iter(records.values()))["url"].rsplit("/",1)[0]+"/"+METADATA.name,
                        bytes=METADATA.stat().st_size,sha256=sha(METADATA)))
    pd.DataFrame(summary).to_csv(OUT/"source_inputs.csv",index=False)


def metadata(records):
    expected = next(iter(records.values()))["metadata_file_sha256"]
    if sha(METADATA) != expected:
        raise ValueError("Complete author metadata differs from the GEO audit")
    frame = pd.read_csv(METADATA, index_col=0)
    frame.index.name = "source_barcode"
    identity = pd.read_csv(GEO / "GSE164522_metadata_identity_columns.csv.gz", index_col=0)
    if (not frame.index.is_unique or not frame.index.equals(identity.index) or
            not np.array_equal(frame[["patient", "tissue"]].to_numpy(), identity[["patient_code", "tissue"]].to_numpy())):
        raise ValueError("Cell/donor/tissue mapping differs from the archived original identity table")
    tissues = {r["metadata_tissue"]: ORGANS[code] for code, r in records.items()}
    frame = frame.loc[frame.tissue.isin(tissues)].copy()
    frame["author_tissue"] = frame.tissue
    frame["tissue"] = frame.author_tissue.map(tissues)
    frame["source_barcode"] = frame.index
    frame["expression_column"] = frame.index.str.replace("-", ".", regex=False)
    frame["lineage"] = frame.celltype_major.map(MAJOR)
    if not frame.expression_column.is_unique or frame.lineage.isna().any():
        raise ValueError("Unsupported major label or normalized barcode collision")
    frame["donor_uid"] = GSE + ":" + frame.patient
    frame["sample_uid"] = frame.donor_uid + ":" + frame["sample"].str.replace(" ", "_", regex=False)
    frame["lesion_uid"] = frame.sample_uid
    return frame.set_index("expression_column", drop=False)


def header(stream):
    return next(csv.reader([stream.readline().decode("utf-8-sig").rstrip("\r\n")]))[1:]


def write_identity(frame, records):
    cov = OUT / "coverage"
    cov.mkdir(parents=True, exist_ok=True)
    checks=[]
    for code, record in records.items():
        raw=RAW/record["url"].rsplit("/",1)[1]
        path=raw if raw.is_file() else GEO/f"GSE164522_{code}_expression_header.csv.gz"
        with gzip.open(path,"rb") as stream:
            line=stream.readline()
        if hashlib.sha256(line).hexdigest()!=record["header_sha256"]:
            raise ValueError(f"Archived and selected expression header differ for {code}")
        cols=pd.Index(next(csv.reader([line.decode("utf-8-sig").rstrip("\r\n")]))[1:])
        meta = frame.index[frame.tissue.eq(ORGANS[code])]
        if not cols.is_unique or len(meta) != len(cols) or not cols.isin(meta).all():
            raise ValueError(f"Header mismatch for {code}")
        checks.append(dict(code=code,tissue=ORGANS[code],n_header_cells=len(cols),n_matched=len(meta),
                           header_source=str(path.relative_to(ROOT)),header_sha256=record["header_sha256"]))
    pd.DataFrame(checks).to_csv(cov/"header_audit.csv",index=False)
    keys = ["donor_uid", "sample_uid", "lesion_uid", "patient", "tissue", "author_tissue", "sample", "ID"]
    manifest = frame.groupby(keys, observed=True).size().rename("input_cells").reset_index()
    if manifest.sample_uid.duplicated().any() or manifest.ID.duplicated().any():
        raise ValueError("Original author ID and patient/specimen mapping disagree")
    manifest["dataset"] = GSE
    manifest["role"] = "external"
    manifest["input_scale"] = "normalized"
    manifest["patient_mapping_verified"] = True
    manifest["treatment"] = "treatment-naive at cohort level in B02"
    manifest["episode"] = "source patient and specimen; exact collection dates unavailable"
    manifest.to_csv(cov / "sample_manifest.csv", index=False)
    frame[["source_barcode", "expression_column", *keys, "lineage", "celltype_global", "celltype_major", "celltype_sub"]].to_csv(
        cov / "cell_identity.csv.gz", index=False, compression={"method":"gzip", "mtime":0})
    frame[["celltype_global", "celltype_major", "lineage", "celltype_sub"]].drop_duplicates().sort_values(
        ["lineage","celltype_major","celltype_sub"]).to_csv(cov / "author_lineage_mapping.csv", index=False)
    return manifest


def parse_published(path, cells):
    """Read one numeric row at a time; only nonzero values enter CSR storage."""
    genes, blocks, positions, values, offsets = [], [], [], [], [0]
    def flush(n_cells):
        if values:
            block = sparse.csr_matrix((np.concatenate(values), np.concatenate(positions), np.asarray(offsets)),
                                      shape=(len(values), n_cells))
            blocks.append(block)
            values.clear()
            positions.clear()
            offsets[:] = [0]
    with gzip.open(path, "rb") as stream:
        actual = pd.Index(header(stream))
        if not actual.is_unique or len(cells) != len(actual) or not actual.isin(cells.index).all():
            raise ValueError(f"{path.name}: author cells and expression columns disagree")
        for line in stream:
            symbol, sep, numeric = line.partition(b",")
            if not sep:
                raise ValueError(f"Missing gene row in {path.name}")
            genes.append(symbol.strip(b'"').decode())
            vector = np.fromstring(numeric, dtype=np.float64, sep=",")
            if len(vector) != len(actual) or not np.isfinite(vector).all() or vector.min() < 0:
                raise ValueError(f"{path.name}: invalid published numeric row {len(genes)}")
            keep = np.flatnonzero(vector).astype(np.int32)
            positions.append(keep)
            values.append(vector[keep])
            offsets.append(offsets[-1] + len(keep))
            if len(values) == 256:
                flush(len(actual))
        flush(len(actual))
    if not genes or "" in genes:
        raise ValueError(f"{path.name}: missing published gene axis")
    return sparse.vstack(blocks, format="csr").T.tocsr(), genes, actual


def convert(code, frame, record):
    path = RAW / record["url"].rsplit("/", 1)[1]
    if path.stat().st_size != int(record["http_content_length"]):
        raise ValueError(f"{path.name}: incomplete GEO file")
    compressed_sha = sha(path)
    dest = OBJECTS / f"{code}.h5ad"
    if dest.exists() and (OBJECTS / f"{code}.json").exists():
        saved = ad.read_h5ad(dest, backed="r")
        same = saved.uns["expression_sha256"] == compressed_sha and saved.uns["metadata_sha256"] == record["metadata_file_sha256"]
        saved.file.close()
        if same:
            return
    start = time.perf_counter()
    subset = frame.loc[frame.tissue.eq(ORGANS[code])]
    matrix, genes, columns = parse_published(path, subset)
    obs = subset.loc[columns].copy()
    obs.index = pd.Index(GSE + ":" + obs.source_barcode.to_numpy(str), name="cell_uid")
    gene_index = pd.Index(genes, name="source_gene")
    var = pd.DataFrame({"gene_symbol": genes, "unambiguous_symbol": ~gene_index.duplicated(keep=False)},
                       index=gene_index)
    obj = ad.AnnData(X=matrix, obs=obs, var=var, uns={
        "expression_scale": "author processed continuous; published value in X",
        "expression_sha256": compressed_sha, "metadata_sha256": record["metadata_file_sha256"],
        "source_url": record["url"], "tissue_code": code})
    obj.var_names_make_unique(join="::duplicate-row-")
    scales = np.empty(len(obs))
    for first in range(0,len(obs),512):
        block=matrix[first:first+512].copy()
        block.data=np.expm1(block.data)
        scales[first:first+len(block.indptr)-1]=np.asarray(block.sum(axis=1)).ravel()
    check=pd.DataFrame({"sample_uid":obs.sample_uid.to_numpy(), "inverse_sum_diagnostic":scales,
                        "within_one_of_10000":np.abs(scales-10000)<=1,
                        "within_9500_10500":(scales>=9500)&(scales<=10500)})
    check=check.loc[np.diff(matrix.indptr)>0]
    scale_summary=check.groupby("sample_uid").agg(n_nonempty=("inverse_sum_diagnostic","size"),
        min_nonempty_inverse_sum=("inverse_sum_diagnostic","min"), median_nonempty_inverse_sum=("inverse_sum_diagnostic","median"),
        max_nonempty_inverse_sum=("inverse_sum_diagnostic","max"),fraction_within_one_of_10000=("within_one_of_10000","mean"),
        fraction_within_five_percent_of_10000=("within_9500_10500","mean")).reindex(sorted(obs.sample_uid.unique()))
    scale_summary["n_nonempty"]=scale_summary.n_nonempty.fillna(0).astype(int)
    scale_summary.to_csv(OUT/"coverage"/f"{code}_scale_check.csv")
    out = dest.with_suffix(".partial.h5ad")
    obj.write_h5ad(out, compression="lzf")
    reread = ad.read_h5ad(out, backed="r")
    checked=0
    try:
        if not reread.obs_names.equals(obj.obs_names) or not reread.var_names.equals(obj.var_names) or reread.X.dtype != np.float64:
            raise ValueError(f"{code}: output axes or dtype differ")
        for row in (0, len(obs)//2, len(obs)-1):
            if (matrix[row] != reread.X[row]).nnz:
                raise ValueError(f"{code}: readback scores differ")
        with gzip.open(path,"rt",newline="") as stream:
            table=csv.reader(stream)
            source_columns=next(table)[1:]
            for gene_no in range(min(8,len(genes))):
                gene,*v=next(table)
                if reread.var.gene_symbol.iloc[gene_no]!=gene:
                    raise ValueError(f"{code}: source feature mismatch")
                positive=next((j for j,value in enumerate(v) if float(value)>0), 0)
                for j in (0,len(obs)//2,len(obs)-1,positive):
                    if (source_columns[j]!=reread.obs.expression_column.iloc[j]
                            or float(v[j])!=reread.X[j,gene_no]):
                        raise ValueError(f"{code}: published value/column mismatch")
                    checked+=1
    finally:
        reread.file.close()
    out.replace(dest)
    entry = dict(code=code, tissue=ORGANS[code], input_bytes=path.stat().st_size, source_sha256=compressed_sha,
                 genes=len(genes), cells=len(obs), nnz=matrix.nnz, source_entries_checked=checked,
                 ambiguous_symbols=sorted(set(gene_index[gene_index.duplicated(keep=False)])),
                 conversion_seconds=round(time.perf_counter()-start,2),
                 h5ad_bytes=dest.stat().st_size, h5ad_sha256=sha(dest), completed_utc=now())
    (OBJECTS / f"{code}.json").write_text(json.dumps(entry, indent=2)+"\n")
    print(f"Converted {code}: {len(genes)} x {len(obs)}", flush=True)


def frozen_definitions():
    reg = pd.read_csv(FREEZE / "frozen_candidate_registry.csv")
    genes = pd.read_csv(FREEZE / "gene_sets/program_genes.csv")
    axes = []
    for code in ORGANS:
        obj = ad.read_h5ad(OBJECTS / f"{code}.h5ad", backed="r")
        axes.append(set(obj.var.loc[obj.var.unambiguous_symbol, "gene_symbol"]))
        obj.file.close()
    axis = sorted(set.intersection(*axes))
    if not axis:
        raise ValueError("No common measured gene axis")
    cov = OUT / "coverage"
    (cov / "external_rank_background.txt").write_text("\n".join(axis)+"\n")
    background = set(axis)
    model_genes = pd.read_csv(FREEZE / "gene_sets/model_program_genes.csv")
    programs, coverage = {}, []
    candidate_keys = set(zip(reg.lineage,reg.program))
    model_keys = set(zip(model_genes.lineage,model_genes.program))
    for lineage, program in sorted(candidate_keys | model_keys):
        source = genes if (lineage,program) in candidate_keys else model_genes
        g = source.loc[source.program.eq(program) & source.lineage.eq(lineage), "gene_symbol"]
        if not len(g) or g.duplicated().any():
            raise ValueError(f"Invalid frozen members: {program}")
        if (lineage,program) in model_keys:
            model_set = model_genes.loc[model_genes.program.eq(program) & model_genes.lineage.eq(lineage), "gene_symbol"]
            if set(g) != set(model_set):
                raise ValueError(f"Frozen program differs between biology and model: {program}")
        measured = [s for s in g if s in background]
        sufficient = len(measured)>=5 and len(measured)/len(g)>=.8
        programs[(lineage, program)] = measured if sufficient else []
        coverage.append(dict(lineage=lineage, program=program, primary_candidate=(lineage,program) in candidate_keys,
                             model_input=(lineage,program) in model_keys, frozen_members=len(g), measured_members=len(measured),
                             measured_fraction=len(measured)/len(g), eligible=sufficient,
                             absent_members="|".join(sorted(set(g)-set(measured)))))
    pd.DataFrame(coverage).to_csv(cov / "candidate_measurement.csv", index=False)
    return axis, programs, candidate_keys, model_keys


def gene_variants(programs):
    variants = {}
    for (lineage,p), genes in programs.items():
        if not genes:
            continue
        variants[(lineage,p)] = {p:genes}
        filtered = [g for g in genes if g not in LIVER_RBC and g not in ("HLA-DRB3","HLA-DRB4","HLA-DRB5")]
        if len(filtered)>=5:
            variants[(lineage,p)][p+"::drop_liver_Hb_DRB345"] = filtered
    return variants


def measure(axis, programs, candidates):
    values = []
    variants = gene_variants(programs)
    for code in ORGANS:
        obj = ad.read_h5ad(OBJECTS / f"{code}.h5ad")
        positions = obj.var_names.get_indexer(axis)
        for lineage in sorted({key[0] for key in variants}):
            names = {}
            for key, definitions in variants.items():
                if key[0] != lineage:
                    continue
                if key in candidates:
                    names.update(definitions)
                elif code in ("PT","MT"):
                    names[key[1]]=definitions[key[1]]
            if not names:
                continue
            rows = np.flatnonzero(obj.obs.lineage.eq(lineage))
            measured = obj.X[rows][:,positions]
            scores = rank_scores(measured, axis, names, max_rank=1500)
            data = obj.obs.iloc[rows][["source_barcode","donor_uid","sample_uid","tissue","lineage",
                                       "celltype_sub","celltype_major","celltype_global"]].reset_index(drop=True).copy()
            data["cell_uid"] = GSE + ":" + data.source_barcode
            data["n_genes_common"] = np.diff(measured.indptr)
            for key in variants:
                p=key[1]
                if key[0]==lineage and p in scores:
                    values.append(data.assign(program=p,score=scores[p],
                                              drop_liver_Hb_DRB345=scores.get(p+"::drop_liver_Hb_DRB345",np.nan)))
    if not values:
        return pd.DataFrame(columns=["cell_uid","source_barcode","donor_uid","sample_uid","tissue","lineage","celltype_sub",
                                     "celltype_major","celltype_global","program","n_genes_common","score","drop_liver_Hb_DRB345"])
    return pd.concat(values, ignore_index=True)


def model_inputs(cells, keys, support):
    rows = []
    for (lineage, program) in sorted(keys):
        values = cells.loc[cells.lineage.eq(lineage) & cells.program.eq(program)]
        for sample_uid, group in values.groupby("sample_uid"):
            if group.tissue.iloc[0] not in ("PT","LM") or len(group)<50:
                continue
            bag = group.sort_values("cell_uid")
            if len(bag)>256:
                rng = np.random.default_rng(cell_seed(f"{sample_uid}|{lineage}"))
                bag=bag.iloc[rng.choice(len(bag),256,replace=False)]
            rows.append(dict(sample_uid=sample_uid,lineage=lineage,program=program,score=bag.score.mean(),
                             n_cells=len(bag),n_population=len(group),donor_uid=group.donor_uid.iloc[0],
                             tissue=group.tissue.iloc[0],n_captured=support[sample_uid]))
    target=OUT/"model_predictions"
    target.mkdir(exist_ok=True)
    pd.DataFrame(rows,columns=["sample_uid","lineage","program","score","n_cells","n_population",
                               "donor_uid","tissue","n_captured"]).to_csv(target/"external_model_inputs.csv",index=False)


def summarize(manifest, meta, cells, programs):
    effects = OUT / "patient_effects"
    effects.mkdir(parents=True, exist_ok=True)
    totals = meta.groupby(["sample_uid","lineage"]).size()
    sizes = meta.groupby("sample_uid").size()
    coverage = []
    for sample in manifest.sample_uid:
        for lineage, p in programs:
            n = int(totals.get((sample,lineage), 0))
            entry=dict(sample_uid=sample,lineage=lineage,program=p,n_cells=n,
                       n_captured=int(sizes[sample]),captured_fraction=n/sizes[sample])
            for minimum in (20,50,100):
                entry[f"status_at_{minimum}"]=("insufficient_gene_coverage" if not programs[(lineage,p)] else
                                               "no_lineage_cells" if not n else
                                               "below_specimen_cell_threshold" if n<minimum else "eligible")
            coverage.append(entry)
    support = pd.DataFrame(coverage)
    support.to_csv(OUT/"coverage/sample_lineage_support.csv",index=False)
    sensitive = []
    policies = [("main",50),("n20",20),("n100",100),("exclude_CD3E_myeloid",50),
                ("major_global_concordance",50),("drop_liver_Hb_DRB345",50)]
    for label, threshold in policies:
        subset = cells
        if label=="exclude_CD3E_myeloid":
            subset=cells.loc[~(cells.lineage.eq("Mono/Mac/DC") & cells.celltype_sub.str.contains("CD3E", regex=False))]
        elif label=="major_global_concordance":
            subset=cells.loc[cells.celltype_global.eq(cells.lineage.map(EXPECTED_GLOBAL))]
        score = label if label=="drop_liver_Hb_DRB345" else "score"
        agg = subset.groupby(["sample_uid","lineage","program"])[score].agg(["mean","size"]).reset_index()
        table = support.merge(agg, how="left",on=["sample_uid","lineage","program"], validate="one_to_one")
        if label in ("exclude_CD3E_myeloid","major_global_concordance"):
            table["n_cells"] = table["size"].fillna(0).astype(int)
        table["score"] = table["mean"].where(table.n_cells>=threshold)
        patients, lesions, means, _ = calculate_patient_contrasts(manifest, table[SCORE_COLUMNS])
        summaries = summarize_patients(patients, seed=SEED)
        summaries.insert(0, "analysis", label)
        sensitive.append(summaries)
        if label == "main":
            patients.to_csv(effects/"patient_contrasts.csv", index=False)
            means.to_csv(effects/"patient_tissue_means.csv", index=False)
            lesions.to_csv(effects/"lesion_means.csv", index=False)
            table[SCORE_COLUMNS].to_csv(effects/"sample_program_scores.csv", index=False)
            summaries.drop(columns="analysis").to_csv(effects/"effects_summary.csv",index=False)
            leave_one_out(patients).to_csv(effects/"leave_one_patient_out.csv",index=False)
        else:
            patients.to_csv(effects/f"patients_{label}.csv",index=False)
    pd.concat(sensitive,ignore_index=True).to_csv(effects/"sensitivity_summary.csv",index=False)
    proportions = support.assign(score=support.captured_fraction)
    comp, _, _, _ = calculate_patient_contrasts(manifest, proportions[SCORE_COLUMNS])
    comp.to_csv(effects/"patient_lineage_composition.csv", index=False)
    def strata(frame, labels, filename):
        frame=frame.assign(subtype=labels)
        sub=frame.groupby(["sample_uid","lineage","program","subtype"]).score.agg(["size","mean"]).reset_index()
        sub=sub.rename(columns={"size":"n_subtype","mean":"mean_score"}).merge(
            support[["sample_uid","lineage","program","n_cells"]],on=["sample_uid","lineage","program"])
        sub.to_csv(effects/(filename+"_sample.csv"),index=False)
        decompose(manifest,sub).to_csv(effects/(filename+".csv"),index=False)
    strata(cells,cells.celltype_sub,"patient_subtype_decomposition")
    b=cells.loc[cells.lineage.eq("B/Plasma")]
    strata(b,np.where(b.celltype_sub.str.contains("plasmaB",regex=False),"Plasma","B"),
           "patient_B_plasma_decomposition")
    bins=np.digitize(cells.n_genes_common.to_numpy(),[500,1000,2000])
    strata(cells,bins.astype(str),"patient_detection_depth_decomposition")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("step", choices=("fetch","identity","convert","score"))
    parser.add_argument("--tissue-code", choices=tuple(ORGANS), nargs="+")
    parser.add_argument("--freeze-commit", help="Git revision containing the frozen M5 registry and tests")
    args=parser.parse_args()
    OUT.mkdir(parents=True,exist_ok=True)
    ad.settings.allow_write_nullable_strings = True
    records=sources()
    if args.step=="fetch":
        fetch(records)
        return
    frame=metadata(records)
    manifest=write_identity(frame,records)
    if args.step=="identity":
        print(f"Author identity: {len(frame)} cells, {len(manifest)} specimens, {manifest.donor_uid.nunique()} donors")
        return
    if not args.freeze_commit:
        parser.error("--freeze-commit is required for expression processing")
    frozen_files=["artifacts/m5/frozen_candidate_registry.csv","artifacts/m5/gene_sets/program_genes.csv",
                  "artifacts/m5/registered_tests.csv","artifacts/m5/gene_sets/model_program_genes.csv",
                  "artifacts/m5/program_scores.py","artifacts/m5/validation_rules.md"]
    for name in frozen_files:
        frozen=subprocess.check_output(["git","show",f"{args.freeze_commit}:{name}"], cwd=ROOT)
        if frozen!=(ROOT/name).read_bytes():
            raise ValueError(f"Working and frozen version disagree for {name}")
    if sha(FREEZE/"program_scores.py")!=sha(ROOT/"src/crlm/program_scores.py"):
        raise ValueError("Current rank function differs from the registered revision")
    OBJECTS.mkdir(parents=True,exist_ok=True)
    if args.step=="convert":
        for code in args.tissue_code or ORGANS:
            convert(code,frame,records[code])
        return
    start=time.perf_counter()
    axis, programs, candidate_keys, model_keys=frozen_definitions()
    cells=measure(axis,programs,candidate_keys)
    cells.to_csv(OBJECTS/"frozen_cell_scores.csv.gz",index=False,compression={"method":"gzip","mtime":0})
    model_inputs(cells,model_keys,frame.groupby("sample_uid").size())
    is_candidate=pd.MultiIndex.from_frame(cells[["lineage","program"]]).isin(list(candidate_keys))
    summarize(manifest,frame,cells.loc[is_candidate],{k:v for k,v in programs.items() if k in candidate_keys})
    tests=pd.read_csv(FREEZE/"registered_tests.csv", keep_default_na=False)
    plan=tests.loc[tests.cohort.eq(GSE)]
    summary=pd.read_csv(OUT/"patient_effects/effects_summary.csv")
    adjusted, extra, families, audit=adjust_registered(plan,summary)
    target=OUT/"registered_tests";target.mkdir(exist_ok=True)
    for name,table in (("registered_results",adjusted),("family_summary",families),("unplanned_results",extra)):
        table.to_csv(target/f"{name}.csv",index=False)
    (target/"audit.json").write_text(json.dumps(audit,indent=2)+"\n")
    profile=dict(completed_utc=now(), command=shlex.join([sys.executable,*sys.argv]), freeze_commit=args.freeze_commit,
                 seconds=time.perf_counter()-start, peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,
                 code_sha256=sha(Path(__file__)),rank_code_sha256=sha(ROOT/"src/crlm/program_scores.py"),
                 n_rank_genes=len(axis),datasets=[json.loads((OBJECTS/f"{c}.json").read_text()) for c in ORGANS])
    (OUT/"run_profile.json").write_text(json.dumps(profile,indent=2)+"\n")


if __name__=="__main__":
    main()
