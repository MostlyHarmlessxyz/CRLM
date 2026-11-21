"""Patient-held-out classification of collected tumor specimens, separately by lineage."""

import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import spearmanr
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src"), str(ROOT / "starter")]
from crlm.model_bags import LINEAGES, SEED, cell_seed, sample_bags
from crlm.program_scores import rank_scores
from crlm.benchmark_statistics import summarize
from linear_attention_mil import (
    build_design,
    fit_model,
    predict_bags,
    save_model,
    unpack,
    bag_forward,
)

OUT = ROOT / "artifacts/m4"
MODELS = ROOT / "models/m4"
SPLIT = ROOT / "data/manifests/m0/donor_folds_pre_qc.json"
MANIFEST = ROOT / "data/manifests/m0/sample_manifest_validated.csv"
PROGRAMS = ROOT / "artifacts/m3/program_definitions_draft/program_genes.csv"


def library_inputs(write_selection=True):
    rows = pd.read_csv(MANIFEST).query(
        "role == 'development' and tissue in ['PT','LM']"
    )
    genes = pd.read_csv(PROGRAMS)
    genes = genes.loc[genes.on_discovery_axis & genes.role.eq("complete_program")]
    groups = {lineage: {} for lineage in LINEAGES}
    samples = []
    cells = []
    for row in rows.itertuples():
        gsm = row.sample_uid
        bags = sample_bags(ROOT / f"data/processed/m2/{gsm}.h5ad")
        for lineage in LINEAGES:
            bag = bags.get(lineage)
            samples.append(
                dict(
                    sample_uid=gsm,
                    donor_uid=row.donor_uid,
                    tissue=row.tissue,
                    lineage=lineage,
                    eligible=bag is not None,
                    reason="" if bag else "fewer_than_50_retained",
                )
            )
            if bag is None:
                continue
            defs = (
                genes[genes.lineage == lineage]
                .groupby("program")
                .gene_symbol.apply(list)
                .to_dict()
            )
            scores = rank_scores(bag.counts, bag.genes, defs)
            X = bag.counts.astype(np.float32)
            X.data *= np.repeat(
                10000 / np.asarray(X.sum(axis=1)).ravel(), np.diff(X.indptr)
            )
            np.log1p(X.data, out=X.data)
            cell = bag.obs.copy()
            cell["population"] = bag.population
            cells.append(cell)
            groups[lineage][gsm] = (bag, X, scores)
    if write_selection:
        pd.DataFrame(samples).to_csv(OUT / "eligible_bags.csv", index=False)
        pd.concat(cells).to_csv(
            OUT / "selected_cells.csv.gz", compression={"method": "gzip", "mtime": 0}
        )
    return rows, groups


def pca_rep(group, train, dest, remove_liver=False):
    all_genes = next(iter(group.values()))[0].genes
    if not all(all_genes.equals(v[0].genes) for v in group.values()):
        raise ValueError("Tumor gene axes differ")
    X = sparse.vstack([group[s][1] for s in train], format="csr")
    var = (
        np.asarray(X.power(2).mean(axis=0)).ravel()
        - np.asarray(X.mean(axis=0)).ravel() ** 2
    )
    if remove_liver:
        contamination = [
            "ALB",
            "APOA1",
            "APOA2",
            "APOC3",
            "APOH",
            "TTR",
            "FGA",
            "FGB",
            "FGG",
            "HP",
            "HPR",
            "HNF4A",
            "ASGR1",
            "ASGR2",
            "CPS1",
            "FTCD",
            "CYP2E1",
            "EPCAM",
            "KRT8",
            "KRT18",
            "KRT19",
        ]
        var[all_genes.isin(contamination)] = -np.inf
    idx = np.sort(np.argsort(-var, kind="stable")[:512])
    model = PCA(n_components=32, svd_solver="randomized", random_state=SEED)
    model.fit(X[:, idx].toarray())
    np.savez_compressed(
        dest,
        genes=np.asarray(all_genes[idx], dtype=str),
        center=model.mean_,
        components=model.components_,
        training_samples=np.asarray(train, dtype=str),
    )
    return {
        s: model.transform(v[1][:, idx].toarray()).astype(float)
        for s, v in group.items()
    }


def elastic(group, rows, test, kind, dest):
    design = build_design(rows)
    sample_names = rows.sample_uid.tolist() + test
    strata = sorted(
        {str(t) for s in rows.sample_uid for t in group[s][0].obs.subtype_m2}
    )
    vectors = []
    for s in sample_names:
        bag, _, scores = group[s]
        if kind == "program_elastic":
            values = scores.mean(axis=0)
        elif kind == "composition_elastic":
            values = bag.obs.subtype_m2.value_counts(normalize=True).reindex(
                strata, fill_value=0
            )
        else:
            q = bag.obs[["n_counts_visible", "n_genes_visible"]].quantile(
                [0.25, 0.5, 0.75]
            )
            values = pd.Series(
                np.log1p(q.to_numpy()).ravel(),
                index=[f"{c}_q{v}" for v in q.index for c in q.columns],
            )
        vectors.append(np.asarray(values, float))
    x = np.stack(vectors)
    n = len(rows)
    weights = design.bce_weights * n
    scaler = StandardScaler().fit(x[:n], sample_weight=weights)
    model = LogisticRegression(
        solver="saga", l1_ratio=0.5, C=1, random_state=SEED, tol=1e-6, max_iter=10000
    )
    z = scaler.transform(x)
    model.fit(z[:n], design.labels, sample_weight=weights)
    converged = int(model.n_iter_.max()) < model.max_iter
    np.savez_compressed(
        dest,
        coef=model.coef_,
        bias=model.intercept_,
        center=scaler.mean_,
        scale=scaler.scale_,
        columns=np.asarray(values.index, dtype=str),
        training_samples=np.asarray(design.samples),
    )
    return (
        model.decision_function(z[n:]) if test else np.array([]),
        converged,
        int(model.n_iter_.max()),
        x.shape[1] + 1,
    )


def balanced_indices(bag, lineage):
    subtype = bag.obs.subtype_m2.astype(str).to_numpy()
    first, second = {
        "T/NK": (["T cell"], ["NK/cytotoxic lymphocyte"]),
        "B/Plasma": (["B cell"], ["plasma cell"]),
        "Mono/Mac/DC": (
            ["monocyte/macrophage"],
            ["conventional DC", "plasmacytoid DC", "mature DC"],
        ),
    }[lineage]
    u, v = (
        np.flatnonzero(np.isin(subtype, first)),
        np.flatnonzero(np.isin(subtype, second)),
    )
    n = min(len(u), len(v))
    if n < 3:
        return np.array([], int)
    rng = np.random.default_rng(cell_seed(str(bag.obs.index[0]) + "|balanced"))
    return np.sort(
        np.r_[rng.choice(u, n, replace=False), rng.choice(v, n, replace=False)]
    )


def scgpt_rep(group):
    out = {}
    for gsm, (bag, _, _) in group.items():
        path = ROOT / f"embeddings/m4/scgpt/{gsm}.npz"
        with np.load(path, allow_pickle=False) as f:
            subset = pd.Index(f["cell_uid"]).get_indexer(bag.obs.index)
            if min(subset) < 0 or len(set(subset)) != len(subset):
                raise ValueError(f"{gsm} scGPT cell identity mismatch")
            out[gsm] = f["embedding"][subset].astype(float)
    return out


def cell_attribution(model, bag, values, scores, pooling, gsm, seed):
    z = (values - model["center"]) / model["scale"]
    w, bias, v = unpack(model["theta"], z.shape[1], pooling)
    logit, _, _, a = bag_forward(z, w, bias, v)
    signed = a * (z @ w)
    if not np.isclose(bias + signed.sum(), logit):
        raise ValueError("Cell contributions do not reconstruct the bag logit")
    k = max(1, len(a) // 10)
    top = np.argsort(-np.abs(signed), kind="stable")[:k]
    rng = np.random.default_rng(cell_seed(f"{gsm}|{seed}|remove"))
    subtype = bag.obs.subtype_m2.astype(str).to_numpy()
    random = np.concatenate(
        [
            rng.choice(
                np.flatnonzero(subtype == t), np.sum(subtype[top] == t), replace=False
            )
            for t in np.unique(subtype[top])
        ]
    )

    def drop(indices):
        return logit - bag_forward(np.delete(z, indices, axis=0), w, bias, v)[0]

    frame = pd.DataFrame(
        {"cell_uid": bag.obs.index, "attention": a, "signed_contribution": signed}
    )
    corr = {
        f"rho_{c}": spearmanr(signed, scores[c]).statistic
        if scores[c].nunique() > 1
        else np.nan
        for c in scores
    }
    return frame, dict(
        delta_remove_top_abs=drop(top),
        delta_remove_random_subtype_matched=drop(random),
        n_removed=k,
        **corr,
    )


def run(args):
    OUT.mkdir(parents=True, exist_ok=True)
    MODELS.mkdir(parents=True, exist_ok=True)
    outdir = OUT / args.representation
    outdir.mkdir(exist_ok=True)
    folds = json.loads(SPLIT.read_text())["folds"]
    rows, groups = library_inputs()
    predictions, runs, perturb = [], [], []
    attribution = []
    for lineage, group in groups.items():
        eligible = rows[rows.sample_uid.isin(group)]
        for fold in folds:
            f = fold["fold"]
            if set(fold["train_donors"]) & set(fold["test_donors"]):
                raise ValueError("Donor leakage")
            training = eligible[
                eligible.donor_uid.isin(fold["train_donors"])
            ].reset_index(drop=True)
            testing = eligible[
                eligible.donor_uid.isin(fold["test_donors"])
            ].reset_index(drop=True)
            dest = MODELS / lineage.replace("/", "_") / f"fold{f}"
            dest.mkdir(parents=True, exist_ok=True)

            def record(method, seed, subset, logits):
                frame = subset[
                    ["sample_uid", "donor_uid", "lesion_uid", "tissue"]
                ].copy()
                frame = frame.assign(
                    method=method, lineage=lineage, fold=f, seed=seed, logit=logits
                )
                predictions.append(frame)

            if not len(testing) or training.tissue.nunique() != 2:
                runs.append(
                    dict(
                        lineage=lineage,
                        fold=f,
                        method=args.representation,
                        status="insufficient support",
                    )
                )
                continue
            train = training.sample_uid.tolist()
            if args.representation == "PCA":
                for kind in ("program_elastic", "composition_elastic", "QC_elastic"):
                    logits, ok, iterations, nparams = elastic(
                        group,
                        training,
                        testing.sample_uid.tolist(),
                        kind,
                        dest / f"{kind}.npz",
                    )
                    runs.append(
                        dict(
                            method=kind,
                            lineage=lineage,
                            fold=f,
                            seed=SEED,
                            n_parameters=nparams,
                            iterations=iterations,
                            status="converged" if ok else "iteration limit",
                        )
                    )
                    if ok:
                        record(kind, SEED, testing, logits)
                vectors = pca_rep(group, train, dest / "pca.npz")
            else:
                vectors = scgpt_rep(group)
            settings = [
                ("mean", 0.1, SEED, "mean"),
                ("attention", 0.1, SEED, "attention"),
                ("attention", 0.0, SEED, "no_pair"),
                ("attention", 0.1, SEED + 1, "attention"),
                ("attention", 0.1, SEED + 2, "attention"),
            ]
            if args.representation == "PCA":
                settings += [
                    ("mean", 0.1, SEED, "subtype_balanced"),
                    ("mean", 0.1, SEED, "permuted"),
                    ("mean", 0.1, SEED, "remove_organ_markers"),
                ]
            for pooling, pair, seed, variant in settings:
                method = f"{args.representation}_{variant}"
                tr, te = training.copy(), testing.copy()
                h = vectors
                if variant == "subtype_balanced":
                    h = {s: h[s][balanced_indices(group[s][0], lineage)] for s in h}
                    supported = [s for s, v in h.items() if len(v)]
                    tr = tr[tr.sample_uid.isin(supported)].reset_index(drop=True)
                    te = te[te.sample_uid.isin(supported)].reset_index(drop=True)
                if variant == "remove_organ_markers":
                    h = pca_rep(
                        group, train, dest / "pca_remove_organ.npz", remove_liver=True
                    )
                if variant == "permuted":
                    flip = tr.donor_uid.map(lambda d: cell_seed(f"{d}|null") % 2 == 1)
                    tr.loc[flip, "tissue"] = tr.loc[flip, "tissue"].map(
                        {"PT": "LM", "LM": "PT"}
                    )
                info = dict(
                    method=method,
                    lineage=lineage,
                    fold=f,
                    seed=seed,
                    lambda_pair=pair,
                    n_train_bags=len(tr),
                    n_test_bags=len(te),
                    n_train_donors=tr.donor_uid.nunique(),
                )
                if tr.tissue.nunique() != 2 or not len(te):
                    runs.append(dict(info, status="insufficient support"))
                    continue
                design = build_design(tr)
                if pair and not design.pair_donors:
                    runs.append(dict(info, status="no training pair"))
                    continue
                start = time.perf_counter()
                model, fit = fit_model(
                    [h[s] for s in design.samples],
                    design,
                    pooling=pooling,
                    lambda_pair=pair,
                    lambda_l2=0.01,
                    seed=seed,
                    standardize=True,
                    maxiter=500,
                )
                runs.append(
                    dict(
                        info,
                        n_parameters=fit["n_parameters"],
                        iterations=fit["iterations"],
                        seconds=time.perf_counter() - start,
                        status=fit["message"] if not fit["converged"] else "converged",
                    )
                )
                if not fit["converged"]:
                    continue
                meta = dict(
                    lineage=lineage,
                    fit=fit,
                    training_samples=design.samples,
                    training_donors=sorted(tr.donor_uid.unique()),
                    test_donors=sorted(te.donor_uid.unique()),
                    representation=args.representation,
                    variant=variant,
                )
                save_model(dest / f"{method}_{seed}.npz", model, meta)
                ids = te.sample_uid.tolist()
                pred, _ = predict_bags(model, meta, [h[s] for s in ids], ids)
                record(method, seed, te, pred.logit.to_numpy())
                if variant in ("mean", "attention") and seed == SEED:
                    for s in ids:
                        frame, stress = cell_attribution(
                            model, group[s][0], h[s], group[s][2], pooling, s, seed
                        )
                        attribution.append(
                            frame.assign(
                                method=method, fold=f, lineage=lineage, sample_uid=s
                            )
                        )
                        perturb.append(
                            dict(
                                method=method,
                                fold=f,
                                lineage=lineage,
                                sample_uid=s,
                                **stress,
                            )
                        )
            print(lineage, f, training.donor_uid.nunique(), len(testing), flush=True)
    pd.DataFrame(runs).to_csv(outdir / "benchmark_runs.csv", index=False)
    pd.concat(predictions).to_csv(outdir / "oof_predictions.csv", index=False)
    pd.DataFrame(perturb).to_csv(outdir / "attribution_stress.csv", index=False)
    pd.concat(attribution).to_csv(
        outdir / "cell_contributions.csv.gz",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )
    available = [
        p for p in (OUT / "PCA", OUT / "scGPT") if (p / "oof_predictions.csv").exists()
    ]
    oof = pd.concat([pd.read_csv(p / "oof_predictions.csv") for p in available])
    oof.to_csv(OUT / "oof_predictions.csv", index=False)
    pd.concat([pd.read_csv(p / "benchmark_runs.csv") for p in available]).to_csv(
        OUT / "benchmark_runs.csv", index=False
    )
    patient, metrics, contrasts = summarize(oof, SEED)
    patient.to_csv(OUT / "patient_site_logits.csv", index=False)
    metrics.to_csv(OUT / "patient_metrics.csv", index=False)
    compare = OUT / "paired_model_comparisons"
    compare.mkdir(exist_ok=True)
    contrasts.to_csv(compare / "auc_differences.csv", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--representation", choices=["PCA", "scGPT"], required=True)
    with threadpool_limits(limits=4):
        run(parser.parse_args())
