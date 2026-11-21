#!/usr/bin/env python3
"""Small NumPy/SciPy MIL core for ALREADY FROZEN cell representations.

attention: a=softmax(H @ v), z=a @ H, logit=z @ w+b.
mean: z=mean(H), logit=z @ w+b.

This intentionally omits the learned tanh projection of standard ABMIL. It is
a transparent implementation/gradient reference, not a claim of a new model.
No gene filtering, scGPT encoding, count reconstruction, batch integration,
candidate discovery, model selection, or biological validation occurs here.
"""
from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import expit, softmax

from _common import bool_column, cli_error, read_csv, sha256_file, write_json
from validate_manifest import validate_frame


@dataclass
class Design:
    samples: list[str]
    labels: np.ndarray
    bce_weights: np.ndarray
    pair_matrix: np.ndarray
    pair_donors: list[str]
    weight_table: pd.DataFrame


def select_training_rows(frame, split=None, fold=None):
    """Validate the COMPLETE role manifest before selecting training donors."""
    samples, _, warnings = validate_frame(frame)
    dev = samples.loc[samples.role.eq("development")].copy()
    if dev.empty or (~dev.patient_mapping_verified).any():
        raise ValueError("Training requires verified development donors")
    dev_donors = set(dev.donor_uid)
    if split is not None:
        if fold is None:
            raise ValueError("--splits requires --fold")
        candidates = [item for item in split.get("folds", []) if item.get("fold") == fold]
        if len(candidates) != 1:
            raise ValueError("Requested fold must occur exactly once in the split JSON")
        item = candidates[0]
        train, test = set(item["train_donors"]), set(item["test_donors"])
        if not train or not test or train & test:
            raise ValueError("Training/test donors must be nonempty and disjoint")
        if train | test != dev_donors:
            raise ValueError("Fold donors must partition exactly the current development donors; external donors cannot enter training")
        for key, donors in [("train_sample_uids", train), ("test_sample_uids", test)]:
            expected = set(dev.loc[dev.donor_uid.isin(donors), "sample_uid"])
            if set(item[key]) != expected:
                raise ValueError(f"{key}: split and manifest disagree; all tissues/lesions of each donor must move together")
        dev = dev.loc[dev.donor_uid.isin(train)].copy()
    elif fold is not None:
        raise ValueError("--fold requires --splits")
    planned_train_donors = sorted(dev.donor_uid.unique())
    if "embedding_eligible" in dev:
        eligible = bool_column(dev.embedding_eligible, "embedding_eligible")
        if (~eligible).any():
            if "embedding_exclusion_reason" not in dev or dev.loc[~eligible, "embedding_exclusion_reason"].eq("").any():
                raise ValueError("Explicit embedding exclusions require embedding_exclusion_reason")
            warnings.append(f"Predeclared embedding exclusions: {sorted(dev.loc[~eligible, 'sample_uid'])}")
        dev = dev.loc[eligible].copy()
    # Normal tissues have no tumor-site binary training label.
    dev = dev.loc[dev.tissue.isin(["PT", "LM"])].sort_values("sample_uid").reset_index(drop=True)
    if dev.empty or set(dev.tissue) != {"PT", "LM"}:
        raise ValueError("At least one eligible PT and one eligible LM development bag are required")
    omitted = sorted(set(planned_train_donors) - set(dev.donor_uid))
    if omitted:
        warnings.append(f"Training donors without eligible PT/LM bags: {omitted}")
    return dev, planned_train_donors, warnings


def build_design(rows):
    """BCE patient -> present tumor site -> lesion -> sample equal weights.

    A paired donor's difference is the LM lesion/sample mean logit minus the
    PT lesion/sample mean logit. Cells never define the statistical weights.
    """
    rows = rows.reset_index(drop=True)
    if rows.empty or rows.sample_uid.duplicated().any():
        raise ValueError("Design requires nonempty, unique sample rows")
    if not set(rows.tissue).issubset({"PT", "LM"}):
        raise ValueError("Only PT/LM rows can receive BCE training labels")
    n_donors = rows.donor_uid.nunique()
    weights = np.zeros(len(rows), dtype=float)
    pairs, pair_donors = [], []
    for donor, donor_rows in rows.groupby("donor_uid", sort=True):
        coefficients = {}
        n_sites = donor_rows.tissue.nunique()
        for site, site_rows in donor_rows.groupby("tissue", sort=True):
            coeff = np.zeros(len(rows), dtype=float)
            n_lesions = site_rows.lesion_uid.nunique()
            for _, lesion_rows in site_rows.groupby("lesion_uid", sort=True):
                idx = lesion_rows.index.to_numpy()
                coeff[idx] = 1.0 / (n_lesions * len(idx))
            weights += coeff / (n_donors * n_sites)
            coefficients[site] = coeff
        if set(coefficients) == {"PT", "LM"}:
            pairs.append(coefficients["LM"] - coefficients["PT"])
            pair_donors.append(str(donor))
    if not np.isclose(weights.sum(), 1.0) or (weights <= 0).any():
        raise ValueError("Internal weighting error")
    table = rows[["sample_uid", "donor_uid", "tissue", "lesion_uid"]].copy()
    table["label_lm"] = rows.tissue.eq("LM").astype(int)
    table["bce_weight"] = weights
    pair_matrix = np.vstack(pairs) if pairs else np.zeros((0, len(rows)))
    return Design(rows.sample_uid.tolist(), table.label_lm.to_numpy(float), weights,
                  pair_matrix, pair_donors, table)


def load_bags(path, samples, expected_dim=None):
    """Load only the requested arrays. External arrays need not even be valid."""
    if not samples or len(samples) != len(set(samples)):
        raise ValueError("Require a nonempty list of unique sample identifiers")
    bags, dimension = [], expected_dim
    with np.load(path, allow_pickle=False) as archive:
        missing = set(samples) - set(archive.files)
        if missing:
            raise ValueError(f"Missing requested embedding bags: {sorted(missing)}; declare justified exclusions in the manifest, do not silently drop patients")
        for sample in samples:
            h = np.asarray(archive[sample])
            if h.ndim != 2 or min(h.shape) < 1 or h.dtype.kind not in "fiu":
                raise ValueError(f"{sample}: require numeric cells-by-features matrix with at least one cell/feature")
            h = np.asarray(h, dtype=np.float64)
            if not np.isfinite(h).all():
                raise ValueError(f"{sample}: non-finite embeddings")
            if dimension is None:
                dimension = h.shape[1]
            if h.shape[1] != dimension:
                raise ValueError(f"{sample}: embedding dimension mismatch")
            bags.append(h)
    return bags


def validate_provenance(info, lineage, planned_train_donors, dimension):
    required = {"lineage", "encoder_id", "representation_frozen", "preprocessing_fit_scope", "gene_universe_policy", "gene_coverage_audited"}
    if not required.issubset(info):
        raise ValueError(f"Embedding provenance missing fields: {sorted(required - set(info))}")
    if info["lineage"] != lineage or info["representation_frozen"] is not True or info["gene_coverage_audited"] is not True:
        raise ValueError("Require matching lineage and literal true for representation_frozen/gene_coverage_audited")
    if not str(info["encoder_id"]).strip() or not str(info["gene_universe_policy"]).strip():
        raise ValueError("Encoder and gene-universe policy must be recorded")
    scopes = {"external_pretrained", "train_fold_only", "cell_intrinsic", "synthetic"}
    if info["preprocessing_fit_scope"] not in scopes:
        raise ValueError(f"Unsupported/unknown preprocessing_fit_scope; use {sorted(scopes)} and audit honestly")
    if info["preprocessing_fit_scope"] == "train_fold_only":
        if set(info.get("preprocessing_training_donors", [])) != set(planned_train_donors):
            raise ValueError("Fold-fitted HVG/PCA/scaling provenance must name exactly the planned training donors")
    if "embedding_dimension" in info and info["embedding_dimension"] != dimension:
        raise ValueError("Provenance embedding_dimension does not match arrays")


def fit_scaler(bags, weights, standardize=False):
    d = bags[0].shape[1]
    if not standardize:
        return np.zeros(d), np.ones(d)
    # Equal patient/site/lesion/sample importance, not pooled-cell weighting.
    center = sum(q * h.mean(axis=0) for q, h in zip(weights, bags))
    variance = sum(q * ((h - center) ** 2).mean(axis=0) for q, h in zip(weights, bags))
    scale = np.sqrt(np.maximum(variance, 0.0))
    scale[scale < 1e-12] = 1.0
    return center, scale


def unpack(theta, d, pooling):
    expected = d + 1 if pooling == "mean" else 2 * d + 1
    if pooling not in {"mean", "attention"} or len(theta) != expected:
        raise ValueError("Invalid pooling or parameter dimension")
    return theta[:d], theta[d], None if pooling == "mean" else theta[d + 1:]


def bag_forward(h, w, bias, v=None):
    a = np.full(len(h), 1.0 / len(h)) if v is None else softmax(h @ v)
    z = a @ h
    score = float(z @ w + bias)
    # d(score)/dv = H^T [a * (H w - w^T z)].
    derivative_v = None if v is None else h.T @ (a * (h @ w - z @ w))
    return score, z, derivative_v, a


def loss_and_gradient(theta, bags, design, pooling, lambda_pair=1.0, lambda_l2=0.01):
    if not np.isfinite([lambda_pair, lambda_l2]).all() or lambda_pair < 0 or lambda_l2 < 0:
        raise ValueError("Loss coefficients must be nonnegative")
    d = bags[0].shape[1]
    w, bias, v = unpack(theta, d, pooling)
    outputs = [bag_forward(h, w, bias, v) for h in bags]
    logits = np.array([x[0] for x in outputs])
    bce = np.dot(design.bce_weights, np.logaddexp(0.0, logits) - design.labels * logits)
    d_score = design.bce_weights * (expit(logits) - design.labels)
    pair_loss = 0.0
    if len(design.pair_donors):
        delta = design.pair_matrix @ logits
        pair_loss = np.logaddexp(0.0, -delta).mean()
        d_score += lambda_pair * (design.pair_matrix.T @ (-expit(-delta) / len(delta)))
    regularized = theta.copy()
    regularized[d] = 0.0  # Never regularize the intercept.
    objective = float(bce + lambda_pair * pair_loss + 0.5 * lambda_l2 * (regularized @ regularized))
    gradient = lambda_l2 * regularized
    gradient[:d] += np.stack([x[1] for x in outputs]).T @ d_score
    gradient[d] += d_score.sum()
    if v is not None:
        gradient[d + 1:] += np.stack([x[2] for x in outputs]).T @ d_score
    return objective, gradient


def fit_model(bags, design, pooling="attention", lambda_pair=1.0, lambda_l2=0.01,
              seed=20261121, standardize=False, maxiter=500):
    if maxiter < 1:
        raise ValueError("maxiter must be positive")
    if lambda_pair > 0 and not design.pair_donors:
        raise ValueError("Positive paired loss requires at least one complete PT/LM patient pair")
    if len(bags) != len(design.samples):
        raise ValueError("Bag/design length mismatch")
    center, scale = fit_scaler(bags, design.bce_weights, standardize)
    transformed = [(h - center) / scale for h in bags]
    d = transformed[0].shape[1]
    n_parameters = d + 1 if pooling == "mean" else 2 * d + 1
    rng = np.random.default_rng(seed)
    theta0 = rng.normal(0, 0.01, n_parameters)
    theta0[d] = 0.0
    objective = lambda theta: loss_and_gradient(theta, transformed, design, pooling, lambda_pair, lambda_l2)
    fitted = minimize(objective, theta0, method="L-BFGS-B", jac=True,
                      options={"maxiter": int(maxiter), "gtol": 1e-7, "ftol": 1e-12, "maxls": 50})
    report = {"converged": bool(fitted.success), "message": str(fitted.message),
              "iterations": int(fitted.nit), "final_training_objective": float(fitted.fun),
              "gradient_infinity_norm": float(np.abs(fitted.jac).max()),
              "n_parameters": n_parameters, "n_features": d,
              "n_bags": len(bags), "n_paired_donors": len(design.pair_donors),
              "seed": int(seed), "lambda_pair": float(lambda_pair), "lambda_l2": float(lambda_l2),
              "pair_margin": 0.0, "standardize": bool(standardize), "pooling": pooling,
              "training_objective_is_not_validation": True}
    return {"theta": fitted.x, "center": center, "scale": scale}, report


def save_model(path, model, metadata):
    path = Path(path)
    if path.suffix != ".npz":
        raise ValueError("Model filename must end in .npz")
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, theta=model["theta"], center=model["center"], scale=model["scale"],
                        metadata=np.array(json.dumps(metadata, ensure_ascii=False, allow_nan=False)))


def load_model(path):
    with np.load(path, allow_pickle=False) as archive:
        model = {key: np.asarray(archive[key], dtype=float) for key in ["theta", "center", "scale"]}
        metadata = json.loads(str(archive["metadata"].item()))
    if not metadata.get("fit", {}).get("converged", False):
        raise ValueError("Saved model does not declare a converged fit")
    d = len(model["center"])
    if model["scale"].shape != (d,) or not np.isfinite(model["theta"]).all() or not np.isfinite(model["center"]).all() or not np.isfinite(model["scale"]).all() or (model["scale"] <= 0).any():
        raise ValueError("Invalid/non-finite model parameters or scaling")
    unpack(model["theta"], d, metadata["fit"]["pooling"])
    return model, metadata


def predict_bags(model, metadata, bags, samples):
    if not bags or len(bags) != len(samples) or len(samples) != len(set(samples)):
        raise ValueError("Inference requires matching nonempty bags and unique sample identifiers")
    d = len(model["center"])
    w, bias, v = unpack(model["theta"], d, metadata["fit"]["pooling"])
    records, attentions = [], {}
    for sample, h in zip(samples, bags):
        if h.shape[1] != d or not np.isfinite(h).all():
            raise ValueError(f"{sample}: invalid inference embedding")
        transformed = (h - model["center"]) / model["scale"]
        score, _, _, a = bag_forward(transformed, w, bias, v)
        records.append({"sample_uid": sample, "lineage": metadata["lineage"], "logit": score,
                        "probability_lm_label": float(expit(score)), "n_cells": len(h),
                        "was_training_sample": sample in metadata.get("training_samples", [])})
        attentions[sample] = a
    return pd.DataFrame(records), attentions


def train_cli(args):
    split = json.loads(Path(args.splits).read_text()) if args.splits else None
    rows, planned_donors, warnings = select_training_rows(read_csv(args.manifest), split, args.fold)
    design = build_design(rows)
    bags = load_bags(args.embeddings, design.samples)
    provenance = json.loads(Path(args.provenance).read_text())
    validate_provenance(provenance, args.lineage, planned_donors, bags[0].shape[1])
    model, fit = fit_model(bags, design, args.pooling, args.lambda_pair, args.lambda_l2,
                           args.seed, args.standardize, args.maxiter)
    metadata = {"format_version": 1, "lineage": args.lineage, "fit": fit,
                "training_donors": sorted(rows.donor_uid.unique()), "planned_training_donors": planned_donors,
                "training_samples": design.samples, "paired_donors": design.pair_donors,
                "manifest_sha256": sha256_file(args.manifest), "embeddings_sha256": sha256_file(args.embeddings),
                "provenance_sha256": sha256_file(args.provenance), "embedding_provenance": provenance,
                "split_sha256": sha256_file(args.splits) if args.splits else None,
                "fold": args.fold, "warnings": warnings,
                "limitations": "Already frozen representations only; no automatic encoding/QC/selection. Simple linear attention, not full tanh ABMIL. No clinical risk or causal interpretation. Confirm canonical donor identity and clinical sampling window upstream."}
    output = Path(args.output)
    write_json(output.with_suffix(".training.json"), metadata)
    if not fit["converged"]:
        raise ValueError(f"Optimizer did not converge; diagnostic saved to {output.with_suffix('.training.json')}; no model written")
    save_model(output, model, metadata)
    design.weight_table.to_csv(output.with_suffix(".weights.csv"), index=False)
    print(json.dumps({"model": str(output), "fit": fit, "warnings": warnings}, ensure_ascii=False))


def predict_cli(args):
    model, metadata = load_model(args.model)
    if args.sample:
        samples = args.sample
        if len(samples) != len(set(samples)):
            raise ValueError("Repeated --sample identifiers")
    else:
        with np.load(args.embeddings, allow_pickle=False) as archive:
            samples = sorted(archive.files)
    bags = load_bags(args.embeddings, samples, len(model["center"]))
    predictions, attentions = predict_bags(model, metadata, bags, samples)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    predictions.to_csv(args.output, index=False)
    if args.attention_output:
        destination = Path(args.attention_output)
        if destination.suffix != ".npz":
            raise ValueError("Attention filename must end in .npz")
        destination.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(destination, **attentions)
    print(f"Scored {len(predictions)} bags independently; no fitting or cross-sample normalization was performed.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    train = commands.add_parser("train", help="Fit only verified development PT/LM bags")
    for name in ["manifest", "embeddings", "provenance", "lineage", "output"]:
        train.add_argument("--" + name, required=True)
    train.add_argument("--pooling", choices=["mean", "attention"], default="attention")
    train.add_argument("--splits")
    train.add_argument("--fold", type=int)
    train.add_argument("--lambda-pair", type=float, default=1.0)
    train.add_argument("--lambda-l2", type=float, default=0.01)
    train.add_argument("--seed", type=int, default=20261121)
    train.add_argument("--maxiter", type=int, default=500)
    train.add_argument("--standardize", action="store_true", help="Fit weighted scaling using selected training bags only")
    infer = commands.add_parser("predict", help="Score independently using stored parameters/scaling")
    for name in ["model", "embeddings", "output"]:
        infer.add_argument("--" + name, required=True)
    infer.add_argument("--sample", action="append", help="Optional exact NPZ key; repeat to select multiple bags")
    infer.add_argument("--attention-output")
    args = parser.parse_args()
    cli_error(parser, lambda: train_cli(args) if args.command == "train" else predict_cli(args))


if __name__ == "__main__":
    main()
