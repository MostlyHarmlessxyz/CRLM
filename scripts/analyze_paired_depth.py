"""Development sensitivity of program effects to library-depth distributions and quantitative
comparison of `ucell_debiased` vs V1 `legacy_positive` depth confounding in GSE245552."""

from __future__ import annotations

import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from crlm.hgnc_alias_map import get_hgnc_mapper
from crlm.program_scores import (
    controls_from_dataframe,
    debiased_ucell_scores,
    legacy_positive_rank_scores,
)


def _canonical_indexer(obj: ad.AnnData, axis: list[str]) -> np.ndarray:
    if "canonical_symbol" in obj.var.columns:
        canon = pd.Index(obj.var["canonical_symbol"].astype(str))
    else:
        canon = pd.Index(get_hgnc_mapper().normalize_gene_axis(obj.var["gene_symbol"])["canonical_symbol"])
    cols = canon.get_indexer(axis)
    if (cols < 0).any():
        raise ValueError("Missing canonical discovery background features")
    return cols


def run() -> None:
    out = ROOT / "artifacts/m3"
    defs_dir = out / "program_definitions_draft"
    members = pd.read_csv(defs_dir / "program_genes.csv")
    members = members[(members.role == "complete_program") & members.on_discovery_axis]
    axis = (defs_dir / "ranking_background.txt").read_text(encoding="utf-8").splitlines()
    ctrl_table = pd.read_csv(defs_dir / "matched_controls.csv")
    controls_by_lin = {
        str(lin): controls_from_dataframe(sub)
        for lin, sub in ctrl_table.groupby("lineage", sort=False)
    }

    manifest = pd.read_csv(ROOT / "data/manifests/m0/sample_manifest_validated.csv")
    if manifest.lesion_uid.duplicated().any():
        raise ValueError("Aggregate technical specimens to lesions before standardizing")

    summary_rows = []
    cell_corr_rows = []

    for r in manifest[manifest.tissue.isin(["PT", "LM"])].itertuples(index=False):
        d = ad.read_h5ad(ROOT / f"data/processed/m2/{r.sample_uid}.h5ad")
        cols = _canonical_indexer(d, axis)
        for lineage, defs in members.groupby("lineage"):
            ix = np.flatnonzero(d.obs.qc_retained.to_numpy(dtype=bool) & (d.obs.lineage_m2.astype(str).to_numpy() == lineage))
            if len(ix) < 50:
                continue
            programs = defs.groupby("program").gene_symbol.apply(list).to_dict()
            mat = d.X[ix, :][:, cols]
            deb_scores = debiased_ucell_scores(
                mat, axis, programs, matched_controls=controls_by_lin[str(lineage)]
            )
            leg_scores = legacy_positive_rank_scores(mat, axis, programs)

            counts_vec = d.obs.n_counts_visible.to_numpy(dtype=np.float64)[ix]
            genes_vec = d.obs.n_genes_visible.to_numpy(dtype=np.float64)[ix]
            bucket = np.digitize(counts_vec, [1000, 3000, 10000])

            for mode_name, score_df in (("ucell_debiased", deb_scores), ("legacy_positive", leg_scores)):
                for name in programs:
                    vec = score_df[name].to_numpy(dtype=np.float64)
                    r_cnt = (
                        float(spearmanr(counts_vec, vec).statistic)
                        if (np.unique(counts_vec).size > 1 and np.unique(vec).size > 1)
                        else np.nan
                    )
                    r_gen = (
                        float(spearmanr(genes_vec, vec).statistic)
                        if (np.unique(genes_vec).size > 1 and np.unique(vec).size > 1)
                        else np.nan
                    )
                    cell_corr_rows.append(
                        {
                            "donor_uid": r.donor_uid,
                            "tissue": r.tissue,
                            "sample_uid": r.sample_uid,
                            "lineage": lineage,
                            "program": name,
                            "scoring_mode": mode_name,
                            "n_cells": len(ix),
                            "median_counts": float(np.median(counts_vec)),
                            "median_genes": float(np.median(genes_vec)),
                            "cell_spearman_counts": r_cnt,
                            "cell_spearman_genes": r_gen,
                        }
                    )
                    for group in range(4):
                        subset = bucket == group
                        summary_rows.append(
                            (
                                mode_name,
                                r.donor_uid,
                                r.tissue,
                                r.lesion_uid,
                                r.sample_uid,
                                lineage,
                                name,
                                group,
                                int(subset.sum()),
                                float(subset.mean()),
                                float(vec[subset].sum() / len(ix)),
                            )
                        )

    cells = pd.DataFrame(
        summary_rows,
        columns=[
            "scoring_mode", "donor_uid", "tissue", "lesion_uid", "sample_uid",
            "lineage", "program", "depth_bin", "n_cells", "proportion", "joint_score",
        ],
    )
    patient = (
        cells.groupby(["scoring_mode", "donor_uid", "tissue", "lineage", "program", "depth_bin"])
        .agg(n_cells=("n_cells", "sum"), p=("proportion", "mean"), joint=("joint_score", "mean"))
        .reset_index()
    )

    rows = []
    for (mode_name, donor, lineage, program), g in patient.groupby(
        ["scoring_mode", "donor_uid", "lineage", "program"]
    ):
        if set(g.tissue) != {"PT", "LM"}:
            continue
        pt = g[g.tissue == "PT"].set_index("depth_bin")
        lm = g[g.tissue == "LM"].set_index("depth_bin")
        supported = (pt.n_cells >= 10) & (lm.n_cells >= 10)
        overlap = np.minimum(pt.p, lm.p).where(supported, 0.0)
        mass = float(overlap.sum())
        pt_cond = (pt.joint / pt.p.where(pt.p > 0)).fillna(0.0)
        lm_cond = (lm.joint / lm.p.where(lm.p > 0)).fillna(0.0)
        differences = lm_cond - pt_cond
        contrast = float((differences[supported] * overlap[supported]).sum() / mass) if mass > 0 else np.nan
        rows.append(
            (
                mode_name,
                donor,
                lineage,
                program,
                float(lm.joint.sum() - pt.joint.sum()),
                contrast,
                mass,
                int(supported.sum()),
            )
        )

    result = pd.DataFrame(
        rows,
        columns=[
            "scoring_mode", "donor_uid", "lineage", "program",
            "observed_D_T", "depth_standardized_D_T", "overlap_mass", "shared_bins",
        ],
    )

    # Verify exact numerical identity against primary ucell_debiased and legacy_positive patient contrasts
    for mode_name, ref_path in (
        ("ucell_debiased", out / "patient_effects/patient_contrasts.csv"),
        ("legacy_positive", out / "patient_effects/legacy_positive/patient_contrasts.csv"),
    ):
        primary = pd.read_csv(ref_path)
        required = primary[primary.program.isin(members.program) & primary.D_T.notna()]
        sub_res = result.loc[result["scoring_mode"].eq(mode_name)]
        proof = sub_res.merge(required, on=["donor_uid", "lineage", "program"], validate="one_to_one")
        if (
            len(proof) != len(required)
            or len(proof) != len(sub_res)
            or not np.allclose(proof.observed_D_T, proof.D_T, rtol=0, atol=1e-12)
        ):
            raise ValueError(f"Paired support or observed effect differs from {mode_name} analysis")

    result.to_csv(out / "patient_depth_standardized.csv", index=False)

    # Build quantitative depth-debiasing comparison table across legacy_positive vs ucell_debiased
    cell_df = pd.DataFrame(cell_corr_rows)
    depth_piv = (
        cell_df.groupby(["donor_uid", "lineage", "program", "tissue"])[["median_counts", "median_genes"]]
        .mean()
        .unstack("tissue")
    )
    delta_counts = (depth_piv["median_counts"]["LM"] - depth_piv["median_counts"]["PT"]).rename("delta_median_counts")
    delta_genes = (depth_piv["median_genes"]["LM"] - depth_piv["median_genes"]["PT"]).rename("delta_median_genes")
    depth_deltas = pd.concat([delta_counts, delta_genes], axis=1).reset_index()

    comp_rows = []
    for (lineage, program), grp in result.groupby(["lineage", "program"]):
        entry: dict[str, object] = {"lineage": lineage, "program": program}
        for mode_name in ("legacy_positive", "ucell_debiased"):
            sub_m = grp.loc[grp["scoring_mode"].eq(mode_name)].merge(
                depth_deltas, on=["donor_uid", "lineage", "program"], how="left"
            )
            sub_c = cell_df.loc[
                cell_df["lineage"].eq(lineage)
                & cell_df["program"].eq(program)
                & cell_df["scoring_mode"].eq(mode_name)
            ]
            obs = sub_m["observed_D_T"].dropna()
            std = sub_m["depth_standardized_D_T"].dropna()
            d_sign = int(np.sign(obs.mean())) if len(obs) else 0
            obs_conc = (
                int((obs > 0).sum() if d_sign >= 0 else (obs < 0).sum()) / max(1, int((obs != 0).sum()))
                if len(obs)
                else np.nan
            )
            std_conc = (
                int((std > 0).sum() if d_sign >= 0 else (std < 0).sum()) / max(1, int((std != 0).sum()))
                if len(std)
                else np.nan
            )
            rho_dt_genes = (
                float(spearmanr(sub_m["delta_median_genes"], sub_m["observed_D_T"]).statistic)
                if sub_m["delta_median_genes"].nunique() > 1 and sub_m["observed_D_T"].nunique() > 1
                else np.nan
            )
            rho_dt_counts = (
                float(spearmanr(sub_m["delta_median_counts"], sub_m["observed_D_T"]).statistic)
                if sub_m["delta_median_counts"].nunique() > 1 and sub_m["observed_D_T"].nunique() > 1
                else np.nan
            )
            entry[f"{mode_name}_n_pairs"] = len(obs)
            entry[f"{mode_name}_mean_D_T"] = float(obs.mean()) if len(obs) else np.nan
            entry[f"{mode_name}_mean_depth_std_D_T"] = float(std.mean()) if len(std) else np.nan
            entry[f"{mode_name}_obs_concordance"] = obs_conc
            entry[f"{mode_name}_depth_std_concordance"] = std_conc
            entry[f"{mode_name}_mean_cell_rho_genes"] = float(sub_c["cell_spearman_genes"].mean())
            entry[f"{mode_name}_mean_cell_rho_counts"] = float(sub_c["cell_spearman_counts"].mean())
            entry[f"{mode_name}_paired_rho_DT_delta_genes"] = rho_dt_genes
            entry[f"{mode_name}_paired_rho_DT_delta_counts"] = rho_dt_counts
        comp_rows.append(entry)

    comp_df = pd.DataFrame(comp_rows)
    comp_df.to_csv(out / "depth_debiasing_proof_summary.csv", index=False)
    print(
        comp_df[
            [
                "lineage",
                "program",
                "legacy_positive_mean_cell_rho_genes",
                "ucell_debiased_mean_cell_rho_genes",
                "legacy_positive_paired_rho_DT_delta_genes",
                "ucell_debiased_paired_rho_DT_delta_genes",
            ]
        ].to_string(index=False)
    )


if __name__ == "__main__":
    run()
