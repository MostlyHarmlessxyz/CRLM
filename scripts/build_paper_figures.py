"""Build final patient, classification and spatial figures from the stage tables."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
A = ROOT / "artifacts"
OUT = A / "m8"
FIG = OUT / "figures"
DATA = OUT / "source_data"
PROGRAMS = {
    "B19_S4_Mac_SPP1": "Myeloid · Mac_SPP1",
    "HALLMARK_UNFOLDED_PROTEIN_RESPONSE": "B / Plasma · UPR",
}
COLORS = ["#175b75", "#b14925", "#577644"]


def save(fig, stem):
    for ext in ("png", "pdf"):
        fig.savefig(FIG / f"{stem}.{ext}", dpi=185, bbox_inches="tight")
    plt.close(fig)


def patient_figure():
    patients = pd.concat(
        [
            pd.read_csv(A / "m3/patient_effects/patient_contrasts.csv"),
            pd.read_csv(A / "m6/patient_contrasts.csv"),
        ],
        ignore_index=True,
    )
    patients = patients[patients.program.isin(PROGRAMS)]
    patients.to_csv(DATA / "candidate_patient_effects.csv", index=False)
    cohorts = list(dict.fromkeys(patients.dataset))
    fig, axes = plt.subplots(2, 2, figsize=(10, 6.8), layout="constrained")
    for col, (key, title) in enumerate(PROGRAMS.items()):
        pp = patients[patients.program.eq(key)]
        for row, contrasts in enumerate((cohorts, ["D_T_Q", "D_N_Q", "I"])):
            ax = axes[row, col]
            ax.axhline(0, color="#aaa", lw=0.8)
            labels = []
            for x, group in enumerate(contrasts):
                if row == 0:
                    values = pp.loc[pp.dataset.eq(group), "D_T"].dropna()
                    label = group.replace("GSE", "")
                else:
                    q = pp[pp.dataset.eq("GSE164522") & pp.I.notna()]
                    values = q[group].dropna()
                    label = {
                        "D_T_Q": "LM − PT",
                        "D_N_Q": "NL − NC",
                        "I": "Interaction",
                    }[group]
                n = len(values)
                offset = np.linspace(-0.12, 0.12, n) if n > 1 else np.zeros(n)
                ax.scatter(
                    x + offset,
                    values,
                    s=26,
                    alpha=0.70,
                    color=COLORS[col],
                    edgecolors="none",
                )
                if n:
                    ax.plot(
                        [x - 0.23, x + 0.23], [values.mean()] * 2, color="black", lw=2
                    )
                labels.append(f"{label}\nN = {n}")
            ax.set_xticks(range(len(labels)), labels)
            ax.set_ylabel("Within-patient rank-score difference")
            ax.set_title(
                title if row == 0 else f"{title} · GSE164522 quartet",
                loc="left",
                fontsize=11,
                weight="bold",
            )
            if row == 0:
                ax.set_xlabel("Patient-paired LM − PT by source")
            else:
                ax.set_xlabel("Same eligible four-tissue patients")
    save(fig, "01_patient_programs")


def model_figure():
    metrics = pd.read_csv(A / "m4/patient_metrics.csv")
    methods = [
        "program_elastic",
        "composition_elastic",
        "PCA_mean",
        "PCA_attention",
        "PCA32_CTA",
        "PCA32_HA",
        "PCA32_MixMIL",
        "scGPT_mean",
        "scGPT_attention",
        "scGPT512_CTA",
        "scGPT512_HA",
    ]
    labels = [
        "Program",
        "Composition",
        "PCA mean",
        "PCA attention",
        "PCA CTA",
        "PCA HA",
        "PCA MixMIL",
        "scGPT mean",
        "scGPT attention",
        "scGPT CTA",
        "scGPT HA",
    ]
    table = metrics[metrics.method.isin(methods) & metrics.seed.eq(20261121)]
    table.to_csv(DATA / "development_classifier_metrics.csv", index=False)
    fig, ax = plt.subplots(figsize=(7.3, 5.3), layout="constrained")
    for i, lineage in enumerate(["T/NK", "Mono/Mac/DC", "B/Plasma"]):
        values = table[table.lineage.eq(lineage)].set_index("method").loc[methods]
        ax.scatter(
            values.auroc,
            np.arange(len(methods)) + 0.19 * (i - 1),
            s=29,
            color=COLORS[i],
            label=lineage,
        )
    ax.set_yticks(range(len(labels)), labels)
    ax.invert_yaxis()
    ax.axvline(0.5, lw=0.9, ls="--", color="#999")
    ax.set_xlim(0, 1.025)
    ax.set_xlabel("AUROC on held-out patients' tissue sites")
    ax.set_title(
        "GSE245552 · fixed patient splits", loc="left", weight="bold", fontsize=12
    )
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.20), ncol=3, frameon=False)
    save(fig, "02_development_classifiers")


def spatial_figure():
    all_rows = pd.read_csv(A / "m7/spatial_measurement_audit/local_distances.csv")
    all_rows.to_csv(DATA / "spatial_local_distances.csv", index=False)
    rows = all_rows[
        all_rows.tile_um.eq(200) & all_rows.spp1_min_transcripts.eq(1)
    ].sort_values(["tissue", "donor_uid"])
    fig, ax = plt.subplots(figsize=(7.2, 3.5), layout="constrained")
    for i, r in enumerate(rows.itertuples()):
        ax.plot(
            [r.null_lo, r.null_hi],
            [i, i],
            color="#999",
            lw=4,
            solid_capstyle="round",
            label="Within-tile permutation reference (95%)" if i == 0 else None,
        )
        ax.scatter(
            r.delta_distance_um,
            i,
            s=42,
            color=COLORS[0],
            zorder=3,
            label="Measured difference" if i == 0 else None,
        )
    ax.axvline(0, color="#444", lw=0.8)
    ax.set_yticks(
        range(len(rows)),
        [f"{r.donor_uid.split(':')[-1]} · {r.tissue}" for r in rows.itertuples()],
    )
    ax.set_xlabel("Nearest CD3/CD8 cell distance difference (SPP1+ minus SPP1−, µm)")
    ax.set_title(
        "GSE335552 · measured myeloid neighbourhoods",
        fontsize=11,
        loc="left",
        weight="bold",
    )
    ax.legend(
        loc="lower center", bbox_to_anchor=(0.5, -0.54), frameon=False, fontsize=8
    )
    save(fig, "03_spatial_distances")


def transfer_figure():
    metrics = pd.read_csv(A / "m6/GSE164522/model_predictions/patient_metrics.csv")
    metrics.to_csv(DATA / "external_classifier_metrics.csv", index=False)
    table = metrics[metrics.scope.eq("common_predicting_methods")]
    methods = [
        "program_elastic",
        "PCA_mean",
        "PCA_attention",
        "scGPT_mean",
        "scGPT_attention",
    ]
    labels = ["Program", "PCA mean", "PCA attention", "scGPT mean", "scGPT attention"]
    fig, axes = plt.subplots(
        1, 3, figsize=(10.4, 3.5), sharey=True, layout="constrained"
    )
    for ax, lineage, color in zip(axes, ["T/NK", "Mono/Mac/DC", "B/Plasma"], COLORS):
        rows = table[table.lineage.eq(lineage)].set_index("method").loc[methods]
        for i, (_, row) in enumerate(rows.iterrows()):
            if pd.notna(row.auroc):
                ax.plot(
                    [row.auroc_lower95, row.auroc_upper95],
                    [i, i],
                    color=color,
                    alpha=0.6,
                )
                ax.scatter(row.auroc, i, color=color, s=26)
            else:
                ax.text(0.5, i, "not estimable", ha="center", fontsize=8)
        ax.axvline(0.5, ls="--", color="#aaa", lw=0.8)
        ax.set_xlim(-0.03, 1.03)
        ax.set_title(
            f"{lineage} · {int(rows.n_donors.max())} donors", loc="left", fontsize=10
        )
        ax.set_xlabel("External AUROC")
    axes[0].set_yticks(range(len(labels)), labels)
    axes[0].invert_yaxis()
    fig.suptitle(
        "GSE164522 · shared patients among measurable classifiers",
        weight="bold",
        fontsize=12,
    )
    save(fig, "04_external_classifiers")


if __name__ == "__main__":
    FIG.mkdir(parents=True, exist_ok=True)
    DATA.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update(
        {
            "axes.spines.top": False,
            "axes.spines.right": False,
            "font.size": 9,
            "pdf.fonttype": 42,
        }
    )
    patient_figure()
    model_figure()
    spatial_figure()
    transfer_figure()
