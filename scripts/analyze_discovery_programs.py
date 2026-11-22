"""Discovery-cohort (GSE245552) debiased program scoring, full H1/H2/B-Plasma hypothesis screen,
fine-subtype symmetric decomposition, and donor-consensus NMF audit.

Upgrades over V1 (strictly restricted to GSE245552 development data):
1. Harmonizes all gene symbols via deterministic `HGNCMapper` onto the 29,569-gene
   `GSE245552` canonical nuclear axis (`245_tumor_canonical_nuclear`, free of all 353
   normal-only MT/ribosomal genes) and cross-checks `245_164_203608_shared_nuclear` (18,733).
2. Computes lineage-specific development `reference_gene_stats.csv` and locks deterministic
   `n_bins=20, n_ctrl=10` detection-matched background control sets (`matched_controls.csv`
   drawn from `245_164_203608_shared_nuclear` non-members for 100% external measurability).
3. Uses depth-debiased Mann-Whitney UCell (`debiased_ucell_scores`) as the primary metric,
   alongside `ucell_raw`, `log1p_cp10k_debiased`, and `legacy_positive` (`positive_midrank_1500`).
4. Evaluates the full `CRLM_Research_Plan.md` §5.1-5.4 & §7.4 candidate landscape (19 complete
   programs across Myeloid R/C/A/Kupffer/Bill-S4, T/NK CXCL13/Cytotoxic/Memory/Treg/B07-Th17-TWEAK/IFNG,
   and B/Plasma UPR/Secretion/APC-TLS) plus multi-seed Donor-Consensus NMF (`denovo_nmf_audit.csv`).
5. Applies a 5-gate pre-specified discovery shortlist rule (N>=8, >=75% concordance & 95% bootstrap CI
   excluding 0, dual-scale `log1p_cp10k_debiased` & QC-floor robustness, within-subtype state + pure
   focal fine-subtype concordance eliminating pure composition artifacts such as `B/Plasma::UPR`,
   and collinearity/mechanism deduplication capping primary candidates at 2-3).
"""

from __future__ import annotations

import argparse
import collections
import hashlib
import itertools
import json
import resource
import sys
import time
import warnings
from datetime import datetime, timezone
from pathlib import Path
from typing import Mapping, Sequence

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import spearmanr
from sklearn.cluster import AgglomerativeClustering
from sklearn.decomposition import NMF, non_negative_factorization

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src"), str(ROOT / "starter")]

from _common import stable_seed
from composition_decomposition import decompose_pair
from crlm.annotation_m1 import HELD_OUT_CORROBORATION_GENES, PANELS
from crlm.hgnc_alias_map import get_hgnc_mapper, load_shared_nuclear_universe
from crlm.program_scores import (
    build_detection_matched_controls,
    compute_reference_gene_stats,
    controls_to_dataframe,
    debiased_ucell_scores,
    legacy_positive_rank_scores,
    log1p_cp10k_scores,
)
from paired_effects import calculate_patient_contrasts, summarize_patients

OUT = ROOT / "artifacts/m3"
DEFS = OUT / "program_definitions_draft"
MANIFEST = ROOT / "data/manifests/m0/sample_manifest_validated.csv"
M2 = ROOT / "data/processed/m2"
M1_RAW_AXIS = ROOT / "starter/audit_data/GSE245552_feature_intersection.txt"
SEEDS = ROOT / "starter/literature_seeds/literature_program_seeds.csv"
S4_TXT = ROOT / "starter/literature_seeds/literature_program_sources/NIHMS1944744-supplement-TABLE_S4.pdf.txt"
BACK = ROOT / "starter/literature_seeds/literature_program_sources/B06_S6_background_sets.csv"
HALLMARK = DEFS / "h.all.v2025.1.Hs.symbols.gmt"

SEED = 20261121
N_BINS_CTRL = 20
N_CTRL_REPS = 10
PARENTS = ("T/NK", "Mono/Mac/DC", "B/Plasma")
RESOLVED_IMMUNE = ("T/NK", "Mono/Mac/DC", "B/Plasma", "Neutrophil", "Mast")
LABELS = RESOLVED_IMMUNE + ("Mixed", "Unknown", "Nonimmune")
SCORE_COLUMNS = ["sample_uid", "lineage", "program", "score", "n_cells"]

LIVER_RBC = {
    "ALB", "APOA1", "APOA2", "APOC3", "APOH", "TTR",
    "FGA", "FGB", "FGG", "HP", "HPR", "HBA1", "HBA2", "HBB",
}
EXPLICIT_LYMPHOID_EPITHELIAL_CONTAM = {
    "GZMB", "GNLY", "LAT", "IL2RG", "PRF1",
    "EPCAM", "KRT8", "KRT18", "KRT19", "ASGR1", "ASGR2", "HNF4A",
    "PECAM1", "VWF", "CLDN5", "FLT1",
}

FINE_POPULATIONS: dict[str, tuple[str, ...]] = {
    "T/NK": ("CD8 T cell", "CD4 T cell", "Double_neg_or_unsplit_T", "NK cell"),
    "Mono/Mac/DC": ("Macrophage", "Monocyte", "conventional DC", "plasmacytoid DC", "mature DC"),
    "B/Plasma": ("B cell", "plasma cell"),
}

# Pre-specified complete literature & hypothesis-driven programs (CRLM_Research_Plan.md §5.1-5.3)
CURATED_PROGRAMS: tuple[dict[str, object], ...] = (
    # --- H1 Myeloid axis (Mono/Mac/DC) ---
    {
        "lineage": "Mono/Mac/DC",
        "program": "Core_Lipid_SAM_R",
        "genes": [
            "SPP1", "TREM2", "APOE", "APOC1", "GPNMB", "FABP5", "LGALS3",
            "CTSB", "CTSD", "CTSL", "MMP9", "FTL", "LIPA", "PLD3", "CD9", "CD63",
        ],
        "role": "complete_program",
        "candidate_class": "candidate_program",
        "mechanism_cluster": "Myeloid_R_Lipid_SAM",
        "cleanliness_tier": 1,
        "focal_subtype": "Macrophage",
        "source_doi": "10.1158/1078-0432.CCR-22-2041;10.1038/s41586-019-1631-3;10.1158/2159-8290.CD-23-1300",
        "source_location": "CRLM_Research_Plan §5.1 H1-R (Sathe 2023 B06 / Ramachandran SAM / Matusiak 2024 B20)",
    },
    {
        "lineage": "Mono/Mac/DC",
        "program": "Core_Chemokine_IFN_C",
        "genes": [
            "CXCL9", "CXCL10", "CXCL11", "STAT1", "IRF1", "GBP1", "GBP4",
            "GBP5", "IDO1", "ISG15", "IFIT1", "IFIT2", "IFIT3", "WARS1",
        ],
        "role": "complete_program",
        "candidate_class": "candidate_program",
        "mechanism_cluster": "Myeloid_C_Chemokine_IFN",
        "cleanliness_tier": 1,
        "focal_subtype": "Macrophage",
        "source_doi": "10.1126/science.ade2292;10.1158/2159-8290.CD-23-1300;10.1136/jitc-2025-013697",
        "source_location": "CRLM_Research_Plan §5.1 H1-C (Bill 2023 B19 S5 / Matusiak 2024 B20 / Wang 2026 B36)",
    },
    {
        "lineage": "Mono/Mac/DC",
        "program": "Core_MHCII_Presentation_A",
        "genes": [
            "CD74", "CIITA", "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1",
            "HLA-DQA1", "HLA-DQB1", "HLA-DMA", "HLA-DMB", "LGMN", "IFI30",
        ],
        "role": "complete_program",
        "candidate_class": "candidate_program",
        "mechanism_cluster": "Myeloid_A_MHCII",
        "cleanliness_tier": 1,
        "focal_subtype": "Macrophage",
        "source_doi": "10.1158/1078-0432.CCR-22-2041",
        "source_location": "CRLM_Research_Plan §5.1 H1-A (Sathe 2023 B06 S5 / haplotype-invariant MHC-II core)",
    },
    {
        "lineage": "Mono/Mac/DC",
        "program": "Resident_Kupffer",
        "genes": [
            "MARCO", "VSIG4", "TIMD4", "CD5L", "FOLR2", "CETP", "SLC40A1",
            "HMOX1", "VCAM1", "LYVE1", "CD163", "C1QA", "C1QB", "C1QC",
        ],
        "role": "complete_program",
        "candidate_class": "organ_resident_control",
        "mechanism_cluster": "Myeloid_Resident_Kupffer",
        "cleanliness_tier": 1,
        "focal_subtype": "Macrophage",
        "source_doi": "10.1016/j.cell.2021.12.018;10.1038/s41467-018-06318-7",
        "source_location": "CRLM_Research_Plan §5.1 & B07 S3 Macro_c12_CD5L hepatic resident Kupffer control",
    },
    # --- H2 T/NK axis ---
    {
        "lineage": "T/NK",
        "program": "CD8_CXCL13_Dysfunction",
        "genes": [
            "CXCL13", "PDCD1", "HAVCR2", "LAG3", "TIGIT", "TOX", "ENTPD1",
            "LAYN", "TNFRSF9", "ITGAE", "RBPJ", "PHLDA1", "HOPX", "SNX9",
        ],
        "role": "complete_program",
        "candidate_class": "candidate_program",
        "mechanism_cluster": "TNK_CD8_CXCL13_Exhaustion",
        "cleanliness_tier": 1,
        "focal_subtype": "CD8 T cell",
        "source_doi": "10.1016/j.ccell.2021.12.009;10.1126/science.abe6474",
        "source_location": "CRLM_Research_Plan §5.2 H2 & B07 S3 CD8_CXCL13_Tex / Liu 2022 B02 / Zheng 2021 B08",
    },
    {
        "lineage": "T/NK",
        "program": "CD8_Cytotoxic_Effector",
        "genes": [
            "GZMA", "GZMB", "GZMH", "GZMK", "GZMM", "PRF1", "GNLY", "NKG7",
            "CST7", "CTSW", "FGFBP2", "KLRG1", "CX3CR1", "FCGR3A", "CCL4", "CCL5", "IFNG",
        ],
        "role": "complete_program",
        "candidate_class": "candidate_program",
        "mechanism_cluster": "TNK_Cytotoxic_Effector",
        "cleanliness_tier": 1,
        "focal_subtype": "CD8 T cell",
        "source_doi": "10.1016/j.ccell.2021.12.009",
        "source_location": "CRLM_Research_Plan §5.2 H2 & B07 S3 CD8_GZMK_Tem / TYROBP_NK&NKT effector axis",
    },
    {
        "lineage": "T/NK",
        "program": "T_Progenitor_Memory",
        "genes": [
            "TCF7", "IL7R", "CCR7", "SELL", "LEF1", "KLF2", "LTB",
            "GPR183", "CXCR4", "ZFP36L2", "TSC22D3", "FOXP1", "SATB1",
        ],
        "role": "complete_program",
        "candidate_class": "candidate_program",
        "mechanism_cluster": "TNK_Progenitor_Memory",
        "cleanliness_tier": 1,
        "focal_subtype": "CD4 T cell|CD8 T cell",
        "source_doi": "10.1126/science.abe6474",
        "source_location": "CRLM_Research_Plan §5.2 H2 & B07 S3 CD4_CCR7_Tn-like1 progenitor/memory axis",
    },
    {
        "lineage": "T/NK",
        "program": "CD4_Treg_Suppressive",
        "genes": [
            "FOXP3", "IL2RA", "CTLA4", "TIGIT", "IKZF2", "CCR8",
            "TNFRSF4", "TNFRSF18", "BATF", "ICOS", "LAIR2",
        ],
        "role": "complete_program",
        "candidate_class": "candidate_program",
        "mechanism_cluster": "TNK_CD4_Treg",
        "cleanliness_tier": 1,
        "focal_subtype": "CD4 T cell",
        "source_doi": "10.1016/j.ccell.2021.12.009;10.1126/science.abe6474",
        "source_location": "CRLM_Research_Plan §5.2 H2 & B07 S3 CD4_FOXP3_Treg suppressive axis",
    },
    {
        "lineage": "T/NK",
        "program": "B07_Th17_TWEAK",
        "genes": [
            "IL17A", "IL17F", "RORC", "RORA", "CCR6", "CCL20", "IL23R",
            "KLRB1", "TNFSF12", "IL22", "IL26", "CTSH", "CXCR6", "BHLHE40",
        ],
        "role": "complete_program",
        "candidate_class": "candidate_program",
        "mechanism_cluster": "TNK_CD4_Th17_TWEAK",
        "cleanliness_tier": 1,
        "focal_subtype": "CD4 T cell",
        "source_doi": "10.1016/j.cell.2024.03.029",
        "source_location": "CRLM_Research_Plan §5.3 GSE245552 original study B07 Th17-TWEAK reference",
    },
    # --- B/Plasma axis ---
    {
        "lineage": "B/Plasma",
        "program": "Plasma_Antibody_Secretion",
        "genes": [
            "MZB1", "JCHAIN", "SDC1", "TNFRSF17", "DERL3", "XBP1", "PRDM1", "IRF4",
            "SSR4", "SEC11C", "FKBP11", "SPCS2", "IGHA1", "IGHA2", "IGHG1", "IGHG2", "IGHG3", "IGHM",
        ],
        "role": "complete_program",
        "candidate_class": "candidate_program",
        "mechanism_cluster": "BPlasma_ER_Secretory_UPR",
        "cleanliness_tier": 1,
        "focal_subtype": "plasma cell",
        "source_doi": "10.1016/j.cell.2024.03.029",
        "source_location": "CRLM_Research_Plan §5.3 & B07 S3 Plasma_c1_IGHA2 constant-chain/ER secretion",
    },
    {
        "lineage": "B/Plasma",
        "program": "B_Antigen_Presentation_TLS",
        "genes": [
            "MS4A1", "CD19", "CD79A", "CD79B", "BANK1", "CD74", "HLA-DRA", "HLA-DRB1",
            "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1", "CXCR4", "CXCR5", "CCR7",
            "SELL", "LTB", "CD83", "CD40", "FCRLA",
        ],
        "role": "complete_program",
        "candidate_class": "candidate_program",
        "mechanism_cluster": "BPlasma_Bcell_APC_TLS",
        "cleanliness_tier": 1,
        "focal_subtype": "B cell",
        "source_doi": "10.1016/j.cell.2024.03.029",
        "source_location": "CRLM_Research_Plan §5.3 & B07 S3 B_c2_CXCR4 follicular/MHC-II antigen presentation",
    },
)

B19_S4_BLOCKS: tuple[tuple[str, str, str, int, str], ...] = (
    ("Mac_SPP1", "B19_S4_Mac_SPP1_clean", "Myeloid_R_Lipid_SAM", 2, "Macrophage"),
    ("Mac_APOE", "B19_S4_Mac_APOE", "Myeloid_APOE_Complement", 2, "Macrophage"),
    ("Mac_CCL18", "B19_S4_Mac_CCL18", "Myeloid_R_Lipid_SAM", 3, "Macrophage"),
    ("Mac_CXCL10", "B19_S4_Mac_CXCL10", "Myeloid_C_Chemokine_IFN", 2, "Macrophage"),
    ("Mono_CCL3L1", "B19_S4_Mono_CCL3L1", "Myeloid_Inflammatory_Mono", 2, "Monocyte"),
    ("Mac_F13A1", "B19_S4_Mac_F13A1", "Myeloid_F13A1_TissueMac", 2, "Macrophage"),
)

HALLMARK_SPECS: dict[str, tuple[str, str, int, str]] = {
    "HALLMARK_INTERFERON_GAMMA_RESPONSE": ("T/NK", "TNK_IFNG_Response", 2, "CD8 T cell|CD4 T cell"),
    "HALLMARK_UNFOLDED_PROTEIN_RESPONSE": ("B/Plasma", "BPlasma_ER_Secretory_UPR", 2, "plasma cell"),
}


def sha(path: Path) -> str:
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def _build_myeloid_contaminant_set() -> set[str]:
    """Deterministic non-myeloid cross-lineage marker filter for broad pan-cancer myeloid blocks."""
    myeloid_allowed: set[str] = set()
    non_myeloid: set[str] = set(EXPLICIT_LYMPHOID_EPITHELIAL_CONTAM) | set(LIVER_RBC)
    for _, (lin, anc, sup) in PANELS.items():
        toks = set((anc + " " + sup).split())
        if lin == "Mono/Mac/DC":
            myeloid_allowed.update(toks)
        elif lin in ("T/NK", "B/Plasma", "Mast", "Nonimmune"):
            non_myeloid.update(toks)
    for lin, gs in HELD_OUT_CORROBORATION_GENES.items():
        if lin == "Mono/Mac/DC":
            myeloid_allowed.update(gs)
        elif lin in ("T/NK", "B/Plasma", "Mast", "Nonimmune"):
            non_myeloid.update(gs)
    return non_myeloid - myeloid_allowed


def _is_myeloid_contaminant(canon_sym: str, contam_set: set[str]) -> bool:
    if canon_sym in contam_set:
        return True
    if canon_sym.startswith(("IGH", "IGK", "IGL", "TRAV", "TRBV", "TRDV", "TRGV")):
        return True
    return False


def _load_b19_s4_raw_blocks() -> dict[str, list[str]]:
    blocks: dict[str, list[str]] = collections.defaultdict(list)
    for line in S4_TXT.read_text(encoding="utf-8").replace("\f", "\n").splitlines():
        bits = line.split()
        if len(bits) == 2 and bits[0] != "cg":
            blocks[bits[0]].append(bits[1])
    return blocks


def build_program_definitions(*, write_files: bool = True) -> tuple[pd.DataFrame, list[str], list[str]]:
    """Construct canonicalized program catalog across H1, H2, B/Plasma, and background controls."""
    mapper = get_hgnc_mapper()
    axis_245 = list(load_shared_nuclear_universe("245_tumor_canonical_nuclear"))
    axis_164_set = set(load_shared_nuclear_universe("245_164_shared_nuclear"))
    axis_all = list(load_shared_nuclear_universe("245_164_203608_shared_nuclear"))
    axis_245_set = set(axis_245)
    axis_all_set = set(axis_all)

    raw_245 = M1_RAW_AXIS.read_text(encoding="utf-8").splitlines()
    raw_245_norm = mapper.normalize_gene_axis(raw_245)
    canon_to_raw245 = dict(
        zip(
            raw_245_norm.loc[raw_245_norm["is_unique_measurable"], "canonical_symbol"],
            raw_245_norm.loc[raw_245_norm["is_unique_measurable"], "source_symbol"],
        )
    )

    myeloid_contam = _build_myeloid_contaminant_set()
    rows: list[dict[str, object]] = []

    def append_program(
        lineage: str,
        program: str,
        raw_genes: Sequence[str],
        *,
        role: str,
        candidate_class: str,
        mechanism_cluster: str,
        cleanliness_tier: int,
        focal_subtype: str,
        source_doi: str,
        source_location: str,
        strip_myeloid_contam: bool = False,
    ) -> None:
        _, norm_df = mapper.normalize_gene_set(raw_genes, deduplicate=True)
        kept_df = norm_df.loc[~norm_df["deduplicated_drop"]]
        for r in kept_df.itertuples(index=False):
            if strip_myeloid_contam and _is_myeloid_contaminant(r.canonical_symbol, myeloid_contam):
                continue
            csym = r.canonical_symbol
            rows.append(
                {
                    "lineage": lineage,
                    "program": program,
                    "gene_symbol": csym,
                    "gene_symbol_original": r.source_symbol,
                    "gene_symbol_245_raw": canon_to_raw245.get(csym, ""),
                    "hgnc_id": r.hgnc_id,
                    "hgnc_match_type": r.match_type,
                    "role": role,
                    "candidate_class": candidate_class,
                    "mechanism_cluster": mechanism_cluster,
                    "cleanliness_tier": int(cleanliness_tier),
                    "focal_subtype": focal_subtype,
                    "source_doi": source_doi,
                    "source_location": source_location,
                    "on_discovery_axis": csym in axis_245_set,
                    "on_shared_nuclear_164": csym in axis_164_set,
                    "on_shared_nuclear_all": csym in axis_all_set,
                }
            )

    # 1. Curated H1 / H2 / B-Plasma literature programs
    for spec in CURATED_PROGRAMS:
        append_program(
            str(spec["lineage"]),
            str(spec["program"]),
            list(spec["genes"]),  # type: ignore[arg-type]
            role=str(spec["role"]),
            candidate_class=str(spec["candidate_class"]),
            mechanism_cluster=str(spec["mechanism_cluster"]),
            cleanliness_tier=int(spec["cleanliness_tier"]),  # type: ignore[arg-type]
            focal_subtype=str(spec["focal_subtype"]),
            source_doi=str(spec["source_doi"]),
            source_location=str(spec["source_location"]),
            strip_myeloid_contam=False,
        )

    # 2. Bill et al., Science 2023 Table S4 complete myeloid state blocks (cross-lineage cleaned)
    s4_blocks = _load_b19_s4_raw_blocks()
    for raw_block, prog_name, cluster, tier, focal in B19_S4_BLOCKS:
        append_program(
            "Mono/Mac/DC",
            prog_name,
            s4_blocks[raw_block],
            role="complete_program",
            candidate_class="candidate_program",
            mechanism_cluster=cluster,
            cleanliness_tier=tier,
            focal_subtype=focal,
            source_doi="10.1126/science.ade2292",
            source_location=f"Bill Science 2023 Table S4 {raw_block} (cross-lineage contaminant stripped)",
            strip_myeloid_contam=True,
        )
    # Also keep unstripped B19_S4_Mac_SPP1 strictly as sensitivity reference against V1_F1
    append_program(
        "Mono/Mac/DC",
        "B19_S4_Mac_SPP1",
        s4_blocks["Mac_SPP1"],
        role="sensitivity_reference",
        candidate_class="sensitivity_reference",
        mechanism_cluster="Myeloid_R_Lipid_SAM",
        cleanliness_tier=9,
        focal_subtype="Macrophage",
        source_doi="10.1126/science.ade2292",
        source_location="Bill Science 2023 Table S4 Mac_SPP1 unstripped V1 reference",
        strip_myeloid_contam=False,
    )

    # 3. MSigDB Hallmark 2025.1.Hs complete sets
    for line in HALLMARK.read_text(encoding="utf-8").splitlines():
        prog, url, *members = line.rstrip().split("\t")
        if prog in HALLMARK_SPECS:
            lin, cluster, tier, focal = HALLMARK_SPECS[prog]
            append_program(
                lin,
                prog,
                members,
                role="complete_program",
                candidate_class="candidate_program",
                mechanism_cluster=cluster,
                cleanliness_tier=tier,
                focal_subtype=focal,
                source_doi="10.1016/j.cels.2015.12.004",
                source_location=f"MSigDB 2025.1.Hs {url}",
                strip_myeloid_contam=False,
            )

    # 4. Sathe 2023 B06 S5 pathway overlaps & S6 non-cancer background sets
    src_seeds = pd.read_csv(SEEDS)
    for sid, grp in src_seeds[src_seeds["seed_set_id"].str.startswith("B06_S5_")].groupby("seed_set_id", sort=True):
        append_program(
            "Mono/Mac/DC",
            str(sid),
            grp["gene_symbol_original"].tolist(),
            role="pathway_overlap",
            candidate_class="pathway_overlap",
            mechanism_cluster="Myeloid_B06_Overlap",
            cleanliness_tier=9,
            focal_subtype="Macrophage",
            source_doi=str(grp["source_doi"].iloc[0]),
            source_location=str(grp["source_location"].iloc[0]),
        )
    back_df = pd.read_csv(BACK)
    for sid, grp in back_df.groupby("set_id", sort=True):
        append_program(
            "Mono/Mac/DC",
            str(sid),
            grp["gene_symbol_original"].tolist(),
            role="noncancer_background",
            candidate_class="noncancer_background",
            mechanism_cluster="Myeloid_Noncancer_Background",
            cleanliness_tier=9,
            focal_subtype="Macrophage",
            source_doi=str(grp["source_doi"].iloc[0]),
            source_location=str(grp["source_location"].iloc[0]),
        )

    genes = pd.DataFrame(rows)
    if genes.duplicated(["program", "gene_symbol"]).any():
        dups = genes.loc[genes.duplicated(["program", "gene_symbol"]), ["program", "gene_symbol"]].to_dict("records")
        raise ValueError(f"Duplicated canonical symbols inside program definitions: {dups[:5]}")

    if write_files:
        DEFS.mkdir(parents=True, exist_ok=True)
        genes.to_csv(DEFS / "program_genes.csv", index=False)
        (DEFS / "ranking_background.txt").write_text("\n".join(axis_245) + "\n", encoding="utf-8")
        (DEFS / "shared_nuclear_ranking_background.txt").write_text("\n".join(axis_all) + "\n", encoding="utf-8")
    return genes, axis_245, axis_all


def build_program_variants(genes: pd.DataFrame, *, write_files: bool = True) -> tuple[pd.DataFrame, dict[str, list[str]]]:
    """Build original + unique/shared/organ-cleaned variants and pairwise Jaccard overlap table."""
    info: list[dict[str, object]] = []
    sets: dict[str, list[str]] = {}
    overlaps: list[dict[str, object]] = []

    measured = genes.loc[genes["on_discovery_axis"]]
    base = {k: set(v["gene_symbol"]) for k, v in measured.groupby("program", sort=True)}
    meta_by_prog = (
        genes.groupby("program", sort=True)
        .agg(
            lineage=("lineage", "first"),
            role=("role", "first"),
            candidate_class=("candidate_class", "first"),
            mechanism_cluster=("mechanism_cluster", "first"),
            cleanliness_tier=("cleanliness_tier", "first"),
            focal_subtype=("focal_subtype", "first"),
            n_total=("gene_symbol", "size"),
        )
        .to_dict("index")
    )

    for a, b in itertools.combinations(sorted(base), 2):
        shared = base[a] & base[b]
        union = base[a] | base[b]
        min_len = min(len(base[a]), len(base[b]))
        overlaps.append(
            {
                "program_a": a,
                "program_b": b,
                "n_a": len(base[a]),
                "n_b": len(base[b]),
                "n_shared": len(shared),
                "jaccard": len(shared) / len(union) if union else 0.0,
                "overlap_coefficient": len(shared) / min_len if min_len else 0.0,
                "shared_genes": " ".join(sorted(shared)),
            }
        )

    for name in sorted(meta_by_prog):
        meta = meta_by_prog[name]
        members = base.get(name, set())
        role = str(meta["role"])
        min_req = 5 if role in ("complete_program", "sensitivity_reference") else 3
        available = len(members) >= min_req and len(members) >= 0.8 * int(meta["n_total"])

        def add_variant(variant: str, selected: set[str]) -> None:
            key = name if variant == "original" else f"{name}__{variant}"
            elig = bool(selected) and available and len(selected) >= 3
            info.append(
                {
                    "score_id": key,
                    "program": name,
                    "lineage": meta["lineage"],
                    "role": role,
                    "candidate_class": meta["candidate_class"],
                    "mechanism_cluster": meta["mechanism_cluster"],
                    "cleanliness_tier": int(meta["cleanliness_tier"]),
                    "focal_subtype": meta["focal_subtype"],
                    "variant": variant,
                    "n_requested_genes": int(meta["n_total"]),
                    "n_scored_genes": len(selected),
                    "eligible": elig,
                }
            )
            if elig:
                sets[key] = sorted(selected)

        add_variant("original", members)
        if role not in ("complete_program", "sensitivity_reference"):
            continue
        same_lin_others = set().union(
            *(
                g
                for p, g in base.items()
                if p != name
                and meta_by_prog[p]["lineage"] == meta["lineage"]
                and meta_by_prog[p]["role"] == "complete_program"
            ),
            set(),
        )
        add_variant("unique", members - same_lin_others)
        add_variant("remove_liver_rbc_HLA_DRB345", members - LIVER_RBC - {"HLA-DRB3", "HLA-DRB4", "HLA-DRB5"})

    if write_files:
        pd.DataFrame(overlaps).to_csv(DEFS / "gene_overlap.csv", index=False)
    return pd.DataFrame(info), sets


def _extract_canonical_columns(obj: ad.AnnData, axis: Sequence[str]) -> np.ndarray:
    if "canonical_symbol" in obj.var.columns:
        canon = pd.Index(obj.var["canonical_symbol"].astype(str))
    else:
        norm = get_hgnc_mapper().normalize_gene_axis(obj.var["gene_symbol"].tolist())
        canon = pd.Index(norm["canonical_symbol"].astype(str))
    cols = canon.get_indexer(axis)
    if (cols < 0).any():
        missing = [axis[i] for i in np.flatnonzero(cols < 0)[:5]]
        raise ValueError(f"AnnData missing canonical nuclear axis features: {missing}")
    return cols


def _fine_population_labels(obs: pd.DataFrame) -> np.ndarray:
    broad = obs["lineage_m2"].astype(str).to_numpy()
    sub = obs["subtype_m2"].astype(str).to_numpy()
    if "myeloid_fine" in obs.columns:
        mf = obs["myeloid_fine"].fillna("").astype(str).to_numpy()
        return np.where((broad == "Mono/Mac/DC") & (sub == "monocyte/macrophage") & (mf != ""), mf, sub)
    return sub


def build_frozen_reference_and_controls(
    objects_by_gsm: Mapping[str, ad.AnnData] | None,
    manifest: pd.DataFrame,
    axis_245: Sequence[str],
    axis_shared_all: Sequence[str],
    info: pd.DataFrame,
    sets: Mapping[str, list[str]],
    *,
    write_files: bool = True,
) -> tuple[pd.DataFrame, dict[str, dict[str, dict[str, list[str]]]]]:
    """Compute lineage-specific development reference stats on PT/LM and lock matched controls."""
    tumor_gsms = manifest.loc[manifest["tissue"].isin(["PT", "LM"]), "sample_uid"].tolist()

    def iter_lineage_matrices(lineage: str):
        for gsm in tumor_gsms:
            d = objects_by_gsm[gsm] if objects_by_gsm is not None else ad.read_h5ad(M2 / f"{gsm}.h5ad")
            cols = _extract_canonical_columns(d, axis_245)
            keep = d.obs["qc_retained"].to_numpy(dtype=bool) & d.obs["lineage_m2"].astype(str).eq(lineage).to_numpy()
            if keep.any():
                yield d.X[np.flatnonzero(keep), :][:, cols]

    stats_frames: list[pd.DataFrame] = []
    controls_by_lineage: dict[str, dict[str, dict[str, list[str]]]] = {}
    control_dfs: list[pd.DataFrame] = []
    legacy_ctrl_rows: list[tuple[str, str]] = []
    shared_set = set(axis_shared_all)

    for lineage in PARENTS:
        lin_stats = compute_reference_gene_stats(
            iter_lineage_matrices(lineage),
            axis_245,
            input_scale="counts",
            n_bins=N_BINS_CTRL,
        )
        lin_stats.insert(0, "lineage", lineage)
        stats_frames.append(lin_stats)

        elig_specs = info.loc[(info["lineage"] == lineage) & info["eligible"]]
        lin_progs = {r.score_id: sets[r.score_id] for r in elig_specs.itertuples(index=False)}
        all_prog_members = {g for gs in lin_progs.values() for g in gs}
        cand_univ = sorted(shared_set | all_prog_members)
        lin_ctrls = build_detection_matched_controls(
            lin_progs,
            lin_stats,
            candidate_universe=cand_univ,
            n_bins=N_BINS_CTRL,
            n_ctrl=N_CTRL_REPS,
            seed=SEED,
        )
        controls_by_lineage[lineage] = lin_ctrls
        cdf = controls_to_dataframe(lin_ctrls, lin_stats)
        cdf.insert(0, "lineage", lineage)
        control_dfs.append(cdf)

        # Legacy flat control sets for original complete_program entries
        for r in elig_specs.loc[
            elig_specs["role"].eq("complete_program") & elig_specs["variant"].eq("original")
        ].itertuples(index=False):
            gmap = lin_ctrls[r.score_id]
            for rep in range(N_CTRL_REPS):
                cname = f"{r.program}__matched{rep + 1:02d}"
                rep_genes = sorted({gmap[tg][rep] for tg in sets[r.score_id]})
                legacy_ctrl_rows.extend((cname, cg) for cg in rep_genes)

    ref_stats_all = pd.concat(stats_frames, ignore_index=True)
    controls_all_df = pd.concat(control_dfs, ignore_index=True)
    if write_files:
        ref_stats_all.to_csv(DEFS / "reference_gene_stats.csv", index=False)
        controls_all_df.to_csv(DEFS / "matched_controls.csv", index=False)
        pd.DataFrame(legacy_ctrl_rows, columns=["score_id", "gene_symbol"]).to_csv(
            DEFS / "matched_control_genes.csv", index=False
        )
    return ref_stats_all, controls_by_lineage


def score_libraries(
    manifest: pd.DataFrame,
    axis_245: Sequence[str],
    axis_shared_all: Sequence[str],
    info: pd.DataFrame,
    sets: Mapping[str, list[str]],
    controls_by_lineage: Mapping[str, Mapping[str, Mapping[str, list[str]]]],
    *,
    objects_by_gsm: Mapping[str, ad.AnnData] | None = None,
    verbose: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, list[tuple[str, str, str, np.ndarray]]]]:
    """Score all 39 GSE245552 libraries across primary debiased UCell and 3 comparison modes."""
    sample_rows: list[dict[str, object]] = []
    comp_rows: list[dict[str, object]] = []
    subtype_rows: list[dict[str, object]] = []
    fine_sample_rows: list[dict[str, object]] = []
    tech_rows: list[dict[str, object]] = []
    precision_rows: list[dict[str, object]] = []
    nmf_micro_pools: dict[str, list[tuple[str, str, str, np.ndarray]]] = {lin: [] for lin in PARENTS}

    shared_idx_in_245 = pd.Index(axis_245).get_indexer(axis_shared_all)
    shared_set = set(axis_shared_all)

    for row in manifest.itertuples(index=False):
        gsm = row.sample_uid
        d = objects_by_gsm[gsm] if objects_by_gsm is not None else ad.read_h5ad(M2 / f"{gsm}.h5ad")
        if set(d.obs.sample_uid) != {gsm} or not set(d.obs.lineage_m2).issubset(LABELS):
            raise ValueError(f"Unexpected M2 identity or annotation in {gsm}")

        keep = d.obs.qc_retained.to_numpy(dtype=bool)
        pass_g800 = (
            d.obs["qc_pass_g800"].to_numpy(dtype=bool)
            if "qc_pass_g800" in d.obs.columns
            else keep & (d.obs.n_genes_visible.to_numpy() >= 800) & (d.obs.n_counts_visible.to_numpy() >= 800)
        )
        labels = d.obs.lineage_m2.astype(str).to_numpy()
        fine_all = _fine_population_labels(d.obs)
        cols_245 = _extract_canonical_columns(d, axis_245)

        n_ret = int(keep.sum())
        n_res_imm = int((keep & np.isin(labels, RESOLVED_IMMUNE)).sum())
        n_core3_imm = int((keep & np.isin(labels, PARENTS)).sum())
        for lineage in LABELS:
            n_l = int((keep & (labels == lineage)).sum())
            comp_rows.append(
                {
                    "sample_uid": gsm,
                    "lineage": lineage,
                    "n_cells": n_l,
                    "n_retained": n_ret,
                    "n_total_resolved_immune": n_res_imm,
                    "n_core_three_immune": n_core3_imm,
                    "proportion": (n_l / n_ret) if n_ret > 0 else np.nan,
                    "prop_total_resolved_immune": (
                        (n_l / n_res_imm) if (lineage in RESOLVED_IMMUNE and n_res_imm > 0) else np.nan
                    ),
                    "prop_core_three_immune": (
                        (n_l / n_core3_imm) if (lineage in PARENTS and n_core3_imm > 0) else np.nan
                    ),
                }
            )

        for lineage in PARENTS:
            idx = np.flatnonzero((labels == lineage) & (d.obs.n_counts_visible.to_numpy() > 0))
            specs = info.loc[info["lineage"] == lineage]
            elig_ids = specs.loc[specs["eligible"], "score_id"].tolist()
            selected = {p: sets[p] for p in elig_ids}
            lin_ctrls = controls_by_lineage[lineage]

            if len(idx) > 0 and selected:
                mat_245 = d.X[idx, :][:, cols_245]
                deb_df, raw_u_df, _ = debiased_ucell_scores(
                    mat_245,
                    axis_245,
                    selected,
                    matched_controls=lin_ctrls,
                    return_components=True,
                )
                log_deb_df = log1p_cp10k_scores(
                    mat_245,
                    axis_245,
                    selected,
                    input_scale="counts",
                    matched_controls=lin_ctrls,
                )
                leg_df = legacy_positive_rank_scores(mat_245, axis_245, selected)

                # Shared 18,733-gene axis sensitivity on original complete programs
                orig_comp_ids = specs.loc[
                    specs["eligible"] & specs["variant"].eq("original") & specs["role"].eq("complete_program"),
                    "score_id",
                ].tolist()
                shared_progs = {
                    p: [g for g in sets[p] if g in shared_set]
                    for p in orig_comp_ids
                    if len([g for g in sets[p] if g in shared_set]) >= 5
                }
                if shared_progs:
                    mat_sh = mat_245[:, shared_idx_in_245]
                    sh_df = debiased_ucell_scores(
                        mat_sh, axis_shared_all, shared_progs, matched_controls=lin_ctrls
                    )
                else:
                    sh_df = pd.DataFrame(index=deb_df.index)
            else:
                deb_df = pd.DataFrame(columns=elig_ids, dtype=np.float64)
                raw_u_df = deb_df.copy()
                log_deb_df = deb_df.copy()
                leg_df = deb_df.copy()
                sh_df = pd.DataFrame(dtype=np.float64)

            kept = keep[idx]
            g800_m = pass_g800[idx]
            strict = kept & ~d.obs.qc_joint_low.to_numpy(dtype=bool)[idx]
            fine_kept = fine_all[idx][kept]
            cells_kept = d.obs.iloc[idx].loc[kept]
            n_kept = int(kept.sum())

            # Collect donor-balanced micro-pools (up to 4 equal pools per PT/LM lesion with >=50 cells) for consensus NMF
            if row.tissue in ("PT", "LM") and n_kept >= 50:
                kept_rows_in_idx = np.flatnonzero(kept)
                rng_pool = np.random.default_rng(stable_seed(SEED, "nmf_pool", gsm, lineage))
                perm = rng_pool.permutation(kept_rows_in_idx)
                n_pools = min(4, max(1, n_kept // 25))
                mat_kept_245 = d.X[idx[perm], :][:, cols_245]
                splits = np.array_split(np.arange(len(perm)), n_pools)
                for p_i, s_ix in enumerate(splits):
                    sub_m = mat_kept_245[s_ix].tocsr()
                    tot_c = float(sub_m.sum())
                    if tot_c > 0:
                        cp10k_vec = np.asarray(sub_m.sum(axis=0)).ravel().astype(np.float64) * (10000.0 / tot_c)
                        nmf_micro_pools[lineage].append((row.donor_uid, row.tissue, f"{gsm}_p{p_i}", np.log1p(cp10k_vec)))

            for item in specs.itertuples(index=False):
                sid = item.score_id
                has_sid = sid in deb_df.columns
                mean_deb = float(deb_df.loc[kept, sid].mean()) if (has_sid and n_kept) else np.nan
                mean_ucell_raw = float(raw_u_df.loc[kept, sid].mean()) if (has_sid and n_kept) else np.nan
                mean_log_deb = float(log_deb_df.loc[kept, sid].mean()) if (has_sid and n_kept) else np.nan
                mean_legacy = float(leg_df.loc[kept, sid].mean()) if (has_sid and n_kept) else np.nan
                mean_sh = float(sh_df.loc[kept, sid].mean()) if (sid in sh_df.columns and n_kept) else np.nan
                nonempty_mean = float(deb_df[sid].mean()) if (has_sid and len(deb_df)) else np.nan
                strict_mean = float(deb_df.loc[strict, sid].mean()) if (has_sid and strict.any()) else np.nan
                g800_mean = float(deb_df.loc[g800_m, sid].mean()) if (has_sid and g800_m.any()) else np.nan

                sample_rows.append(
                    {
                        "sample_uid": gsm,
                        "lineage": lineage,
                        "program": sid,
                        "unfiltered_score": mean_deb,
                        "ucell_debiased_score": mean_deb,
                        "ucell_raw_score": mean_ucell_raw,
                        "log1p_cp10k_debiased_score": mean_log_deb,
                        "legacy_positive_score": mean_legacy,
                        "shared_nuclear_18733_score": mean_sh,
                        "n_cells": n_kept,
                        "nonempty_score": nonempty_mean,
                        "n_nonempty": len(deb_df),
                        "strict_score": strict_mean,
                        "n_strict": int(strict.sum()),
                        "g800_score": g800_mean,
                        "n_g800": int(g800_m.sum()),
                    }
                )

                if item.variant != "original" or item.role != "complete_program" or not has_sid or n_kept == 0:
                    continue

                obs_deb = deb_df.loc[kept, sid].to_numpy(dtype=np.float64)
                obs_raw = raw_u_df.loc[kept, sid].to_numpy(dtype=np.float64)
                obs_leg = leg_df.loc[kept, sid].to_numpy(dtype=np.float64)

                if n_kept >= 2:
                    order = np.random.default_rng(stable_seed(SEED, gsm, lineage)).permutation(n_kept)
                    half = n_kept // 2
                    first, second = float(obs_deb[order[:half]].mean()), float(obs_deb[order[half:]].mean())
                    precision_rows.append(
                        {
                            "sample_uid": gsm,
                            "lineage": lineage,
                            "program": item.program,
                            "n_cells": n_kept,
                            "half_a": first,
                            "half_b": second,
                            "abs_gap": abs(first - second),
                        }
                    )

                for fine_label in FINE_POPULATIONS[lineage]:
                    f_mask = fine_kept == fine_label
                    n_f = int(f_mask.sum())
                    f_mean = float(obs_deb[f_mask].mean()) if n_f > 0 else np.nan
                    if n_f > 0:
                        subtype_rows.append(
                            {
                                "sample_uid": gsm,
                                "lineage": lineage,
                                "program": item.program,
                                "subtype": fine_label,
                                "n_subtype": n_f,
                                "n_cells": n_kept,
                                "mean_score": f_mean,
                            }
                        )
                    fine_sample_rows.append(
                        {
                            "sample_uid": gsm,
                            "lineage": fine_label,
                            "parent_lineage": lineage,
                            "program": item.program,
                            "score": f_mean if n_f >= 20 else np.nan,
                            "raw_mean_score": f_mean,
                            "n_cells": n_f,
                        }
                    )

                for var in ("n_counts_visible", "n_genes_visible"):
                    vals = cells_kept[var].to_numpy(dtype=np.float64)
                    for mode_label, vec in (
                        ("ucell_debiased", obs_deb),
                        ("ucell_raw", obs_raw),
                        ("legacy_positive", obs_leg),
                    ):
                        corr = (
                            float(spearmanr(vals, vec).statistic)
                            if (np.unique(vals).size > 1 and np.unique(vec).size > 1)
                            else np.nan
                        )
                        tech_rows.append(
                            {
                                "sample_uid": gsm,
                                "lineage": lineage,
                                "program": item.program,
                                "scoring_mode": mode_label,
                                "qc_metric": var,
                                "median_qc": float(np.median(vals)),
                                "cell_spearman": corr,
                                "n_cells": len(vals),
                            }
                        )
        if verbose:
            print(f"{gsm} {n_ret} QC-retained cells; completed multi-mode debiased program scoring", flush=True)

    return (
        pd.DataFrame(sample_rows),
        pd.DataFrame(comp_rows),
        pd.DataFrame(subtype_rows),
        pd.DataFrame(fine_sample_rows),
        pd.DataFrame(tech_rows),
        pd.DataFrame(precision_rows),
        nmf_micro_pools,
    )


def run_donor_consensus_nmf(
    nmf_micro_pools: Mapping[str, Sequence[tuple[str, str, str, np.ndarray]]],
    axis_245: Sequence[str],
    genes_df: pd.DataFrame,
    *,
    n_components: int = 6,
    n_runs: int = 10,
    n_hvg: int = 1500,
    top_k_genes: int = 50,
) -> pd.DataFrame:
    """Multi-seed donor-lesion balanced NMF consensus extraction & overlap audit against fixed programs."""
    measured_lit = genes_df.loc[genes_df["on_discovery_axis"] & genes_df["role"].eq("complete_program")]
    lit_sets_by_lin: dict[str, dict[str, set[str]]] = {
        lin: {p: set(g["gene_symbol"]) for p, g in sub.groupby("program")}
        for lin, sub in measured_lit.groupby("lineage")
    }

    axis_arr = np.asarray(axis_245, dtype=str)
    valid_hvg_mask = np.array(
        [
            (g not in LIVER_RBC)
            and not g.startswith(("IGHV", "IGKV", "IGLV", "TRAV", "TRBV", "TRDV", "TRGV"))
            for g in axis_arr
        ],
        dtype=bool,
    )
    valid_idx = np.flatnonzero(valid_hvg_mask)
    audit_rows: list[dict[str, object]] = []

    for lineage in PARENTS:
        pools = list(nmf_micro_pools.get(lineage, ()))
        if len(pools) < max(8, n_components + 2):
            continue
        donors = [p[0] for p in pools]
        tissues = [p[1] for p in pools]
        X_full = np.vstack([p[3] for p in pools])[:, valid_idx]
        syms_sub = axis_arr[valid_idx]

        var = X_full.var(axis=0)
        pos_var = np.flatnonzero(var > 1e-8)
        if len(pos_var) < n_components * 2:
            continue
        k_hvg = min(n_hvg, len(pos_var))
        top_var_idx = pos_var[np.argsort(-var[pos_var], kind="mergesort")[:k_hvg]]
        X_hvg = X_full[:, top_var_idx]
        hvg_syms = syms_sub[top_var_idx]
        scale = np.maximum(X_hvg.std(axis=0, keepdims=True), 1e-6)
        X_scaled = X_hvg / scale

        basis_list = []
        k_eff = min(n_components, X_scaled.shape[0] - 1, X_scaled.shape[1] - 1)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            for r_idx in range(n_runs):
                r_seed = stable_seed(SEED, "nmf_run", lineage, r_idx) % (2**31 - 1)
                model = NMF(
                    n_components=k_eff,
                    init="random",
                    random_state=int(r_seed),
                    max_iter=1000,
                    tol=1e-4,
                )
                model.fit(X_scaled)
                H = model.components_
                norms = np.maximum(np.linalg.norm(H, axis=1, keepdims=True), 1e-12)
                basis_list.append(H / norms)

        H_all = np.vstack(basis_list)
        sim = np.clip(H_all @ H_all.T, 0.0, 1.0)
        np.fill_diagonal(sim, 0.0)
        k_nn = min(max(2, n_runs // 2), sim.shape[0] - 1)
        nn_sim = np.sort(sim, axis=1)[:, -k_nn:].mean(axis=1)
        keep_comp = nn_sim >= 0.45
        if keep_comp.sum() < k_eff:
            keep_comp = np.ones(len(H_all), dtype=bool)
        H_clean = H_all[keep_comp]

        dist = np.clip(1.0 - (H_clean @ H_clean.T), 0.0, 2.0)
        np.fill_diagonal(dist, 0.0)
        n_clust = min(k_eff, len(H_clean))
        if n_clust == 1:
            labels = np.zeros(len(H_clean), dtype=int)
        else:
            clust = AgglomerativeClustering(n_clusters=n_clust, metric="precomputed", linkage="average")
            labels = clust.fit_predict(dist)

        # Build consensus H matrix and refit non-negative usages W + Kotliar-2019 Z-score gene spectra
        H_cons = np.zeros((n_clust, X_scaled.shape[1]), dtype=np.float64)
        stabilities = np.zeros(n_clust, dtype=np.float64)
        cluster_sizes = np.zeros(n_clust, dtype=np.int32)
        for c_id in range(n_clust):
            members_h = H_clean[labels == c_id]
            cluster_sizes[c_id] = len(members_h)
            c_vec = np.maximum(np.median(members_h, axis=0), 0.0)
            norm_c = float(np.linalg.norm(c_vec))
            if norm_c > 0:
                c_vec = c_vec / norm_c
            H_cons[c_id] = c_vec
            stabilities[c_id] = float((members_h @ c_vec).mean())

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            W_cons, _, _ = non_negative_factorization(
                X_scaled,
                H=H_cons.copy(),
                update_H=False,
                n_components=n_clust,
                max_iter=1000,
                tol=1e-4,
            )
        W_prob = W_cons / np.maximum(W_cons.sum(axis=1, keepdims=True), 1e-12)
        Z_expr = (X_hvg - X_hvg.mean(axis=0, keepdims=True)) / scale
        gene_z_spectra, _, _, _ = np.linalg.lstsq(W_prob, Z_expr, rcond=None)

        pool_meta = pd.DataFrame({"donor_uid": donors, "tissue": tissues})
        lit_map = lit_sets_by_lin.get(lineage, {})

        for c_id in range(n_clust):
            stability = float(stabilities[c_id])
            spec_vec = gene_z_spectra[c_id]
            k_top = min(top_k_genes, len(hvg_syms))
            top_idx = np.argsort(-spec_vec, kind="mergesort")[:k_top]
            top_genes = hvg_syms[top_idx].tolist()
            top_set = set(top_genes)

            # Paired donor D_T across micro-pools using consensus factor usage W_prob
            scores_pool = W_prob[:, c_id]
            df_sc = pool_meta.assign(score=scores_pool)
            dt_piv = df_sc.groupby(["donor_uid", "tissue"])["score"].mean().unstack("tissue")
            if {"PT", "LM"}.issubset(dt_piv.columns):
                dt_vals = (dt_piv["LM"] - dt_piv["PT"]).dropna()
                n_p = len(dt_vals)
                mean_dt = float(dt_vals.mean()) if n_p else np.nan
                pos_p = int((dt_vals > 0).sum())
                neg_p = int((dt_vals < 0).sum())
                conc = max(pos_p, neg_p) / max(1, pos_p + neg_p) if (pos_p + neg_p) else 0.0
            else:
                n_p, mean_dt, pos_p, neg_p, conc = 0, np.nan, 0, 0, 0.0

            best_prog = "none"
            best_jacc = 0.0
            best_ov = 0.0
            shared_preview = ""
            for p_name, p_set in sorted(lit_map.items()):
                inter = top_set & p_set
                if not inter:
                    continue
                jacc = len(inter) / len(top_set | p_set)
                ov = len(inter) / min(len(top_set), len(p_set))
                if (ov, jacc) > (best_ov, best_jacc):
                    best_prog = p_name
                    best_jacc = jacc
                    best_ov = ov
                    shared_preview = " ".join(sorted(inter)[:12])

            cycling_overlap = len(top_set & {"MKI67", "TOP2A", "STMN1", "PCNA", "TYMS", "PCLAF", "NUSAP1"})
            stress_overlap = len(top_set & {"HSPA1A", "HSPA1B", "DNAJB1", "FOS", "JUN", "HSPH1", "IER2"})
            if best_ov >= 0.25 or best_jacc >= 0.10:
                disposition = "converges_to_known_literature_program"
            elif cycling_overlap >= 3:
                disposition = "cell_cycle_proliferation_module"
            elif stress_overlap >= 3:
                disposition = "dissociation_or_heat_shock_stress_module"
            elif stability >= 0.75 and n_p >= 8 and conc >= 0.75:
                disposition = "unmatched_consistent_module_review"
            else:
                disposition = "low_concordance_or_donor_specific_not_promoted"

            audit_rows.append(
                {
                    "lineage": lineage,
                    "factor_id": f"{lineage.replace('/', '_')}__cNMF_{c_id + 1:02d}",
                    "n_runs_clustered": int(cluster_sizes[c_id]),
                    "consensus_cosine_stability": round(stability, 4),
                    "n_paired_donors": n_p,
                    "mean_pool_D_T": mean_dt,
                    "n_positive": pos_p,
                    "n_negative": neg_p,
                    "direction_concordance": round(conc, 4),
                    "best_matched_literature_program": best_prog,
                    "best_jaccard": round(best_jacc, 4),
                    "best_overlap_coefficient": round(best_ov, 4),
                    "shared_literature_genes": shared_preview,
                    "top_15_genes": " ".join(top_genes[:15]),
                    "top_50_genes": ";".join(top_genes),
                    "disposition": disposition,
                }
            )
    return pd.DataFrame(audit_rows)


def patient_outputs(
    manifest: pd.DataFrame,
    scores: pd.DataFrame,
    where: Path | None,
    *,
    summary: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame | None]:
    patients, lesions, tissue, _ = calculate_patient_contrasts(manifest, scores[SCORE_COLUMNS])
    result = summarize_patients(patients, seed=SEED) if summary else None
    if where is not None:
        where.mkdir(parents=True, exist_ok=True)
        patients.to_csv(where / "patient_contrasts.csv", index=False)
        tissue.to_csv(where / "patient_tissue_means.csv", index=False)
        if summary and result is not None:
            lesions.to_csv(where / "lesion_means.csv", index=False)
            result.to_csv(where / "effects_summary.csv", index=False)
    return patients, result


def leave_one_out(patients: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (lineage, program), group in patients.groupby(["lineage", "program"]):
        available = group.dropna(subset=["D_T"])
        for patient in available.donor_uid:
            others = available[available.donor_uid != patient].D_T
            rows.append((lineage, program, patient, len(others), float(others.mean())))
    return pd.DataFrame(rows, columns=["lineage", "program", "omitted_patient", "n_remaining", "mean_D_T"])


def decompose(manifest: pd.DataFrame, subtypes: pd.DataFrame) -> pd.DataFrame:
    """Exact lesion-equal-weighted symmetric subtype composition/state decomposition (0% unresolved)."""
    usable = subtypes[subtypes.n_cells >= 50].merge(
        manifest[["sample_uid", "donor_uid", "lesion_uid", "tissue"]],
        on="sample_uid",
        validate="many_to_one",
    )
    usable = usable[usable.tissue.isin(["PT", "LM"])].copy()
    if manifest.duplicated("lesion_uid").any():
        raise ValueError("Separate technical captures must first be averaged within each lesion")
    usable["proportion"] = usable.n_subtype / usable.n_cells
    usable["score_mass"] = usable.proportion * usable.mean_score
    rows = []
    keys = ["donor_uid", "lineage", "program"]
    for grouping, group in usable.groupby(keys):
        frame = []
        observed = {}
        names = sorted(group.subtype.unique())
        for tissue, tissue_rows in group.groupby("tissue"):
            n_lesions = tissue_rows.lesion_uid.nunique()
            aggregated = (
                tissue_rows.groupby("subtype")[["n_subtype", "proportion", "score_mass"]]
                .sum()
                .reindex(names, fill_value=0)
            )
            p = aggregated.proportion / n_lesions
            mu = (aggregated.score_mass / aggregated.proportion).where(p > 0)
            observed[tissue] = float(aggregated.score_mass.sum() / n_lesions)
            frame.extend(
                dict(
                    tissue=tissue,
                    subtype=s,
                    proportion=float(p.at[s]),
                    mean_score=float(mu.at[s]) if pd.notna(mu.at[s]) else np.nan,
                    n_subtype=int(aggregated.n_subtype.at[s]),
                )
                for s in names
            )
        values = pd.DataFrame(frame)
        for support in (1, 20):
            checked = values.copy()
            checked.loc[checked.n_subtype < support, "mean_score"] = np.nan
            pair = decompose_pair(checked)
            common = pair.get("common_support") or {}
            full = pair.get("full_decomposition") or {}
            cond = common.get("renormalized_components", {})
            actual = observed.get("LM", np.nan) - observed.get("PT", np.nan)
            if full and (
                abs(full["identity_residual"]) > 1e-10
                or not np.isclose(actual, full["total_difference"], rtol=0, atol=1e-10)
            ):
                raise ValueError(f"Marginal mean decomposition failed for {grouping}")
            rows.append(
                dict(zip(keys, grouping))
                | {
                    "min_subtype_cells_per_tissue": support,
                    "status": pair["status"],
                    "observed_total_D_T": actual,
                    "full_composition": full.get("composition"),
                    "full_within_subtype_state": full.get("within_type_state"),
                    "common_mass_PT": common.get("retained_mass_PT"),
                    "common_mass_LM": common.get("retained_mass_LM"),
                    "common_total_D_T": cond.get("total_difference"),
                    "common_composition": cond.get("composition"),
                    "common_within_subtype_state": cond.get("within_type_state"),
                    "common_subtypes": " | ".join(common.get("subtypes", [])),
                }
            )
    return pd.DataFrame(rows)


def state_summaries(
    manifest: pd.DataFrame,
    raw: pd.DataFrame,
    fine_samples: pd.DataFrame,
    info: pd.DataFrame,
    where: Path | None,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Compute patient contrasts across primary debiased UCell, alternative scales, QC tiers, and fine subtypes."""
    rows = []
    base_ids = set(info.loc[info.variant == "original", "score_id"])
    branches = (
        ("main", 50, "ucell_debiased_score", "n_cells", raw.program.isin(base_ids)),
        ("ucell_raw", 50, "ucell_raw_score", "n_cells", raw.program.isin(base_ids)),
        ("log1p_cp10k_debiased", 50, "log1p_cp10k_debiased_score", "n_cells", raw.program.isin(base_ids)),
        ("legacy_positive", 50, "legacy_positive_score", "n_cells", raw.program.isin(base_ids)),
        ("shared_nuclear_18733", 50, "shared_nuclear_18733_score", "n_cells", raw.program.isin(base_ids)),
        ("min_genes_800_align164", 50, "g800_score", "n_g800", raw.program.isin(base_ids)),
        ("n20", 20, "ucell_debiased_score", "n_cells", raw.program.isin(base_ids)),
        ("n100", 100, "ucell_debiased_score", "n_cells", raw.program.isin(base_ids)),
        ("before_lower_tail_qc", 50, "nonempty_score", "n_nonempty", raw.program.isin(base_ids)),
        ("strict_library_qc", 50, "strict_score", "n_strict", raw.program.isin(base_ids)),
        ("gene_sensitivity", 50, "ucell_debiased_score", "n_cells", ~raw.program.isin(base_ids)),
    )

    main_patients = pd.DataFrame()
    main_effects = pd.DataFrame()
    loo = pd.DataFrame()

    for name, threshold, value, ncol, mask in branches:
        current = raw.loc[mask].copy()
        if current.empty:
            continue
        current["score"] = current[value].where(current[ncol] >= threshold)
        current["n_cells"] = current[ncol]
        if name == "main" and where is not None:
            current[SCORE_COLUMNS].to_csv(OUT / "sample_program_scores.csv", index=False)
        dest = None if where is None else (where if name == "main" else where / name)
        need_summary = name in ("main", "log1p_cp10k_debiased", "legacy_positive", "ucell_raw")
        patients, summary_df = patient_outputs(manifest, current, dest, summary=need_summary)
        for (lineage, program), group in patients.groupby(["lineage", "program"]):
            val = group.D_T.dropna()
            rows.append(
                (
                    name,
                    lineage,
                    program,
                    len(val),
                    float(val.mean()) if len(val) else np.nan,
                    float(val.median()) if len(val) else np.nan,
                    int((val > 0).sum()),
                    int((val < 0).sum()),
                )
            )
        if name == "main":
            main_patients = patients
            main_effects = summary_df if summary_df is not None else pd.DataFrame()
            loo = leave_one_out(patients)
            if where is not None:
                loo.to_csv(where / "leave_one_patient_out.csv", index=False)

    fine_dest = None if where is None else where / "fine_subtypes"
    fine_patients, fine_effects = patient_outputs(
        manifest, fine_samples[SCORE_COLUMNS], fine_dest, summary=True
    )
    assert fine_effects is not None

    sens_df = pd.DataFrame(
        rows,
        columns=["analysis", "lineage", "program", "n_pairs", "mean_D_T", "median_D_T", "positive", "negative"],
    )
    return sens_df, loo, main_patients, main_effects, fine_effects


def shortlist_candidates(
    info: pd.DataFrame,
    genes_df: pd.DataFrame,
    sensitivity: pd.DataFrame,
    loo: pd.DataFrame,
    main_patients: pd.DataFrame,
    main_effects: pd.DataFrame,
    decomp_df: pd.DataFrame,
    fine_effects: pd.DataFrame,
    raw_samples: pd.DataFrame,
    manifest: pd.DataFrame,
) -> pd.DataFrame:
    """Apply the 5 pre-specified discovery shortlist gates and cap primary shortlist at 2-3."""
    dt_eff = main_effects.loc[main_effects["contrast"].eq("D_T")].set_index(["lineage", "program"])
    fine_dt = fine_effects.loc[fine_effects["contrast"].eq("D_T")].set_index(["lineage", "program"])
    decomp_s1 = decomp_df.loc[
        decomp_df["min_subtype_cells_per_tissue"].eq(1) & decomp_df["observed_total_D_T"].notna()
    ]

    # Unpaired normal NL - NC organ benchmark in GSE245552 (6 normal libraries across donors)
    norm_samples = raw_samples.merge(
        manifest[["sample_uid", "donor_uid", "tissue"]], on="sample_uid", validate="many_to_one"
    )
    norm_samples = norm_samples.loc[
        norm_samples["tissue"].isin(["NC", "NL"]) & (norm_samples["n_cells"] >= 20)
    ]

    measured_sets = {
        p: set(g["gene_symbol"])
        for p, g in genes_df.loc[genes_df["on_discovery_axis"]].groupby("program")
    }

    outcomes: list[dict[str, object]] = []
    complete_specs = info.loc[(info["role"] == "complete_program") & (info["variant"] == "original")]

    for r in complete_specs.itertuples(index=False):
        metrics = sensitivity.loc[
            (sensitivity.program == r.program) & (sensitivity.lineage == r.lineage)
        ].set_index("analysis")
        main = metrics.loc["main"]
        d = int(np.sign(main.mean_D_T)) if pd.notna(main.mean_D_T) else 0
        concordant = int(main.positive if d == 1 else main.negative if d == -1 else 0)
        n_nonzero = int(main.positive + main.negative)
        share = concordant / n_nonzero if n_nonzero else 0.0

        eff_row = dt_eff.loc[(r.lineage, r.program)] if (r.lineage, r.program) in dt_eff.index else None
        ci_lo = (
            float(eff_row["bootstrap_mean_ci_low"])
            if (eff_row is not None and pd.notna(eff_row["bootstrap_mean_ci_low"]))
            else np.nan
        )
        ci_hi = (
            float(eff_row["bootstrap_mean_ci_high"])
            if (eff_row is not None and pd.notna(eff_row["bootstrap_mean_ci_high"]))
            else np.nan
        )
        raw_sign_p = (
            float(eff_row["sign_p_two_sided"])
            if (eff_row is not None and pd.notna(eff_row["sign_p_two_sided"]))
            else np.nan
        )
        ci_excludes_zero = bool((d == 1 and ci_lo > 0.0) or (d == -1 and ci_hi < 0.0))

        influence = loo.loc[(loo.program == r.program) & (loo.lineage == r.lineage), "mean_D_T"]

        # Log1p CP10k debiased dual-scale metrics
        log_m = metrics.loc["log1p_cp10k_debiased"]
        log_conc = int(log_m.positive if d == 1 else log_m.negative if d == -1 else 0)
        log_nz = int(log_m.positive + log_m.negative)
        log_share = log_conc / log_nz if log_nz else 0.0

        # Within-subtype symmetric decomposition
        dec_sub = decomp_s1.loc[decomp_s1["lineage"].eq(r.lineage) & decomp_s1["program"].eq(r.program)]
        mean_obs_dt = float(dec_sub["observed_total_D_T"].mean()) if len(dec_sub) else np.nan
        mean_full_comp = float(dec_sub["full_composition"].mean()) if len(dec_sub) else np.nan
        mean_full_state = float(dec_sub["full_within_subtype_state"].mean()) if len(dec_sub) else np.nan
        mean_common_state = float(dec_sub["common_within_subtype_state"].mean()) if len(dec_sub) else np.nan
        within_state_share = (
            abs(mean_full_state) / max(1e-12, abs(mean_obs_dt))
            if (pd.notna(mean_full_state) and pd.notna(mean_obs_dt) and abs(mean_obs_dt) > 1e-12)
            else 0.0
        )

        # Focal pure fine subtype metrics
        focal_candidates = [s.strip() for s in str(r.focal_subtype).split("|") if s.strip()]
        best_focal_label = focal_candidates[0]
        focal_n = 0
        focal_mean = np.nan
        focal_share = 0.0
        for f_lab in focal_candidates:
            if (f_lab, r.program) not in fine_dt.index:
                continue
            fr = fine_dt.loc[(f_lab, r.program)]
            fnz = int(fr["n_positive"] + fr["n_negative"])
            fc = int(fr["n_positive"] if d == 1 else fr["n_negative"] if d == -1 else 0)
            fsh = fc / fnz if fnz else 0.0
            fmean = float(fr["mean_difference"]) if pd.notna(fr["mean_difference"]) else np.nan
            if (int(fr["n_patients"]) >= 5 and pd.notna(fmean) and np.sign(fmean) == d and fsh > focal_share) or focal_n == 0:
                best_focal_label = f_lab
                focal_n = int(fr["n_patients"])
                focal_mean = fmean
                focal_share = fsh

        # GSE245552 normal NL vs NC donor mean shift
        ns = norm_samples.loc[
            norm_samples["lineage"].eq(r.lineage) & norm_samples["program"].eq(r.program)
        ]
        if not ns.empty:
            by_dt = ns.groupby(["tissue", "donor_uid"])["ucell_debiased_score"].mean()
            nl_mean = float(by_dt.loc["NL"].mean()) if "NL" in by_dt.index.get_level_values(0) else np.nan
            nc_mean = float(by_dt.loc["NC"].mean()) if "NC" in by_dt.index.get_level_values(0) else np.nan
            normal_nl_minus_nc = nl_mean - nc_mean
        else:
            normal_nl_minus_nc = np.nan

        reasons: list[str] = []
        if main.n_pairs < 8:
            reasons.append("fewer_than_8_pairs")
        if not n_nonzero or share < 0.75:
            reasons.append("direction_support_below_75pct")
        if not ci_excludes_zero:
            reasons.append("bootstrap_ci_crosses_zero")
        if d == 0 or len(influence) < 2 or not (np.sign(influence) == d).all():
            reasons.append("leave_one_out_direction")
        if pd.isna(log_m.mean_D_T) or np.sign(log_m.mean_D_T) != d or log_share < 0.70:
            reasons.append("log1p_scale_discordant")
        for subset in ("n20", "n100", "min_genes_800_align164"):
            if subset in metrics.index and metrics.at[subset, "n_pairs"] >= 4:
                if np.sign(metrics.at[subset, "mean_D_T"]) != d:
                    reasons.append(f"{subset}_direction")
        if pd.isna(mean_common_state) or np.sign(mean_common_state) != d or np.sign(mean_full_state) != d:
            reasons.append("within_subtype_state_discordant")
        elif within_state_share < 0.35:
            reasons.append("composition_dominated_below_35pct_state")
        if focal_n < 5 or pd.isna(focal_mean) or np.sign(focal_mean) != d or focal_share < 0.70:
            reasons.append("focal_subtype_not_concordant")
        if str(r.candidate_class) == "organ_resident_control":
            reasons.append("organ_resident_control_benchmark")

        outcomes.append(
            {
                "lineage": r.lineage,
                "program": r.program,
                "candidate_class": r.candidate_class,
                "mechanism_cluster": r.mechanism_cluster,
                "cleanliness_tier": int(r.cleanliness_tier),
                "n_pairs": int(main.n_pairs),
                "mean_D_T": float(main.mean_D_T),
                "median_D_T": float(main.median_D_T),
                "bootstrap_ci_low": ci_lo,
                "bootstrap_ci_high": ci_hi,
                "discovery_sign_p": raw_sign_p,
                "direction": d,
                "n_nonzero": n_nonzero,
                "direction_fraction": share,
                "log1p_mean_D_T": float(log_m.mean_D_T) if pd.notna(log_m.mean_D_T) else np.nan,
                "log1p_direction_fraction": log_share,
                "full_composition_D_T": mean_full_comp,
                "full_within_subtype_state_D_T": mean_full_state,
                "common_within_subtype_state_D_T": mean_common_state,
                "within_state_share": within_state_share,
                "focal_subtype": best_focal_label,
                "focal_n_pairs": focal_n,
                "focal_mean_D_T": focal_mean,
                "focal_direction_fraction": focal_share,
                "normal_unpaired_NL_minus_NC": normal_nl_minus_nc,
                "loo_min": float(influence.min()) if len(influence) else np.nan,
                "loo_max": float(influence.max()) if len(influence) else np.nan,
                "gates_1_to_4_pass": len(reasons) == 0,
                "screen_pass": len(reasons) == 0,
                "dedup_status": "not_passed",
                "shortlisted": False,
                "screen_reason": " | ".join(reasons),
            }
        )

    res_df = pd.DataFrame(outcomes)

    # Gate 5: Within-lineage mechanism & collinearity deduplication, then cap total primary candidates at <= 3
    dt_piv = main_patients.pivot(index="donor_uid", columns=["lineage", "program"], values="D_T")
    passing_idx = res_df.index[res_df["gates_1_to_4_pass"]].tolist()
    # Order passing candidates within each mechanism cluster by cleanliness_tier ASC, direction_fraction DESC
    passing_sorted = sorted(
        passing_idx,
        key=lambda ix: (
            int(res_df.at[ix, "cleanliness_tier"]),
            -float(res_df.at[ix, "direction_fraction"]),
            -float(res_df.at[ix, "focal_direction_fraction"]),
            -int(res_df.at[ix, "n_pairs"]),
            str(res_df.at[ix, "program"]),
        ),
    )

    kept_reps: list[int] = []
    for ix in passing_sorted:
        lin = str(res_df.at[ix, "lineage"])
        prog = str(res_df.at[ix, "program"])
        clust = str(res_df.at[ix, "mechanism_cluster"])
        gset = measured_sets.get(prog, set())
        duplicate_of: str | None = None

        for prev_ix in kept_reps:
            if str(res_df.at[prev_ix, "lineage"]) != lin:
                continue
            prev_prog = str(res_df.at[prev_ix, "program"])
            prev_clust = str(res_df.at[prev_ix, "mechanism_cluster"])
            prev_gset = measured_sets.get(prev_prog, set())
            ov = (
                len(gset & prev_gset) / max(1, min(len(gset), len(prev_gset)))
                if (gset and prev_gset)
                else 0.0
            )
            rho = 0.0
            if (lin, prog) in dt_piv.columns and (lin, prev_prog) in dt_piv.columns:
                pair_vals = dt_piv[[(lin, prog), (lin, prev_prog)]].dropna()
                if len(pair_vals) >= 6:
                    r_val = spearmanr(pair_vals.iloc[:, 0], pair_vals.iloc[:, 1]).statistic
                    rho = float(abs(r_val)) if pd.notna(r_val) else 0.0
            if clust == prev_clust or ov >= 0.30 or rho >= 0.75:
                duplicate_of = prev_prog
                break

        if duplicate_of is not None:
            res_df.at[ix, "dedup_status"] = f"secondary_collinear_with:{duplicate_of}"
            res_df.at[ix, "screen_reason"] = f"collinear_or_same_mechanism_as:{duplicate_of}"
        else:
            kept_reps.append(ix)
            res_df.at[ix, "dedup_status"] = "primary_cluster_representative"

    # Rank independent cluster representatives by statistical concordance & precision; cap at <= 3
    reps_ranked = sorted(
        kept_reps,
        key=lambda ix: (
            -float(res_df.at[ix, "direction_fraction"]),
            -float(res_df.at[ix, "focal_direction_fraction"]),
            -int(res_df.at[ix, "n_pairs"]),
            int(res_df.at[ix, "cleanliness_tier"]),
            str(res_df.at[ix, "program"]),
        ),
    )
    for rank_pos, ix in enumerate(reps_ranked):
        if rank_pos < 3:
            res_df.at[ix, "shortlisted"] = True
        else:
            res_df.at[ix, "dedup_status"] = "passed_exceeded_top3_cap"
            res_df.at[ix, "screen_reason"] = "exceeded_prespecified_max_3_primary_candidates"

    res_df = res_df.sort_values(
        ["shortlisted", "gates_1_to_4_pass", "direction_fraction", "focal_direction_fraction", "n_pairs", "program"],
        ascending=[False, False, False, False, False, True],
        kind="mergesort",
    ).reset_index(drop=True)
    return res_df


def write_evidence_markdown(result: pd.DataFrame, nmf_df: pd.DataFrame) -> None:
    report = [
        "# GSE245552 发现集去偏完整候选筛选草案（V1.1）",
        "",
        "本表严格在 `GSE245552` 开发集内基于 HGNC 规范化 `29,569` 肿瘤核基因轴与检测率分箱背景扣除 `debiased_ucell_scores` 生成（主候选上限 2–3 个）。",
        "",
        "| 谱系 | 程序 | 配对 N | 去偏均值 D_T (95% CI) | 同向率 | log1p 同向率 | 亚型内状态分量 | 纯焦点细类 (同向率) | 筛选结论 |",
        "|---|---|---:|---|---:|---:|---:|---|---|",
    ]
    for r in result.itertuples(index=False):
        n_conc = round(r.n_nonzero * r.direction_fraction)
        f_conc = f"{r.focal_subtype}: {r.focal_direction_fraction * 100:.0f}% (N={r.focal_n_pairs})"
        verdict = (
            "**入选主短名单**"
            if r.shortlisted
            else (f"同机制次要代表 ({r.dedup_status.split(':', 1)[-1]})" if r.gates_1_to_4_pass else r.screen_reason)
        )
        report.append(
            f"| {r.lineage} | `{r.program}` | {r.n_pairs} | "
            f"{r.mean_D_T:+.5f} ([{r.bootstrap_ci_low:+.5f}, {r.bootstrap_ci_high:+.5f}]) | "
            f"{n_conc}/{r.n_nonzero} ({r.direction_fraction * 100:.1f}%) | "
            f"{r.log1p_direction_fraction * 100:.1f}% | {r.full_within_subtype_state_D_T:+.5f} | "
            f"{f_conc} | {verdict} |"
        )
    report.extend(
        [
            "",
            f"- **无监督 Donor-Consensus NMF 审计**：共提取 `{len(nmf_df)}` 个共识模块（详见 [denovo_nmf_audit.csv](denovo_nmf_audit.csv)），与固定文献集合重合的模块均按已知程序收敛印证归档，不另造新名。",
            "- **双尺度与组成敏感性**：完整多尺度对照见 [sensitivity_summary.csv](patient_effects/sensitivity_summary.csv)，零 `unresolved` 对称细类分解见 [patient_subtype_decomposition.csv](composition_results/patient_subtype_decomposition.csv)，纯细类患者效应见 [fine_subtypes/effects_summary.csv](patient_effects/fine_subtypes/effects_summary.csv)。",
        ]
    )
    (OUT / "candidate_evidence_draft.md").write_text("\n".join(report) + "\n", encoding="utf-8")


def run() -> None:
    start = time.monotonic()
    started = datetime.now(timezone.utc).isoformat()
    manifest = pd.read_csv(MANIFEST)
    if len(manifest) != 39 or set(manifest.dataset) != {"GSE245552"} or set(manifest.role) != {"development"}:
        raise ValueError("Expected the validated 39-sample GSE245552 discovery manifest")
    if not (ROOT / "artifacts/m2/pseudobulk_manifest.csv").exists():
        raise FileNotFoundError("M2 has not published the complete pseudobulk manifest")

    genes, axis_245, axis_shared_all = build_program_definitions(write_files=True)
    info, sets = build_program_variants(genes, write_files=True)
    _, controls_by_lineage = build_frozen_reference_and_controls(
        None, manifest, axis_245, axis_shared_all, info, sets, write_files=True
    )
    info.to_csv(DEFS / "scored_definitions.csv", index=False)

    raw, comp, subtype, fine_samples, tech, precision, nmf_pools = score_libraries(
        manifest, axis_245, axis_shared_all, info, sets, controls_by_lineage
    )

    composition = OUT / "composition_results"
    composition.mkdir(parents=True, exist_ok=True)
    raw.to_csv(OUT / "sample_program_measurements.csv.gz", index=False)
    fine_samples.to_csv(OUT / "fine_subtype_program_scores.csv", index=False)
    comp.to_csv(composition / "sample_proportions.csv", index=False)
    tech.to_csv(OUT / "depth_sensitivity.csv", index=False)
    precision.to_csv(OUT / "split_half_precision.csv", index=False)
    subtype.to_csv(composition / "subtype_contributions.csv", index=False)

    decomp_df = decompose(manifest, subtype)
    decomp_df.to_csv(composition / "patient_subtype_decomposition.csv", index=False)

    # Immune composition contrasts on resolved-immune denominator (primary) + core-3 + all-retained
    res_imm_comp = comp.loc[comp["lineage"].isin(RESOLVED_IMMUNE)].assign(
        program="resolved_immune_fraction", score=comp["prop_total_resolved_immune"]
    )
    patient_outputs(manifest, res_imm_comp, composition)
    core3_comp = comp.loc[comp["lineage"].isin(PARENTS)].assign(
        program="core_three_immune_fraction", score=comp["prop_core_three_immune"]
    )
    patient_outputs(manifest, core3_comp, composition / "core_three_immune")
    all_ret_comp = comp.assign(program="captured_cell_fraction", score=comp["proportion"])
    patient_outputs(manifest, all_ret_comp, composition / "all_retained")

    where = OUT / "patient_effects"
    sensitivity, loo, main_patients, main_effects, fine_effects = state_summaries(
        manifest, raw, fine_samples, info, where
    )
    sensitivity.to_csv(where / "sensitivity_summary.csv", index=False)

    nmf_df = run_donor_consensus_nmf(nmf_pools, axis_245, genes)
    nmf_df.to_csv(OUT / "denovo_nmf_audit.csv", index=False)

    result = shortlist_candidates(
        info, genes, sensitivity, loo, main_patients, main_effects, decomp_df, fine_effects, raw, manifest
    )
    result.to_csv(OUT / "candidate_screen.csv", index=False)

    # Paired patient-level depth delta vs program delta on primary ucell_debiased vs legacy_positive
    depth_rows = []
    tech_deb = tech.loc[(tech["n_cells"] >= 50) & tech["scoring_mode"].eq("ucell_debiased")]
    for metric, statistical in tech_deb.groupby("qc_metric"):
        qc = statistical.merge(
            manifest[["sample_uid", "donor_uid", "tissue", "lesion_uid"]],
            on="sample_uid",
            validate="many_to_one",
        )
        qc_piv = qc.groupby(["donor_uid", "lineage", "program", "tissue"])["median_qc"].mean().unstack("tissue")
        delta = (qc_piv.LM - qc_piv.PT).rename("depth_D_T").reset_index()
        joined = main_patients.merge(delta, on=["donor_uid", "lineage", "program"]).dropna(subset=["depth_D_T", "D_T"])
        for (lineage, program), g in joined.groupby(["lineage", "program"]):
            rho = (
                float(spearmanr(g.depth_D_T, g.D_T).statistic)
                if min(g.depth_D_T.nunique(), g.D_T.nunique()) > 1
                else np.nan
            )
            depth_rows.append(
                {
                    "qc_metric": metric,
                    "lineage": lineage,
                    "program": program,
                    "n_pairs": len(g),
                    "spearman_depth_delta_program_delta": rho,
                }
            )
    pd.DataFrame(depth_rows).to_csv(OUT / "paired_depth_sensitivity.csv", index=False)

    write_evidence_markdown(result, nmf_df)

    inputs = [
        Path(__file__),
        ROOT / "src/crlm/program_scores.py",
        ROOT / "src/crlm/hgnc_alias_map.py",
        MANIFEST,
        SEEDS,
        S4_TXT,
        BACK,
        HALLMARK,
        M1_RAW_AXIS,
        ROOT / "starter/paired_effects.py",
        ROOT / "starter/composition_decomposition.py",
    ]
    for gsm in manifest.sample_uid:
        inputs.append(M2 / f"{gsm}.h5ad")
    profile = {
        "started_utc": started,
        "ended_utc": datetime.now(timezone.utc).isoformat(),
        "seed": SEED,
        "scoring_method_primary": "ucell_debiased_20bin_10ctrl",
        "ranking_universe": "245_tumor_canonical_nuclear",
        "n_ranking_genes": len(axis_245),
        "n_shared_nuclear_genes": len(axis_shared_all),
        "n_complete_programs_screened": int(len(result)),
        "command": " ".join([sys.executable, *sys.argv]),
        "inputs_sha256": {str(p.relative_to(ROOT)): sha(p) for p in inputs if p.is_file()},
        "peak_rss_bytes": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024,
        "elapsed_seconds": round(time.monotonic() - start, 3),
        "discovery_samples": len(manifest),
        "retained_cells": int(comp.loc[comp["lineage"].eq("T/NK"), "n_retained"].sum()),
        "development_patients": int(manifest.donor_uid.nunique()),
        "shortlisted": result.loc[result.shortlisted, ["lineage", "program", "direction", "mean_D_T"]].to_dict("records"),
    }
    (OUT / "resource_profile.json").write_text(json.dumps(profile, indent=2) + "\n", encoding="utf-8")
    print(
        result[
            [
                "lineage", "program", "n_pairs", "mean_D_T", "direction_fraction",
                "within_state_share", "focal_subtype", "focal_direction_fraction", "shortlisted",
            ]
        ].to_string(index=False),
        flush=True,
    )


def self_check() -> None:
    """Pure in-memory verification of catalog, debiased scoring, Gate 4 UPR rejection, Gate 5 dedup, and NMF."""
    genes, axis_245, axis_shared_all = build_program_definitions(write_files=False)
    info, sets = build_program_variants(genes, write_files=False)
    complete_orig = info.loc[(info["role"] == "complete_program") & (info["variant"] == "original")]
    assert len(complete_orig) == 19, f"Expected 19 complete programs, got {len(complete_orig)}"
    assert (complete_orig["eligible"]).all(), "All 19 complete programs must be eligible on 245 nuclear axis"
    # Verify cross-lineage T/NK contaminants stripped from B19_S4_Mac_SPP1_clean
    spp1_clean_set = set(sets["B19_S4_Mac_SPP1_clean"])
    assert not (spp1_clean_set & {"GZMB", "GNLY", "LAT", "IL2RG"})

    # Construct compact 600-gene canonical toy universe covering all 19 programs + background genes
    prog_union = sorted({g for p in complete_orig["score_id"] for g in sets[p]})
    bg_extra = [g for g in axis_shared_all if g not in set(prog_union)][:400]
    toy_axis = sorted(set(prog_union) | set(bg_extra))
    toy_info = complete_orig.copy()
    toy_sets = {r.score_id: sets[r.score_id] for r in toy_info.itertuples(index=False)}

    rng = np.random.default_rng(SEED)
    n_donors = 9
    manifest_rows = []
    objs: dict[str, ad.AnnData] = {}
    var_df = pd.DataFrame({"gene_symbol": toy_axis, "canonical_symbol": toy_axis}, index=toy_axis)

    sam_idx = pd.Index(toy_axis).get_indexer(sets["Core_Lipid_SAM_R"])
    spp1_idx = pd.Index(toy_axis).get_indexer(sets["B19_S4_Mac_SPP1_clean"])
    upr_idx = pd.Index(toy_axis).get_indexer(sets["HALLMARK_UNFOLDED_PROTEIN_RESPONSE"])

    for d_i in range(n_donors):
        donor = f"GSE245552:ToyPt{d_i + 1:02d}"
        for tissue in ("PT", "LM"):
            gsm = f"GSM_TOY_{d_i + 1:02d}_{tissue}"
            manifest_rows.append(
                {
                    "dataset": "GSE245552",
                    "donor_uid": donor,
                    "sample_uid": gsm,
                    "lesion_uid": f"{donor}_{tissue}_L1",
                    "tissue": tissue,
                    "role": "development",
                    "input_scale": "raw_counts",
                    "patient_mapping_verified": True,
                }
            )
            # 60 T/NK (30 CD8, 30 CD4), 60 Mono/Mac/DC (45 Macrophage, 15 Monocyte), 60 B/Plasma
            # Make B/Plasma shift purely compositional: PT = 50 plasma + 10 B; LM = 10 plasma + 50 B,
            # with identical within-plasma and within-B UPR expression across PT and LM!
            n_pl = 50 if tissue == "PT" else 10
            n_bc = 60 - n_pl
            lin_vec = ["T/NK"] * 60 + ["Mono/Mac/DC"] * 60 + ["B/Plasma"] * 60
            sub_vec = (
                ["CD8 T cell"] * 30 + ["CD4 T cell"] * 30
                + ["monocyte/macrophage"] * 60
                + ["plasma cell"] * n_pl + ["B cell"] * n_bc
            )
            mf_vec = [""] * 60 + ["Macrophage"] * 45 + ["Monocyte"] * 15 + [""] * 60
            mat = rng.poisson(0.4, size=(180, len(toy_axis))).astype(np.int32)
            # Plasma cells constitutively express higher UPR than B cells in BOTH tissues (pure composition shift!)
            mat[120 : 120 + n_pl][:, upr_idx] += rng.poisson(3.0, size=(n_pl, len(upr_idx)))
            # True within-Macrophage LM upregulation of Core_Lipid_SAM_R & collinear B19_S4_Mac_SPP1_clean
            if tissue == "LM":
                mat[60:105][:, sam_idx] += rng.poisson(4.0, size=(45, len(sam_idx)))
                mat[60:105][:, spp1_idx[:25]] += rng.poisson(4.0, size=(45, 25))
            obs = pd.DataFrame(
                {
                    "sample_uid": gsm,
                    "donor_uid": donor,
                    "tissue": tissue,
                    "lineage_m2": lin_vec,
                    "subtype_m2": sub_vec,
                    "myeloid_fine": mf_vec,
                    "qc_retained": True,
                    "qc_joint_low": False,
                    "qc_pass_g800": True,
                    "n_counts_visible": np.maximum(1, mat.sum(axis=1)),
                    "n_genes_visible": np.maximum(1, (mat > 0).sum(axis=1)),
                },
                index=[f"{gsm}:c{k}" for k in range(180)],
            )
            objs[gsm] = ad.AnnData(X=sparse.csr_matrix(mat), obs=obs, var=var_df.copy())

    toy_manifest = pd.DataFrame(manifest_rows)
    _, ctrls = build_frozen_reference_and_controls(
        objs, toy_manifest, toy_axis, toy_axis, toy_info, toy_sets, write_files=False
    )
    raw, comp, subtype, fine_samples, tech, precision, nmf_pools = score_libraries(
        toy_manifest, toy_axis, toy_axis, toy_info, toy_sets, ctrls, objects_by_gsm=objs, verbose=False
    )
    decomp_df = decompose(toy_manifest, subtype)
    sens, loo, main_pts, main_eff, fine_eff = state_summaries(toy_manifest, raw, fine_samples, toy_info, None)
    screen = shortlist_candidates(
        toy_info, genes.loc[genes["gene_symbol"].isin(toy_axis)], sens, loo, main_pts, main_eff, decomp_df, fine_eff, raw, toy_manifest
    )
    nmf_df = run_donor_consensus_nmf(nmf_pools, toy_axis, genes.loc[genes["gene_symbol"].isin(toy_axis)], n_components=3, n_runs=4)

    by_prog = screen.set_index("program")
    # 1) Pure compositional UPR shift MUST be rejected by Gate 4!
    assert not bool(by_prog.at["HALLMARK_UNFOLDED_PROTEIN_RESPONSE", "shortlisted"])
    upr_reason = str(by_prog.at["HALLMARK_UNFOLDED_PROTEIN_RESPONSE", "screen_reason"])
    assert ("composition_dominated" in upr_reason) or ("within_subtype_state_discordant" in upr_reason) or ("focal_subtype" in upr_reason), upr_reason
    # 2) True within-Macrophage Core_Lipid_SAM_R MUST pass Gates 1-4 and be shortlisted, while collinear B19_S4_Mac_SPP1_clean is deduplicated by Gate 5!
    assert bool(by_prog.at["Core_Lipid_SAM_R", "shortlisted"])
    if bool(by_prog.at["B19_S4_Mac_SPP1_clean", "gates_1_to_4_pass"]):
        assert not bool(by_prog.at["B19_S4_Mac_SPP1_clean", "shortlisted"])
        assert "secondary_collinear_with:Core_Lipid_SAM_R" == str(by_prog.at["B19_S4_Mac_SPP1_clean", "dedup_status"])
    # 3) Consensus NMF recovers Core_Lipid_SAM_R convergence!
    assert (nmf_df["best_matched_literature_program"] == "Core_Lipid_SAM_R").any()
    print(
        f"analyze_discovery_programs --self-check PASSED: 19 complete programs verified, "
        f"UPR composition artifact rejected ({upr_reason}), "
        f"Core_Lipid_SAM_R shortlisted & collinear Mac_SPP1 deduplicated, "
        f"NMF factors={len(nmf_df)} (zero disk writes)."
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--definitions-only", action="store_true", help="Write draft program definitions only.")
    parser.add_argument("--self-check", action="store_true", help="Run pure in-memory verification with zero disk writes.")
    args = parser.parse_args()
    if args.self_check:
        self_check()
    elif args.definitions_only:
        build_program_variants(build_program_definitions(write_files=True)[0], write_files=True)
    else:
        run()
