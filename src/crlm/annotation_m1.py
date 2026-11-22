"""Two-stage hierarchical cell-type annotation and doublet stripping."""

from __future__ import annotations

from pathlib import Path
from typing import Mapping

import numpy as np
import pandas as pd
from scipy import sparse

LEVELS = (
    "T/NK",
    "Mono/Mac/DC",
    "B/Plasma",
    "Neutrophil",
    "Mast",
    "Nonimmune",
    "Mixed",
    "Unknown",
)

# Stage-1 broad-lineage panels: (broad_lineage, identity_anchors, compatible_support).
# Each panel requires >=1 detected anchor gene and >=2 detected panel genes total.
# CSF3R is intentionally excluded from neutrophil anchors because classical CD14+
# monocytes and inflammatory macrophages physiologically co-express CSF3R + S100A8/9.
PANELS: dict[str, tuple[str, str, str]] = {
    "T cell": ("T/NK", "CD3D CD3E CD3G TRAC TRBC1 TRBC2", "CD2 CD7 LCK IL7R CD4 CD8A CD8B"),
    "NK cell": ("T/NK", "NCR1 KLRF1 KIR2DL1 KIR2DL3 KLRD1 KLRK1", "GNLY NKG7 PRF1 CTSW CD247 FCGR3A"),
    "monocyte/macrophage": (
        "Mono/Mac/DC",
        "CSF1R MS4A7 C1QA C1QB C1QC CD163 CD14",
        "FCGR3A FCN1 LST1 AIF1 CD68 S100A8 S100A9 CSF3R",
    ),
    "conventional DC": ("Mono/Mac/DC", "CLEC9A XCR1 CD1C CLEC10A", "FLT3 FCER1A CST3"),
    "plasmacytoid DC": ("Mono/Mac/DC", "CLEC4C LILRA4", "IL3RA PLD4 TCF4"),
    "mature DC": ("Mono/Mac/DC", "LAMP3", "FSCN1 CCL22 CD86 CCR7"),
    "B cell": ("B/Plasma", "MS4A1 CD19 CD79A CD79B", "BANK1 VPREB3 FCRLA CD22"),
    "plasma cell": ("B/Plasma", "MZB1 TNFRSF17 DERL3 SDC1", "JCHAIN IGKC"),
    "neutrophil": ("Neutrophil", "FCGR3B CEACAM8 CXCR2 CMTM2", "G0S2 CXCL8 PROK2"),
    "mast cell": ("Mast", "TPSAB1 TPSB2 CPA3", "KIT HDC MS4A2"),
    "epithelial": ("Nonimmune", "EPCAM CLDN3 CLDN4 KRT19", "KRT8 KRT18 ELF3"),
    "hepatocyte": ("Nonimmune", "ASGR1 ASGR2 HNF4A", "CYP2E1 CPS1 FTCD"),
    "endothelial": ("Nonimmune", "VWF CLDN5 EMCN PLVAP", "PECAM1 KDR RAMP2"),
    "fibroblast": ("Nonimmune", "COL1A1 COL1A2 COL3A1 PDGFRA", "DCN LUM COL6A2"),
    "mural": ("Nonimmune", "RGS5 CSPG4 MYH11", "ACTA2 TAGLN PDGFRB"),
}

# Auxiliary markers used strictly for cross-lineage doublet checks or Stage-2 fine splits.
AUXILIARY_MARKERS = ("PTPRC", "FOXP3", "APOE", "MRC1", "MARCO", "VCAN", "S100A12")

HIERARCHICAL_SUBTYPES: Mapping[str, tuple[str, ...]] = {
    "T/NK": ("CD8 T cell", "CD4 T cell", "Double_neg_or_unsplit_T", "NK cell"),
    "Mono/Mac/DC": ("monocyte/macrophage", "conventional DC", "plasmacytoid DC", "mature DC"),
    "B/Plasma": ("B cell", "plasma cell"),
    "Neutrophil": ("neutrophil",),
    "Mast": ("mast cell",),
    "Nonimmune": ("epithelial", "hepatocyte", "endothelial", "fibroblast", "mural"),
    "Mixed": ("Cross_lineage_doublet",),
    "Unknown": ("Unknown",),
}

# Independent corroborating markers completely disjoint from PANELS and AUXILIARY_MARKERS.
HELD_OUT_CORROBORATION_GENES: Mapping[str, tuple[str, ...]] = {
    "T/NK": ("CD6", "CCL5", "GZMA", "ITK", "SIT1"),
    "Mono/Mac/DC": ("LYZ", "MS4A6A", "MAFB", "SPI1", "TYROBP"),
    "B/Plasma": ("IGHA1", "IGHA2", "IGHG1", "CD72", "POU2AF1"),
    "Neutrophil": ("PI3", "MMP9", "SELL"),
    "Mast": ("TPSD1", "GATA2", "LTC4S"),
    "Nonimmune": ("TFF3", "CEACAM5", "PHGR1", "SPARCL1", "IGFBP7", "MGP"),
}


def _s3_lineages(marker_table: Path | str | None) -> dict[str, set[str]]:
    if marker_table is None:
        return {}
    path = Path(marker_table)
    if not path.is_file():
        return {}
    rows = pd.read_csv(path, encoding="utf-8-sig")
    by_table = {1: "Nonimmune", 2: "T/NK", 3: "Mono/Mac/DC", 4: "B/Plasma", 5: "Nonimmune"}
    exceptions = {"Mast_cell": "Mast", "Neutrophil": "Neutrophil"}
    result: dict[str, set[str]] = {}
    for row in rows.itertuples(index=False):
        lineage = exceptions.get(row.author_group, by_table[row.table_index])
        result.setdefault(row.gene, set()).add(lineage)
    return result


def _format_value(val: float | int) -> str:
    if isinstance(val, (int, np.integer)) or float(val).is_integer():
        return str(int(val))
    return f"{float(val):.4g}"


def classify(counts, symbols, marker_table: Path | str | None = None):
    """Return hierarchical broad labels and row-aligned subtype/marker evidence.

    Stage 1 resolves broad lineages (`T/NK`, `Mono/Mac/DC`, `B/Plasma`,
    `Neutrophil`, `Mast`, `Nonimmune`, `Mixed`, `Unknown`), decoupling classical
    `CSF3R+` monocytes/macrophages from `FCGR3B+/CEACAM8+` neutrophils while
    routing explicit cross-lineage doublets (`CD3+CD79A/MS4A1+`, `PTPRC+EPCAM+`,
    `Immune+EPCAM/ASGR1+`, `CD3+multi-myeloid+`, `B+multi-myeloid+`) to `Mixed`
    (`subtype_hierarchical = "Cross_lineage_doublet"`).

    Stage 2 deterministically disambiguates within-lineage subtypes without
    emitting artificial `"* unresolved"` buckets.
    """
    panel_genes = {
        gene
        for _, anchors, support in PANELS.values()
        for gene in (anchors + " " + support).split()
    }
    all_genes = sorted(panel_genes | set(AUXILIARY_MARKERS))
    symbols_arr = np.asarray(list(symbols), dtype=str)

    if not sparse.issparse(counts):
        counts_csr = sparse.csr_matrix(counts)
    else:
        counts_csr = counts.tocsr()

    if counts_csr.shape[1] != len(symbols_arr):
        raise ValueError(f"counts columns ({counts_csr.shape[1]}) != symbols ({len(symbols_arr)})")

    selected = np.flatnonzero(np.isin(symbols_arr, all_genes))
    genes = symbols_arr[selected]
    gene_to_col = {g: idx for idx, g in enumerate(genes)}

    matrix = counts_csr[:, selected].tocsr()
    matrix.eliminate_zeros()
    matrix.sort_indices()
    dense = matrix.toarray()
    if (dense < 0).any() or not np.isfinite(dense).all():
        raise ValueError("Expression matrix for annotation must be finite and nonnegative")
    pos = dense > 0
    n_cells = matrix.shape[0]

    def any_pos(glist: list[str] | tuple[str, ...]) -> np.ndarray:
        cols = [gene_to_col[g] for g in glist if g in gene_to_col]
        return pos[:, cols].any(axis=1) if cols else np.zeros(n_cells, dtype=bool)

    def count_pos(glist: list[str] | tuple[str, ...]) -> np.ndarray:
        cols = [gene_to_col[g] for g in glist if g in gene_to_col]
        return pos[:, cols].sum(axis=1) if cols else np.zeros(n_cells, dtype=np.int32)

    def sum_log_expr(glist: list[str] | tuple[str, ...]) -> np.ndarray:
        cols = [gene_to_col[g] for g in glist if g in gene_to_col]
        return np.log1p(dense[:, cols]).sum(axis=1) if cols else np.zeros(n_cells, dtype=np.float64)

    panel_hits: dict[str, np.ndarray] = {}
    groups: list[list[str]] = [[] for _ in range(n_cells)]
    lineages: list[set[str]] = [set() for _ in range(n_cells)]

    for panel, (lineage, anchors, support) in PANELS.items():
        has_anchor = any_pos(anchors.split())
        sufficient = count_pos((anchors + " " + support).split()) >= 2
        hit = has_anchor & sufficient
        panel_hits[panel] = hit
        for cell in np.flatnonzero(hit):
            groups[cell].append(panel)
            lineages[cell].add(lineage)

    # Stage-1 decoupling & cross-lineage doublet masks
    mono_dc_core = any_pos(
        ("C1QA", "C1QB", "C1QC", "CD163", "MS4A7", "CSF1R", "CLEC9A", "XCR1", "CD1C", "CLEC10A", "LAMP3", "CLEC4C", "LILRA4")
    )
    hard_gran = any_pos(("FCGR3B", "CEACAM8"))
    cd14_pos = any_pos(("CD14",))
    mye_panel_score = count_pos(
        ("CSF1R", "MS4A7", "C1QA", "C1QB", "C1QC", "CD163", "CD14", "FCGR3A", "FCN1", "LST1", "AIF1", "CD68", "CSF3R")
    )
    neu_panel_score = count_pos(("FCGR3B", "CEACAM8", "CXCR2", "CMTM2", "G0S2", "CXCL8", "PROK2"))

    cd3_tcr = any_pos(("CD3D", "CD3E", "CD3G", "TRAC"))
    trbc_pos = any_pos(("TRBC1", "TRBC2"))
    b_core_anchor = any_pos(("MS4A1", "CD19", "CD79A", "TNFRSF17")) | (
        any_pos(("MZB1",)) & any_pos(("JCHAIN", "SDC1", "DERL3"))
    )
    mye_multi = (
        count_pos(
            ("CSF1R", "MS4A7", "C1QA", "C1QB", "C1QC", "CD163", "CD14", "FCN1", "LST1", "CD68", "AIF1", "CLEC9A", "CD1C", "LILRA4")
        )
        >= 2
    )
    ptprc_pos = any_pos(("PTPRC",))
    # Note: ASGR1/ASGR2 alone are endocytic C-type lectins natively expressed by
    # human monocytes/macrophages; true hepatocyte contamination requires HNF4A
    # or >=2 hepatocyte panel markers (captured by panel_hits["hepatocyte"]).
    epi_hep_anchor = (
        any_pos(("EPCAM", "HNF4A"))
        | panel_hits["hepatocyte"]
        | (any_pos(("ASGR1", "ASGR2")) & any_pos(("CYP2E1", "CPS1", "FTCD")))
    )
    krt19_epi = any_pos(("KRT19",)) & any_pos(("CLDN3", "CLDN4", "KRT8", "KRT18", "ELF3"))

    # Stage-2 within-lineage disambiguation vectors
    cd8_pos = any_pos(("CD8A", "CD8B"))
    cd4_gene = any_pos(("CD4",))
    foxp3_pos = any_pos(("FOXP3",))
    il7r_pos = any_pos(("IL7R",))
    cyto_pos = any_pos(("GNLY", "PRF1", "NKG7"))
    cd8_expr = sum_log_expr(("CD8A", "CD8B"))
    cd4_expr = sum_log_expr(("CD4", "FOXP3"))

    b_core_count = count_pos(("MS4A1", "CD19", "BANK1", "CD22", "FCRLA", "VPREB3"))
    plasma_core_count = count_pos(("MZB1", "JCHAIN", "SDC1", "TNFRSF17", "DERL3"))
    b_core_expr = sum_log_expr(("MS4A1", "CD19", "BANK1", "CD22", "FCRLA", "VPREB3"))
    plasma_core_expr = sum_log_expr(("MZB1", "JCHAIN", "SDC1", "TNFRSF17", "DERL3"))
    ms4a1_pos = any_pos(("MS4A1",))

    pdc_n = count_pos(("CLEC4C", "LILRA4", "IL3RA", "PLD4", "TCF4"))
    pdc_anc = any_pos(("CLEC4C", "LILRA4"))
    mdc_n = count_pos(("LAMP3", "FSCN1", "CCL22", "CD86", "CCR7"))
    mdc_anc = any_pos(("LAMP3",)) & any_pos(("FSCN1", "CCL22", "CCR7", "CD86"))
    cdc1_anc = any_pos(("CLEC9A", "XCR1"))
    cdc2_anc = any_pos(("CD1C", "CLEC10A")) & any_pos(("FLT3", "FCER1A", "CLEC10A", "CD1C"))
    mac_core_pos = any_pos(("C1QA", "C1QB", "C1QC", "CD163", "MRC1", "MARCO", "APOE"))
    mono_sig_pos = any_pos(("FCN1", "S100A8", "S100A9", "S100A12", "VCAN"))

    nonimm_panels = ("epithelial", "hepatocyte", "endothelial", "fibroblast", "mural")
    nonimm_scores = np.column_stack(
        [
            count_pos((PANELS[p][1] + " " + PANELS[p][2]).split()) * 10.0
            + sum_log_expr((PANELS[p][1] + " " + PANELS[p][2]).split())
            for p in nonimm_panels
        ]
    )

    s3 = _s3_lineages(marker_table)
    panel_gene_mask = np.isin(genes, sorted(panel_genes))
    immune_set = {"T/NK", "Mono/Mac/DC", "B/Plasma", "Neutrophil", "Mast"}

    labels: list[str] = []
    subtypes: list[str] = []
    myeloid_fine: list[str] = []
    evidence: list[str] = []
    s3_evidence: list[str] = []

    for cell in range(n_cells):
        sup = set(lineages[cell])
        if "Mono/Mac/DC" in sup and hard_gran[cell]:
            if not mono_dc_core[cell]:
                sup.discard("Mono/Mac/DC")
                sup.add("Neutrophil")
            else:
                sup.add("Neutrophil")
        elif sup == {"Mono/Mac/DC", "Neutrophil"}:
            if not hard_gran[cell] and (mono_dc_core[cell] or mye_panel_score[cell] >= neu_panel_score[cell]):
                sup = {"Mono/Mac/DC"}

        is_imm_cand = bool(sup & immune_set)
        cross_doublet = (
            (cd3_tcr[cell] and b_core_anchor[cell])
            or (cd3_tcr[cell] and mye_multi[cell])
            or (b_core_anchor[cell] and mye_multi[cell])
            or (ptprc_pos[cell] and epi_hep_anchor[cell])
            or (is_imm_cand and (epi_hep_anchor[cell] or krt19_epi[cell]))
        )

        mf = ""
        if cross_doublet or len(sup) > 1:
            label = "Mixed"
            sub = "Cross_lineage_doublet"
        elif not sup:
            label = "Unknown"
            sub = "Unknown"
        else:
            label = next(iter(sup))
            if label == "T/NK":
                is_t = cd3_tcr[cell] or (trbc_pos[cell] and not panel_hits["NK cell"][cell])
                if is_t:
                    if cd8_pos[cell] and not cd4_gene[cell]:
                        sub = "CD8 T cell"
                    elif cd4_gene[cell] and not cd8_pos[cell]:
                        sub = "CD4 T cell"
                    elif cd8_pos[cell] and cd4_gene[cell]:
                        if cd8_expr[cell] > cd4_expr[cell]:
                            sub = "CD8 T cell"
                        elif cd4_expr[cell] > cd8_expr[cell]:
                            sub = "CD4 T cell"
                        else:
                            sub = "Double_neg_or_unsplit_T"
                    elif foxp3_pos[cell] or (il7r_pos[cell] and not cyto_pos[cell]):
                        sub = "CD4 T cell"
                    else:
                        sub = "Double_neg_or_unsplit_T"
                else:
                    sub = "NK cell"
            elif label == "B/Plasma":
                if not ms4a1_pos[cell] and plasma_core_count[cell] >= 2:
                    sub = "plasma cell"
                elif plasma_core_count[cell] > b_core_count[cell] or (
                    plasma_core_count[cell] == b_core_count[cell] and plasma_core_expr[cell] > b_core_expr[cell]
                ):
                    sub = "plasma cell"
                else:
                    sub = "B cell"
            elif label == "Mono/Mac/DC":
                if pdc_anc[cell] and pdc_n[cell] >= 2 and not mac_core_pos[cell]:
                    sub = mf = "plasmacytoid DC"
                elif mdc_anc[cell] and mdc_n[cell] >= 2 and not mac_core_pos[cell]:
                    sub = mf = "mature DC"
                elif cdc1_anc[cell] or (cdc2_anc[cell] and not mac_core_pos[cell] and not cd14_pos[cell]):
                    sub = mf = "conventional DC"
                else:
                    sub = "monocyte/macrophage"
                    mf = "Macrophage" if mac_core_pos[cell] or not mono_sig_pos[cell] else "Monocyte"
            elif label == "Neutrophil":
                sub = "neutrophil"
            elif label == "Mast":
                sub = "mast cell"
            elif label == "Nonimmune":
                hit_indices = [i for i, p in enumerate(nonimm_panels) if panel_hits[p][cell]]
                best_idx = max(hit_indices, key=lambda i: nonimm_scores[cell, i])
                sub = nonimm_panels[best_idx]
            else:
                sub = "Unknown"

        labels.append(label)
        subtypes.append(sub)
        myeloid_fine.append(mf)

        start, stop = matrix.indptr[cell : cell + 2]
        cols_c = matrix.indices[start:stop]
        vals_c = matrix.data[start:stop]
        expressed = genes[cols_c]
        evidence.append(" ".join(f"{gene}={_format_value(val)}" for gene, val in zip(expressed, vals_c)))
        if s3:
            matched: dict[str, list[str]] = {}
            for gene, is_p in zip(expressed, panel_gene_mask[cols_c]):
                if not is_p:
                    continue
                for s3_lin in s3.get(gene, ()):
                    matched.setdefault(s3_lin, []).append(gene)
            s3_evidence.append(" | ".join(f"{key}: {' '.join(matched[key])}" for key in sorted(matched)))
        else:
            s3_evidence.append("")

    details = pd.DataFrame(
        {
            "subtype_hierarchical": subtypes,
            "myeloid_fine": myeloid_fine,
            "m1_marker_groups": [" | ".join(names) for names in groups],
            "m1_marker_evidence": evidence,
            "m1_s3_marker_evidence": s3_evidence,
        }
    )
    return np.asarray(labels, dtype=str), details
