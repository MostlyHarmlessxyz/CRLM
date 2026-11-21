"""Provisional broad lineage calls on the deposited GSE245552 counts."""

import numpy as np
import pandas as pd


# Each panel requires two detected genes from its two columns combined, including
# at least one in the first column. Separate panels keep, for example, endothelial
# and epithelial transcripts from combining into a single nonimmune call.
# panel: (broad lineage, identity markers, additional compatible markers)
PANELS = {
    "T cell": ("T/NK", "CD3D CD3E CD3G TRAC TRBC1 TRBC2", "CD2 CD7 LCK"),
    "NK/cytotoxic lymphocyte": ("T/NK", "KLRD1 GNLY NCR1 KLRK1", "NKG7 PRF1 CTSW CD247"),
    "monocyte/macrophage": ("Mono/Mac/DC", "CSF1R MS4A7 C1QA C1QB C1QC CD163 CD14",
                            "FCN1 LST1 AIF1 CD68"),
    "conventional DC": ("Mono/Mac/DC", "CLEC9A XCR1 CD1C CLEC10A", "FLT3 FCER1A CST3"),
    "plasmacytoid DC": ("Mono/Mac/DC", "CLEC4C LILRA4", "IL3RA PLD4 TCF4"),
    "mature DC": ("Mono/Mac/DC", "LAMP3", "FSCN1 CCL22 CD86"),
    "B cell": ("B/Plasma", "MS4A1 CD19 CD79A CD79B", "BANK1 VPREB3 FCRLA"),
    "plasma cell": ("B/Plasma", "MZB1 TNFRSF17 DERL3", "JCHAIN SDC1 IGKC"),
    "neutrophil": ("Neutrophil", "FCGR3B CSF3R CEACAM8 CXCR2", "FPR1 S100A8 S100A9"),
    "mast cell": ("Mast", "TPSAB1 TPSB2 CPA3", "KIT HDC MS4A2"),
    "epithelial": ("Nonimmune", "EPCAM CLDN3 CLDN4 KRT19", "KRT8 KRT18 ELF3"),
    "hepatocyte": ("Nonimmune", "ASGR1 ASGR2 HNF4A", "CYP2E1 CPS1 FTCD"),
    "endothelial": ("Nonimmune", "VWF CLDN5 EMCN PLVAP", "PECAM1 KDR RAMP2"),
    "fibroblast": ("Nonimmune", "COL1A1 COL1A2 COL3A1 PDGFRA", "DCN LUM COL6A2"),
    "mural": ("Nonimmune", "RGS5 CSPG4 MYH11", "ACTA2 TAGLN PDGFRB"),
}


def _s3_lineages(marker_table):
    rows = pd.read_csv(marker_table, encoding="utf-8-sig")
    by_table = {1: "Nonimmune", 2: "T/NK", 3: "Mono/Mac/DC", 4: "B/Plasma", 5: "Nonimmune"}
    exceptions = {"Mast_cell": "Mast", "Neutrophil": "Neutrophil"}
    result = {}
    for row in rows.itertuples(index=False):
        lineage = exceptions.get(row.author_group, by_table[row.table_index])
        result.setdefault(row.gene, set()).add(lineage)
    return result


def classify(counts, symbols, marker_table):
    """Return broad labels and row-aligned evidence for each input cell.

    Input is a cells-by-genes CSR of deposited integer counts with unique symbols.
    A cell supporting several panels in one lineage receives the lineage label.
    Competing lineage calls remain Mixed; cells without a qualifying panel are
    Unknown. Mixed is a marker conflict, not an inferred physical doublet.
    S3 overlap is a cross-reference to the author's gene column, separate from
    the call itself. All evidence uses visible counts with no normalization.
    """
    marker_genes = {gene for _, anchors, support in PANELS.values()
                    for gene in (anchors + " " + support).split()}
    symbols = np.asarray(symbols, dtype=str)
    selected = np.flatnonzero(np.isin(symbols, sorted(marker_genes)))
    genes = symbols[selected]
    matrix = counts[:, selected].tocsr()
    matrix.eliminate_zeros()
    matrix.sort_indices()
    groups = [[] for _ in range(matrix.shape[0])]
    lineages = [set() for _ in range(matrix.shape[0])]

    for panel, (lineage, anchors, support) in PANELS.items():
        anchor_cols = np.flatnonzero(np.isin(genes, anchors.split()))
        panel_cols = np.flatnonzero(np.isin(genes, (anchors + " " + support).split()))
        has_anchor = np.diff(matrix[:, anchor_cols].indptr) > 0
        sufficient = np.diff(matrix[:, panel_cols].indptr) >= 2
        for cell in np.flatnonzero(has_anchor & sufficient):
            groups[cell].append(panel)
            lineages[cell].add(lineage)

    s3 = _s3_lineages(marker_table)
    labels, evidence, s3_evidence = [], [], []
    for cell, support in enumerate(lineages):
        label = "Unknown"
        if len(support) == 1:
            label = next(iter(support))
        elif support:
            label = "Mixed"
        labels.append(label)
        start, stop = matrix.indptr[cell:cell + 2]
        expressed = genes[matrix.indices[start:stop]]
        values = matrix.data[start:stop]
        evidence.append(" ".join(f"{gene}={int(value)}" for gene, value in zip(expressed, values)))
        matched = {}
        for gene in expressed:
            for lineage in s3.get(gene, ()):
                matched.setdefault(lineage, []).append(gene)
        s3_evidence.append(" | ".join(f"{key}: {' '.join(matched[key])}" for key in sorted(matched)))

    details = pd.DataFrame({
        "m1_marker_groups": [" | ".join(names) for names in groups],
        "m1_marker_evidence": evidence,
        "m1_s3_marker_evidence": s3_evidence,
    })
    return np.asarray(labels, dtype=str), details
