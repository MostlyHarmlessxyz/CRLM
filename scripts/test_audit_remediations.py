#!/usr/bin/env python3
"""Automated regression test suite for CRLM audit remediations (Part 1 & extensible Part 2).

Verifies:
1. Deterministic HGNC symbol harmonization, 1-to-1 injectivity, >=98.8% external M4 PCA
   HVG recovery on GSE164522 & GSE203608, and strict disjointness of shared nuclear axes
   from the 353 tumor-deleted mitochondrial/ribosomal genes.
2. Two-stage hierarchical cell annotation (`src/crlm/annotation_m1.py`) & M2 materialized
   cohort (`CD14+CSF3R+FCGR3B=0 -> Mono/Mac/DC`, `CD3D+CD8A+GNLY+ -> CD8 T cell`,
   `MZB1+JCHAIN+CD79A+MS4A1=0 -> plasma cell`, cross-lineage doublets -> `Mixed`,
   `0` unresolved cells across all 39 M2 libraries, dual immune denominators summing to 1).
3. Exact full-axis Mann-Whitney UCell & detection-matched depth-debiased scoring
   (`src/crlm/program_scores.py`): synthetic Poisson depth independence (`max |rho| < 0.08`),
   anti-leakage `ValueError` guard when called without frozen controls/reference stats,
   sub-axis `log1p(CP10k)` renormalization invariance, and bit-exact V1 legacy compatibility.
4. Direct prior MIL bugfixes (`src/crlm/prior_mil.py`: in-repo `references/code/`, removal
   of `[::2]` half-bag drop, `1/sqrt(Q)` feature scaling, signed-logit additivity) plus
   M6/M7/Track-1 remediation contracts.
5. Full cryptographic SHA-256 verification across `starter/`, `artifacts/archive/CRLM_V1_245_F1/`,
   `artifacts/m0/`, `artifacts/m1/`, `artifacts/m2/`, `artifacts/m7/`, and
   `artifacts/techfix_v1_f1/{programs,models}/`.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import rankdata, spearmanr

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))
if str(ROOT / "starter") not in sys.path:
    sys.path.insert(1, str(ROOT / "starter"))
if str(ROOT / "scripts") not in sys.path:
    sys.path.insert(2, str(ROOT / "scripts"))

from crlm.annotation_m1 import AUXILIARY_MARKERS, HIERARCHICAL_SUBTYPES, LEVELS, PANELS, classify  # noqa: E402
from crlm.hgnc_alias_map import (  # noqa: E402
    HGNC_PROVENANCE_PATH,
    HGNC_TSV_PATH,
    SHARED_NUCLEAR_UNIVERSE,
    align_axis_to_canonical_targets,
    canonical_to_target_vocab,
    get_hgnc_mapper,
    load_shared_nuclear_universe,
    normalize_gene_axis,
)
from crlm.program_scores import (  # noqa: E402
    build_detection_matched_controls,
    compute_reference_gene_stats,
    controls_from_dataframe,
    controls_to_dataframe,
    debiased_ucell_scores,
    legacy_positive_rank_scores,
    log1p_cp10k_scores,
    rank_scores,
    ucell_scores,
)


def _sha256(path: Path) -> str:
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def _verify_sha_manifest(manifest_path: Path, base_dir: Path) -> int:
    if not manifest_path.is_file():
        raise FileNotFoundError(f"Missing SHA256 manifest: {manifest_path}")
    lines = [ln.strip() for ln in manifest_path.read_text(encoding="utf-8").splitlines() if ln.strip()]
    if not lines:
        raise AssertionError(f"Empty SHA256 manifest: {manifest_path}")
    for line in lines:
        parts = line.split(None, 1)
        if len(parts) != 2:
            raise ValueError(f"Malformed line in {manifest_path}: {line!r}")
        expected_hex, rel_str = parts[0], parts[1].lstrip("*").strip()
        target = base_dir / rel_str
        if not target.is_file():
            raise FileNotFoundError(f"Missing hashed file {target} referenced by {manifest_path}")
        actual_hex = _sha256(target)
        if actual_hex != expected_hex:
            raise AssertionError(
                f"SHA256 mismatch in {manifest_path} for {rel_str}: expected {expected_hex}, got {actual_hex}"
            )
    return len(lines)


def test_hgnc_harmonization_and_universes() -> dict[str, object]:
    mapper = get_hgnc_mapper()
    prov = json.loads(HGNC_PROVENANCE_PATH.read_text(encoding="utf-8"))
    assert _sha256(HGNC_TSV_PATH) == prov["sha256"], "HGNC TSV SHA mismatch vs PROVENANCE.json"

    u245 = set(load_shared_nuclear_universe("245_tumor_canonical_nuclear"))
    u164 = set(load_shared_nuclear_universe("245_164_shared_nuclear"))
    u203 = set(load_shared_nuclear_universe("245_164_203608_shared_nuclear"))
    excl353 = set(load_shared_nuclear_universe("excluded_353_canonical"))

    assert len(u245) == 29569, f"Expected 29,569 245 nuclear genes, got {len(u245)}"
    assert len(u164) == 18815, f"Expected 18,815 shared 245-164 nuclear genes, got {len(u164)}"
    assert len(u203) == 18733, f"Expected 18,733 shared 3-cohort nuclear genes, got {len(u203)}"
    assert len(excl353) == 353, f"Expected 353 excluded normal-only genes, got {len(excl353)}"
    assert len(u245 & excl353) == 0, "245 nuclear universe overlaps excluded_353"
    assert len(u164 & excl353) == 0, "245_164 shared nuclear universe overlaps excluded_353"
    assert len(u203 & excl353) == 0, "245_164_203608 shared nuclear universe overlaps excluded_353"
    assert u203.issubset(u164) and u164.issubset(u245), "Shared nuclear universes must nest monotonically"

    # Idempotence & collision safety
    raw_164 = (
        ROOT / "artifacts/archive/CRLM_V1_245_F1/artifacts/m6/GSE164522/coverage/external_rank_background.txt"
    ).read_text(encoding="utf-8").splitlines()
    norm_164 = normalize_gene_axis(raw_164, excluded_canonical=excl353)
    assert norm_164["canonical_symbol"].is_unique, "Duplicate canonical symbols on normalized 164 axis"
    meas_164 = norm_164.loc[norm_164["is_unique_measurable"], "canonical_symbol"].tolist()
    norm_164_twice = normalize_gene_axis(meas_164)
    assert (norm_164_twice["canonical_symbol"] == norm_164_twice["source_symbol"]).all(), "HGNC axis not idempotent"

    # Forced collision test: Approved symbol SEPTIN7 alongside legacy prev_symbol SEPT7
    coll_df = normalize_gene_axis(["SEPT7", "SEPTIN7", "DUP_RAW", "DUP_RAW"])
    assert coll_df["canonical_symbol"].is_unique, "Forced collision produced duplicate canonical_symbol"
    assert coll_df.loc[1, "canonical_symbol"] == "SEPTIN7" and bool(coll_df.loc[1, "is_unique_measurable"])
    assert not bool(coll_df.loc[2, "is_unique_measurable"]) and not bool(coll_df.loc[3, "is_unique_measurable"])

    # GSE203608 ENSG disambiguation & M4 512 PCA HVG external coverage >= 98.8%
    g203 = pd.read_csv(ROOT / "data/raw/GSE203608/GSE203608_gene_annotation.csv.gz", dtype=str, keep_default_na=False)
    norm_203 = normalize_gene_axis(g203["gene_short_name"].tolist(), g203["gene_id"].tolist())
    assert norm_203.loc[norm_203["canonical_symbol"].ne(""), "canonical_symbol"].is_unique
    assert int(norm_203["is_unique_measurable"].sum()) >= 33658

    vocab_scgpt = json.loads((ROOT / "models/scGPT-whole-human/vocab.json").read_text(encoding="utf-8"))
    gene_vocab = {k: v for k, v in vocab_scgpt.items() if not k.startswith("<")}
    bridge_164 = canonical_to_target_vocab(raw_164, gene_vocab)
    matched_idx = bridge_164.loc[bridge_164["matched"], "target_index"]
    assert matched_idx.is_unique, "canonical_to_target_vocab violated 1-to-1 injectivity on scGPT vocab"

    pca_cov = {}
    v1_archive_pca_cov = {}
    for lin in ("T_NK", "Mono_Mac_DC", "B_Plasma"):
        hvgs_v11 = [str(g) for g in np.load(ROOT / f"models/m4/full/{lin}/pca.npz", allow_pickle=False)["genes"]]
        c164 = float(align_axis_to_canonical_targets(raw_164, hvgs_v11)["present"].mean())
        c203 = float(
            align_axis_to_canonical_targets(
                g203["gene_short_name"].tolist(), hvgs_v11, g203["gene_id"].tolist()
            )["present"].mean()
        )
        csc = float(canonical_to_target_vocab(hvgs_v11, gene_vocab)["matched"].mean())
        assert c164 == 1.0, f"{lin} V1.1 164 PCA HVG coverage {c164:.4f} != 1.0"
        assert c203 >= 0.994, f"{lin} V1.1 203608 PCA HVG coverage {c203:.4f} < 0.994"
        assert csc == 1.0, f"{lin} V1.1 scGPT HVG coverage {csc:.4f} < 1.0"
        pca_cov[lin] = {"164": round(c164, 4), "203608": round(c203, 4), "scgpt": round(csc, 4)}

        arch_pca = ROOT / f"models/archive/CRLM_V1_245_F1/m4/full/{lin}/pca.npz"
        if arch_pca.is_file():
            hvgs_arch = [str(g) for g in np.load(arch_pca, allow_pickle=False)["genes"]]
            a164 = float(align_axis_to_canonical_targets(raw_164, hvgs_arch)["present"].mean())
            a203 = float(
                align_axis_to_canonical_targets(
                    g203["gene_short_name"].tolist(), hvgs_arch, g203["gene_id"].tolist()
                )["present"].mean()
            )
            assert a164 >= 0.988 and a203 >= 0.988, f"{lin} archived V1 PCA HVG coverage < 0.988"
            v1_archive_pca_cov[lin] = {"164": round(a164, 4), "203608": round(a203, 4)}

    return {
        "v1_1_pca_512_coverage": pca_cov,
        "v1_archive_pca_512_coverage": v1_archive_pca_cov,
        "shared_nuclear_245_164": len(u164),
        "shared_nuclear_all": len(u203),
    }


def test_hierarchical_cell_annotation_and_m2() -> dict[str, object]:
    panel_genes = sorted(
        {g for _, a, s in PANELS.values() for g in (a + " " + s).split()} | set(AUXILIARY_MARKERS)
    )
    g2i = {g: i for i, g in enumerate(panel_genes)}

    def make_row(expr: dict[str, int]) -> np.ndarray:
        row = np.zeros(len(panel_genes), dtype=np.int32)
        for g, val in expr.items():
            row[g2i[g]] = val
        return row

    synthetic_cells = np.vstack(
        [
            # 0: CD14+ CSF3R+ FCGR3B=0 macrophage -> Mono/Mac/DC (monocyte/macrophage, Macrophage)
            make_row({"CD14": 6, "C1QA": 8, "CSF1R": 3, "CSF3R": 4, "S100A8": 5, "FCGR3B": 0, "CEACAM8": 0}),
            # 1: CD14+ CSF3R+ FCN1+ S100A8/9+ FCGR3B=0 monocyte -> Mono/Mac/DC (Monocyte)
            make_row({"CD14": 7, "FCN1": 5, "CSF3R": 4, "S100A8": 10, "S100A9": 9, "FCGR3B": 0}),
            # 2: True neutrophil FCGR3B+ CEACAM8+ CXCR2+ CSF3R+ -> Neutrophil
            make_row({"FCGR3B": 8, "CEACAM8": 3, "CXCR2": 5, "G0S2": 6}),
            # 3: Cytotoxic CD8 T cell co-expressing GNLY/KLRD1/NKG7 -> T/NK (CD8 T cell, NOT unresolved!)
            make_row({"CD3D": 5, "CD3E": 4, "CD8A": 6, "GNLY": 9, "KLRD1": 3, "NKG7": 8}),
            # 4: Hepatic NK cell with sterile TRBC1+ but CD3D/E/G/TRAC=0 -> T/NK (NK cell)
            make_row({"NCR1": 4, "KLRD1": 5, "GNLY": 10, "NKG7": 7, "TRBC1": 2}),
            # 5: Plasma cell co-expressing CD79A/CD79B with MZB1+JCHAIN+SDC1+ and MS4A1=0 -> B/Plasma (plasma cell)
            make_row({"MZB1": 12, "JCHAIN": 15, "SDC1": 4, "CD79A": 5, "CD79B": 2, "MS4A1": 0}),
            # 6: Follicular B cell MS4A1+ CD19+ CD79A+ -> B/Plasma (B cell)
            make_row({"MS4A1": 8, "CD19": 4, "CD79A": 6, "BANK1": 3}),
            # 7: Cross-lineage T-B doublet -> Mixed (Cross_lineage_doublet)
            make_row({"CD3D": 6, "CD3E": 5, "MS4A1": 7, "CD79A": 5}),
            # 8: Cross-lineage Immune + EPCAM+ doublet -> Mixed (Cross_lineage_doublet)
            make_row({"CD3D": 5, "CD8A": 4, "PTPRC": 4, "EPCAM": 9, "KRT19": 6}),
        ]
    )
    labels, details = classify(sparse.csr_matrix(synthetic_cells), panel_genes)
    subs = details["subtype_hierarchical"].tolist()
    mf = details["myeloid_fine"].tolist()

    assert labels[0] == "Mono/Mac/DC" and subs[0] == "monocyte/macrophage" and mf[0] == "Macrophage"
    assert labels[1] == "Mono/Mac/DC" and subs[1] == "monocyte/macrophage" and mf[1] == "Monocyte"
    assert labels[2] == "Neutrophil" and subs[2] == "neutrophil"
    assert labels[3] == "T/NK" and subs[3] == "CD8 T cell"
    assert labels[4] == "T/NK" and subs[4] == "NK cell"
    assert labels[5] == "B/Plasma" and subs[5] == "plasma cell"
    assert labels[6] == "B/Plasma" and subs[6] == "B cell"
    assert labels[7] == "Mixed" and subs[7] == "Cross_lineage_doublet"
    assert labels[8] == "Mixed" and subs[8] == "Cross_lineage_doublet"
    assert not details["subtype_hierarchical"].str.contains("unresolved", case=False).any()

    # Check materialized M2 cohort invariants across all 39 libraries
    qc_m2 = pd.read_csv(ROOT / "artifacts/m2/qc_and_exclusions.csv")
    pb_m2 = pd.read_csv(ROOT / "artifacts/m2/pseudobulk_manifest.csv")
    assert len(qc_m2) == 39, f"Expected 39 M2 libraries, got {len(qc_m2)}"
    assert int(qc_m2["input_cells"].sum()) == 190190
    assert int(qc_m2["retained_cells"].sum()) == 173450
    assert int(qc_m2["excluded_total"].sum()) == 16740
    assert (qc_m2["high_genes"] <= qc_m2["native_genes"]).all(), "high_genes exceeded native_genes in M2"
    gsm850 = qc_m2.loc[qc_m2["sample_uid"].eq("GSM7844850")].iloc[0]
    assert gsm850["pct_mt_scope"] == "visible_axis_depleted_source_anomaly" and not bool(gsm850["pct_mt_usable_for_qc"])

    assert not pb_m2["group"].astype(str).str.contains("unresolved", case=False).any(), "Found unresolved group in M2"
    broad_pb = pb_m2.loc[pb_m2["level"].eq("lineage")]
    for gsm, grp in broad_pb.groupby("sample_uid"):
        assert np.isclose(grp["prop_total_resolved_immune"].dropna().sum(), 1.0, atol=1e-12), f"{gsm} total_imm != 1"
        assert np.isclose(grp["prop_core_three_lineage_immune"].dropna().sum(), 1.0, atol=1e-12), f"{gsm} core3 != 1"

    # Spot-check an actual materialized M2 H5AD (`GSM7844812.h5ad`) for 0 unresolved obs labels
    ad_spot = ad.read_h5ad(ROOT / "data/processed/m2/GSM7844812.h5ad", backed="r")
    try:
        assert not ad_spot.obs["subtype_m2"].astype(str).str.contains("unresolved", case=False).any()
    finally:
        ad_spot.file.close()

    return {
        "m2_retained_cells": int(qc_m2["retained_cells"].sum()),
        "m2_unresolved_cells": 0,
    }


def _bruteforce_ucell(dense: np.ndarray, sig_idx: list[int], max_rank: int) -> np.ndarray:
    m = len(sig_idx)
    out = np.empty(dense.shape[0], dtype=np.float64)
    for i in range(dense.shape[0]):
        r = rankdata(-dense[i], method="average")
        r_cap = np.where(r <= max_rank, r, max_rank + 1.0)
        u_stat = r_cap[sig_idx].sum() - 0.5 * m * (m + 1)
        out[i] = max(0.0, 1.0 - u_stat / (m * max_rank))
    return out


def test_program_scores_math_and_debiasing() -> dict[str, object]:
    rng = np.random.default_rng(20261121)

    # 1. Exact Mann-Whitney UCell vs brute-force scipy.stats.rankdata on G=309 (< max_rank=1500)
    g_small = 309
    syms_small = [f"G{i:04d}" for i in range(g_small)]
    dense_small = rng.poisson(0.4, size=(40, g_small)).astype(np.float64)
    dense_small[0, :] = 0.0  # All-zero cell (K_c = 0)
    dense_small[1, :] = rng.integers(1, 15, size=g_small)  # Fully positive cell (K_c = G)
    dense_small[2, :10] = np.arange(100, 90, -1)  # Top-10 saturation -> score == 1.0
    progs_small = {"TOP10": syms_small[:10], "MIX25": syms_small[5:30]}
    for r_max in (1500, 200, 80):
        got = ucell_scores(sparse.csr_matrix(dense_small), syms_small, progs_small, max_rank=r_max)
        for pname, pmembers in progs_small.items():
            idx = [syms_small.index(g) for g in pmembers]
            exp = _bruteforce_ucell(dense_small, idx, r_max)
            err = float(np.max(np.abs(got[pname].to_numpy() - exp)))
            assert err < 1e-12, f"UCell mismatch at max_rank={r_max} {pname}: {err}"
        assert np.isclose(got.loc[2, "TOP10"], 1.0, atol=1e-15), "Top-m rank saturation must equal 1.0"

    # 2. Bit-exact legacy_positive compatibility vs archived V1_F1 program_scores.py
    spec = importlib.util.spec_from_file_location(
        "archived_v1_program_scores",
        ROOT / "artifacts/archive/CRLM_V1_245_F1/artifacts/m5/program_scores.py",
    )
    archived_mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(archived_mod)
    old_df = archived_mod.rank_scores(sparse.csr_matrix(dense_small), syms_small, progs_small, max_rank=1500)
    new_df = rank_scores(sparse.csr_matrix(dense_small), syms_small, progs_small, max_rank=1500)
    leg_df = legacy_positive_rank_scores(sparse.csr_matrix(dense_small), syms_small, progs_small, max_rank=1500)
    assert float(np.max(np.abs(old_df.to_numpy() - new_df.to_numpy()))) == 0.0
    assert float(np.max(np.abs(old_df.to_numpy() - leg_df.to_numpy()))) == 0.0

    # 3. Anti-leakage guard: unparameterized external call MUST raise ValueError
    for fn in (
        lambda: debiased_ucell_scores(sparse.csr_matrix(dense_small), syms_small, progs_small),
        lambda: rank_scores(sparse.csr_matrix(dense_small), syms_small, progs_small, mode="ucell_debiased"),
        lambda: log1p_cp10k_scores(sparse.csr_matrix(dense_small), syms_small, progs_small, debiased=True),
    ):
        raised = False
        try:
            fn()
        except ValueError:
            raised = True
        assert raised, "Expected ValueError when calling debiased scorer without frozen controls/reference stats"

    # 4. Synthetic Poisson sparse depth-independence test (K_c spanning [180, 2300], zero bio signal)
    n_genes = 3000
    n_cells = 3000
    n_ref = 1500
    syms = [f"GENE_{i:04d}" for i in range(n_genes)]
    base_rates = rng.lognormal(mean=-2.2, sigma=1.1, size=n_genes)
    base_rates /= base_rates.sum()
    depths = np.exp(rng.uniform(np.log(260), np.log(7500), size=n_cells))
    sim_counts = sparse.csr_matrix(rng.poisson(depths[:, None] * base_rates[None, :]).astype(np.int32))
    kc = np.diff(sim_counts.indptr)
    assert (kc < 1500).mean() > 0.50, "Synthetic test must exercise K_c < 1500 regime"

    ref_stats = compute_reference_gene_stats(sim_counts[:n_ref], syms, input_scale="counts", n_bins=20)
    rand_progs = {
        f"RAND_{p:02d}": rng.choice(syms, size=80, replace=False).tolist()
        for p in range(10)
    }
    ctrls = build_detection_matched_controls(rand_progs, ref_stats, n_bins=20, n_ctrl=10, seed=20261121)
    # Round-trip serialization check
    rt_ctrls = controls_from_dataframe(controls_to_dataframe(ctrls, ref_stats))
    assert rt_ctrls == ctrls, "controls_to_dataframe <-> controls_from_dataframe roundtrip mismatch"

    deb_eval, raw_eval, _ = debiased_ucell_scores(
        sim_counts[n_ref:], syms, rand_progs, matched_controls=ctrls, max_rank=1500, return_components=True
    )
    kc_eval = kc[n_ref:]
    raw_rhos = [abs(float(spearmanr(kc_eval, raw_eval[c]).statistic)) for c in rand_progs]
    deb_rhos = [abs(float(spearmanr(kc_eval, deb_eval[c]).statistic)) for c in rand_progs]
    max_deb_rho = max(deb_rhos)
    assert max(raw_rhos) > 0.45, "Raw UCell should exhibit expected depth correlation on K_c < 1500"
    assert max_deb_rho < 0.08, f"Debiased UCell max |rho| with K_c = {max_deb_rho:.4f} >= 0.08"

    # 5. Sub-axis log1p(CP10k) renormalization invariance to extra non-nuclear (MT/Ribo) counts
    extra_ribo = rng.poisson(depths[:100, None] * 0.25, size=(100, 20)).astype(np.float64)
    full_counts = np.hstack([sim_counts[:100].toarray().astype(np.float64), extra_ribo])
    row_tot = full_counts.sum(axis=1, keepdims=True)
    published_log1p = np.log1p(10000.0 * full_counts / row_tot)
    nuclear_slice_log1p = sparse.csr_matrix(published_log1p[:, :n_genes])
    s_from_counts = log1p_cp10k_scores(
        sim_counts[:100], syms, {"RAND_00": rand_progs["RAND_00"]}, matched_controls=ctrls, input_scale="counts"
    )
    s_from_log1p = log1p_cp10k_scores(
        nuclear_slice_log1p,
        syms,
        {"RAND_00": rand_progs["RAND_00"]},
        matched_controls=ctrls,
        input_scale="log1p",
        renormalize_subaxis=True,
    )
    assert np.allclose(s_from_counts.to_numpy(), s_from_log1p.to_numpy(), atol=1e-12)

    return {"max_raw_rho": round(max(raw_rhos), 4), "max_debiased_rho": round(max_deb_rho, 4)}


def _run_prior_mil_torch_assertions() -> dict[str, object]:
    import math
    import torch
    from crlm.prior_mil import (
        BagBatch,
        DEFAULT_HIERMIL_DIR,
        DEFAULT_MIXMIL_DIR,
        HierarchicalMIL,
        MixedMIL,
        mean_glmm_init,
        official_hiermil_model,
        official_posterior,
        resolve_prior_repo,
    )

    # 1. In-repo prior repositories (zero ~/.cache dependency)
    assert resolve_prior_repo("mixmil") == DEFAULT_MIXMIL_DIR and DEFAULT_MIXMIL_DIR.is_dir()
    assert resolve_prior_repo("hier-mil") == DEFAULT_HIERMIL_DIR and (DEFAULT_HIERMIL_DIR / "model.py").is_file()
    post_cls = official_posterior()

    # 2. Verify removal of [::2] odd-bag drop & 1/sqrt(Q) scaling in MixedMIL
    rng = np.random.default_rng(20261121)
    q_dim = 16
    bags_np = []
    labels_np = np.array([0, 1] * 8, dtype=float)
    for idx in range(16):
        scale = 0.2 if idx % 2 == 0 else 4.5  # Make odd-indexed bags have much larger variance than even bags
        shift = 0.8 if labels_np[idx] > 0 else -0.8
        bags_np.append((rng.normal(shift, scale, size=(20, q_dim))).astype(np.float32))

    m_full, _, _ = mean_glmm_init(bags_np, labels_np, legacy_even_slice=False)
    m_even, _, _ = mean_glmm_init(bags_np, labels_np, legacy_even_slice=True)
    assert not torch.allclose(m_full, m_even), "mean_glmm_init full vs legacy [::2] slice should differ on asymmetric bags"

    torch.manual_seed(20261121)
    mix_model = MixedMIL(bags_np, labels_np, post_cls)
    assert mix_model.scale_by_sqrt_q is True
    assert math.isclose(float(mix_model.feature_scale.item()), 1.0 / math.sqrt(q_dim), rel_tol=1e-6)
    batch = BagBatch(bags_np, [np.array(["subA"] * 10 + ["subB"] * 10) for _ in range(16)])
    mix_model.bind_training_scale(batch)
    with torch.inference_mode():
        logits, att, signed, bias = mix_model(batch)
        owner = batch.owner.numpy()
        assert np.allclose(np.bincount(owner, weights=att.numpy()), np.ones(16), atol=1e-5)
        assert np.allclose(np.bincount(owner, weights=signed.numpy()) + float(bias.item()), logits.numpy(), atol=1e-5)

    hm = HierarchicalMIL(q_dim, cell_attention=True, dropout=0.25)
    hm.eval()
    with torch.inference_mode():
        h_logits, h_att, h_signed, h_bias = hm(batch)
        assert np.allclose(np.bincount(owner, weights=h_att.numpy()), np.ones(16), atol=1e-5)
        assert np.allclose(np.bincount(owner, weights=h_signed.numpy()) + float(h_bias.item()), h_logits.numpy(), atol=1e-5)

    return {"scale_by_sqrt_q": True, "feature_scale_q16": round(float(mix_model.feature_scale.item()), 6)}


def test_prior_mil_and_stage_contracts() -> dict[str, object]:
    if importlib.util.find_spec("torch") is not None:
        torch_res = _run_prior_mil_torch_assertions()
    else:
        scgpt_py = ROOT / ".venv-scgpt/bin/python"
        if not scgpt_py.is_file():
            raise RuntimeError("PyTorch not installed in active interpreter and .venv-scgpt/bin/python is missing")
        env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1")
        out = subprocess.check_output(
            [
                str(scgpt_py),
                "-B",
                "-c",
                (
                    "import json, sys; "
                    f"sys.path[:0] = [{str(ROOT / 'scripts')!r}, {str(ROOT / 'src')!r}, {str(ROOT / 'starter')!r}]; "
                    "import test_audit_remediations as t; "
                    "print(json.dumps(t._run_prior_mil_torch_assertions()))"
                ),
            ],
            cwd=ROOT,
            env=env,
            text=True,
        )
        torch_res = json.loads(out.strip().splitlines()[-1])

    # 3. Verify M6 external script contracts (missing matched_controls.csv guard)
    import validate_gse164522 as v164
    import validate_gse203608 as v203

    with tempfile.TemporaryDirectory() as tmpdir:
        reg_dummy = pd.DataFrame([{"lineage": "Mono/Mac/DC", "program": "P1", "scoring_method": "ucell_debiased"}])
        for mod in (v164, v203):
            raised = False
            try:
                mod.load_frozen_matched_controls(Path(tmpdir), [("Mono/Mac/DC", "P1")], reg_dummy)
            except FileNotFoundError:
                raised = True
            assert raised, f"{mod.__name__} failed to raise FileNotFoundError when matched_controls.csv is absent"

    # 4. Verify M7 spatial & experimental design remediations
    cov_m7 = pd.read_csv(ROOT / "artifacts/m7/spatial_measurement_audit/measured_coverage.csv")
    decontam = cov_m7.loc[cov_m7["program"].eq("B19_S4_Mac_SPP1_decontaminated"), "measured_symbols"]
    assert len(decontam) == 4
    for s in decontam:
        toks = set(str(s).split())
        assert toks.isdisjoint({"GZMB", "GNLY", "LAT", "IL2RG"}), f"T/NK contaminant leaked into M7: {toks}"

    supp_m7 = pd.read_csv(ROOT / "artifacts/m7/spatial_measurement_audit/sample_support.csv")
    m397 = supp_m7.loc[supp_m7["sample"].str.contains("M397T")].iloc[0]
    assert int(m397["n_cd8_contact_rescued"]) > 3000, "Expected >3,000 contact-rescued CD8 T cells in M397T"
    dist_m7 = pd.read_csv(ROOT / "artifacts/m7/spatial_measurement_audit/local_distances.csv")
    assert set(dist_m7["tile_um"].unique()) == {50, 100, 200, 400}
    assert (dist_m7["n_permutations"] == 1000).all()

    design_md = (ROOT / "artifacts/m7/candidate_validation_design.md").read_text(encoding="utf-8")
    assert "Sun" in design_md and "Immunity" in design_md, "candidate_validation_design.md missing Sun Immunity 2026"
    assert "CXCL9" in design_md, "candidate_validation_design.md missing CXCL9/10 neutralization control"

    # 5. Verify Track 1 techfix external models (15/15 predicted, 0 rejected) & Step 0 program parity
    tf_models = pd.read_csv(ROOT / "artifacts/techfix_v1_f1/models/techfix_external_classifier_metrics.csv")
    step_b = tf_models.loc[(tf_models["step"] == "step_b_full_techfix") & (tf_models["scope"] == "method_available")]
    assert len(step_b) == 15 and step_b["auroc"].notna().all(), "Expected all 15 frozen models predicted in Track 1"
    pca_tnk_auc = float(step_b.loc[(step_b["lineage"] == "T/NK") & (step_b["method"] == "PCA_mean"), "auroc"].iloc[0])
    assert pca_tnk_auc >= 0.95, f"Expected T/NK PCA_mean external AUROC >= 0.95 after HGNC fix, got {pca_tnk_auc}"

    return {
        "techfix_models_predicted": int(len(step_b)),
        "m397_rescued_cd8": int(m397["n_cd8_contact_rescued"]),
        **torch_res,
    }


def test_m4_m5_m6_v1_1_contracts() -> dict[str, object]:
    """Part 2A: End-to-end verification of V1.1 M4, M5, and M6 empirical artifacts & contracts."""
    # 1. M4: 67 hashed outputs, zero-starvation deconfounding, 200-perm null, XPU profile, 15/15 transfer heads
    m4_metrics = pd.read_csv(ROOT / "artifacts/m4/patient_metrics.csv")
    m4_s = m4_metrics.loc[m4_metrics["seed"].eq(20261121)].set_index(["lineage", "method"])

    tnk_sb = m4_s.loc[("T/NK", "PCA_subtype_balanced")]
    assert int(tnk_sb["patients"]) == 16 and int(tnk_sb["n_patient_sites"]) == 28 and int(tnk_sb["paired_n"]) == 12
    assert float(tnk_sb["auroc"]) >= 0.98, f"Expected T/NK PCA_subtype_balanced >= 0.98, got {tnk_sb['auroc']}"
    assert float(m4_s.loc[("T/NK", "PCA_pure_primary"), "auroc"]) >= 0.95
    assert float(m4_s.loc[("T/NK", "PCA_pure_secondary"), "auroc"]) >= 0.98
    assert float(m4_s.loc[("T/NK", "PCA_remove_organ_markers"), "auroc"]) >= 0.98
    assert float(m4_s.loc[("Mono/Mac/DC", "program_elastic"), "auroc"]) >= 0.90
    assert float(m4_s.loc[("B/Plasma", "program_elastic"), "auroc"]) > 0.50, "B/Plasma sign inversion recurred"
    assert np.isclose(float(m4_s.loc[("All3Lineages", "PCA32_MixMIL"), "auroc"]), 1.0, atol=1e-6)
    assert float(m4_s.loc[("All3Lineages", "PCA_mean"), "auroc"]) >= 0.98

    perm_null = pd.read_csv(ROOT / "artifacts/m4/PCA/permutation_null_summary.csv")
    assert len(perm_null) == 3 and (perm_null["n_permutations"] == 200).all()
    assert ((perm_null["null_auroc_mean"] >= 0.45) & (perm_null["null_auroc_mean"] <= 0.56)).all()
    assert (perm_null["empirical_p_one_sided"] <= 0.02).all()

    scgpt_prof = json.loads((ROOT / "artifacts/m4/resource_profiles/scgpt_profile.json").read_text(encoding="utf-8"))
    xpu_info = scgpt_prof["devices"]["xpu"]
    assert bool(xpu_info["accepted"]) and float(xpu_info["max_abs_diff_vs_cpu"]) < 1e-5
    assert float(xpu_info["min_cosine_vs_cpu"]) > 0.999999

    full_models = pd.read_csv(ROOT / "artifacts/m4/full_data_models.csv")
    assert len(full_models) == 15 and full_models["accepted"].all()
    for r in full_models.itertuples(index=False):
        assert _sha256(ROOT / r.path) == r.sha256, f"M4 full model checksum mismatch: {r.path}"

    # 2. M5: 72 frozen inputs, 19 M5 outputs, config parity, matched_controls lock, exposure disclosure
    assert _sha256(ROOT / "configs/project.json") == _sha256(ROOT / "artifacts/m5/development_config.json")
    proj_cfg = json.loads((ROOT / "configs/project.json").read_text(encoding="utf-8"))
    assert proj_cfg["research_version"] == "V1.1_AUDIT_REMEDIATED"
    assert proj_cfg["candidate_registration"]["freeze_version"] == "CRLM_V1.1_245_F1"

    reg_m5 = pd.read_csv(ROOT / "artifacts/m5/frozen_candidate_registry.csv")
    plan_m5 = pd.read_csv(ROOT / "artifacts/m5/registered_tests.csv")
    assert len(reg_m5) == 2 and (reg_m5["status"] == "frozen").all()
    assert set(zip(reg_m5["lineage"], reg_m5["program"])) == {
        ("Mono/Mac/DC", "Core_Lipid_SAM_R"),
        ("T/NK", "B07_Th17_TWEAK"),
    }
    assert len(plan_m5) == 8 and plan_m5["family"].value_counts().eq(2).all()
    mc_sha = _sha256(ROOT / "artifacts/m3/program_definitions_draft/matched_controls.csv")
    assert _sha256(ROOT / "artifacts/m5/matched_controls.csv") == mc_sha
    assert _sha256(ROOT / "artifacts/m5/gene_sets/matched_controls.csv") == mc_sha

    # 3. M6: Clean cohorts + rule_annotation_m1 branch + GSE178318 terms_pending + 15/15 zero-rejection transfer
    reg_res = pd.read_csv(ROOT / "artifacts/m6/registered_tests/registered_test_results.csv")
    assert len(reg_res) == 8
    th17_i = reg_res.loc[
        reg_res["cohort"].eq("GSE164522")
        & reg_res["program"].eq("B07_Th17_TWEAK")
        & reg_res["contrast"].eq("I")
    ].iloc[0]
    assert int(th17_i["n_patients"]) == 10 and int(th17_i["n_positive"]) == 1 and int(th17_i["n_negative"]) == 9
    assert float(th17_i["holm_adjusted_p"]) <= 0.05 and bool(th17_i["holm_reject"])

    g178 = reg_res.loc[reg_res["cohort"].eq("GSE178318")]
    assert len(g178) == 2 and g178["planned_test_missing"].all()
    assert g178["raw_p"].isna().all() and (g178["adjustment_only_p"] == 1.0).all() and (g178["holm_adjusted_p"] == 1.0).all()

    unif_m6 = pd.read_csv(ROOT / "artifacts/m6/uniform_annotation_effects_summary.csv")
    assert set(unif_m6["dataset"].unique()) == {"GSE164522", "GSE203608"}

    man_164 = pd.read_csv(ROOT / "artifacts/m6/GSE164522/coverage/sample_manifest.csv")
    assert int(man_164["clean_cells"].sum()) == 156511 and int(man_164["excluded_cd3e_doublets"].sum()) == 6836
    qc_203 = pd.read_csv(ROOT / "artifacts/m6/GSE203608/coverage/cell_qc_audit.csv")
    assert int(qc_203["clean_author_cells"].sum()) == 47525 and int(qc_203["excluded_total_author"].sum()) == 12136

    pca_ext = pd.read_csv(ROOT / "artifacts/m6/GSE164522/model_predictions/pca_measurement.csv")
    assert len(pca_ext) == 57 and (pca_ext["n_mapped_hvg"] == 512).all() and (pca_ext["mapped_fraction"] == 1.0).all()
    assert pca_ext["numeric_scale_pass"].all()

    pred_cov = pd.read_csv(ROOT / "artifacts/m6/GSE164522/model_predictions/prediction_coverage.csv")
    assert set(pred_cov["status"].unique()).issubset({"predicted", "fewer_than_50_cells"})
    assert int(pred_cov.loc[pred_cov["status"].eq("predicted"), "n_bags"].sum()) == 57 * 5

    ext_metrics = pd.read_csv(ROOT / "artifacts/m6/GSE164522/model_predictions/patient_metrics.csv")
    ext_common = ext_metrics.loc[ext_metrics["scope"].eq("common_predicting_methods")]
    assert len(ext_common) == 15 and ext_common["auroc"].notna().all()
    b_prog_auc = float(
        ext_common.loc[ext_common["lineage"].eq("B/Plasma") & ext_common["method"].eq("program_elastic"), "auroc"].iloc[0]
    )
    assert b_prog_auc >= 0.85, f"Expected B/Plasma external program_elastic >= 0.85, got {b_prog_auc}"

    return {
        "m4_tnk_subtype_balanced_auroc": round(float(tnk_sb["auroc"]), 4),
        "m6_th17_I_holm_p": round(float(th17_i["holm_adjusted_p"]), 5),
        "m6_external_models_active": 15,
        "m6_external_clean_bags": 57,
    }


def test_m8_figures_source_data_and_manifest() -> dict[str, object]:
    """Part 2B: Verify M8 PNG/PDF figures, exact source_data parity vs upstream CSVs, and 24-row manifest."""
    fig_dir = ROOT / "artifacts/m8/figures"
    src_dir = ROOT / "artifacts/m8/source_data"
    stems = (
        "01_patient_programs",
        "02_development_classifiers",
        "03_spatial_distances",
        "04_external_classifiers",
    )
    for stem in stems:
        png_p = fig_dir / f"{stem}.png"
        pdf_p = fig_dir / f"{stem}.pdf"
        assert png_p.is_file() and png_p.stat().st_size > 15_000, f"Invalid/missing PNG {png_p}"
        assert pdf_p.is_file() and pdf_p.stat().st_size > 8_000, f"Invalid/missing PDF {pdf_p}"
        assert png_p.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n", f"Corrupt PNG header in {png_p}"
        assert pdf_p.read_bytes()[:5] == b"%PDF-", f"Corrupt PDF header in {pdf_p}"

    # Element-by-element equality of all 7 M8 source_data CSVs against upstream M3/M4/M5/M6/M7 CSVs
    reg_m5 = pd.read_csv(ROOT / "artifacts/m5/frozen_candidate_registry.csv")
    target_keys = set(zip(reg_m5["lineage"].astype(str), reg_m5["program"].astype(str))) | {
        ("B/Plasma", "HALLMARK_UNFOLDED_PROTEIN_RESPONSE")
    }
    up_pat = pd.concat(
        [
            pd.read_csv(ROOT / "artifacts/m3/patient_effects/patient_contrasts.csv"),
            pd.read_csv(ROOT / "artifacts/m6/patient_contrasts.csv"),
        ],
        ignore_index=True,
    )
    exp_cand_pat = up_pat.loc[
        [(str(l), str(p)) in target_keys for l, p in zip(up_pat["lineage"], up_pat["program"])]
    ].reset_index(drop=True)
    pd.testing.assert_frame_equal(
        pd.read_csv(src_dir / "candidate_patient_effects.csv"), exp_cand_pat, check_exact=False, rtol=1e-14, atol=1e-15
    )

    m3_dec = pd.read_csv(ROOT / "artifacts/m3/composition_results/patient_subtype_decomposition.csv")
    d245 = m3_dec.loc[
        m3_dec["program"].eq("HALLMARK_UNFOLDED_PROTEIN_RESPONSE")
        & m3_dec["min_subtype_cells_per_tissue"].eq(1)
        & m3_dec["status"].eq("complete_support"),
        ["donor_uid", "observed_total_D_T", "full_composition", "full_within_subtype_state"],
    ].rename(
        columns={
            "observed_total_D_T": "total_D_T",
            "full_composition": "composition",
            "full_within_subtype_state": "within_subtype_state",
        }
    )
    d245.insert(0, "dataset", "GSE245552")

    m6_dec = pd.read_csv(ROOT / "artifacts/m6/GSE164522/patient_effects/patient_B_plasma_decomposition.csv")
    d164 = m6_dec.loc[
        m6_dec["program"].eq("HALLMARK_UNFOLDED_PROTEIN_RESPONSE")
        & m6_dec["min_subtype_cells"].eq(1)
        & m6_dec["status"].eq("complete_support"),
        ["donor_uid", "observed_total_difference", "full_composition", "full_within_type_state"],
    ].rename(
        columns={
            "observed_total_difference": "total_D_T",
            "full_composition": "composition",
            "full_within_type_state": "within_subtype_state",
        }
    )
    d164.insert(0, "dataset", "GSE164522")
    exp_upr = pd.concat([d245, d164], ignore_index=True)
    pd.testing.assert_frame_equal(
        pd.read_csv(src_dir / "upr_decomposition_reference.csv"), exp_upr, check_exact=False, rtol=1e-14, atol=1e-15
    )

    exp_dev_cls = (
        pd.read_csv(ROOT / "artifacts/m4/patient_metrics.csv")
        .query("seed == 20261121")
        .reset_index(drop=True)
    )
    pd.testing.assert_frame_equal(
        pd.read_csv(src_dir / "development_classifier_metrics.csv"), exp_dev_cls, check_exact=False, rtol=1e-14, atol=1e-15
    )
    pd.testing.assert_frame_equal(
        pd.read_csv(src_dir / "development_permutation_null.csv"),
        pd.read_csv(ROOT / "artifacts/m4/PCA/permutation_null_summary.csv"),
        check_exact=False,
        rtol=1e-14,
        atol=1e-15,
    )
    pd.testing.assert_frame_equal(
        pd.read_csv(src_dir / "spatial_local_distances.csv"),
        pd.read_csv(ROOT / "artifacts/m7/spatial_measurement_audit/local_distances.csv"),
        check_exact=False,
        rtol=1e-14,
        atol=1e-15,
    )
    pd.testing.assert_frame_equal(
        pd.read_csv(src_dir / "spatial_sample_support.csv"),
        pd.read_csv(ROOT / "artifacts/m7/spatial_measurement_audit/sample_support.csv"),
        check_exact=False,
        rtol=1e-14,
        atol=1e-15,
    )
    pd.testing.assert_frame_equal(
        pd.read_csv(src_dir / "external_classifier_metrics.csv"),
        pd.read_csv(ROOT / "artifacts/m6/GSE164522/model_predictions/patient_metrics.csv"),
        check_exact=False,
        rtol=1e-14,
        atol=1e-15,
    )

    # Verify 24-row reproducibility_manifest.csv: all 120 file/venv paths exist & executables have +x
    repro = pd.read_csv(ROOT / "artifacts/m8/reproducibility_manifest.csv")
    assert len(repro) == 24, f"Expected 24 rows in reproducibility_manifest.csv, got {len(repro)}"
    path_cols = ("entrypoint", "principal_table", "verification_or_run", "environment", "package_lock")
    for idx_r, row in repro.iterrows():
        for col in path_cols:
            target = ROOT / str(row[col]).strip()
            assert target.is_file(), f"reproducibility_manifest row {idx_r} missing {col}: {target}"
        env_bin = ROOT / str(row["environment"]).strip()
        assert os.access(env_bin, os.X_OK), f"Non-executable environment binary: {env_bin}"

    # Verify relative Markdown links across root & M8 documents
    link_re = re.compile(r"\[[^\]]+\]\(([^)#]+)(?:#[^)]*)?\)")
    md_checked = 0
    for md_path in (
        ROOT / "README.md",
        ROOT / "STATUS.md",
        ROOT / "artifacts/m8/final_report.md",
        ROOT / "artifacts/m8/limitations_and_deviations.md",
        ROOT / "artifacts/m8/reproduce.md",
        ROOT / "artifacts/m8/stage_report.md",
    ):
        txt = md_path.read_text(encoding="utf-8")
        for m in link_re.finditer(txt):
            href = m.group(1).strip()
            if not href or "://" in href or href.startswith("mailto:"):
                continue
            resolved = (md_path.parent / href).resolve()
            assert resolved.exists(), f"Broken relative link in {md_path.relative_to(ROOT)}: {href} -> {resolved}"
            md_checked += 1

    # Isolated verification of starter 106 synthetic checks without mutating starter/
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_starter = Path(tmpdir) / "starter"
        shutil.copytree(ROOT / "starter", tmp_starter)
        env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1")
        starter_counts = []
        for script_name in ("selfcheck.py", "mil_selfcheck.py", "registered_tests_selfcheck.py"):
            out = subprocess.check_output(
                [sys.executable, "-B", str(tmp_starter / script_name)],
                cwd=tmp_starter,
                env=env,
                text=True,
            )
            payload = json.loads(out)
            assert payload["status"] == "passed"
            starter_counts.append(int(payload["n_checks"]))
        assert sum(starter_counts) == 106, f"Expected 106 starter checks, got {starter_counts}"

    return {
        "m8_figures_png_pdf": 8,
        "m8_source_tables_exact": 7,
        "reproducibility_manifest_paths": len(repro) * len(path_cols),
        "markdown_links_verified": md_checked,
        "isolated_starter_checks": 106,
    }


def test_sha256_manifests() -> dict[str, int]:
    root_manifest = ROOT / "SHA256SUMS.txt"
    if root_manifest.is_file():
        for ln in root_manifest.read_text(encoding="utf-8").splitlines():
            if ln.strip():
                assert ln.split(None, 1)[1].strip() != "SHA256SUMS.txt", "SHA256SUMS.txt must not hash itself"
    manifests = [
        ("root_sha256sums", root_manifest, ROOT),
        ("starter", ROOT / "starter/SHA256SUMS.txt", ROOT / "starter"),
        ("archive_v1_f1", ROOT / "artifacts/archive/CRLM_V1_245_F1/archive_sha256.txt", ROOT / "artifacts/archive/CRLM_V1_245_F1"),
        ("archive_models_v1_f1", ROOT / "models/archive/CRLM_V1_245_F1/ARCHIVE_SHA256.txt", ROOT),
        ("archive_embeddings_v1_f1", ROOT / "embeddings/archive/CRLM_V1_245_F1/ARCHIVE_SHA256.txt", ROOT),
        ("m0", ROOT / "artifacts/m0/output_sha256.txt", ROOT),
        ("m1", ROOT / "artifacts/m1/output_sha256.txt", ROOT),
        ("m2", ROOT / "artifacts/m2/output_sha256.txt", ROOT),
        ("m3", ROOT / "artifacts/m3/output_sha256.txt", ROOT),
        ("m4_direct_prior", ROOT / "artifacts/m4/direct_prior/output_sha256.txt", ROOT),
        ("m4", ROOT / "artifacts/m4/output_sha256.txt", ROOT),
        ("m5_frozen_inputs", ROOT / "artifacts/m5/frozen_inputs_sha256.txt", ROOT),
        ("m5", ROOT / "artifacts/m5/output_sha256.txt", ROOT),
        ("m6_164_models", ROOT / "artifacts/m6/GSE164522/model_predictions/output_sha256.txt", ROOT),
        ("m6_164", ROOT / "artifacts/m6/GSE164522/output_sha256.txt", ROOT),
        ("m6_203608", ROOT / "artifacts/m6/GSE203608/output_sha256.txt", ROOT),
        ("m6", ROOT / "artifacts/m6/output_sha256.txt", ROOT),
        ("m7", ROOT / "artifacts/m7/output_sha256.txt", ROOT),
        ("m8", ROOT / "artifacts/m8/output_sha256.txt", ROOT),
        ("techfix_programs", ROOT / "artifacts/techfix_v1_f1/programs/output_sha256.txt", ROOT),
        ("techfix_models", ROOT / "artifacts/techfix_v1_f1/models/output_sha256.txt", ROOT),
    ]
    counts: dict[str, int] = {}
    for label, m_path, base in manifests:
        counts[label] = _verify_sha_manifest(m_path, base)
    return counts


def main() -> None:
    t0 = time.perf_counter()
    suites = [
        ("hgnc_and_universes", test_hgnc_harmonization_and_universes),
        ("hierarchical_annotation_and_m2", test_hierarchical_cell_annotation_and_m2),
        ("program_scores_math_and_debiasing", test_program_scores_math_and_debiasing),
        ("prior_mil_and_stage_contracts", test_prior_mil_and_stage_contracts),
        ("m4_m5_m6_v1_1_contracts", test_m4_m5_m6_v1_1_contracts),
        ("m8_figures_source_data_and_manifest", test_m8_figures_source_data_and_manifest),
        ("sha256_manifests", test_sha256_manifests),
    ]
    results = {}
    for name, fn in suites:
        s0 = time.perf_counter()
        out = fn()
        results[name] = {"status": "passed", "seconds": round(time.perf_counter() - s0, 3), "details": out}
        print(f"[PASS] {name} ({results[name]['seconds']}s): {json.dumps(out, ensure_ascii=False)}", flush=True)

    summary = {
        "ok": True,
        "suites_passed": len(suites),
        "total_seconds": round(time.perf_counter() - t0, 3),
        "results": results,
    }
    print("\nALL AUDIT REMEDIATION TESTS PASSED:", json.dumps({"ok": True, "suites": len(suites), "seconds": summary["total_seconds"]}))


if __name__ == "__main__":
    main()
