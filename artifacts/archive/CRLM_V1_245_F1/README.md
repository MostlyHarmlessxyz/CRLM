# CRLM_V1_245_F1 — Codex Original Baseline Snapshot Archive

## 1. Snapshot Provenance

- **Archived Analysis Version**: `CRLM_V1_245_F1`
- **Source Git Commit (HEAD at Archive Time)**: `91a589f5ec58736753e5fd71138e0733e3bc8ab1` (`91a589f`)
- **Local Git Tag**: `CRLM_V1_245_F1_codex_snapshot` -> `91a589f`
- **Original M5 Freeze Commit**: `72ee86bb6a9891609077a08f8350a3e84e102ffd` (`72ee86b`, `freeze_timestamp_utc = 2026-11-21T12:44:41+00:00`)
- **Archive Role**: Read-only historical snapshot of Codex's original M3–M8 delivery prior to Wave 1–Wave 7 audit remediation. Preserved verbatim so downstream three-column comparisons (`CRLM_V1_245_F1` vs `CRLM_V1_245_F1_techfix` vs `CRLM_V1.1_AUDIT_REMEDIATED`) can programmatically read exact original numbers without checking out historical commits.

## 2. Archived Contents

Relative paths inside `artifacts/archive/CRLM_V1_245_F1/` mirror the repository root layout:

- `STATUS.md` — Original Codex status ledger at `91a589f`
- `configs/project.json` — Original frozen `CRLM_V1_245_F1` project configuration
- `artifacts/m3/` — Complete M3 candidate-program discovery & paired-counts outputs
- `artifacts/m4/` — M4 linear-attention MIL, PCA/scGPT ablations, and `HA/CTA/MixMIL` benchmark outputs (excluding `artifacts/m4/selected_cells.csv.gz`; see Section 3)
- `artifacts/m5/` — Frozen `CRLM_V1_245_F1` registration bundle (`candidate_programs.csv`, `registered_tests_template.csv`, `gene_sets/`, `model_transfer/`)
- `artifacts/m6/` — External validation outputs on `GSE164522` and `GSE203608`, `GSE178318` `terms_pending` status record, and `registered_tests_results.csv`
- `artifacts/m7/` — `GSE335552` Xenium spatial measurement audit and functional validation protocol draft
- `artifacts/m8/` — Original figures, `source_data/`, `final_report.md`, `limitations_and_deviations.md`, and `reproducibility_manifest.csv`

## 3. Omitted Cell-Selection Table (Recorded by Original SHA256)

Per archive policy, `artifacts/m4/selected_cells.csv.gz` was not duplicated on disk; its exact original SHA256 and byte count at `91a589f` (also recoverable via `git show CRLM_V1_245_F1_codex_snapshot:artifacts/m4/selected_cells.csv.gz`) are recorded in `omitted_large_files_sha256.txt`:

| Original Repo Path | Size (Bytes) | Original SHA256 |
|---|---:|---|
| `artifacts/m4/selected_cells.csv.gz` | `269198` | `be50de6317c63dbb3535d474e43f598ffecf6ce516fe465474c3ef5d6e65cb62` |

## 4. Integrity Verification

Every file inside this archive (including this `README.md` and `omitted_large_files_sha256.txt`) is checksummed in `archive_sha256.txt`:

```bash
(cd artifacts/archive/CRLM_V1_245_F1 && sha256sum -c archive_sha256.txt)
```
