#!/usr/bin/env python3
"""Check the handoff structure only; never mark real research as completed."""
from __future__ import annotations

import csv
import json
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
REQUIRED = [
    "README.md", "AGENTS.md", "STATUS.md", "configs/project.json",
    "docs/CRLM_Research_Plan.md", "docs/IMPLEMENTATION_PLAN.md", "docs/COMPUTE_GUIDE.md",
    "docs/PAPERS_AND_DATA_PREPARATION.md", "docs/source_brief/说明.md",
    "prompts/PLAN_MODE_PROMPT.md", "prompts/EXECUTE_M0_M1_PROMPT.md",
    "references/bibliography.json", "references/reading_manifest.csv",
    "starter/README.txt", "starter/selfcheck.py", "starter/mil_selfcheck.py",
    "starter/registered_tests_selfcheck.py", "scripts/preflight.py",
]


def main():
    errors = [f"Missing required file: {rel}" for rel in REQUIRED if not (ROOT / rel).is_file()]
    stats = {}
    config_path = ROOT / "configs/project.json"
    if config_path.exists():
        try:
            config = json.loads(config_path.read_text(encoding="utf-8"))
            for key, relative in config.get("paths", {}).items():
                resolved = (ROOT / relative).resolve()
                if Path(relative).is_absolute() or not resolved.is_relative_to(ROOT):
                    errors.append(f"Nonportable or out-of-project path: {key}")
                elif not resolved.is_dir():
                    errors.append(f"Missing configured directory: {relative}")
            dataset_rows = config.get("datasets", [])
            accessions = [row["accession"] for row in dataset_rows]
            if len(accessions) != len(set(accessions)):
                errors.append("Duplicate dataset accession; conflicting roles cannot be silently collapsed")
            roles = {row["accession"]: row["role"] for row in dataset_rows}
            if config.get("research_version") in ("V1_DRAFT", "V1_FROZEN", "V1.1_AUDIT_REMEDIATED"):
                developers = {row["accession"] for row in dataset_rows if row["role"] == "development"}
                if developers != {"GSE245552"}:
                    errors.append("V1 development sources must be exactly {GSE245552}")
                expected = {"GSE245552": "development", "GSE164522": "external", "GSE178318": "external", "GSE225857": "excluded"}
                for accession, role in expected.items():
                    if roles.get(accession) != role:
                        errors.append(f"V1 role differs from dataset contract: {accession}")
                if roles.get("GSE203608") not in ("external", "excluded", None):
                    errors.append("V1 GSE203608 must remain external, explicitly excluded, or not selected")
            stats["research_version"] = config.get("research_version")
            stats["milestone_states"] = config.get("milestones")
        except (ValueError, KeyError, TypeError) as error:
            errors.append(f"Invalid project config: {error}")
    refs = ROOT / "references/bibliography.json"
    if refs.exists():
        try:
            records = json.loads(refs.read_text(encoding="utf-8"))
            ids = [item["id"] for item in records]
            stats["bibliography_records"] = len(records)
            stats["unique_nonempty_dois"] = len({str(item.get("doi", "")).lower().strip() for item in records if item.get("doi")})
            if len(ids) != len(set(ids)):
                errors.append("Duplicate bibliography record IDs")
            manifest_path = ROOT / "references/reading_manifest.csv"
            if manifest_path.exists():
                with manifest_path.open(encoding="utf-8-sig", newline="") as handle:
                    rows = list(csv.DictReader(handle))
                id_column = next((name for name in ("id", "reference_id", "ref_id") if rows and name in rows[0]), None)
                if id_column is None:
                    errors.append("Reading manifest has no recognized reference ID column")
                elif sorted(row[id_column] for row in rows) != sorted(ids):
                    errors.append("Reading manifest and bibliography have different record IDs")
                stats["reading_manifest_records"] = len(rows)
        except (ValueError, KeyError, TypeError) as error:
            errors.append(f"Invalid reference inventory: {error}")
    report = ROOT / "docs/CRLM_Research_Plan.md"
    if report.exists():
        text = report.read_text(encoding="utf-8")
        figures = re.findall(r"!\[[^\]]*\]\((report_assets/[^)]+)\)", text)
        for relative in figures:
            if not (report.parent / relative).is_file():
                errors.append(f"Missing report image: {relative}")
        stats["report_figures"] = len(figures)
        if "!FIG:" in text:
            errors.append("Unconverted report figure directive")
    print(json.dumps({"ok": not errors, "scope": "Handoff structure and consistency only, not a scientific validation or M0 acceptance.", "errors": errors, **stats}, ensure_ascii=False, indent=2))
    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())
