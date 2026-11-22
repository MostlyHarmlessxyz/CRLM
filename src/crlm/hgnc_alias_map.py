"""Deterministic HGNC gene symbol harmonization and shared nuclear gene universes.

Provides unambiguous, idempotent mapping across heterogeneous cohort feature axes
(GSE245552, GSE164522, GSE203608, GSE335552), literature/MSigDB gene sets, M4 PCA
features, and external model vocabularies (scGPT whole-human vocab.json), backed
by a SHA256-locked local copy of the official HGNC complete set.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from functools import lru_cache
import hashlib
import json
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import h5py
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
HGNC_DIR = ROOT / "data/manifests/hgnc"
HGNC_TSV_PATH = HGNC_DIR / "hgnc_complete_set.txt"
HGNC_PROVENANCE_PATH = HGNC_DIR / "PROVENANCE.json"

TUMOR_NUCLEAR_245_PATH = HGNC_DIR / "tumor_canonical_nuclear_245.txt"
SHARED_NUCLEAR_245_164_PATH = HGNC_DIR / "shared_nuclear_245_164.txt"
SHARED_NUCLEAR_245_164_203608_PATH = HGNC_DIR / "shared_nuclear_245_164_203608.txt"
EXCLUDED_353_CANONICAL_PATH = HGNC_DIR / "excluded_353_canonical.txt"
HARMONIZATION_AUDIT_CSV_PATH = HGNC_DIR / "cross_cohort_gene_harmonization_audit.csv"

UNIVERSE_FILES = {
    "245_tumor_canonical_nuclear": TUMOR_NUCLEAR_245_PATH,
    "245_164_shared_nuclear": SHARED_NUCLEAR_245_164_PATH,
    "245_164_203608_shared_nuclear": SHARED_NUCLEAR_245_164_203608_PATH,
    "excluded_353_canonical": EXCLUDED_353_CANONICAL_PATH,
}


@dataclass(frozen=True)
class SymbolResolution:
    """Single-feature resolution proposal prior to axis-wide collision checks."""

    source_symbol: str
    ensembl_gene_id: str
    canonical_symbol: str
    hgnc_id: str
    match_type: str
    priority_tier: int
    note: str = ""


def _sha256(path: Path) -> str:
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def _clean_ensembl_id(raw: object) -> str:
    if raw is None:
        return ""
    token = str(raw).strip()
    if not token or token == "NA":
        return ""
    base = token.split(".", 1)[0].strip()
    return base if base.startswith("ENSG") else ""


class HGNCMapper:
    """Deterministic, collision-safe HGNC symbol canonicalizer and vocabulary bridge."""

    def __init__(self, hgnc_tsv: Path = HGNC_TSV_PATH, provenance_json: Path = HGNC_PROVENANCE_PATH, *, verify_sha: bool = True):
        if not hgnc_tsv.is_file():
            raise FileNotFoundError(f"Missing official HGNC table: {hgnc_tsv}")
        self.hgnc_tsv = hgnc_tsv
        self.sha256 = _sha256(hgnc_tsv)
        if verify_sha and provenance_json.is_file():
            expected = json.loads(provenance_json.read_text(encoding="utf-8")).get("sha256")
            if expected and expected != self.sha256:
                raise ValueError(f"HGNC TSV SHA256 mismatch: expected {expected}, got {self.sha256}")

        table = pd.read_csv(hgnc_tsv, sep="\t", dtype=str, keep_default_na=False)
        approved = table.loc[table["status"].eq("Approved")].copy()
        if approved.empty or approved["symbol"].duplicated().any() or approved["hgnc_id"].duplicated().any():
            raise ValueError("Invalid HGNC table: expected unique Approved symbols and hgnc_ids")

        self.symbol_to_hgnc: dict[str, str] = dict(zip(approved["symbol"], approved["hgnc_id"]))
        self.hgnc_to_symbol: dict[str, str] = dict(zip(approved["hgnc_id"], approved["symbol"]))

        ensg_map: dict[str, set[str]] = {}
        prev_map: dict[str, set[str]] = {}
        alias_map: dict[str, set[str]] = {}
        hgnc_prev: dict[str, list[str]] = {}
        hgnc_alias: dict[str, list[str]] = {}

        for row in approved.itertuples(index=False):
            hid = row.hgnc_id
            if row.ensembl_gene_id:
                for item in row.ensembl_gene_id.split("|"):
                    clean = _clean_ensembl_id(item)
                    if clean:
                        ensg_map.setdefault(clean, set()).add(hid)
            p_list: list[str] = []
            if row.prev_symbol:
                for item in row.prev_symbol.split("|"):
                    sym = item.strip().strip('"')
                    if sym:
                        prev_map.setdefault(sym, set()).add(hid)
                        p_list.append(sym)
            a_list: list[str] = []
            if row.alias_symbol:
                for item in row.alias_symbol.split("|"):
                    sym = item.strip().strip('"')
                    if sym:
                        alias_map.setdefault(sym, set()).add(hid)
                        a_list.append(sym)
            hgnc_prev[hid] = p_list
            hgnc_alias[hid] = a_list

        self.ensg_to_hgnc: dict[str, tuple[str, ...]] = {k: tuple(sorted(v)) for k, v in ensg_map.items()}
        self.prev_to_hgnc: dict[str, tuple[str, ...]] = {k: tuple(sorted(v)) for k, v in prev_map.items()}
        self.alias_to_hgnc: dict[str, tuple[str, ...]] = {k: tuple(sorted(v)) for k, v in alias_map.items()}
        self.hgnc_prev_symbols: dict[str, tuple[str, ...]] = {k: tuple(v) for k, v in hgnc_prev.items()}
        self.hgnc_alias_symbols: dict[str, tuple[str, ...]] = {k: tuple(v) for k, v in hgnc_alias.items()}

    def resolve_symbol(self, symbol: object, ensembl_gene_id: object = None) -> SymbolResolution:
        """Resolve a single feature using ENSG -> Approved -> 1:1 prev_symbol -> 1:1 alias_symbol."""
        sym = "" if symbol is None else str(symbol).strip()
        ensg = _clean_ensembl_id(ensembl_gene_id)
        if not sym and not ensg:
            return SymbolResolution("", "", "", "", "empty_symbol", 99, "empty_input")

        if ensg and ensg in self.ensg_to_hgnc:
            hids = self.ensg_to_hgnc[ensg]
            if len(hids) == 1:
                hid = hids[0]
                return SymbolResolution(sym, ensg, self.hgnc_to_symbol[hid], hid, "ensembl_gene_id", 1)
            if sym in self.symbol_to_hgnc and self.symbol_to_hgnc[sym] in hids:
                hid = self.symbol_to_hgnc[sym]
                return SymbolResolution(sym, ensg, self.hgnc_to_symbol[hid], hid, "ensembl_plus_approved_symbol", 1)
            candidates = "|".join(sorted(self.hgnc_to_symbol[h] for h in hids))
            return SymbolResolution(sym, ensg, sym, "", "ambiguous_ensembl", 9, f"ensembl_maps_to:{candidates}")

        if not sym:
            return SymbolResolution("", ensg, "", "", "empty_symbol", 99, "unmapped_ensembl_with_empty_symbol")

        if sym in self.symbol_to_hgnc:
            hid = self.symbol_to_hgnc[sym]
            return SymbolResolution(sym, ensg, sym, hid, "approved_symbol", 2)

        if sym in self.prev_to_hgnc:
            hids = self.prev_to_hgnc[sym]
            if len(hids) == 1:
                hid = hids[0]
                return SymbolResolution(sym, ensg, self.hgnc_to_symbol[hid], hid, "prev_symbol_unique", 3)
            candidates = "|".join(sorted(self.hgnc_to_symbol[h] for h in hids))
            return SymbolResolution(sym, ensg, sym, "", "ambiguous_prev_symbol", 9, f"prev_maps_to:{candidates}")

        if sym in self.alias_to_hgnc:
            hids = self.alias_to_hgnc[sym]
            if len(hids) == 1:
                hid = hids[0]
                return SymbolResolution(sym, ensg, self.hgnc_to_symbol[hid], hid, "alias_symbol_unique", 4)
            candidates = "|".join(sorted(self.hgnc_to_symbol[h] for h in hids))
            return SymbolResolution(sym, ensg, sym, "", "ambiguous_alias_symbol", 9, f"alias_maps_to:{candidates}")

        return SymbolResolution(sym, ensg, sym, "", "unmapped_preserved", 9)

    def normalize_gene_axis(
        self,
        symbols: Sequence[object],
        ensembl_gene_ids: Sequence[object] | None = None,
        *,
        excluded_canonical: set[str] | None = None,
    ) -> pd.DataFrame:
        """Canonicalize a matrix gene axis with deterministic collision protection.

        Guarantees:
        1. Idempotent on unique measurable canonical axes (`f(f(X)) == f(X)`).
        2. Every non-empty `canonical_symbol` in the returned table is strictly unique.
        3. Higher-priority matches (`ensembl_gene_id` > `approved_symbol` >
           `prev_symbol_unique` > `alias_symbol_unique`) retain canonical ownership;
           lower-priority or tied collisions revert safely without overwriting or
           creating duplicate column names.
        """
        n = len(symbols)
        if ensembl_gene_ids is not None and len(ensembl_gene_ids) != n:
            raise ValueError(f"symbols ({n}) and ensembl_gene_ids ({len(ensembl_gene_ids)}) length mismatch")
        ensg_seq = [None] * n if ensembl_gene_ids is None else list(ensembl_gene_ids)
        proposals = [self.resolve_symbol(s, e) for s, e in zip(symbols, ensg_seq)]

        by_target: dict[str, list[tuple[int, int]]] = {}
        for idx, prop in enumerate(proposals):
            if prop.canonical_symbol:
                by_target.setdefault(prop.canonical_symbol, []).append((prop.priority_tier, idx))

        winners: dict[int, str] = {}
        claimed_targets: set[str] = set()
        collision_notes = [p.note for p in proposals]

        for target, items in by_target.items():
            if len(items) == 1:
                winner_idx = items[0][1]
                winners[winner_idx] = target
                claimed_targets.add(target)
                continue
            best_tier = min(tier for tier, _ in items)
            top = [idx for tier, idx in items if tier == best_tier]
            if len(top) == 1:
                winner_idx = top[0]
                winners[winner_idx] = target
                claimed_targets.add(target)
                winner_raw = proposals[winner_idx].source_symbol
                for _, idx in items:
                    if idx != winner_idx:
                        extra = f"target_already_claimed:{target}_by_{winner_raw}(row_{winner_idx})"
                        collision_notes[idx] = f"{collision_notes[idx]}; {extra}".strip("; ")
            else:
                tied_sources = "|".join(proposals[idx].source_symbol for idx in top[:6])
                for _, idx in items:
                    extra = f"same_tier_collision:{target}_among_{tied_sources}"
                    collision_notes[idx] = f"{collision_notes[idx]}; {extra}".strip("; ")

        non_winner_raw_counts: dict[str, int] = {}
        for idx, prop in enumerate(proposals):
            if idx not in winners and prop.source_symbol:
                non_winner_raw_counts[prop.source_symbol] = non_winner_raw_counts.get(prop.source_symbol, 0) + 1

        final_symbols: list[str] = []
        final_hgnc_ids: list[str] = []
        final_types: list[str] = []
        is_measurable: list[bool] = []
        used_symbols: set[str] = set()

        for idx, prop in enumerate(proposals):
            raw = prop.source_symbol
            if not raw and idx not in winners:
                final_symbols.append("")
                final_hgnc_ids.append("")
                final_types.append(prop.match_type)
                is_measurable.append(False)
                continue

            if idx in winners:
                canon = winners[idx]
                final_symbols.append(canon)
                final_hgnc_ids.append(prop.hgnc_id)
                final_types.append(prop.match_type)
                is_measurable.append(True)
                used_symbols.add(canon)
                continue

            # Non-winner: revert to raw symbol if unambiguous and unclaimed; else disambiguate.
            if raw in claimed_targets or raw in used_symbols or non_winner_raw_counts.get(raw, 0) > 1:
                suffix = prop.ensembl_gene_id if prop.ensembl_gene_id else f"row{idx}"
                disamb = f"{raw}::{suffix}"
                while disamb in used_symbols or disamb in claimed_targets:
                    disamb = f"{disamb}_{idx}"
                final_symbols.append(disamb)
                final_hgnc_ids.append("")
                final_types.append(f"{prop.match_type}_shadowed_disambiguated")
                is_measurable.append(False)
                used_symbols.add(disamb)
            else:
                final_symbols.append(raw)
                final_hgnc_ids.append(self.symbol_to_hgnc.get(raw, ""))
                final_types.append(f"{prop.match_type}_reverted_collision")
                is_measurable.append(True)
                used_symbols.add(raw)

        excl = excluded_canonical or set()
        return pd.DataFrame(
            {
                "source_index": np.arange(n, dtype=np.int32),
                "source_symbol": [str(s).strip() if s is not None else "" for s in symbols],
                "ensembl_gene_id": [p.ensembl_gene_id for p in proposals],
                "canonical_symbol": final_symbols,
                "hgnc_id": final_hgnc_ids,
                "match_type": final_types,
                "collision_note": collision_notes,
                "is_unique_measurable": is_measurable,
                "in_excluded_353": [sym in excl for sym in final_symbols],
            }
        )

    def normalize_gene_set(
        self,
        symbols: Iterable[object],
        *,
        deduplicate: bool = True,
    ) -> tuple[list[str], pd.DataFrame]:
        """Normalize a biological signature / marker list to canonical HGNC symbols."""
        items = [str(s).strip() for s in symbols if s is not None and str(s).strip()]
        rows = []
        canonical_list: list[str] = []
        seen: set[str] = set()
        for s in items:
            res = self.resolve_symbol(s)
            dup = res.canonical_symbol in seen
            if not (deduplicate and dup):
                canonical_list.append(res.canonical_symbol)
                seen.add(res.canonical_symbol)
            rows.append(
                {
                    "source_symbol": s,
                    "canonical_symbol": res.canonical_symbol,
                    "hgnc_id": res.hgnc_id,
                    "match_type": res.match_type,
                    "deduplicated_drop": bool(deduplicate and dup),
                    "note": res.note,
                }
            )
        return canonical_list, pd.DataFrame(rows)

    def canonical_to_target_vocab(
        self,
        symbols: Sequence[object],
        target_vocab: Mapping[str, int] | Sequence[object] | set[str],
        ensembl_gene_ids: Sequence[object] | None = None,
    ) -> pd.DataFrame:
        """Deterministically bridge symbols into a target vocabulary (e.g. scGPT vocab.json).

        Guarantees strict 1-to-1 injectivity onto target_vocab: two distinct input
        entries never map to the same target token/index.
        """
        if isinstance(target_vocab, Mapping):
            target_lookup: dict[str, int] = {str(k): int(v) for k, v in target_vocab.items()}
        elif isinstance(target_vocab, set):
            target_lookup = {str(k): i for i, k in enumerate(sorted(target_vocab))}
        else:
            target_lookup = {}
            for idx, item in enumerate(target_vocab):
                key = str(item).strip()
                if key and key not in target_lookup:
                    target_lookup[key] = idx
        target_set = set(target_lookup)

        n = len(symbols)
        ensg_seq = [None] * n if ensembl_gene_ids is None else list(ensembl_gene_ids)
        proposals: list[tuple[tuple[int, int], str | None, str, SymbolResolution]] = []

        for sym_raw, ensg_raw in zip(symbols, ensg_seq):
            sym = "" if sym_raw is None else str(sym_raw).strip()
            res = self.resolve_symbol(sym, ensg_raw)
            if sym in target_set:
                proposals.append(((0, 0), sym, "exact_target_match", res))
                continue
            if not res.hgnc_id:
                proposals.append(((99, 99), None, f"unmapped:{res.match_type}", res))
                continue

            hid = res.hgnc_id
            in_tier = res.priority_tier
            approved_sym = self.hgnc_to_symbol[hid]
            if approved_sym in target_set:
                proposals.append(((in_tier, 1), approved_sym, f"{res.match_type}->target_approved", res))
                continue

            prev_hits = [
                cand
                for cand in self.hgnc_prev_symbols.get(hid, ())
                if cand in target_set and len(self.prev_to_hgnc.get(cand, ())) == 1
            ]
            if len(prev_hits) == 1:
                proposals.append(((in_tier, 2), prev_hits[0], f"{res.match_type}->target_prev_unique", res))
                continue

            alias_hits = [
                cand
                for cand in self.hgnc_alias_symbols.get(hid, ())
                if cand in target_set and len(self.alias_to_hgnc.get(cand, ())) == 1
            ]
            if len(alias_hits) == 1:
                proposals.append(((in_tier, 3), alias_hits[0], f"{res.match_type}->target_alias_unique", res))
                continue

            proposals.append(((99, 99), None, "no_unambiguous_target_token", res))

        by_token: dict[str, list[tuple[tuple[int, int], int]]] = {}
        for idx, (prio, tok, _, _) in enumerate(proposals):
            if tok is not None:
                by_token.setdefault(tok, []).append((prio, idx))

        winners: dict[int, str] = {}
        for tok, items in by_token.items():
            if len(items) == 1:
                winners[items[0][1]] = tok
            else:
                best = min(prio for prio, _ in items)
                top = [idx for prio, idx in items if prio == best]
                if len(top) == 1:
                    winners[top[0]] = tok

        records = []
        for idx, (sym_raw, (_, tok_candidate, btype, res)) in enumerate(zip(symbols, proposals)):
            sym = "" if sym_raw is None else str(sym_raw).strip()
            if idx in winners:
                tok = winners[idx]
                records.append(
                    {
                        "input_index": idx,
                        "input_symbol": sym,
                        "canonical_symbol": res.canonical_symbol,
                        "hgnc_id": res.hgnc_id,
                        "target_symbol": tok,
                        "target_index": target_lookup[tok],
                        "bridge_type": btype,
                        "matched": True,
                    }
                )
            else:
                reason = btype if tok_candidate is None else f"collision_prevented:{tok_candidate}"
                records.append(
                    {
                        "input_index": idx,
                        "input_symbol": sym,
                        "canonical_symbol": res.canonical_symbol,
                        "hgnc_id": res.hgnc_id,
                        "target_symbol": "",
                        "target_index": -1,
                        "bridge_type": reason,
                        "matched": False,
                    }
                )
        return pd.DataFrame(records)

    def align_axis_to_canonical_targets(
        self,
        source_symbols: Sequence[object],
        canonical_targets: Sequence[str],
        source_ensembl_ids: Sequence[object] | None = None,
    ) -> pd.DataFrame:
        """Align a cohort's raw matrix axis onto an ordered list of canonical target symbols.

        Returns one row per entry in `canonical_targets` with `source_index` (0-based
        column index in the cohort matrix, or -1 if unmeasured) and `present` flag.
        """
        axis_df = self.normalize_gene_axis(source_symbols, source_ensembl_ids)
        measurable = axis_df.loc[axis_df["is_unique_measurable"]]
        canon_to_row = {
            row.canonical_symbol: (int(row.source_index), row.source_symbol, row.match_type)
            for row in measurable.itertuples(index=False)
        }
        # Normalize targets as well in case caller passed pre-canonical symbols (e.g. M4 HVGs)
        rows = []
        for pos, target in enumerate(canonical_targets):
            t_str = str(target).strip()
            t_res = self.resolve_symbol(t_str)
            lookup_key = t_str if t_str in canon_to_row else t_res.canonical_symbol
            if lookup_key in canon_to_row:
                src_idx, src_sym, mtype = canon_to_row[lookup_key]
                rows.append(
                    {
                        "target_position": pos,
                        "requested_symbol": t_str,
                        "canonical_symbol": lookup_key,
                        "source_index": src_idx,
                        "source_symbol": src_sym,
                        "match_type": mtype,
                        "present": True,
                    }
                )
            else:
                rows.append(
                    {
                        "target_position": pos,
                        "requested_symbol": t_str,
                        "canonical_symbol": t_res.canonical_symbol,
                        "source_index": -1,
                        "source_symbol": "",
                        "match_type": "absent",
                        "present": False,
                    }
                )
        return pd.DataFrame(rows)


@lru_cache(maxsize=1)
def get_hgnc_mapper() -> HGNCMapper:
    """Return process-cached singleton HGNCMapper verified against PROVENANCE.json."""
    return HGNCMapper()


def normalize_gene_axis(
    symbols: Sequence[object],
    ensembl_gene_ids: Sequence[object] | None = None,
    *,
    excluded_canonical: set[str] | None = None,
) -> pd.DataFrame:
    """Module-level convenience wrapper for HGNCMapper.normalize_gene_axis."""
    return get_hgnc_mapper().normalize_gene_axis(
        symbols, ensembl_gene_ids, excluded_canonical=excluded_canonical
    )


def normalize_gene_set(symbols: Iterable[object], *, deduplicate: bool = True) -> tuple[list[str], pd.DataFrame]:
    """Module-level convenience wrapper for HGNCMapper.normalize_gene_set."""
    return get_hgnc_mapper().normalize_gene_set(symbols, deduplicate=deduplicate)


def canonicalize_symbols(symbols: Iterable[object], *, deduplicate: bool = True) -> list[str]:
    """Return canonical HGNC Approved Symbols for a gene set."""
    return get_hgnc_mapper().normalize_gene_set(symbols, deduplicate=deduplicate)[0]


def canonical_to_target_vocab(
    symbols: Sequence[object],
    target_vocab: Mapping[str, int] | Sequence[object] | set[str],
    ensembl_gene_ids: Sequence[object] | None = None,
) -> pd.DataFrame:
    """Module-level convenience wrapper for HGNCMapper.canonical_to_target_vocab."""
    return get_hgnc_mapper().canonical_to_target_vocab(symbols, target_vocab, ensembl_gene_ids)


def align_axis_to_canonical_targets(
    source_symbols: Sequence[object],
    canonical_targets: Sequence[str],
    source_ensembl_ids: Sequence[object] | None = None,
) -> pd.DataFrame:
    """Module-level convenience wrapper for HGNCMapper.align_axis_to_canonical_targets."""
    return get_hgnc_mapper().align_axis_to_canonical_targets(
        source_symbols, canonical_targets, source_ensembl_ids
    )


@lru_cache(maxsize=8)
def load_shared_nuclear_universe(name: str = "245_164_shared_nuclear") -> tuple[str, ...]:
    """Load one of the precomputed canonical gene universes from data/manifests/hgnc/."""
    if name not in UNIVERSE_FILES:
        raise KeyError(f"Unknown universe {name!r}; choose from {sorted(UNIVERSE_FILES)}")
    path = UNIVERSE_FILES[name]
    if not path.is_file():
        raise FileNotFoundError(f"Universe file {path} not yet generated; run build_hgnc_manifests()")
    return tuple(line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip())


class _LazyUniverseRegistry(Mapping[str, tuple[str, ...]]):
    """Read-only dictionary view of canonical nuclear universes."""

    def __getitem__(self, key: str) -> tuple[str, ...]:
        return load_shared_nuclear_universe(key)

    def __iter__(self):
        return iter(UNIVERSE_FILES)

    def __len__(self) -> int:
        return len(UNIVERSE_FILES)


SHARED_NUCLEAR_UNIVERSE: Mapping[str, tuple[str, ...]] = _LazyUniverseRegistry()


def build_hgnc_manifests() -> dict[str, object]:
    """Generate shared nuclear universes, audit table, and run full verification."""
    mapper = get_hgnc_mapper()
    HGNC_DIR.mkdir(parents=True, exist_ok=True)

    raw_245_tumor = (ROOT / "starter/audit_data/GSE245552_feature_intersection.txt").read_text(encoding="utf-8").splitlines()
    raw_245_extra353 = (ROOT / "starter/audit_data/GSE245552_normal_only_features_353.txt").read_text(encoding="utf-8").splitlines()
    raw_164 = (ROOT / "artifacts/m6/GSE164522/coverage/external_rank_background.txt").read_text(encoding="utf-8").splitlines()
    g203_df = pd.read_csv(ROOT / "data/raw/GSE203608/GSE203608_gene_annotation.csv.gz", dtype=str, keep_default_na=False)
    vocab_scgpt = json.loads((ROOT / "models/scGPT-whole-human/vocab.json").read_text(encoding="utf-8"))

    df_extra353 = mapper.normalize_gene_axis(raw_245_extra353)
    excluded_353_set = set(df_extra353.loc[df_extra353["is_unique_measurable"], "canonical_symbol"])
    if len(excluded_353_set) != 353:
        raise RuntimeError(f"Expected 353 unique canonical normal-only genes, got {len(excluded_353_set)}")

    df_245 = mapper.normalize_gene_axis(raw_245_tumor, excluded_canonical=excluded_353_set)
    df_164 = mapper.normalize_gene_axis(raw_164, excluded_canonical=excluded_353_set)
    df_203 = mapper.normalize_gene_axis(
        g203_df["gene_short_name"].tolist(),
        g203_df["gene_id"].tolist(),
        excluded_canonical=excluded_353_set,
    )

    u245 = set(df_245.loc[df_245["is_unique_measurable"], "canonical_symbol"])
    u164 = set(df_164.loc[df_164["is_unique_measurable"], "canonical_symbol"])
    u203 = set(df_203.loc[df_203["is_unique_measurable"], "canonical_symbol"])

    if u245 & excluded_353_set:
        raise RuntimeError(f"Unexpected overlap between 245 tumor and excluded 353: {u245 & excluded_353_set}")

    tumor_nuclear_245 = sorted(u245 - excluded_353_set)
    shared_245_164 = sorted((u245 & u164) - excluded_353_set)
    shared_245_164_203 = sorted((u245 & u164 & u203) - excluded_353_set)
    excluded_353_sorted = sorted(excluded_353_set)

    TUMOR_NUCLEAR_245_PATH.write_text("\n".join(tumor_nuclear_245) + "\n", encoding="utf-8")
    SHARED_NUCLEAR_245_164_PATH.write_text("\n".join(shared_245_164) + "\n", encoding="utf-8")
    SHARED_NUCLEAR_245_164_203608_PATH.write_text("\n".join(shared_245_164_203) + "\n", encoding="utf-8")
    EXCLUDED_353_CANONICAL_PATH.write_text("\n".join(excluded_353_sorted) + "\n", encoding="utf-8")
    load_shared_nuclear_universe.cache_clear()

    # Also inspect Xenium GSE335552 union features
    xenium_syms, xenium_ids = [], []
    seen_xen: set[tuple[str, str]] = set()
    for h5_path in sorted((ROOT / "data/raw/GSE335552").glob("*.h5")):
        with h5py.File(h5_path, "r") as handle:
            names = [x.decode("utf-8") for x in handle["matrix/features/name"][:]]
            ids = [x.decode("utf-8") for x in handle["matrix/features/id"][:]]
            ftypes = [x.decode("utf-8") for x in handle["matrix/features/feature_type"][:]]
        for nm, eid, ft in zip(names, ids, ftypes):
            if ft == "Gene Expression" and (nm, eid) not in seen_xen:
                seen_xen.add((nm, eid))
                xenium_syms.append(nm)
                xenium_ids.append(eid)
    df_xen = mapper.normalize_gene_axis(xenium_syms, xenium_ids, excluded_canonical=excluded_353_set)

    # Build comprehensive cross-cohort audit CSV across all cohorts
    cohort_frames = []
    for cohort_name, frame in (
        ("GSE245552_tumor", df_245),
        ("GSE245552_normal_extra353", df_extra353.assign(in_excluded_353=True)),
        ("GSE164522", df_164),
        ("GSE203608", df_203),
        ("GSE335552_Xenium", df_xen),
    ):
        sub = frame.copy()
        sub.insert(0, "cohort", cohort_name)
        cohort_frames.append(sub)

    all_rows = pd.concat(cohort_frames, ignore_index=True)
    shared_164_set = set(shared_245_164)
    shared_all_set = set(shared_245_164_203)
    all_rows["in_shared_nuclear_245_164"] = all_rows["canonical_symbol"].isin(shared_164_set) & all_rows["is_unique_measurable"]
    all_rows["in_shared_nuclear_245_164_203608"] = all_rows["canonical_symbol"].isin(shared_all_set) & all_rows["is_unique_measurable"]

    # Attach scGPT target token mapping for every unique canonical symbol
    unique_canons = sorted(set(all_rows.loc[all_rows["canonical_symbol"].ne(""), "canonical_symbol"]))
    scgpt_bridge = mapper.canonical_to_target_vocab(unique_canons, vocab_scgpt)
    canon_to_scgpt = dict(zip(scgpt_bridge["input_symbol"], scgpt_bridge["target_symbol"]))
    all_rows["scgpt_vocab_token"] = all_rows["canonical_symbol"].map(canon_to_scgpt).fillna("")
    all_rows.to_csv(HARMONIZATION_AUDIT_CSV_PATH, index=False)

    # Verification checks: idempotence, uniqueness, M4 PCA HVGs, Mac_SPP1, UPR, scGPT
    for label, frame in (("245", df_245), ("extra353", df_extra353), ("164", df_164), ("203608", df_203), ("Xenium", df_xen)):
        non_empty = frame.loc[frame["canonical_symbol"].ne(""), "canonical_symbol"]
        if not non_empty.is_unique:
            raise AssertionError(f"{label} canonical axis contains duplicate symbols")
        measurable_syms = frame.loc[frame["is_unique_measurable"], "canonical_symbol"].tolist()
        second = mapper.normalize_gene_axis(measurable_syms)
        if not (second["canonical_symbol"] == second["source_symbol"]).all():
            raise AssertionError(f"{label} normalization is not idempotent")

    pca_report = {}
    for lin in ("T_NK", "Mono_Mac_DC", "B_Plasma"):
        hvgs = [str(x) for x in np.load(ROOT / f"models/m4/full/{lin}/pca.npz", allow_pickle=True)["genes"]]
        a164 = mapper.align_axis_to_canonical_targets(raw_164, hvgs)
        a203 = mapper.align_axis_to_canonical_targets(
            g203_df["gene_short_name"].tolist(), hvgs, g203_df["gene_id"].tolist()
        )
        bsc = mapper.canonical_to_target_vocab(hvgs, vocab_scgpt)
        pca_report[lin] = {
            "n_hvg": len(hvgs),
            "before_164": int(sum(g in set(raw_164) for g in hvgs)),
            "after_164": int(a164["present"].sum()),
            "fraction_164": float(a164["present"].mean()),
            "absent_164": a164.loc[~a164["present"], "requested_symbol"].tolist(),
            "before_203608": int(sum(g in set(g203_df.loc[~g203_df["gene_short_name"].duplicated(keep=False), "gene_short_name"]) for g in hvgs)),
            "after_203608": int(a203["present"].sum()),
            "fraction_203608": float(a203["present"].mean()),
            "absent_203608": a203.loc[~a203["present"], "requested_symbol"].tolist(),
            "before_scgpt": int(sum(g in vocab_scgpt for g in hvgs)),
            "after_scgpt": int(bsc["matched"].sum()),
        }
        if pca_report[lin]["fraction_164"] < 0.988 or pca_report[lin]["fraction_203608"] < 0.988:
            raise AssertionError(f"{lin} PCA HVG coverage below 98.8%: {pca_report[lin]}")

    # Program recovery checks (both frozen M5 and raw source seeds/GMT)
    prog_frozen = pd.read_csv(ROOT / "artifacts/m5/gene_sets/program_genes.csv")
    model_frozen = pd.read_csv(ROOT / "artifacts/m5/gene_sets/model_program_genes.csv")
    frozen_all = pd.concat([prog_frozen, model_frozen]).drop_duplicates(["lineage", "program", "gene_symbol"])
    program_report = {}
    for (lin, prog), grp in frozen_all.groupby(["lineage", "program"]):
        gs = grp["gene_symbol"].tolist()
        al_164 = mapper.align_axis_to_canonical_targets(raw_164, gs)
        al_203 = mapper.align_axis_to_canonical_targets(
            g203_df["gene_short_name"].tolist(), gs, g203_df["gene_id"].tolist()
        )
        program_report[f"{lin}::{prog}"] = {
            "n_frozen_245": len(gs),
            "before_164": int(sum(g in set(raw_164) for g in gs)),
            "after_164": int(al_164["present"].sum()),
            "before_203608": int(sum(g in set(g203_df.loc[~g203_df["gene_short_name"].duplicated(keep=False), "gene_short_name"]) for g in gs)),
            "after_203608": int(al_203["present"].sum()),
        }

    # Raw seed/GMT recovery on 245 canonical nuclear
    hallmark = {}
    for line in (ROOT / "artifacts/m3/program_definitions_draft/h.all.v2025.1.Hs.symbols.gmt").read_text(encoding="utf-8").splitlines():
        parts = line.rstrip().split("\t")
        if parts[0] in ("HALLMARK_UNFOLDED_PROTEIN_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE"):
            hallmark[parts[0]] = parts[2:]
    seeds_df = pd.read_csv(ROOT / "starter/literature_seeds/literature_program_seeds.csv")
    raw_sets = {
        "HALLMARK_UNFOLDED_PROTEIN_RESPONSE": hallmark["HALLMARK_UNFOLDED_PROTEIN_RESPONSE"],
        "HALLMARK_INTERFERON_GAMMA_RESPONSE": hallmark["HALLMARK_INTERFERON_GAMMA_RESPONSE"],
        "B19_S4_Mac_SPP1": seeds_df.loc[seeds_df["seed_set_id"].eq("B19_S4_Mac_SPP1"), "gene_symbol_original"].tolist(),
        "B19_S4_Mac_CXCL10": seeds_df.loc[seeds_df["seed_set_id"].eq("B19_S4_Mac_CXCL10"), "gene_symbol_original"].tolist(),
    }
    raw_program_report = {}
    for pname, raw_members in raw_sets.items():
        canon_members = mapper.normalize_gene_set(raw_members)[0]
        raw_program_report[pname] = {
            "n_raw": len(raw_members),
            "n_canonical_dedup": len(canon_members),
            "before_245": int(sum(g in set(raw_245_tumor) for g in raw_members)),
            "after_245_nuclear": int(sum(g in set(tumor_nuclear_245) for g in canon_members)),
            "in_excluded_353": sorted(set(canon_members) & excluded_353_set),
            "after_shared_245_164": int(sum(g in shared_164_set for g in canon_members)),
            "after_shared_245_164_203608": int(sum(g in shared_all_set for g in canon_members)),
        }

    scgpt_245 = mapper.canonical_to_target_vocab(raw_245_tumor, vocab_scgpt)
    scgpt_164 = mapper.canonical_to_target_vocab(raw_164, vocab_scgpt)
    scgpt_shared = mapper.canonical_to_target_vocab(shared_245_164, vocab_scgpt)

    # Update PROVENANCE.json with derived manifest SHA256s
    prov = json.loads(HGNC_PROVENANCE_PATH.read_text(encoding="utf-8"))
    prov["derived_manifests"] = {
        path.name: {
            "path": str(path.relative_to(ROOT)),
            "bytes": path.stat().st_size,
            "lines": sum(1 for _ in path.open("rb")),
            "sha256": _sha256(path),
        }
        for path in (
            TUMOR_NUCLEAR_245_PATH,
            SHARED_NUCLEAR_245_164_PATH,
            SHARED_NUCLEAR_245_164_203608_PATH,
            EXCLUDED_353_CANONICAL_PATH,
            HARMONIZATION_AUDIT_CSV_PATH,
        )
    }
    prov["summary_metrics"] = {
        "n_245_tumor_canonical_nuclear": len(tumor_nuclear_245),
        "n_245_164_shared_nuclear": len(shared_245_164),
        "n_245_164_203608_shared_nuclear": len(shared_245_164_203),
        "n_excluded_353_canonical": len(excluded_353_sorted),
        "n_164_in_excluded_353": int((u164 & excluded_353_set).__len__()),
        "n_203608_in_excluded_353": int((u203 & excluded_353_set).__len__()),
        "n_203608_unique_measurable_before": int((~g203_df["gene_short_name"].duplicated(keep=False) & g203_df["gene_short_name"].ne("")).sum()),
        "n_203608_unique_measurable_after": int(df_203["is_unique_measurable"].sum()),
        "scgpt_245_tumor_matched_before": int(sum(g in vocab_scgpt for g in raw_245_tumor)),
        "scgpt_245_tumor_matched_after": int(scgpt_245["matched"].sum()),
        "scgpt_164_matched_before": int(sum(g in vocab_scgpt for g in raw_164)),
        "scgpt_164_matched_after": int(scgpt_164["matched"].sum()),
        "scgpt_shared_245_164_matched": int(scgpt_shared["matched"].sum()),
        "m4_pca_512_coverage": pca_report,
        "frozen_program_coverage": program_report,
        "raw_program_coverage": raw_program_report,
    }
    HGNC_PROVENANCE_PATH.write_text(json.dumps(prov, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return prov["summary_metrics"]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-manifests", action="store_true", help="Build data/manifests/hgnc/ universe & audit files and run self-checks.")
    args = parser.parse_args()
    if args.build_manifests:
        summary = build_hgnc_manifests()
        print(json.dumps(summary, indent=2, ensure_ascii=False))
    else:
        mapper = get_hgnc_mapper()
        print(f"HGNCMapper ready: {len(mapper.symbol_to_hgnc)} approved symbols, sha256={mapper.sha256}")


if __name__ == "__main__":
    main()
