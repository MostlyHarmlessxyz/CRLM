"""Shared, small I/O helpers. No biological inference is implemented here."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


def read_csv(path):
    # Keep identifiers such as NA, 001, and 1e3 as literal strings.
    return pd.read_csv(path, dtype=str, keep_default_na=False).apply(
        lambda x: x.str.strip() if x.dtype == object else x
    )


def require_columns(frame, columns, name="input"):
    missing = sorted(set(columns) - set(frame.columns))
    if missing:
        raise ValueError(f"{name}: missing columns: {missing}")
    if frame.empty:
        raise ValueError(f"{name}: no rows")


def require_nonempty(frame, columns, name="input"):
    bad = frame[list(columns)].eq("").any(axis=1)
    if bad.any():
        raise ValueError(f"{name}: empty required value in CSV row(s) {(np.flatnonzero(bad) + 2).tolist()[:10]}")


def bool_column(series, name):
    mapping = {"true": True, "false": False, "1": True, "0": False}
    normalized = series.astype(str).str.lower().str.strip()
    bad = ~normalized.isin(mapping)
    if bad.any():
        raise ValueError(f"{name}: use true/false or 1/0, found {series[bad].tolist()[:5]}")
    return normalized.map(mapping).astype(bool)


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def stable_seed(seed, *parts):
    payload = json.dumps([int(seed), *map(str, parts)], ensure_ascii=False).encode()
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], "big")


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def cli_error(parser, callback):
    try:
        callback()
    except (ValueError, FileNotFoundError, KeyError, json.JSONDecodeError) as error:
        parser.exit(2, f"ERROR: {error}\n")
