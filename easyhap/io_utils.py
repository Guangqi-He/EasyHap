from __future__ import annotations

from typing import Dict, Optional
import os
import re

import pandas as pd


def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def _iter_tsv_lines(path: str):
    with open(path, "r", encoding="utf-8") as fh:
        for line_no, raw in enumerate(fh, 1):
            line = raw.rstrip("\n\r")
            if not line or line.startswith("#"):
                continue
            yield line_no, line


def _require_tabular_line(path: str, line_no: int, line: str, min_cols: int, desc: str) -> list[str]:
    if "\t" not in line:
        raise ValueError(
            f"{path}:{line_no} must be TAB-delimited; spaces are not accepted. Expected columns: {desc}"
        )
    parts = line.split("\t")
    if len(parts) < min_cols:
        raise ValueError(f"{path}:{line_no} should have at least {min_cols} TAB-delimited columns: {desc}")
    if any(part == "" for part in parts[:min_cols]):
        raise ValueError(f"{path}:{line_no} contains an empty required field. Expected columns: {desc}")
    return parts


def read_group_file(path: str) -> Dict[str, str]:
    groups: Dict[str, str] = {}
    for line_no, line in _iter_tsv_lines(path):
        parts = _require_tabular_line(path, line_no, line, 2, "sample<TAB>group")
        groups[parts[0]] = parts[1]
    if not groups:
        raise ValueError(f"No samples found in group file: {path}")
    return groups


def read_trait_file(path: Optional[str]) -> Optional[pd.DataFrame]:
    if not path:
        return None
    with open(path, "r", encoding="utf-8") as fh:
        for line_no, line in enumerate(fh, 1):
            if not line.strip() or line.startswith("#"):
                continue
            if "\t" not in line:
                raise ValueError(
                    f"{path}:{line_no} must be TAB-delimited; spaces and comma-delimited CSV are not accepted."
                )
            break
    df = pd.read_csv(path, sep="\t", comment="#")
    if df.shape[1] < 2:
        raise ValueError("Trait file should contain a sample column and at least one trait column")
    first = df.columns[0]
    if first != "Accession":
        df = df.rename(columns={first: "Accession"})
    return df


def sanitize_filename(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", text)
