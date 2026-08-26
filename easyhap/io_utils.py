from __future__ import annotations

from typing import Dict, Optional, Tuple
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


def read_group_metadata(path: Optional[str]) -> pd.DataFrame:
    """Read sample groups plus optional geographic coordinates.

    Backward-compatible accepted layouts (no header required):
      sample<TAB>group
      sample<TAB>group<TAB>latitude<TAB>longitude
      sample<TAB>group<TAB>latitude<TAB>longitude<TAB>location

    Lines beginning with ``#`` are ignored, so a commented column description may be
    included in example files without changing the historical headerless format.
    """
    cols = ["Accession", "Type", "Latitude", "Longitude", "Location"]
    if not path:
        return pd.DataFrame(columns=cols)
    rows = []
    seen = set()
    for line_no, line in _iter_tsv_lines(path):
        parts = _require_tabular_line(
            path, line_no, line, 2,
            "sample<TAB>group[<TAB>latitude<TAB>longitude[<TAB>location]]",
        )
        if len(parts) == 3:
            raise ValueError(
                f"{path}:{line_no} provides latitude without longitude. "
                "Use either 2 columns or at least 4 columns."
            )
        if parts[0].strip().lower() in {"sample", "accession"} and parts[1].strip().lower() in {"group", "type"}:
            continue
        sample, group = parts[0], parts[1]
        if sample in seen:
            raise ValueError(f"Duplicate sample in group file {path}:{line_no}: {sample}")
        seen.add(sample)
        lat = lon = pd.NA
        location = pd.NA
        if len(parts) >= 4:
            try:
                lat = float(parts[2])
                lon = float(parts[3])
            except ValueError as exc:
                raise ValueError(
                    f"{path}:{line_no} latitude/longitude must be numeric: {parts[2]!r}, {parts[3]!r}"
                ) from exc
            if not (-90 <= lat <= 90):
                raise ValueError(f"{path}:{line_no} latitude must be between -90 and 90")
            if not (-180 <= lon <= 180):
                raise ValueError(f"{path}:{line_no} longitude must be between -180 and 180")
            if len(parts) >= 5 and parts[4].strip():
                location = parts[4].strip()
        rows.append({
            "Accession": sample, "Type": group,
            "Latitude": lat, "Longitude": lon, "Location": location,
        })
    if not rows:
        raise ValueError(f"No samples found in group file: {path}")
    return pd.DataFrame(rows, columns=cols)


def read_group_file(path: Optional[str]) -> Dict[str, str]:
    meta = read_group_metadata(path)
    if meta.empty:
        return {}
    return dict(zip(meta["Accession"].astype(str), meta["Type"].astype(str)))


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
