from __future__ import annotations

from typing import Dict, Optional, Tuple
import gzip
import os
import re

import pandas as pd


def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def _open_tsv_file(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def _iter_tsv_lines(path: str):
    with _open_tsv_file(path) as fh:
        for line_no, raw in enumerate(fh, 1):
            line = raw.rstrip("\n\r")
            if not line or line.startswith("#"):
                continue
            yield line_no, line


def _require_tabular_line(
    path: str, line_no: int, line: str, min_cols: int, desc: str
) -> list[str]:
    if "\t" not in line:
        raise ValueError(
            f"{path}:{line_no} must be TAB-delimited; spaces are not accepted. Expected columns: {desc}"
        )
    parts = line.split("\t")
    if len(parts) < min_cols:
        raise ValueError(
            f"{path}:{line_no} should have at least {min_cols} TAB-delimited columns: {desc}"
        )
    if any(part == "" for part in parts[:min_cols]):
        raise ValueError(
            f"{path}:{line_no} contains an empty required field. Expected columns: {desc}"
        )
    return parts


def read_group_metadata(path: Optional[str]) -> pd.DataFrame:
    """Read sample groups plus optional geographic coordinates.

    Backward-compatible accepted layouts (no header required):
      sample<TAB>group
      sample<TAB>group<TAB>latitude<TAB>longitude
      sample<TAB>group<TAB>latitude<TAB>longitude<TAB>location

    Named latitude/longitude headers may appear in either order. Missing coordinate
    values retain the sample's group but omit that sample from geographic plots.
    Lines beginning with ``#`` are ignored, so a commented column description may be
    included in example files without changing the historical headerless format.
    Supports both uncompressed and .gz compressed files.
    """
    cols = ["Accession", "Type", "Latitude", "Longitude", "Location"]
    if not path:
        return pd.DataFrame(columns=cols)
    rows = []
    seen = set()
    coordinate_indices = (2, 3)
    location_index = 4
    has_header = False
    for line_no, line in _iter_tsv_lines(path):
        parts = _require_tabular_line(
            path,
            line_no,
            line,
            2,
            "sample<TAB>group[<TAB>latitude<TAB>longitude[<TAB>location]]",
        )
        if parts[0].strip().lower() in {"sample", "accession"} and parts[
            1
        ].strip().lower() in {"group", "type"}:
            header = [p.strip().lower() for p in parts]
            latitude = [i for i, h in enumerate(header) if h in {"latitude", "lat", "纬度"}]
            longitude = [i for i, h in enumerate(header) if h in {"longitude", "lon", "lng", "经度"}]
            if len(latitude) > 1 or len(longitude) > 1 or bool(latitude) != bool(longitude):
                raise ValueError(f"{path}:{line_no} requires one latitude and one longitude column")
            coordinate_indices = (latitude[0], longitude[0]) if latitude else None
            location_index = next((i for i, h in enumerate(header) if h in {"location", "地点"}), None)
            has_header = True
            continue
        if not has_header and len(parts) == 3:
            raise ValueError(
                f"{path}:{line_no} provides latitude without longitude. "
                "Use either 2 columns or at least 4 columns."
            )
        sample, group = parts[0].strip(), parts[1].strip()
        if not sample:
            raise ValueError(f"Empty sample ID in {path}:{line_no}")
        if sample in seen:
            raise ValueError(
                f"Duplicate sample in group file {path}:{line_no}: {sample}"
            )
        seen.add(sample)
        # Handle NA in group: ensure it is a string "NA" rather than dropped or corrupted
        if not group or group.strip() in {"", "NA", "NaN", "nan", "None", "<NA>"}:
            group = "NA"
        else:
            group = group.strip()

        lat = lon = pd.NA
        location = pd.NA
        if coordinate_indices is not None:
            lat_index, lon_index = coordinate_indices
            lat_text = parts[lat_index].strip() if lat_index < len(parts) else ""
            lon_text = parts[lon_index].strip() if lon_index < len(parts) else ""
            missing = {"", ".", "na", "nan", "none", "null", "<na>"}
            values = []
            for value, name, limit in [(lat_text, "latitude", 90), (lon_text, "longitude", 180)]:
                if value.lower() in missing:
                    values.append(pd.NA)
                    continue
                try:
                    number = float(value)
                except ValueError as exc:
                    raise ValueError(
                        f"{path}:{line_no} {name} must be numeric: {value!r}"
                    ) from exc
                if not (-limit <= number <= limit):
                    raise ValueError(f"{path}:{line_no} {name} must be between {-limit} and {limit}")
                values.append(number)
            lat, lon = values
        if location_index is not None and location_index < len(parts) and parts[location_index].strip():
            location = parts[location_index].strip()
        rows.append(
            {
                "Accession": sample,
                "Type": group,
                "Latitude": lat,
                "Longitude": lon,
                "Location": location,
            }
        )
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
    with _open_tsv_file(path) as fh:
        for line_no, line in enumerate(fh, 1):
            if not line.strip() or line.startswith("#"):
                continue
            if "\t" not in line:
                raise ValueError(
                    f"{path}:{line_no} must be TAB-delimited; spaces and comma-delimited CSV are not accepted."
                )
            header = [x.strip() for x in line.rstrip("\n\r").split("\t")]
            if len(header) != len(set(header)) or any(not x for x in header):
                raise ValueError(f"Trait file {path} has empty or duplicate headers")
            break
    compression = "gzip" if path.endswith(".gz") else "infer"
    df = pd.read_csv(
        path,
        sep="\t",
        comment="#",
        compression=compression,
        dtype=str,
        keep_default_na=False,
    )
    if df.shape[1] < 2:
        raise ValueError(
            "Trait file should contain a sample column and at least one trait column"
        )
    df.columns = [str(c).strip() for c in df.columns]
    first = df.columns[0]
    if first != "Accession":
        df = df.rename(columns={first: "Accession"})
    df["Accession"] = df["Accession"].astype(str).str.strip()
    if df["Accession"].eq("").any() or df["Accession"].duplicated().any():
        raise ValueError(f"Trait file {path} has empty or duplicate sample IDs")
    if len(set(df.columns)) != len(df.columns):
        raise ValueError(f"Trait file {path} has duplicate column names")
    return df


def sanitize_filename(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", text)
