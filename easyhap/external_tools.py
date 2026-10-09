"""Resolve EasyHap's subprocess tools without depending on the working directory.

Windows builds intentionally never run an unrelated executable from PATH: they
use the vendor binaries embedded by PyInstaller (or the local vendor tree when
running from source).  Linux uses its installed command-line tools first.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Optional, Tuple


def is_windows() -> bool:
    return os.name == "nt" or sys.platform.startswith("win")


def creationflags() -> int:
    return int(getattr(subprocess, "CREATE_NO_WINDOW", 0)) if is_windows() else 0


def _vendor_root() -> Path:
    # In a one-file EXE, external executables are unpacked under _MEIPASS.
    # Do not fall back to the user's CWD/PATH or a side-loaded executable.
    frozen_root = getattr(sys, "_MEIPASS", None)
    if frozen_root:
        return Path(frozen_root) / "vendor"
    return Path(__file__).resolve().parent.parent / "vendor"


def _vendor_exe(folder: str, names: Tuple[str, ...]) -> Optional[str]:
    base = _vendor_root() / folder
    if not base.is_dir():
        return None
    # Accept layouts such as vendor/bcftools/bin/bcftools.exe while preserving
    # the DLLs in their original folders.
    for path in (base, *sorted(p for p in base.rglob('*') if p.is_dir())):
        for name in names:
            item = path / name
            if item.is_file():
                return str(item.resolve())
    return None


def find_bcftools() -> Optional[str]:
    if is_windows():
        return _vendor_exe("bcftools", ("bcftools.exe",))
    # This override is intentionally Linux-only. Windows must use the bundle.
    override = os.environ.get("EASYHAP_BCFTOOLS", "")
    if override and os.path.isfile(override):
        return os.path.abspath(override)
    return shutil.which("bcftools") or _vendor_exe("bcftools", ("bcftools",))


def find_plink(kind: str) -> Tuple[Optional[str], Optional[str]]:
    """Return (executable, plink flavor); PGEN requires the PLINK 2 flavor."""
    if kind not in {"bfile", "pfile"}:
        raise ValueError(f"Invalid PLINK input type: {kind}")
    if is_windows():
        return _vendor_exe("plink2", ("plink2.exe",)), "plink2"
    plink2 = shutil.which("plink2") or _vendor_exe("plink2", ("plink2",))
    if plink2:
        return plink2, "plink2"
    if kind == "bfile":
        plink1 = shutil.which("plink")
        if plink1:
            return plink1, "plink1"
    return None, None
