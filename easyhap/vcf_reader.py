from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Tuple
import bisect
import gzip
import os
import re
import shutil
import subprocess
import sys
import tempfile

MISSING_ALLELE = "N"
ABSENCE_ALLELE = "ABS"


@dataclass
class Region:
    chrom: str
    start: int
    end: int

    @property
    def label(self) -> str:
        return f"{self.chrom}_{self.start}_{self.end}"

    @property
    def vcf_label(self) -> str:
        return f"{self.chrom}:{self.start}-{self.end}"


@dataclass
class VariantCall:
    chrom: str
    pos: int
    vid: str
    ref: str
    alts: Tuple[str, ...]
    allele_tokens: Tuple[str, ...]
    genotypes: Dict[str, Tuple[Optional[int], ...]]
    phased: Dict[str, bool]
    info: Dict[str, str]

    @property
    def variant_id(self) -> str:
        return self.vid if self.vid and self.vid != "." else f"{self.chrom}:{self.pos}"

    @property
    def alt_tokens(self) -> Tuple[str, ...]:
        return self.allele_tokens[1:]


def parse_region_string(text: str) -> Region:
    m = re.match(r"^([^:]+):(\d+)-(\d+)$", text.strip())
    if not m:
        raise ValueError(f"Invalid region string: {text!r}; expected Chr10:1-500")
    chrom, start, end = m.group(1), int(m.group(2)), int(m.group(3))
    if start < 1 or end < start:
        raise ValueError(f"Invalid region coordinate: {text!r}")
    return Region(chrom, start, end)


def read_regions(region: Optional[str] = None, region_file: Optional[str] = None) -> List[Region]:
    regions: List[Region] = []
    if region:
        regions.append(parse_region_string(region))
    if region_file:
        opener = gzip.open if region_file.endswith(".gz") else open
        with opener(region_file, "rt") as fh:
            for line_no, raw in enumerate(fh, 1):
                line = raw.rstrip("\n\r")
                if not line or line.startswith("#"):
                    continue
                if "\t" not in line:
                    raise ValueError(
                        f"{region_file}:{line_no} must be TAB-delimited; spaces are not accepted. "
                        "Expected columns: chr<TAB>start<TAB>end"
                    )
                parts = line.split("\t")
                if len(parts) < 3 or not parts[0] or not parts[1] or not parts[2]:
                    raise ValueError(
                        f"{region_file}:{line_no} should have at least 3 TAB-delimited columns: chr<TAB>start<TAB>end"
                    )
                regions.append(Region(parts[0], int(parts[1]), int(parts[2])))
    if not regions:
        raise ValueError("Please provide --region or --region-file")
    return regions


def _parse_info(info: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    if not info or info == ".":
        return out
    for item in info.split(";"):
        if not item:
            continue
        if "=" in item:
            k, v = item.split("=", 1)
            out[k] = v
        else:
            out[item] = "True"
    return out


def _info_get_int(info: Dict[str, str], key: str) -> Optional[int]:
    val = info.get(key)
    if val is None:
        return None
    if isinstance(val, (list, tuple)):
        val = val[0] if val else None
    if val is None:
        return None
    try:
        text = str(val).split(",")[0]
        m = re.search(r"[-+]?\d+", text)
        return int(m.group(0)) if m else None
    except Exception:
        return None


def _symbolic_token(allele: str, info: Optional[Dict[str, str]] = None, allele_index: Optional[int] = None) -> str:
    """Return a compact but allele-identity-preserving symbolic-variant token."""
    name = allele.strip("<>").upper() or "SV"
    svlen = _info_get_int(info or {}, "SVLEN")
    suffix = f":{svlen:+d}" if svlen not in (None, 0) else ""
    prefix = f"ALT{allele_index}:" if allele_index is not None else ""
    return f"{prefix}{name}{suffix}"


def encode_alleles(ref: str, alts: Sequence[str], info: Optional[Dict[str, str]] = None, max_literal_len: int = 20) -> Tuple[str, ...]:
    """Encode alleles while preserving every distinct VCF allele state.

    The first token is the REF state.  ALT tokens include their VCF allele index
    for indels/SVs so equal-length but sequence-distinct ALT alleles can never
    collapse to the same EasyHap state.  The VCF spanning-deletion allele `*`
    is structural absence (ABS), not an ordinary missing genotype.
    """
    ref = (ref or "N").upper()
    regular_alt_lengths = [len(a) for a in alts if a and a not in {".", "*"} and not (a.startswith("<") and a.endswith(">"))]
    has_indel_or_symbolic = any(length != len(ref) for length in regular_alt_lengths) or any(
        bool(a) and (a == "*" or (a.startswith("<") and a.endswith(">"))) for a in alts
    )
    if has_indel_or_symbolic and (len(ref) > 1 or ref == "N"):
        ref_token = "REF" if ref == "N" else f"REF{len(ref)}"
    elif len(ref) <= max_literal_len and all(c in "ACGTN" for c in ref):
        ref_token = ref
    else:
        ref_token = f"REF{len(ref)}"
    tokens: List[str] = [ref_token]

    for idx, alt0 in enumerate(alts, 1):
        alt = alt0 or "."
        if alt == ".":
            tokens.append(MISSING_ALLELE)
            continue
        if alt == "*":
            tokens.append(ABSENCE_ALLELE)
            continue
        if alt.startswith("<") and alt.endswith(">"):
            tokens.append(_symbolic_token(alt, info, idx))
            continue
        alt_u = alt.upper()
        diff = len(alt_u) - len(ref)
        if diff == 0 and len(alt_u) <= max_literal_len and all(c in "ACGTN" for c in alt_u):
            tokens.append(alt_u)
            continue
        if diff > 0:
            kind = f"INS+{diff}"
        elif diff < 0:
            kind = f"DEL-{abs(diff)}"
        else:
            kind = f"SEQ{len(alt_u)}"
        literal = alt_u if len(alt_u) <= max_literal_len else alt_u[:8] + "..." + alt_u[-8:]
        tokens.append(f"ALT{idx}:{kind}:{literal}")
    return tuple(tokens)


def deletion_allele_indices(call: VariantCall) -> Tuple[int, ...]:
    """VCF allele indices representing deletion/absence events at this record."""
    out: List[int] = []
    svtype = str(call.info.get("SVTYPE", "")).upper()
    for idx, alt in enumerate(call.alts, 1):
        a = (alt or "").upper()
        if a == "*" or "<DEL" in a or "<CN0" in a or "<ABS" in a:
            out.append(idx)
        elif not a.startswith("<") and a not in {"", "."} and len(a) < len(call.ref):
            out.append(idx)
        elif svtype == "DEL" and idx == 1:
            out.append(idx)
    return tuple(sorted(set(out)))


def deletion_end(call: VariantCall) -> int:
    end = _info_get_int(call.info, "END")
    if end is not None and end >= call.pos:
        return end
    return call.pos + max(1, len(call.ref)) - 1

def token_for_gt(allele_tokens: Sequence[str], allele_index: Optional[int]) -> str:
    if allele_index is None or allele_index < 0:
        return MISSING_ALLELE
    if allele_index >= len(allele_tokens):
        return MISSING_ALLELE
    return allele_tokens[allele_index]


class VCFReader:
    """VCF/BCF reader with platform-aware backend selection.

    Backend policy
    --------------
    * Windows + BGZF-compressed VCF/BCF: prefer bundled/system bcftools.
      Missing regional indexes are created automatically with ``bcftools index -c``.
    * Windows + small uncompressed VCF (<= 100 MiB by default): use the plain
      parser to avoid unnecessary preprocessing.
    * Windows + large uncompressed VCF: fail fast with an actionable message
      instead of silently scanning a multi-GB file for every region.
    * Non-Windows ``auto`` mode preserves the original cyvcf2 -> pysam -> plain
      selection logic.

    The explicit ``prefer`` values are: ``auto``, ``bcftools``, ``cyvcf2``,
    ``pysam``, and ``plain``.
    """

    WINDOWS_PLAIN_MAX_BYTES = 100 * 1024 * 1024  # 100 MiB

    def __init__(self, path: str, prefer: str = "auto") -> None:
        self.path = os.fspath(path)
        self.prefer = (prefer or "auto").lower()
        if self.prefer not in {"auto", "bcftools", "cyvcf2", "pysam", "plain"}:
            raise ValueError(
                f"Unsupported VCF backend {prefer!r}; expected auto, bcftools, cyvcf2, pysam, or plain"
            )
        if not os.path.isfile(self.path):
            raise FileNotFoundError(f"VCF/BCF file not found: {self.path}")

        self.backend = "plain"
        self.backend_detail = ""
        self._vcf = None
        self._bcftools: Optional[str] = None
        self.samples: List[str] = []
        self.index_created = False
        self.index_path: Optional[str] = self._existing_region_index()
        self._init_backend()

    @staticmethod
    def _is_windows() -> bool:
        return os.name == "nt" or sys.platform.startswith("win")

    def _file_size(self) -> int:
        try:
            return os.path.getsize(self.path)
        except OSError:
            return 0

    def _existing_region_index(self) -> Optional[str]:
        """Return the first adjacent tabix/CSI index, if present."""
        for candidate in (self.path + ".csi", self.path + ".tbi"):
            if os.path.isfile(candidate):
                return candidate
        return None

    def _has_region_index(self) -> bool:
        self.index_path = self._existing_region_index()
        return self.index_path is not None

    def _is_binary_bcf(self) -> bool:
        return self.path.lower().endswith(".bcf")

    def _is_gzip_file(self) -> bool:
        """Detect gzip/BGZF by magic bytes, not only by file extension."""
        if self._is_binary_bcf():
            return False
        try:
            with open(self.path, "rb") as fh:
                return fh.read(2) == b"\x1f\x8b"
        except OSError:
            return False

    def _is_compressed_text_vcf(self) -> bool:
        lower = self.path.lower()
        return (
            lower.endswith((".vcf.gz", ".vcf.bgz", ".vcf.bgzip"))
            or self._is_gzip_file()
        )

    def _is_uncompressed_text_vcf(self) -> bool:
        lower = self.path.lower()
        return lower.endswith(".vcf") and not self._is_gzip_file()

    def _is_text_vcf(self) -> bool:
        return self._is_uncompressed_text_vcf() or self._is_compressed_text_vcf()

    @staticmethod
    def _candidate_bcftools_names() -> Tuple[str, ...]:
        return ("bcftools.exe", "bcftools") if os.name == "nt" else ("bcftools", "bcftools.exe")

    def _find_bcftools(self) -> Optional[str]:
        """Locate bcftools in an EasyHap/PyInstaller bundle or on PATH."""
        env_path = os.environ.get("EASYHAP_BCFTOOLS")
        if env_path and os.path.isfile(env_path):
            return os.path.abspath(env_path)

        roots: List[str] = []
        meipass = getattr(sys, "_MEIPASS", None)
        if meipass:
            roots.append(os.fspath(meipass))
        if getattr(sys, "executable", None):
            roots.append(os.path.dirname(os.path.abspath(sys.executable)))
        roots.append(os.path.dirname(os.path.abspath(__file__)))

        rel_dirs = ("", "bin", "tools", "bcftools", os.path.join("vendor", "bcftools"))
        seen = set()
        for root in roots:
            for rel in rel_dirs:
                for name in self._candidate_bcftools_names():
                    candidate = os.path.abspath(os.path.join(root, rel, name))
                    if candidate in seen:
                        continue
                    seen.add(candidate)
                    if os.path.isfile(candidate):
                        return candidate

        for name in self._candidate_bcftools_names():
            found = shutil.which(name)
            if found:
                return os.path.abspath(found)
        return None

    @staticmethod
    def _subprocess_creationflags() -> int:
        # Avoid flashing a console window when the GUI invokes bundled bcftools.
        return int(getattr(subprocess, "CREATE_NO_WINDOW", 0)) if os.name == "nt" else 0

    def _run_bcftools(self, args: Sequence[str], *, context: str) -> subprocess.CompletedProcess:
        assert self._bcftools is not None
        cmd = [self._bcftools, *map(str, args)]
        proc = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            encoding="utf-8",
            errors="replace",
            creationflags=self._subprocess_creationflags(),
        )
        if proc.returncode != 0:
            stderr = (proc.stderr or "").strip()
            raise RuntimeError(
                f"bcftools failed while {context} (exit code {proc.returncode}).\n"
                f"Command: {' '.join(cmd)}\n"
                f"{stderr or 'No error message was returned by bcftools.'}"
            )
        return proc

    def _ensure_bcftools_index(self) -> None:
        """Create a CSI index for BGZF VCF/BCF when no adjacent index exists."""
        if self._has_region_index():
            return
        if self._is_uncompressed_text_vcf():
            raise RuntimeError(
                f"Regional random access requires a BGZF-compressed VCF or BCF, but {self.path!r} is uncompressed. "
                "Compress it with bgzip and then index it, or use --vcf-backend plain for a small file."
            )

        try:
            self._run_bcftools(["index", "-f", "-c", self.path], context=f"indexing {self.path!r}")
        except RuntimeError as exc:
            text = str(exc)
            if self._is_compressed_text_vcf():
                text += (
                    "\nEasyHap detected a compressed VCF but bcftools could not index it. "
                    "A common cause is ordinary gzip compression rather than BGZF. "
                    "Re-compress the VCF with bgzip and try again."
                )
            raise RuntimeError(text) from exc

        self.index_path = self._existing_region_index()
        if self.index_path is None:
            raise RuntimeError(
                f"bcftools index completed but no .csi/.tbi index was found beside {self.path!r}."
            )
        self.index_created = True

    def _read_samples_bcftools(self) -> List[str]:
        proc = self._run_bcftools(["query", "-l", self.path], context=f"reading samples from {self.path!r}")
        return [line.strip() for line in proc.stdout.splitlines() if line.strip()]

    def _init_bcftools(self, *, create_index: bool = True) -> bool:
        self._bcftools = self._find_bcftools()
        if not self._bcftools:
            return False
        if create_index and (self._is_compressed_text_vcf() or self._is_binary_bcf()):
            self._ensure_bcftools_index()
        self.samples = self._read_samples_bcftools()
        self.backend = "bcftools"
        index_note = self.index_path or "none"
        created_note = " (created by EasyHap)" if self.index_created else ""
        self.backend_detail = f"bcftools={self._bcftools}; index={index_note}{created_note}"
        return True

    def _try_cyvcf2(self) -> bool:
        try:
            from cyvcf2 import VCF  # type: ignore

            self._vcf = VCF(self.path)
            self.samples = list(self._vcf.samples)
            self.backend = "cyvcf2"
            self.backend_detail = "cyvcf2"
            return True
        except Exception:
            self._vcf = None
            return False

    def _try_pysam(self) -> bool:
        try:
            import pysam  # type: ignore

            self._vcf = pysam.VariantFile(self.path)
            self.samples = list(self._vcf.header.samples)
            self.backend = "pysam"
            self.backend_detail = "pysam"
            return True
        except Exception:
            self._vcf = None
            return False

    def _init_backend(self) -> None:
        # Explicit choices always take precedence over automatic platform policy.
        if self.prefer == "plain":
            if self._is_binary_bcf():
                raise RuntimeError("BCF is binary and cannot be read by the plain-text backend.")
            self.samples = self._read_samples_plain()
            self.backend = "plain"
            self.backend_detail = "plain (explicit)"
            return

        if self.prefer == "bcftools":
            if not self._init_bcftools(create_index=True):
                raise RuntimeError(
                    "The bcftools backend was requested but bcftools could not be found. "
                    "For the Windows EXE, bundle bcftools.exe with EasyHap; for source installs, add bcftools to PATH "
                    "or set EASYHAP_BCFTOOLS to its full path."
                )
            return

        if self.prefer == "cyvcf2":
            if not self._try_cyvcf2():
                raise RuntimeError(f"Unable to open {self.path!r} with cyvcf2")
            return

        if self.prefer == "pysam":
            if not self._try_pysam():
                raise RuntimeError(f"Unable to open {self.path!r} with pysam")
            return

        # Windows auto policy: compressed VCF/BCF should use indexed bcftools
        # access. This avoids the very expensive plain scan of a multi-GB VCF.
        if self._is_windows():
            if self._is_compressed_text_vcf() or self._is_binary_bcf():
                if self._init_bcftools(create_index=True):
                    return

                # Source/debug installations may not bundle bcftools. If a valid
                # index already exists, cyvcf2/pysam are still safe random-access
                # fallbacks. Never silently plain-scan a large compressed file.
                if self._has_region_index():
                    if self._try_cyvcf2() or self._try_pysam():
                        return

                raise RuntimeError(
                    "EasyHap is running on Windows with a compressed VCF/BCF, but bcftools could not be found and "
                    "no indexed cyvcf2/pysam backend could be opened. Refusing to sequentially scan a large file. "
                    "Bundle bcftools.exe with the EasyHap executable, add bcftools to PATH, or set EASYHAP_BCFTOOLS."
                )

            if self._is_uncompressed_text_vcf():
                if self._file_size() <= self.WINDOWS_PLAIN_MAX_BYTES:
                    self.samples = self._read_samples_plain()
                    self.backend = "plain"
                    self.backend_detail = (
                        f"plain (Windows uncompressed VCF <= {self.WINDOWS_PLAIN_MAX_BYTES // (1024 * 1024)} MiB)"
                    )
                    return
                raise RuntimeError(
                    f"The uncompressed VCF is {self._file_size() / (1024 * 1024):.1f} MiB, larger than EasyHap's "
                    f"{self.WINDOWS_PLAIN_MAX_BYTES // (1024 * 1024)} MiB plain-reader threshold on Windows. "
                    "To avoid repeatedly scanning a large VCF, bgzip-compress it to .vcf.gz and rerun EasyHap; "
                    "EasyHap will create a CSI index automatically and use bcftools regional access."
                )

        # Non-Windows auto behavior: preserve the original implementation.
        # Unindexed textual VCF remains usable through the sequential parser.
        if self._is_text_vcf() and not self._has_region_index():
            self.samples = self._read_samples_plain()
            self.backend = "plain"
            self.backend_detail = "plain (unindexed text VCF)"
            return

        if self._try_cyvcf2():
            return
        if self._try_pysam():
            return

        if self._is_binary_bcf():
            raise RuntimeError(
                f"Unable to open BCF file {self.path!r}. Install/repair cyvcf2 or pysam, or use the bcftools backend. "
                "For fast regional access, create a CSI index."
            )

        self.samples = self._read_samples_plain()
        self.backend = "plain"
        self.backend_detail = "plain (fallback)"

    def _open_text(self):
        compressed = self._is_compressed_text_vcf()
        return gzip.open(self.path, "rt") if compressed else open(self.path, "rt")

    def _read_samples_plain(self) -> List[str]:
        with self._open_text() as fh:
            for line in fh:
                if line.startswith("#CHROM"):
                    parts = line.rstrip("\n").split("\t")
                    return parts[9:]
        raise ValueError(f"No #CHROM header line found in {self.path}")

    def iter_region(self, region: Region) -> Iterator[VariantCall]:
        if self.backend == "bcftools":
            yield from self._iter_bcftools(region)
        elif self.backend == "cyvcf2":
            yield from self._iter_cyvcf2(region)
        elif self.backend == "pysam":
            yield from self._iter_pysam(region)
        else:
            yield from self._iter_plain(region)

    def _parse_vcf_record_line(self, line: str) -> Optional[VariantCall]:
        """Parse one data line in VCF text format into EasyHap's VariantCall."""
        if not line or line.startswith("#"):
            return None
        parts = line.rstrip("\n\r").split("\t")
        if len(parts) < 8:
            return None
        chrom, pos_s, vid, ref, alt_s, _qual, _filt, info_s = parts[:8]
        try:
            pos = int(pos_s)
        except ValueError:
            return None
        fmt = parts[8].split(":") if len(parts) > 8 else []
        gt_idx = fmt.index("GT") if "GT" in fmt else None
        info = _parse_info(info_s)
        alts = tuple([] if alt_s == "." else alt_s.split(","))
        tokens = encode_alleles(ref, alts, info)
        genotypes: Dict[str, Tuple[Optional[int], ...]] = {}
        phased: Dict[str, bool] = {}
        for sample, sample_field in zip(self.samples, parts[9:]):
            gt_text = "."
            if gt_idx is not None:
                fields = sample_field.split(":")
                if gt_idx < len(fields):
                    gt_text = fields[gt_idx]
            sep = "|" if "|" in gt_text else "/"
            phased[sample] = sep == "|"
            alleles: List[Optional[int]] = []
            for a in re.split(r"[|/]", gt_text):
                if a in {".", ""}:
                    alleles.append(None)
                else:
                    try:
                        alleles.append(int(a))
                    except ValueError:
                        alleles.append(None)
            genotypes[sample] = tuple(alleles)
        return VariantCall(chrom, pos, vid, ref, alts, tokens, genotypes, phased, info)

    def _iter_bcftools(self, region: Region) -> Iterator[VariantCall]:
        assert self._bcftools is not None
        if not self._has_region_index():
            # This should normally have been handled during initialization, but
            # re-check in case files were moved/replaced after reader creation.
            self._ensure_bcftools_index()

        cmd = [self._bcftools, "view", "-H", "-r", region.vcf_label, self.path]
        # stderr goes to a temporary file so a verbose bcftools process can
        # never deadlock while stdout is streamed record-by-record.
        with tempfile.TemporaryFile(mode="w+b") as err_fh:
            proc = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=err_fh,
                text=True,
                encoding="utf-8",
                errors="replace",
                creationflags=self._subprocess_creationflags(),
            )
            assert proc.stdout is not None
            try:
                for line in proc.stdout:
                    call = self._parse_vcf_record_line(line)
                    if call is not None:
                        yield call
            finally:
                proc.stdout.close()
            returncode = proc.wait()
            if returncode != 0:
                err_fh.seek(0)
                stderr = err_fh.read().decode("utf-8", errors="replace").strip()
                raise RuntimeError(
                    f"bcftools failed while reading region {region.vcf_label!r} from {self.path!r} "
                    f"(exit code {returncode}).\n{stderr or 'No error message was returned by bcftools.'}"
                )

    def _iter_cyvcf2(self, region: Region) -> Iterator[VariantCall]:
        assert self._vcf is not None

        # Indexed inputs use true random access. For an explicitly selected
        # cyvcf2 backend on an unindexed input, reopen the file and scan it
        # sequentially for every region so batch analyses remain correct.
        if self._has_region_index():
            iterator: Iterable = self._vcf(region.vcf_label)
        else:
            from cyvcf2 import VCF  # type: ignore
            iterator = VCF(self.path)

        for rec in iterator:
            if rec.CHROM != region.chrom or rec.POS < region.start or rec.POS > region.end:
                continue
            info = {k: str(v) for k, v in dict(rec.INFO).items()}
            alts = tuple(rec.ALT or [])
            tokens = encode_alleles(rec.REF, alts, info)
            genotypes: Dict[str, Tuple[Optional[int], ...]] = {}
            phased: Dict[str, bool] = {}
            for sample, gt in zip(self.samples, rec.genotypes):
                if not gt:
                    genotypes[sample] = tuple()
                    phased[sample] = False
                    continue
                # cyvcf2 appends phased flag as the last element.
                gt_alleles = tuple(None if a is None or int(a) < 0 else int(a) for a in gt[:-1])
                genotypes[sample] = gt_alleles
                phased[sample] = bool(gt[-1])
            yield VariantCall(rec.CHROM, int(rec.POS), rec.ID or ".", rec.REF, alts, tokens, genotypes, phased, info)

    def _iter_pysam(self, region: Region) -> Iterator[VariantCall]:
        assert self._vcf is not None
        if self._has_region_index():
            iterator = self._vcf.fetch(region.chrom, region.start - 1, region.end)
        else:
            # Reopen for each region because a sequential iterator is consumed
            # after one pass and region-file analysis may contain many regions.
            import pysam  # type: ignore
            fresh = pysam.VariantFile(self.path)
            iterator = iter(fresh)
        for rec in iterator:
            pos = int(rec.pos)
            if rec.chrom != region.chrom or pos < region.start or pos > region.end:
                continue
            info = {k: str(v) for k, v in dict(rec.info).items()}
            alts = tuple(rec.alts or [])
            tokens = encode_alleles(rec.ref, alts, info)
            genotypes: Dict[str, Tuple[Optional[int], ...]] = {}
            phased: Dict[str, bool] = {}
            for sample in self.samples:
                call = rec.samples[sample]
                gt = call.get("GT")
                genotypes[sample] = tuple(None if a is None or int(a) < 0 else int(a) for a in (gt or tuple()))
                phased[sample] = bool(getattr(call, "phased", False))
            yield VariantCall(rec.chrom, pos, rec.id or ".", rec.ref, alts, tokens, genotypes, phased, info)

    def _iter_plain(self, region: Region) -> Iterator[VariantCall]:
        with self._open_text() as fh:
            for line in fh:
                if not line or line.startswith("#"):
                    continue
                call = self._parse_vcf_record_line(line)
                if call is None:
                    continue
                if call.chrom != region.chrom or call.pos < region.start or call.pos > region.end:
                    continue
                yield call


class PlinkBedReader:
    """PLINK 1 binary BED/BIM/FAM reader for phase-independent analysis.

    PLINK BED is biallelic and unphased. EasyHap therefore exposes it only to
    genotype/inbred mode. BIM allele 1 is represented as ALT and allele 2 as
    REF to match PLINK's own VCF export convention; the output metadata records
    this explicitly and does not claim that BIM allele 2 is a reference-genome
    allele.

    Region access is indexed from BIM coordinates. Instead of scanning the
    complete BED matrix for every requested interval, ``iter_region`` locates
    only the relevant variant records and seeks directly to their fixed-width
    SNP-major BED blocks. This is especially important for batch analyses of
    large PLINK datasets.
    """

    def __init__(self, prefix: str) -> None:
        self.prefix = self._normalize_prefix(prefix)
        self.bed = self.prefix + ".bed"
        self.bim = self.prefix + ".bim"
        self.fam = self.prefix + ".fam"
        for path in (self.bed, self.bim, self.fam):
            if not os.path.exists(path):
                raise FileNotFoundError(f"Missing PLINK file: {path}")
        self.samples = self._read_fam()
        self.variants = self._read_bim()
        self._variant_index = self._build_variant_index()
        self.backend = "plink-bed"
        self._bytes_per_variant = (len(self.samples) + 3) // 4
        expected = 3 + len(self.variants) * self._bytes_per_variant
        actual = os.path.getsize(self.bed)
        if actual < expected:
            raise ValueError(f"PLINK BED appears truncated: expected at least {expected} bytes, found {actual}")
        with open(self.bed, "rb") as fh:
            magic = fh.read(3)
        if magic != bytes((0x6C, 0x1B, 0x01)):
            raise ValueError("Unsupported PLINK BED: expected current SNP-major magic bytes 6c 1b 01")

    @staticmethod
    def _normalize_prefix(prefix: str) -> str:
        """Accept either a PLINK prefix or one member of a BED/BIM/FAM trio."""
        value = os.fspath(prefix).strip()
        lower = value.lower()
        for suffix in (".bed", ".bim", ".fam"):
            if lower.endswith(suffix):
                return value[:-len(suffix)]
        return value

    def _read_fam(self) -> List[str]:
        samples: List[str] = []
        seen = set()
        with open(self.fam, "rt", encoding="utf-8") as fh:
            for line_no, raw in enumerate(fh, 1):
                if not raw.strip():
                    continue
                parts = raw.split()
                if len(parts) < 2:
                    raise ValueError(f"{self.fam}:{line_no} requires at least FID and IID")
                iid = parts[1]
                sid = iid if iid not in seen else f"{parts[0]}:{iid}"
                if sid in seen:
                    raise ValueError(f"Duplicate PLINK sample identifier after FID/IID resolution: {sid}")
                seen.add(sid)
                samples.append(sid)
        return samples

    def _read_bim(self) -> List[Tuple[str, str, int, str, str]]:
        rows: List[Tuple[str, str, int, str, str]] = []
        with open(self.bim, "rt", encoding="utf-8") as fh:
            for line_no, raw in enumerate(fh, 1):
                if not raw.strip():
                    continue
                parts = raw.split()
                if len(parts) < 6:
                    raise ValueError(f"{self.bim}:{line_no} requires 6 columns")
                chrom, vid, _cm, pos_s, a1, a2 = parts[:6]
                try:
                    pos = int(pos_s)
                except ValueError as exc:
                    raise ValueError(f"{self.bim}:{line_no} invalid bp coordinate: {pos_s}") from exc
                rows.append((chrom, vid, pos, a1, a2))
        return rows

    def _build_variant_index(self) -> Dict[str, Tuple[List[int], List[int]]]:
        """Build chromosome -> (sorted positions, BED variant indices)."""
        by_chrom: Dict[str, List[Tuple[int, int]]] = {}
        for idx, (chrom, _vid, pos, _a1, _a2) in enumerate(self.variants):
            by_chrom.setdefault(chrom, []).append((pos, idx))

        index: Dict[str, Tuple[List[int], List[int]]] = {}
        for chrom, entries in by_chrom.items():
            entries.sort(key=lambda item: (item[0], item[1]))
            index[chrom] = ([p for p, _ in entries], [i for _, i in entries])
        return index

    def _decode_variant(self, variant_index: int, block: bytes) -> VariantCall:
        chrom, vid, pos, a1, a2 = self.variants[variant_index]
        genotypes: Dict[str, Tuple[Optional[int], ...]] = {}
        phased: Dict[str, bool] = {}
        for i, sample in enumerate(self.samples):
            code = (block[i // 4] >> ((i % 4) * 2)) & 0b11
            # PLINK .bed: 00=A1/A1, 01=missing, 10=A1/A2, 11=A2/A2.
            # EasyHap uses A2 as REF index 0 and A1 as ALT index 1.
            if code == 0b00:
                gt = (1, 1)
            elif code == 0b01:
                gt = (None, None)
            elif code == 0b10:
                gt = (0, 1)
            else:
                gt = (0, 0)
            genotypes[sample] = gt
            phased[sample] = False
        info = {"SOURCE": "PLINK1_BED", "PLINK_A1": a1, "PLINK_A2": a2}
        tokens = encode_alleles(a2, (a1,), info)
        return VariantCall(chrom, pos, vid or ".", a2, (a1,), tokens, genotypes, phased, info)

    def iter_region(self, region: Region) -> Iterator[VariantCall]:
        chrom_index = self._variant_index.get(region.chrom)
        if chrom_index is None:
            return

        positions, variant_indices = chrom_index
        left = bisect.bisect_left(positions, region.start)
        right = bisect.bisect_right(positions, region.end)
        if left >= right:
            return

        selected = variant_indices[left:right]
        with open(self.bed, "rb") as fh:
            # In standard BIM files, variants for a chromosome are already in
            # physical/BED order. In that common case this loop performs one
            # initial seek and then sequential reads. If the BIM is unusually
            # unsorted, direct seeks still return variants in genomic order.
            next_file_index: Optional[int] = None
            for variant_index in selected:
                if next_file_index != variant_index:
                    fh.seek(3 + variant_index * self._bytes_per_variant)
                block = fh.read(self._bytes_per_variant)
                if len(block) != self._bytes_per_variant:
                    raise ValueError(
                        f"Unexpected end of PLINK BED while reading variant index {variant_index}"
                    )
                next_file_index = variant_index + 1
                yield self._decode_variant(variant_index, block)

