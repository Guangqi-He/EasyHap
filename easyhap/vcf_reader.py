from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field, replace
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Tuple, Union
import bisect
import gzip
import os
import re
import subprocess
import struct
import tempfile

from .external_tools import creationflags, find_bcftools, find_plink, is_windows

MISSING_ALLELE = "N"
ABSENCE_ALLELE = "ABS"


@dataclass
class Region:
    chrom: str
    start: int
    end: int

    def __post_init__(self):
        if not self.chrom or self.start < 1 or self.end < self.start:
            raise ValueError(f"Invalid region: {self.chrom}:{self.start}-{self.end}")

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
    phase_sets: Dict[str, Optional[str]] = field(default_factory=dict)

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


def read_regions(
    region: Optional[str] = None,
    region_file: Optional[str] = None,
    variant_list: Optional[Union[str, Sequence[str]]] = None,
) -> List[Region]:
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
        if variant_list:
            return [Region("CustomVariants", 1, 1)]
        raise ValueError("Please provide --region, --region-file, or --variant-list")
    from .io_utils import sanitize_filename

    if len({sanitize_filename(r.label) for r in regions}) != len(regions):
        raise ValueError("Duplicate regions would overwrite output files")
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


def _symbolic_token(
    allele: str,
    info: Optional[Dict[str, str]] = None,
    allele_index: Optional[int] = None,
) -> str:
    """Return a compact but allele-identity-preserving symbolic-variant token."""
    name = allele.strip("<>").upper() or "SV"
    svlen = _info_get_int(info or {}, "SVLEN")
    suffix = f":{svlen:+d}" if svlen not in (None, 0) else ""
    prefix = f"ALT{allele_index}:" if allele_index is not None else ""
    return f"{prefix}{name}{suffix}"


def _allele_info_int(info, key, index=1):
    value = info.get(key)
    if isinstance(value, (list, tuple)):
        values = value
    else:
        values = str(value).split(",") if value is not None else []
    value = (
        values[index - 1]
        if len(values) > 1 and index <= len(values)
        else (values[0] if values else None)
    )
    try:
        return int(value)
    except (ValueError, TypeError):
        return None


def encode_alleles(ref, alts, info=None, max_literal_len=20):
    """Readable states with unique identity for every ALT, including symbolic SVs."""
    ref = (ref or "N").upper()
    info = info or {}
    tokens = [ref if len(ref) <= max_literal_len else f"REF{len(ref)}"]
    for idx, alt in enumerate(alts, 1):
        alt = (alt or ".").upper()
        if alt == "*":
            token = ABSENCE_ALLELE
        elif alt == ".":
            token = MISSING_ALLELE
        elif "[" in alt or "]" in alt:
            token = f"ALT{idx}:BND"
        elif alt.startswith("<"):
            name = alt.strip("<>")
            length = _allele_info_int(info, "SVLEN", idx)
            token = f"ALT{idx}:{name}" + (f" {abs(length)}bp" if length else "")
        elif len(alt) == len(ref) and len(alt) <= max_literal_len:
            token = alt
        elif len(alt) < len(ref):
            token = f"DEL {len(ref) - len(alt)}bp"
        elif len(alt) > len(ref):
            token = f"INS {len(alt) - len(ref)}bp"
        else:
            token = f"ALT{idx}:SEQ{len(alt)}"
        tokens.append(token)
    counts = Counter(tokens)
    return tuple(
        t
        if i == 0 or counts[t] == 1 or t in {MISSING_ALLELE, ABSENCE_ALLELE}
        else f"ALT{i}:{t}"
        for i, t in enumerate(tokens)
    )


def deletion_interval(call, allele_index):
    """Inclusive removed reference coordinates; complex replacements are not masked."""
    if not 0 < allele_index <= len(call.alts):
        return None
    alt = call.alts[allele_index - 1].upper()
    if alt == "*":
        return None  # spanning-deletion state, not a new causal event
    if alt.startswith("<"):
        name = alt.strip("<>").split(":")[0]
        if name not in {"DEL", "CN0", "ABS"}:
            return None
        end = _allele_info_int(call.info, "END", allele_index)
        length = _allele_info_int(call.info, "SVLEN", allele_index)
        svlens = call.info.get("SVLEN", "")
        n_lengths = (
            len(svlens)
            if isinstance(svlens, (tuple, list))
            else len(str(svlens).split(","))
        )
        # A site-level END can describe the longest ALT. Use a supplied
        # per-ALT length for multiallelic deletions rather than masking that span.
        if length and (end is None or n_lengths > 1):
            end = call.pos + abs(length)
        return (call.pos + 1, end) if end is not None and end > call.pos else None
    if "[" in alt or "]" in alt or alt == ".":
        return None
    ref = call.ref.upper()
    prefix = 0
    while prefix < min(len(ref), len(alt)) and ref[prefix] == alt[prefix]:
        prefix += 1
    ref_tail, alt_tail = ref[prefix:], alt[prefix:]
    suffix = 0
    while (
        suffix < min(len(ref_tail), len(alt_tail))
        and ref_tail[-suffix - 1] == alt_tail[-suffix - 1]
    ):
        suffix += 1
    r = ref_tail[:-suffix] if suffix else ref_tail
    a = alt_tail[:-suffix] if suffix else alt_tail
    if r and not a:
        return call.pos + prefix, call.pos + prefix + len(r) - 1
    return None


def deletion_allele_indices(call):
    return tuple(
        i
        for i in range(1, len(call.alts) + 1)
        if deletion_interval(call, i) is not None
    )


def deletion_end(call):
    return max(
        [call.pos]
        + [deletion_interval(call, i)[1] for i in deletion_allele_indices(call)]
    )


def token_for_gt(allele_tokens: Sequence[str], allele_index: Optional[int]) -> str:
    if allele_index is None or allele_index < 0:
        return MISSING_ALLELE
    if allele_index >= len(allele_tokens):
        return MISSING_ALLELE
    return allele_tokens[allele_index]


def variant_selected(chrom, pos, vid, allowed):
    return allowed is None or bool(
        {str(pos), f"{chrom}:{pos}", f"{chrom}_{pos}"} & allowed
        or set((vid or "").split(";")) & allowed
    )


def _indexed_variant_sites(allowed, known_contigs=None):
    """Extract indexable CHROM:POS/CHROM_POS targets from a variant list.

    Return None when the selection includes bare POS values or IDs.  Neither
    can be located with a chromosome-coordinate VCF index, so scanning is
    required to avoid silently dropping matching variants.  An empty set is
    distinguishable from None: it selects no variants.

    A token such as rs_123 could also be a variant ID.  Accept a token
    as an indexable coordinate only when its chromosome is confirmed by
    the indexed file, otherwise retain the original ID scan behavior.
    """
    if allowed is None:
        return None
    if not allowed:
        return {}
    if known_contigs is None:
        return None
    known_contigs = set(known_contigs)
    if not known_contigs:
        return None
    positions = {}
    for item in allowed:
        item = str(item)
        # Allow chromosome names containing underscores/colons. Split at the
        # final coordinate separator so Chr_04_123 or scaffold:part:123 works.
        match = re.fullmatch(r"(.+)[:_](\d+)", item)
        if match is None:
            return None
        chrom, pos_string = match.groups()
        if chrom not in known_contigs:
            return None
        pos = int(pos_string)
        if pos < 1:
            return None
        positions.setdefault(chrom, set()).add(pos)
    return positions


def _intersect_variant_sites(sites, region):
    """Clip indexable sites to a requested region without changing IDs."""
    if sites is None or region is None:
        return sites
    if region.chrom not in sites:
        return {}
    return {
        region.chrom: {
            pos for pos in sites[region.chrom]
            if region.start <= pos <= region.end
        }
    }


def _variant_site_windows(sites, max_gap=1024):
    """Coalesce nearby indexed lookups to avoid a fetch for every SNP.

    Each yielded interval is disjoint. Readers check the original POS after
    fetching, since tabix/HTSlib may also return SVs overlapping a window.
    """
    for chrom in sorted(sites):
        positions = sorted(sites[chrom])
        if not positions:
            continue
        start = end = positions[0]
        for pos in positions[1:]:
            if pos - end <= max_gap:
                end = pos
            else:
                yield Region(chrom, start, end)
                start = end = pos
        yield Region(chrom, start, end)


def _deletion_candidate(ref, alts):
    return any(
        a.startswith(("<DEL", "<CN0", "<ABS"))
        or (
            a not in {".", "*", ""}
            and not a.startswith("<")
            and "[" not in a
            and "]" not in a
            and len(a) < len(ref)
        )
        for a in alts
    )


class VCFReader:
    """VCF/BCF reader with platform-aware backend selection.

    Backend policy
    --------------
    * Windows + BGZF-compressed VCF/BCF: prefer bundled bcftools only.
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
        if len(set(self.samples)) != len(self.samples) or any(
            not x.strip() for x in self.samples
        ):
            raise ValueError("VCF sample IDs must be nonempty and unique")

    @staticmethod
    def _is_windows() -> bool:
        return is_windows()

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

    def _find_bcftools(self) -> Optional[str]:
        """Windows always uses the PyInstaller-bundled bcftools.exe."""
        return find_bcftools()

    @staticmethod
    def _subprocess_creationflags() -> int:
        return creationflags()

    def _run_bcftools(
        self, args: Sequence[str], *, context: str
    ) -> subprocess.CompletedProcess:
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
            self._run_bcftools(
                ["index", "-f", "-c", self.path], context=f"indexing {self.path!r}"
            )
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

    def _indexed_contigs_bcftools(self):
        """Get reference names from the CSI/TBI index (lazily, once per VCF)."""
        if hasattr(self, "_cached_index_contigs"):
            return self._cached_index_contigs
        try:
            output = self._run_bcftools(
                ["index", "-s", self.path], context="reading index contigs"
            ).stdout
            contigs = {line.split("\t", 1)[0] for line in output.splitlines() if line.strip()}
        except (OSError, RuntimeError):
            contigs = set()
        # Unavailable index metadata must trigger the safe scanning fallback.
        self._cached_index_contigs = contigs or None
        return self._cached_index_contigs

    def _read_samples_bcftools(self) -> List[str]:
        proc = self._run_bcftools(
            ["query", "-l", self.path], context=f"reading samples from {self.path!r}"
        )
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
        self.backend_detail = (
            f"bcftools={self._bcftools}; index={index_note}{created_note}"
        )
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
                raise RuntimeError(
                    "BCF is binary and cannot be read by the plain-text backend."
                )
            self.samples = self._read_samples_plain()
            self.backend = "plain"
            self.backend_detail = "plain (explicit)"
            return

        if self.prefer == "bcftools":
            if not self._init_bcftools(create_index=True):
                raise RuntimeError(
                    "The bcftools backend was requested but bcftools could not be found. "
                    "For Windows, bundle vendor/bcftools/bcftools.exe with EasyHap; "
                    "for Linux, install bcftools or set EASYHAP_BCFTOOLS."
                )
            return

        # Whenever bcftools is available on the system or bundled with the exe, prefer it for fast C-level extraction and filtering
        if self.prefer == "auto" and self._is_gzip_file():
            with open(self.path, "rb") as fh:
                header = fh.read(64)
            if b"BC\x02\x00" not in header:
                self.samples = self._read_samples_plain()
                self.backend_detail = "plain (ordinary gzip VCF)"
                return

        if self.prefer == "auto" and self._find_bcftools():
            if self._is_compressed_text_vcf() or self._is_binary_bcf():
                if self._init_bcftools(create_index=True):
                    return
            elif self._is_uncompressed_text_vcf() and self._has_region_index():
                if self._init_bcftools(create_index=False):
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

                if self._has_region_index():
                    if self._try_cyvcf2() or self._try_pysam():
                        return

                raise RuntimeError(
                    "EasyHap is running on Windows with a compressed VCF/BCF, but bcftools could not be found and "
                    "no indexed cyvcf2/pysam backend could be opened. Refusing to sequentially scan a large file. "
                    "Bundle vendor/bcftools/bcftools.exe for Windows, or install bcftools on Linux."
                )

            if self._is_uncompressed_text_vcf():
                if self._file_size() <= self.WINDOWS_PLAIN_MAX_BYTES:
                    self.samples = self._read_samples_plain()
                    self.backend = "plain"
                    self.backend_detail = f"plain (Windows uncompressed VCF <= {self.WINDOWS_PLAIN_MAX_BYTES // (1024 * 1024)} MiB)"
                    return
                raise RuntimeError(
                    f"The uncompressed VCF is {self._file_size() / (1024 * 1024):.1f} MiB, larger than EasyHap's "
                    f"{self.WINDOWS_PLAIN_MAX_BYTES // (1024 * 1024)} MiB plain-reader threshold on Windows. "
                    "To avoid repeatedly scanning a large VCF, bgzip-compress it to .vcf.gz and rerun EasyHap; "
                    "EasyHap will create a CSI index automatically and use bcftools regional access."
                )

        # Non-Windows auto behavior:
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

    def iter_region(
        self,
        region=None,
        variant_list=None,
        max_missing=None,
        min_maf=None,
        samples=None,
    ):
        from .core import load_variant_list, filter_variants_by_qc

        allowed = load_variant_list(variant_list) if variant_list else None
        selected = list(self.samples) if samples is None else list(samples)
        if len(set(selected)) != len(selected) or not set(selected).issubset(
            self.samples
        ):
            raise ValueError(
                "Requested sample IDs are duplicate or absent from the VCF"
            )
        if self.backend == "bcftools":
            iterator = self._iter_bcftools(region, samples=selected, allowed=allowed)
        elif self.backend == "cyvcf2":
            iterator = self._iter_cyvcf2(region, allowed=allowed)
        elif self.backend == "pysam":
            iterator = self._iter_pysam(region, allowed=allowed)
        else:
            iterator = self._iter_plain(region, allowed=allowed)
        for call in iterator:
            call = replace(
                call,
                genotypes={s: call.genotypes[s] for s in selected},
                phased={s: call.phased[s] for s in selected},
                phase_sets={s: call.phase_sets.get(s) for s in selected},
            )
            if call.alts and all(a in {"<NON_REF>", "<*>"} for a in call.alts):
                continue  # gVCF reference-confidence blocks are not variant markers
            if not filter_variants_by_qc([call], selected, max_missing, min_maf)[0]:
                continue
            yield call

    def iter_all_variants(
        self, variant_list=None, max_missing=None, min_maf=None, samples=None
    ):
        yield from self.iter_region(None, variant_list, max_missing, min_maf, samples)

    def context_query_region(self, region, flank_bp=0):
        """Bound candidate POS to the requested interval plus explicit flanks."""
        if not isinstance(flank_bp, int) or isinstance(flank_bp, bool) or flank_bp < 0:
            raise ValueError("overlap_flank_bp must be a nonnegative integer (bp)")
        return Region(region.chrom, max(1, region.start - flank_bp), region.end + flank_bp)

    def iter_context_deletions(self, region, samples=None, flank_bp=0):
        # No header-dependent extension, upstream chromosome scan or global cache.
        query_region = self.context_query_region(region, flank_bp)
        for call in self.iter_region(query_region, samples=samples):
            if any(
                start <= region.end and end >= region.start
                for allele in deletion_allele_indices(call)
                for start, end in [deletion_interval(call, allele)]
            ):
                yield call

    def _parse_vcf_record_line(self, line, samples=None):
        if not line or line.startswith("#"):
            return None
        parts = line.rstrip("\n\r").split("\t")
        target = self.samples if samples is None else list(samples)
        if len(parts) != (9 + len(target) if target else 8) and not (
            not target and len(parts) == 9
        ):
            raise ValueError(f"VCF column/sample count mismatch: {parts[0:3]}")
        chrom, pos_s, vid, ref, alt_s, _, _, info_s = parts[:8]
        pos = int(pos_s)
        if pos < 1:
            raise ValueError(f"Invalid VCF position {chrom}:{pos}")
        fmt = parts[8].split(":") if len(parts) > 8 else []
        if target and "GT" not in fmt:
            raise ValueError(f"Missing FORMAT/GT at {chrom}:{pos}")
        info = _parse_info(info_s)
        alts = tuple(alt_s.split(",")) if alt_s != "." else ()
        genotypes, phased, phase_sets = {}, {}, {}
        for sample, field in zip(target, parts[9:]):
            fields = dict(zip(fmt, field.split(":")))
            gt_text = fields.get("GT", ".")
            separators = re.findall(r"[|/]", gt_text)
            phased[sample] = bool(separators) and all(x == "|" for x in separators)
            gt = tuple(None if a == "." else int(a) for a in re.split(r"[|/]", gt_text))
            if any(a is not None and not 0 <= a <= len(alts) for a in gt):
                raise ValueError(
                    f"GT allele index outside REF/ALT at {sample}@{chrom}:{pos}: {gt_text}"
                )
            genotypes[sample] = tuple(
                None
                if a is not None and a > 0 and alts[a - 1] in {"<NON_REF>", "<*>"}
                else a
                for a in gt
            )
            ps = fields.get("PS", ".")
            phase_sets[sample] = None if ps in {".", ""} else ps
        return VariantCall(
            chrom,
            pos,
            vid,
            ref,
            alts,
            encode_alleles(ref, alts, info),
            genotypes,
            phased,
            info,
            phase_sets,
        )

    def _line_selected(self, line, region=None, allowed=None, deletions_only=False):
        parts = line.split("\t", 8)
        if len(parts) < 8:
            raise ValueError("Malformed VCF data row")
        chrom, pos, vid, ref, alt = parts[:5]
        pos = int(pos)
        return (
            (
                region is None
                or chrom == region.chrom
                and region.start <= pos <= region.end
            )
            and variant_selected(chrom, pos, vid, allowed)
            and (not deletions_only or _deletion_candidate(ref, alt.split(",")))
        )

    def _iter_bcftools(
        self, region=None, samples=None, allowed=None, deletions_only=False
    ):
        target = self.samples if samples is None else list(samples)
        args = ["view", "-H"]
        sample_path = None
        regions_path = None
        proc = None
        try:
            if list(target) != self.samples:
                with tempfile.NamedTemporaryFile(
                    mode="w", encoding="utf-8", delete=False
                ) as fh:
                    fh.write("\n".join(target) + "\n")
                    sample_path = fh.name
                args += ["-S", sample_path]

            # bcftools -R uses indexed random access.  A three-column TSV
            # uses inclusive 1-based coordinates, not BED's zero-based ones.
            # Only use it when *every* requested token is indexable; otherwise
            # an ID-only or mixed list must still be scanned by VCF ID.
            sites = (
                _intersect_variant_sites(
                    _indexed_variant_sites(allowed, self._indexed_contigs_bcftools()), region
                )
                if allowed is not None and self._has_region_index()
                else None
            )
            if sites is not None:
                if not any(sites.values()):
                    return
                with tempfile.NamedTemporaryFile(
                    mode="w", encoding="utf-8", suffix=".tsv", delete=False
                ) as fh:
                    regions_path = fh.name
                    for chrom in sorted(sites):
                        for pos in sorted(sites[chrom]):
                            fh.write(f"{chrom}\t{pos}\t{pos}\n")
                args += ["-R", regions_path, "--regions-overlap", "0"]
            elif region is not None:
                args += (
                    ["-r", region.vcf_label, "--regions-overlap", "0"]
                    if self._has_region_index()
                    else ["-t", region.vcf_label, "--targets-overlap", "0"]
                )

            with tempfile.TemporaryFile(mode="w+t", encoding="utf-8") as error_file:
                proc = subprocess.Popen(
                    [self._bcftools, *args, self.path],
                    stdout=subprocess.PIPE,
                    stderr=error_file,
                    text=True,
                    encoding="utf-8",
                    errors="replace",
                    creationflags=self._subprocess_creationflags(),
                )
                for line in proc.stdout:
                    # Keep the strict original Python-level match, e.g. for
                    # multiallelic records sharing one indexed coordinate.
                    if self._line_selected(line, region, allowed, deletions_only):
                        yield self._parse_vcf_record_line(line, target)
                code = proc.wait()
                if code:
                    error_file.seek(0)
                    raise RuntimeError(
                        f"bcftools extraction failed (exit {code}): {error_file.read().strip()}"
                    )
        finally:
            if proc is not None:
                if proc.stdout:
                    proc.stdout.close()
                if proc.poll() is None:
                    proc.terminate()
                    try:
                        proc.wait(timeout=10)
                    except subprocess.TimeoutExpired:
                        proc.kill()
                        proc.wait()
            for path in (sample_path, regions_path):
                if path:
                    os.unlink(path)

    def _iter_cyvcf2(self, region=None, allowed=None, deletions_only=False):
        from cyvcf2 import VCF

        with VCF(self.path) as fresh:
            sites = (
                _intersect_variant_sites(
                    _indexed_variant_sites(allowed, getattr(fresh, "seqnames", None)), region
                )
                if allowed is not None and self._has_region_index()
                else None
            )
            if sites is not None:
                # Reuse one opened VCF and seek across disjoint merged windows.
                query_windows = _variant_site_windows(sites)
            elif region is not None and self._has_region_index():
                query_windows = (region,)
            else:
                query_windows = None

            if query_windows is None:
                iterator = iter(fresh)
            else:
                def fetch():
                    for window in query_windows:
                        for rec in fresh(window.vcf_label):
                            # HTSlib may return spanning records whose POS is
                            # before the queried interval. Exclude them here.
                            if rec.CHROM == window.chrom and window.start <= rec.POS <= window.end:
                                yield rec
                iterator = fetch()

            for rec in iterator:
                if region is not None and (
                    rec.CHROM != region.chrom
                    or not region.start <= rec.POS <= region.end
                ):
                    continue
                if not variant_selected(rec.CHROM, rec.POS, rec.ID, allowed):
                    continue
                if deletions_only and not _deletion_candidate(rec.REF, rec.ALT or []):
                    continue
                yield self._parse_vcf_record_line(str(rec))

    def _iter_pysam(self, region=None, allowed=None, deletions_only=False):
        import pysam

        with pysam.VariantFile(self.path) as fresh:
            sites = (
                _intersect_variant_sites(
                    _indexed_variant_sites(allowed, fresh.header.contigs.keys()), region
                )
                if allowed is not None and self._has_region_index()
                else None
            )
            if sites is not None:
                query_windows = _variant_site_windows(sites)
            elif region is not None and self._has_region_index():
                query_windows = (region,)
            else:
                query_windows = None

            if query_windows is None:
                iterator = iter(fresh)
            else:
                def fetch():
                    for window in query_windows:
                        for rec in fresh.fetch(window.chrom, window.start - 1, window.end):
                            if rec.chrom == window.chrom and window.start <= rec.pos <= window.end:
                                yield rec
                iterator = fetch()

            for rec in iterator:
                if region is not None and (
                    rec.chrom != region.chrom
                    or not region.start <= rec.pos <= region.end
                ):
                    continue
                # Avoid the expensive conversion/GT parsing of unrelated
                # sites; the previous pysam implementation skipped this check.
                if not variant_selected(rec.chrom, rec.pos, rec.id, allowed):
                    continue
                if deletions_only and not _deletion_candidate(rec.ref, rec.alts or ()):
                    continue
                call = self._parse_vcf_record_line(str(rec))
                # pysam exposes END as stop instead of rec.info['END'].
                if rec.stop > rec.pos + len(rec.ref) - 1:
                    call.info.setdefault("END", str(rec.stop))
                yield call

    def _iter_plain(self, region=None, allowed=None, deletions_only=False):
        with self._open_text() as fh:
            for line in fh:
                if not line or line.startswith("#"):
                    continue
                if self._line_selected(line, region, allowed, deletions_only):
                    yield self._parse_vcf_record_line(line)


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
        if not os.path.exists(self.bim) and os.path.exists(self.prefix + ".bim.gz"):
            self.bim = self.prefix + ".bim.gz"
        self.fam = self.prefix + ".fam"
        if not os.path.exists(self.fam) and os.path.exists(self.prefix + ".fam.gz"):
            self.fam = self.prefix + ".fam.gz"
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
        if actual != expected:
            raise ValueError(
                f"PLINK BED appears truncated: expected {expected} bytes, found {actual}"
            )
        with open(self.bed, "rb") as fh:
            magic = fh.read(3)
        if magic != bytes((0x6C, 0x1B, 0x01)):
            raise ValueError(
                "Unsupported PLINK BED: expected current SNP-major magic bytes 6c 1b 01"
            )

    @staticmethod
    def _normalize_prefix(prefix: str) -> str:
        """Accept either a PLINK prefix or one member of a BED/BIM/FAM trio."""
        value = os.fspath(prefix).strip()
        lower = value.lower()
        for suffix in (".bed", ".bim.gz", ".bim", ".fam.gz", ".fam"):
            if lower.endswith(suffix):
                return value[: -len(suffix)]
        return value

    def _read_fam(self) -> List[str]:
        samples: List[str] = []
        seen = set()
        opener = gzip.open if self.fam.endswith(".gz") else open
        with opener(self.fam, "rt", encoding="utf-8") as fh:
            for line_no, raw in enumerate(fh, 1):
                if not raw.strip():
                    continue
                parts = raw.split()
                if len(parts) < 2:
                    raise ValueError(
                        f"{self.fam}:{line_no} requires at least FID and IID"
                    )
                iid = parts[1]
                sid = iid
                if sid in seen:
                    raise ValueError(
                        f"Duplicate PLINK sample identifier after FID/IID resolution: {sid}"
                    )
                seen.add(sid)
                samples.append(sid)
        return samples

    def _read_bim(self) -> List[Tuple[str, str, int, str, str]]:
        rows: List[Tuple[str, str, int, str, str]] = []
        opener = gzip.open if self.bim.endswith(".gz") else open
        with opener(self.bim, "rt", encoding="utf-8") as fh:
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
                    raise ValueError(
                        f"{self.bim}:{line_no} invalid bp coordinate: {pos_s}"
                    ) from exc
                if pos < 1:
                    raise ValueError(f"Invalid BIM coordinate {chrom}:{pos}")
                rows.append((chrom, vid, pos, a1, a2))
        if not rows or any(r[2] < 1 for r in rows):
            raise ValueError(
                "PLINK variant metadata is empty or has invalid coordinates"
            )
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
        return VariantCall(
            chrom, pos, vid or ".", a2, (a1,), tokens, genotypes, phased, info
        )

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

    def iter_variants(
        self, allowed_variants: set[str], region: Optional[Region] = None
    ) -> Iterator[VariantCall]:
        selected = []
        for idx, (chrom, vid, pos, a1, a2) in enumerate(self.variants):
            if (
                vid in allowed_variants
                or str(pos) in allowed_variants
                or f"{chrom}:{pos}" in allowed_variants
                or f"{chrom}_{pos}" in allowed_variants
            ):
                if region is None or (
                    chrom == region.chrom and region.start <= pos <= region.end
                ):
                    selected.append(idx)
        if not selected:
            return
        with open(self.bed, "rb") as fh:
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

    def iter_all_variants(self) -> Iterator[VariantCall]:
        with open(self.bed, "rb") as fh:
            fh.seek(3)
            for variant_index in range(len(self.variants)):
                block = fh.read(self._bytes_per_variant)
                if len(block) != self._bytes_per_variant:
                    raise ValueError(
                        f"Unexpected end of PLINK BED while reading variant index {variant_index}"
                    )
                yield self._decode_variant(variant_index, block)


class Plink2PgenReader:
    """PLINK 2 binary (PGEN/PVAR/PSAM) reader for phase-independent analysis.

    Supports .pgen, .pvar / .pvar.gz, and .psam / .psam.gz.
    Also handles hybrid PLINK 1 / PLINK 2 filesets (.bim, .fam).
    Supports direct gz reading for pvar/psam without manual decompression.
    """

    def __init__(self, prefix: str) -> None:
        self.prefix = self._normalize_prefix(prefix)
        self.pgen = self._find_file([self.prefix + ".pgen", self.prefix + ".bed"])
        self.pvar = self._find_file(
            [
                self.prefix + ".pvar",
                self.prefix + ".pvar.gz",
                self.prefix + ".bim",
                self.prefix + ".bim.gz",
            ]
        )
        self.psam = self._find_file(
            [
                self.prefix + ".psam",
                self.prefix + ".psam.gz",
                self.prefix + ".fam",
                self.prefix + ".fam.gz",
            ]
        )
        if not self.pgen:
            raise FileNotFoundError(
                f"Missing PLINK 2 .pgen/.bed file for prefix: {self.prefix}"
            )
        if not self.pvar:
            raise FileNotFoundError(
                f"Missing PLINK 2 .pvar/.bim file for prefix: {self.prefix}"
            )
        if not self.psam:
            raise FileNotFoundError(
                f"Missing PLINK 2 .psam/.fam file for prefix: {self.prefix}"
            )

        self.samples = self._read_psam()
        self.variants = self._read_pvar()
        self._variant_index = self._build_variant_index()
        self.backend = "plink2-pgen"

        import numpy as np

        try:
            import pgenlib
        except ImportError as exc:
            raise RuntimeError(
                "PGEN requires pgenlib; install EasyHap dependencies or pip install 'pgenlib>=0.95.1'"
            ) from exc
        self._offsets = np.asarray(
            [0]
            + list(
                __import__("itertools").accumulate(
                    max(2, 1 + len(v[4])) for v in self.variants
                )
            ),
            dtype=np.uintp,
        )
        self._pgenlib = pgenlib
        with self._open_pgen() as reader:
            if reader.get_raw_sample_ct() != len(
                self.samples
            ) or reader.get_variant_ct() != len(self.variants):
                raise ValueError("PGEN dimensions disagree with PSAM/PVAR metadata")

    def _open_pgen(self):
        return self._pgenlib.PgenReader(
            os.fsencode(self.pgen),
            raw_sample_ct=len(self.samples),
            variant_ct=len(self.variants),
            allele_idx_offsets=self._offsets,
        )

    @staticmethod
    def _normalize_prefix(prefix: str) -> str:
        value = os.fspath(prefix).strip()
        lower = value.lower()
        for suffix in (
            ".pgen",
            ".pvar.gz",
            ".pvar",
            ".psam.gz",
            ".psam",
            ".bed",
            ".bim.gz",
            ".bim",
            ".fam.gz",
            ".fam",
        ):
            if lower.endswith(suffix):
                return value[: -len(suffix)]
        return value

    @staticmethod
    def _find_file(candidates: List[str]) -> Optional[str]:
        for c in candidates:
            if os.path.exists(c):
                return c
        return None

    @staticmethod
    def _open_text(path: str):
        if path.endswith(".gz"):
            return gzip.open(path, "rt", encoding="utf-8", errors="replace")
        return open(path, "rt", encoding="utf-8", errors="replace")

    def _read_psam(self) -> List[str]:
        samples: List[str] = []
        seen = set()
        iid_idx = 1
        fid_idx = 0
        has_header = False
        with self._open_text(self.psam) as fh:
            for line_no, raw in enumerate(fh, 1):
                line = raw.strip()
                if not line:
                    continue
                if line.startswith("#"):
                    cols = line.lstrip("#").split()
                    cols_upper = [c.upper() for c in cols]
                    if "IID" in cols_upper:
                        iid_idx = cols_upper.index("IID")
                        fid_idx = cols_upper.index("FID") if "FID" in cols_upper else -1
                        has_header = True
                    continue
                parts = line.split()
                if has_header:
                    if iid_idx >= len(parts):
                        raise ValueError(f"Malformed PSAM row {self.psam}:{line_no}")
                    iid = parts[iid_idx]
                    fid = (
                        parts[fid_idx] if fid_idx >= 0 and fid_idx < len(parts) else ""
                    )
                else:
                    if len(parts) >= 2:
                        fid, iid = parts[0], parts[1]
                    else:
                        fid, iid = "", parts[0]
                sid = iid
                if sid in seen:
                    raise ValueError(
                        f"Duplicate PLINK IID in {self.psam}:{line_no}: {iid}; provide unique IIDs"
                    )
                seen.add(sid)
                samples.append(sid)
        if not samples:
            raise ValueError(f"No samples found in {self.psam}")
        return samples

    def _read_pvar(self) -> List[Tuple[str, str, int, str, Tuple[str, ...]]]:
        rows: List[Tuple[str, str, int, str, Tuple[str, ...]]] = []
        is_bim = self.pvar.endswith((".bim", ".bim.gz"))
        with self._open_text(self.pvar) as fh:
            for line_no, raw in enumerate(fh, 1):
                line = raw.strip()
                if not line or line.startswith("##"):
                    continue
                if line.startswith("#"):
                    continue
                parts = line.split()
                if is_bim:
                    if len(parts) < 6:
                        raise ValueError(
                            f"Malformed PVAR/BIM row {self.pvar}:{line_no}"
                        )
                    chrom, vid, _cm, pos_s, a1, a2 = parts[:6]
                    try:
                        pos = int(pos_s)
                    except ValueError:
                        raise ValueError(
                            f"Malformed PVAR/BIM row {self.pvar}:{line_no}"
                        )
                    rows.append((chrom, vid, pos, a2, (a1,)))
                else:
                    if len(parts) < 5:
                        raise ValueError(
                            f"Malformed PVAR/BIM row {self.pvar}:{line_no}"
                        )
                    chrom, pos_s, vid, ref, alt_s = parts[:5]
                    try:
                        pos = int(pos_s)
                    except ValueError:
                        raise ValueError(
                            f"Malformed PVAR/BIM row {self.pvar}:{line_no}"
                        )
                    alts = tuple(alt_s.split(",")) if alt_s and alt_s != "." else ()
                    rows.append((chrom, vid, pos, ref, alts))
        if not rows or any(r[2] < 1 for r in rows):
            raise ValueError(
                "PLINK variant metadata is empty or has invalid coordinates"
            )
        return rows

    def _build_variant_index(self) -> Dict[str, Tuple[List[int], List[int]]]:
        by_chrom: Dict[str, List[Tuple[int, int]]] = {}
        for idx, (chrom, _vid, pos, _ref, _alts) in enumerate(self.variants):
            by_chrom.setdefault(chrom, []).append((pos, idx))
        index: Dict[str, Tuple[List[int], List[int]]] = {}
        for chrom, entries in by_chrom.items():
            entries.sort(key=lambda item: (item[0], item[1]))
            index[chrom] = ([p for p, _ in entries], [i for _, i in entries])
        return index

    def _read_indices(self, indices):
        import numpy as np

        buffer = np.empty(2 * len(self.samples), dtype=np.int32)
        with self._open_pgen() as reader:
            for index in indices:
                reader.read_alleles(index, buffer)
                chrom, vid, pos, ref, alts = self.variants[index]
                gts = {
                    s: tuple(
                        None if int(a) < 0 else int(a)
                        for a in buffer[2 * i : 2 * i + 2]
                    )
                    for i, s in enumerate(self.samples)
                }
                if any(
                    a is not None and a > len(alts) for gt in gts.values() for a in gt
                ):
                    raise ValueError(
                        f"PGEN alleles disagree with PVAR at {chrom}:{pos}"
                    )
                info = {
                    "SOURCE": "PLINK2_PGEN",
                    "CALL_TYPE": "hardcall",
                    "PLOIDY": "diploid",
                }
                yield VariantCall(
                    chrom,
                    pos,
                    vid or ".",
                    ref,
                    alts,
                    encode_alleles(ref, alts, info),
                    gts,
                    {s: False for s in self.samples},
                    info,
                )

    def iter_region(self, region):
        positions, indices = self._variant_index.get(region.chrom, ([], []))
        yield from self._read_indices(
            indices[
                bisect.bisect_left(positions, region.start) : bisect.bisect_right(
                    positions, region.end
                )
            ]
        )

    def iter_variants(self, allowed_variants, region=None):
        indices = [
            i
            for i, (chrom, vid, pos, ref, alts) in enumerate(self.variants)
            if variant_selected(chrom, pos, vid, allowed_variants)
            and (
                region is None
                or chrom == region.chrom
                and region.start <= pos <= region.end
            )
        ]
        yield from self._read_indices(indices)

    def iter_all_variants(self):
        yield from self._read_indices(range(len(self.variants)))


class ExternalPlinkReader:
    """Read BED/PGEN using an installed PLINK (Linux) or bundled PLINK2 (Windows).

    Only each requested region is exported as a temporary VCF.  This avoids a
    whole-genome export for ordinary regional analyses and requires no pgenlib.
    Sample names remain IIDs, matching the original EasyHap PLINK readers.
    """

    def __init__(self, prefix: str, kind: str, executable: str, flavor: str):
        self.kind = kind
        self.executable = executable
        self.flavor = flavor
        self.prefix = (
            PlinkBedReader._normalize_prefix(prefix)
            if kind == "bfile"
            else Plink2PgenReader._normalize_prefix(prefix)
        )
        self.backend = f"{flavor}-{kind}-vcf"
        self.backend_detail = f"PLINK={executable}; input={self.prefix}"

        if kind == "bfile":
            # The existing reader validates BED magic/size and samples, and
            # indexes BIM coordinates without decoding any genotypes.
            meta = PlinkBedReader(self.prefix)
            # Preserve input-file fingerprints in EasyHap analysis metadata.
            self.bed, self.bim, self.fam = meta.bed, meta.bim, meta.fam
            self.samples = meta.samples
            self.variants = meta.variants
            self._variant_index = meta._variant_index
            self.input_args = ["--bfile", self.prefix]
        else:
            if flavor != "plink2":
                raise RuntimeError("PLINK 2 is required for PGEN input")
            self.pgen = Plink2PgenReader._find_file([self.prefix + ".pgen"])
            self.pvar = Plink2PgenReader._find_file(
                [self.prefix + ".pvar", self.prefix + ".pvar.gz",
                 self.prefix + ".bim", self.prefix + ".bim.gz"]
            )
            self.psam = Plink2PgenReader._find_file(
                [self.prefix + ".psam", self.prefix + ".psam.gz",
                 self.prefix + ".fam", self.prefix + ".fam.gz"]
            )
            for extension, path in ((".pgen", self.pgen), (".pvar", self.pvar), (".psam", self.psam)):
                if path is None:
                    raise FileNotFoundError(f"Missing PLINK 2 {extension} file for prefix: {self.prefix}")
            # Reuse original PVAR/PSAM metadata readers without initializing
            # pgenlib (which is unavailable in the Windows build).
            self.samples = Plink2PgenReader._read_psam(self)
            self.variants = Plink2PgenReader._read_pvar(self)
            self._variant_index = Plink2PgenReader._build_variant_index(self)
            is_bim = self.pvar.endswith((".bim", ".bim.gz"))
            is_fam = self.psam.endswith((".fam", ".fam.gz"))
            if is_bim and is_fam and not (
                self.pvar.endswith(".gz") or self.psam.endswith(".gz")
            ):
                self.input_args = ["--bpfile", self.prefix]
            elif not (is_bim or is_fam) and not (
                self.pvar.endswith(".gz") or self.psam.endswith(".gz")
            ):
                self.input_args = ["--pfile", self.prefix]
            else:
                self.input_args = [
                    "--pgen", self.pgen,
                    "--bim" if is_bim else "--pvar", self.pvar,
                    "--fam" if is_fam else "--psam", self.psam,
                ]

    @staticmethod
    def _open_text(path: str):
        return Plink2PgenReader._open_text(path)

    def _has_variants(self, region: Optional[Region]) -> bool:
        if region is None:
            return True
        positions = self._variant_index.get(region.chrom, ([], []))[0]
        return bisect.bisect_left(positions, region.start) < bisect.bisect_right(
            positions, region.end
        )

    def _export_vcf(
        self, region: Optional[Region], allowed_variants=None
    ) -> Iterator[VariantCall]:
        if not self._has_variants(region):
            return
        # For --variant-list without a region, use PLINK 2's coordinate-BED
        # filter to avoid converting a whole chromosome/genome into VCF.
        selected_positions = None
        if region is None and allowed_variants and self.flavor == "plink2":
            selected_positions = {
                (chrom, pos)
                for chrom, vid, pos, *_rest in self.variants
                if variant_selected(chrom, pos, vid, allowed_variants)
            }
            if not selected_positions:
                return
        with tempfile.TemporaryDirectory(prefix="EasyHap_plink_") as tmpdir:
            output_prefix = os.path.join(tmpdir, "subset")
            cmd = [self.executable, *self.input_args, "--allow-extra-chr"]
            if selected_positions is not None:
                regions_path = os.path.join(tmpdir, "selected.bed")
                with open(regions_path, "w", encoding="utf-8") as fh:
                    for chrom, pos in sorted(selected_positions):
                        fh.write(f"{chrom}\t{pos - 1}\t{pos}\n")
                cmd += ["--extract", "bed0", regions_path]
            if region is not None:
                cmd += [
                    "--chr", region.chrom,
                    "--from-bp", str(region.start),
                    "--to-bp", str(region.end),
                ]
            if self.flavor == "plink2":
                cmd += ["--export", "vcf", "id-paste=iid"]
            else:
                cmd += ["--keep-allele-order", "--recode", "vcf-iid"]
            cmd += ["--out", output_prefix]
            try:
                completed = subprocess.run(
                    cmd,
                    cwd=tmpdir,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                    encoding="utf-8",
                    errors="replace",
                    creationflags=creationflags(),
                )
            except OSError as exc:
                raise RuntimeError(f"Cannot start PLINK at {self.executable}: {exc}") from exc
            vcf_path = output_prefix + ".vcf"
            if completed.returncode:
                details = (completed.stderr or completed.stdout or "").strip()
                raise RuntimeError(
                    f"PLINK export failed (exit {completed.returncode}) for {self.prefix}"
                    + (f" / {region.vcf_label}" if region else "")
                    + f"\nCommand: {subprocess.list2cmdline(cmd)}\n{details}"
                )
            if not os.path.isfile(vcf_path):
                raise RuntimeError(f"PLINK reported success without creating {vcf_path}")
            reader = VCFReader(vcf_path, prefer="plain")
            if set(reader.samples) != set(self.samples):
                raise ValueError(
                    "PLINK VCF export changed sample IIDs; check duplicate FID/IID IDs "
                    f"or PLINK version. Expected {len(self.samples)} samples, "
                    f"exported {len(reader.samples)}."
                )
            for call in reader.iter_region(region, samples=self.samples):
                # PLINK files are intentionally analyzed as unphased genotypes,
                # consistent with the original BED/PGEN reader contract.
                info = dict(call.info)
                info["SOURCE"] = "PLINK1_BED" if self.kind == "bfile" else "PLINK2_PGEN"
                info["CALL_TYPE"] = "hardcall"
                yield replace(
                    call,
                    phased={sample: False for sample in call.genotypes},
                    phase_sets={},
                    info=info,
                )

    def iter_region(self, region: Region) -> Iterator[VariantCall]:
        yield from self._export_vcf(region)

    def iter_variants(self, allowed_variants, region=None) -> Iterator[VariantCall]:
        # When region is given, PLINK exports just its variants; otherwise the
        # caller requested a genome-wide variant-list selection.
        for call in self._export_vcf(region, allowed_variants):
            if variant_selected(call.chrom, call.pos, call.vid, allowed_variants):
                yield call

    def iter_all_variants(self) -> Iterator[VariantCall]:
        yield from self._export_vcf(None)


class HapMapReader:
    """HapMap (.hmp, .hmp.txt, .hmp.gz) reader for regional haplotype analysis.
    Supports .gz compressed files directly without manual decompression.
    """

    IUPAC_HET = {
        "R": ("A", "G"),
        "Y": ("C", "T"),
        "S": ("G", "C"),
        "W": ("A", "T"),
        "K": ("G", "T"),
        "M": ("A", "C"),
        "B": ("C", "G"),
        "D": ("A", "G"),
        "H": ("A", "C"),
        "V": ("A", "C"),
    }

    def __init__(self, path: str) -> None:
        self.path = os.fspath(path)
        if not os.path.exists(self.path):
            raise FileNotFoundError(f"HapMap file not found: {self.path}")
        self.samples, self.variants, self._variant_index = self._parse_file()
        self.backend = "hapmap"

    @staticmethod
    def _open_text(path: str):
        if path.endswith(".gz"):
            return gzip.open(path, "rt", encoding="utf-8", errors="replace")
        return open(path, "rt", encoding="utf-8", errors="replace")

    def _parse_file(self):
        samples: List[str] = []
        variants: List[Tuple[str, str, int, str, Tuple[str, ...], List[str]]] = []
        by_chrom: Dict[str, List[Tuple[int, int]]] = {}

        with self._open_text(self.path) as fh:
            header_found = False
            for line_no, raw in enumerate(fh, 1):
                line = raw.rstrip()

                if not line or line.startswith("##"):
                    continue
                parts = line.split("	")
                if not header_found:
                    if len(parts) < 12:
                        raise ValueError(
                            f"HapMap file {self.path} has fewer than 12 columns; expected at least 11 metadata columns and samples"
                        )
                    samples = [p.strip() for p in parts[11:]]
                    if len(samples) != len(set(samples)) or any(not x for x in samples):
                        raise ValueError(
                            "HapMap sample IDs must be nonempty and unique"
                        )
                    header_found = True
                    continue

                if len(parts) != 11 + len(samples):
                    raise ValueError(
                        f"HapMap column count mismatch at {self.path}:{line_no}"
                    )
                vid = parts[0].strip()
                alleles_s = parts[1].strip()
                chrom = parts[2].strip()
                try:
                    pos = int(parts[3].strip())
                except ValueError as exc:
                    raise ValueError(
                        f"Invalid HapMap coordinate at {self.path}:{line_no}"
                    ) from exc
                if pos < 1:
                    raise ValueError(
                        f"Invalid HapMap coordinate at {self.path}:{line_no}"
                    )

                allele_list = [
                    a.strip().upper() for a in alleles_s.split("/") if a.strip()
                ]
                ref = allele_list[0] if allele_list else "N"
                alts = tuple(allele_list[1:]) if len(allele_list) > 1 else ()
                gt_cols = parts[11:]

                v_idx = len(variants)
                variants.append((chrom, vid, pos, ref, alts, gt_cols))
                by_chrom.setdefault(chrom, []).append((pos, v_idx))

        index: Dict[str, Tuple[List[int], List[int]]] = {}
        for chrom, entries in by_chrom.items():
            entries.sort(key=lambda item: (item[0], item[1]))
            index[chrom] = ([p for p, _ in entries], [i for _, i in entries])
        return samples, variants, index

    def _decode_variant(self, variant_tuple) -> VariantCall:
        chrom, vid, pos, ref, alts, gt_cols = variant_tuple
        allele_pool = [ref] + list(alts)

        def get_allele_idx(a: str) -> int:
            if a not in allele_pool:
                raise ValueError(
                    f"HapMap allele {a!r} is not declared at {chrom}:{pos}"
                )
            return allele_pool.index(a)

        genotypes: Dict[str, Tuple[Optional[int], ...]] = {}
        phased: Dict[str, bool] = {}

        for i, sample in enumerate(self.samples):
            val = gt_cols[i].strip().upper() if i < len(gt_cols) else "N"
            if (
                val in allele_pool
                and val not in {"N", "."}
                and any(a not in "ACGT" and len(a) == 1 for a in allele_pool)
            ):
                gt = (get_allele_idx(val), get_allele_idx(val))
                is_phased = False
            elif val in {
                "NN",
                "N",
                "--",
                "-",
                "XX",
                "X",
                ".",
                "./.",
                ".|.",
                "00",
                "0",
                "",
            }:
                gt = (None, None)
                is_phased = False
            elif "|" in val or "/" in val:
                parts = re.split(r"[|/]", val)
                gt = tuple(
                    None
                    if a in {"N", ".", ""} or a == "-" and a not in allele_pool
                    else get_allele_idx(a)
                    for a in parts
                )
                is_phased = "|" in val and "/" not in val
            elif len(val) == 2:
                if any(len(a) > 1 for a in allele_pool):
                    raise ValueError(
                        f"Indel HapMap GT requires explicit allele separators at {sample}@{chrom}:{pos}"
                    )
                gt = tuple(
                    None
                    if a in "N.0" or a == "-" and a not in allele_pool
                    else get_allele_idx(a)
                    for a in val
                )
                is_phased = False
            elif len(val) == 1:
                if val in "ACGT":
                    idx = get_allele_idx(val)
                    gt = (idx, idx)
                elif val in {"B", "D", "H", "V"}:
                    raise ValueError(
                        f"Ambiguous HapMap IUPAC genotype {val} at {sample}@{chrom}:{pos}; provide explicit alleles"
                    )
                elif val in self.IUPAC_HET:
                    b1, b2 = self.IUPAC_HET[val]
                    gt = (get_allele_idx(b1), get_allele_idx(b2))
                else:
                    raise ValueError(
                        f"Unsupported HapMap genotype {val!r} at {sample}@{chrom}:{pos}"
                    )
                is_phased = False
            else:
                raise ValueError(
                    f"Unsupported HapMap genotype {val!r} at {sample}@{chrom}:{pos}; provide explicit alleles"
                )

            genotypes[sample] = gt
            phased[sample] = is_phased

        final_alts = tuple(allele_pool[1:])
        info = {"SOURCE": "HAPMAP"}
        tokens = encode_alleles(ref, final_alts, info)
        return VariantCall(
            chrom, pos, vid or ".", ref, final_alts, tokens, genotypes, phased, info
        )

    def iter_region(self, region: Region) -> Iterator[VariantCall]:
        chrom_index = self._variant_index.get(region.chrom)
        if chrom_index is None:
            return
        positions, variant_indices = chrom_index
        left = bisect.bisect_left(positions, region.start)
        right = bisect.bisect_right(positions, region.end)
        if left >= right:
            return
        for variant_index in variant_indices[left:right]:
            yield self._decode_variant(self.variants[variant_index])

    def iter_variants(
        self, allowed_variants: set[str], region: Optional[Region] = None
    ) -> Iterator[VariantCall]:
        for var_tuple in self.variants:
            chrom, vid, pos, ref, alts, gt_cols = var_tuple
            if (
                vid in allowed_variants
                or str(pos) in allowed_variants
                or f"{chrom}:{pos}" in allowed_variants
                or f"{chrom}_{pos}" in allowed_variants
            ):
                if region is None or (
                    chrom == region.chrom and region.start <= pos <= region.end
                ):
                    yield self._decode_variant(var_tuple)

    def iter_all_variants(self) -> Iterator[VariantCall]:
        for var_tuple in self.variants:
            yield self._decode_variant(var_tuple)


def open_genotype_reader(
    vcf: Optional[str] = None,
    bfile: Optional[str] = None,
    pfile: Optional[str] = None,
    hmp: Optional[str] = None,
    prefer: str = "auto",
):
    """Choose the best genotype backend; Windows PLINK always comes from vendor."""
    if vcf and not (bfile or pfile or hmp):
        path_str = str(vcf).lower()
        if path_str.endswith((".pgen", ".pvar", ".pvar.gz", ".psam", ".psam.gz")):
            pfile, vcf = vcf, None
        elif path_str.endswith((".bed", ".bim", ".bim.gz", ".fam", ".fam.gz")):
            bfile, vcf = vcf, None
        elif path_str.endswith((".hmp", ".hmp.txt", ".hmp.gz")):
            hmp, vcf = vcf, None
    if bfile or pfile:
        kind = "bfile" if bfile else "pfile"
        prefix = bfile or pfile
        executable, flavor = find_plink(kind)
        if executable:
            return ExternalPlinkReader(prefix, kind, executable, flavor)
        if is_windows():
            raise RuntimeError(
                "The Windows version of EasyHap requires bundled PLINK2 for "
                "BED/PGEN input. Put plink2.exe (and its DLLs) under "
                "vendor/plink2/ and rebuild EasyHap.exe."
            )
        # Linux can still use native readers when PLINK is not installed.
        return PlinkBedReader(prefix) if bfile else Plink2PgenReader(prefix)
    if hmp:
        return HapMapReader(hmp)
    if vcf:
        return VCFReader(str(vcf), prefer=prefer)
    raise ValueError(
        "No genotype input provided: specify --vcf, --bfile, --pfile, or --hmp"
    )
