from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Tuple
import bisect
import gzip
import os
import re

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
    """Small wrapper around cyvcf2/pysam with a plain-text fallback.

    Indexed VCF/BCF random access is used when cyvcf2 or pysam is installed and
    the input has a tabix/CSI index. The fallback parser is slower but keeps the
    tool usable for small examples and Windows GUI use.
    """

    def __init__(self, path: str, prefer: str = "auto") -> None:
        self.path = path
        self.prefer = prefer
        self.backend = "plain"
        self._vcf = None
        self.samples: List[str] = []
        self._init_backend()

    def _has_region_index(self) -> bool:
        """Return True when a tabix/CSI index usable for regional access exists."""
        candidates = [
            self.path + ".tbi",
            self.path + ".csi",
        ]
        # Some workflows keep the CSI beside the stem (e.g. sample.bcf.csi is
        # already covered above); keep this helper intentionally conservative.
        return any(os.path.exists(x) for x in candidates)

    def _is_binary_bcf(self) -> bool:
        return self.path.lower().endswith(".bcf")

    def _is_text_vcf(self) -> bool:
        lower = self.path.lower()
        return lower.endswith(".vcf") or lower.endswith(".vcf.gz") or lower.endswith(".vcf.bgz") or lower.endswith(".vcf.bgzip")

    def _init_backend(self) -> None:
        # Auto mode deliberately keeps unindexed textual VCF files on the
        # sequential plain-text parser.  cyvcf2/pysam regional queries are
        # index based and can otherwise raise lazily during iteration.
        if self.prefer == "auto" and self._is_text_vcf() and not self._has_region_index():
            self.samples = self._read_samples_plain()
            self.backend = "plain"
            return

        if self.prefer in {"auto", "cyvcf2"}:
            try:
                from cyvcf2 import VCF  # type: ignore

                self._vcf = VCF(self.path)
                self.samples = list(self._vcf.samples)
                self.backend = "cyvcf2"
                return
            except Exception:
                if self.prefer == "cyvcf2":
                    raise
        if self.prefer in {"auto", "pysam"}:
            try:
                import pysam  # type: ignore

                self._vcf = pysam.VariantFile(self.path)
                self.samples = list(self._vcf.header.samples)
                self.backend = "pysam"
                return
            except Exception:
                if self.prefer == "pysam":
                    raise

        # Binary BCF cannot be parsed by the text fallback.  Give a useful
        # error instead of trying to decode binary bytes as VCF text.
        if self._is_binary_bcf():
            raise RuntimeError(
                f"Unable to open BCF file {self.path!r}. Install/repair cyvcf2 or pysam. "
                "For fast regional access, create a CSI index."
            )

        self.samples = self._read_samples_plain()
        self.backend = "plain"

    def _open_text(self):
        lower = self.path.lower()
        compressed = lower.endswith(".gz") or lower.endswith(".bgz") or lower.endswith(".bgzip")
        return gzip.open(self.path, "rt") if compressed else open(self.path, "rt")

    def _read_samples_plain(self) -> List[str]:
        with self._open_text() as fh:
            for line in fh:
                if line.startswith("#CHROM"):
                    parts = line.rstrip("\n").split("\t")
                    return parts[9:]
        raise ValueError(f"No #CHROM header line found in {self.path}")

    def iter_region(self, region: Region) -> Iterator[VariantCall]:
        if self.backend == "cyvcf2":
            yield from self._iter_cyvcf2(region)
        elif self.backend == "pysam":
            yield from self._iter_pysam(region)
        else:
            yield from self._iter_plain(region)

    def _iter_cyvcf2(self, region: Region) -> Iterator[VariantCall]:
        assert self._vcf is not None

        # Indexed inputs use true random access.  For an explicitly selected
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
            samples = self.samples
            for line in fh:
                if not line or line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 8:
                    continue
                chrom, pos_s, vid, ref, alt_s, _qual, _filt, info_s = parts[:8]
                pos = int(pos_s)
                if chrom != region.chrom or pos < region.start or pos > region.end:
                    continue
                fmt = parts[8].split(":") if len(parts) > 8 else []
                gt_idx = fmt.index("GT") if "GT" in fmt else None
                info = _parse_info(info_s)
                alts = tuple([] if alt_s == "." else alt_s.split(","))
                tokens = encode_alleles(ref, alts, info)
                genotypes: Dict[str, Tuple[Optional[int], ...]] = {}
                phased: Dict[str, bool] = {}
                for sample, sample_field in zip(samples, parts[9:]):
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
                yield VariantCall(chrom, pos, vid, ref, alts, tokens, genotypes, phased, info)


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

