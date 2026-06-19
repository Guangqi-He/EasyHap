from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Tuple
import gzip
import os
import re

MISSING_ALLELE = "N"


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
        return int(str(val).split(",")[0])
    except Exception:
        return None


def _symbolic_token(allele: str, info: Optional[Dict[str, str]] = None) -> str:
    name = allele.strip("<>").upper()
    svlen = _info_get_int(info or {}, "SVLEN")
    if svlen is not None and svlen != 0:
        return f"{svlen:+d}"
    if "DEL" in name or "ABS" in name:
        return "DEL"
    if "INS" in name or "DUP" in name or "PRE" in name:
        return "INS"
    if "PAV" in name:
        return "PAV"
    return name or "SV"


def encode_alleles(ref: str, alts: Sequence[str], info: Optional[Dict[str, str]] = None, max_literal_len: int = 20) -> Tuple[str, ...]:
    """Encode REF/ALT alleles into compact tokens for haplotype tables.

    SNPs and short MNPs are kept as sequence strings. For ALT alleles, length
    differences relative to REF are encoded as +N or -N. Symbolic alleles such
    as <INS>, <DEL>, and <PAV> are converted using SVLEN when available.
    """
    ref = ref or "N"
    tokens: List[str] = []
    regular_alt_lengths = [len(a) for a in alts if a and a not in {".", "*"} and not (a.startswith("<") and a.endswith(">"))]
    has_indel_or_symbolic = any(length != len(ref) for length in regular_alt_lengths) or any(
        bool(a) and a.startswith("<") and a.endswith(">") for a in alts
    )
    # REF: keep SNP/MNP literal for ordinary SNP/MNP sites. For indel/SV/PAV
    # sites, avoid writing a long REF sequence into the haplotype table.
    if has_indel_or_symbolic and len(ref) > 1:
        tokens.append(f"REF{len(ref)}")
    elif len(ref) <= max_literal_len and all(c.upper() in "ACGTN." for c in ref):
        tokens.append(ref.upper().replace(".", "N"))
    else:
        tokens.append(f"REF{len(ref)}")

    for alt in alts:
        alt = alt or "."
        if alt in {".", "*"}:
            tokens.append(MISSING_ALLELE)
            continue
        if alt.startswith("<") and alt.endswith(">"):
            tokens.append(_symbolic_token(alt, info))
            continue
        diff = len(alt) - len(ref)
        if diff > 0:
            tokens.append(f"+{diff}")
        elif diff < 0:
            tokens.append(f"-{abs(diff)}")
        else:
            alt_u = alt.upper()
            if len(alt_u) <= max_literal_len and all(c in "ACGTN" for c in alt_u):
                tokens.append(alt_u)
            else:
                tokens.append(f"SEQ{len(alt)}")
    return tuple(tokens)


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

    def _init_backend(self) -> None:
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
        self.samples = self._read_samples_plain()
        self.backend = "plain"

    def _open_text(self):
        return gzip.open(self.path, "rt") if self.path.endswith(".gz") else open(self.path, "rt")

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
        iterator: Iterable = None  # type: ignore
        try:
            iterator = self._vcf(region.vcf_label)
        except Exception:
            # no index or contig naming mismatch: fall back to whole-file scan.
            iterator = self._vcf
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
        try:
            iterator = self._vcf.fetch(region.chrom, region.start - 1, region.end)
        except Exception:
            iterator = self._vcf.fetch()
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
