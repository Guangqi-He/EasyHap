from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Tuple
import gzip

from .vcf_reader import Region, VariantCall


@dataclass
class AnnotationFeature:
    seqid: str
    start: int
    end: int
    strand: str
    ftype: str
    attrs: Dict[str, str]


def _parse_attrs(text: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for item in text.split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
        elif " " in item:
            key, value = item.split(" ", 1)
            value = value.strip('"')
        else:
            continue
        out[key.strip()] = value.strip()
    return out


def read_annotation_features(path: str, chrom: str, start: int, end: int) -> List[AnnotationFeature]:
    """Read gene-structure features overlapping a genomic interval from GFF3/GTF."""
    opener = gzip.open if path.endswith(".gz") else open
    type_map = {
        "gene": "gene",
        "mrna": "mRNA",
        "transcript": "transcript",
        "exon": "exon",
        "cds": "CDS",
        "utr": "UTR",
        "five_prime_utr": "five_prime_UTR",
        "5utr": "five_prime_UTR",
        "5_prime_utr": "five_prime_UTR",
        "three_prime_utr": "three_prime_UTR",
        "3utr": "three_prime_UTR",
        "3_prime_utr": "three_prime_UTR",
    }
    keep_types = set(type_map.values())
    features: List[AnnotationFeature] = []
    with opener(path, "rt", encoding="utf-8") as fh:  # type: ignore[arg-type]
        for raw in fh:
            if not raw or raw.startswith("#"):
                continue
            parts = raw.rstrip("\n\r").split("\t")
            if len(parts) < 9:
                continue
            seqid, _, raw_type, s, e, _, strand, _, attrs = parts
            if seqid != chrom:
                continue
            try:
                si, ei = int(s), int(e)
            except ValueError:
                continue
            if ei < start or si > end:
                continue
            ftype = type_map.get(raw_type.lower(), raw_type)
            if ftype not in keep_types:
                continue
            features.append(AnnotationFeature(seqid, si, ei, strand, ftype, _parse_attrs(attrs)))
    return features


def _overlap_len(start1: int, end1: int, start2: int, end2: int) -> int:
    return max(0, min(end1, end2) - max(start1, start2) + 1)


def choose_primary_gene(features: Sequence[AnnotationFeature], region: Region) -> Optional[AnnotationFeature]:
    """Choose the gene/transcript with the largest overlap with the requested region."""
    genes = [f for f in features if f.ftype == "gene"]
    if not genes:
        genes = [f for f in features if f.ftype in {"mRNA", "transcript"}]
    if not genes:
        return None
    return max(
        genes,
        key=lambda f: (
            _overlap_len(f.start, f.end, region.start, region.end),
            -(f.end - f.start),
        ),
    )


def _merge_intervals(intervals: Iterable[Tuple[int, int]]) -> List[Tuple[int, int]]:
    vals = sorted((min(a, b), max(a, b)) for a, b in intervals)
    if not vals:
        return []
    merged = [vals[0]]
    for s, e in vals[1:]:
        ps, pe = merged[-1]
        if s <= pe + 1:
            merged[-1] = (ps, max(pe, e))
        else:
            merged.append((s, e))
    return merged


def _subtract_intervals(base: Tuple[int, int], covered: Sequence[Tuple[int, int]]) -> List[Tuple[int, int]]:
    segments = [base]
    for cs, ce in _merge_intervals(covered):
        updated: List[Tuple[int, int]] = []
        for s, e in segments:
            if ce < s or cs > e:
                updated.append((s, e))
                continue
            if cs > s:
                updated.append((s, cs - 1))
            if ce < e:
                updated.append((ce + 1, e))
        segments = updated
    return [(s, e) for s, e in segments if s <= e]


def feature_intervals(
    features: Sequence[AnnotationFeature],
    region: Region,
    feature: str,
) -> Tuple[List[Tuple[int, int]], Optional[AnnotationFeature]]:
    """Return intervals for exon/intron/CDS/UTR from the primary overlapping gene."""
    requested = str(feature).strip().lower()
    if requested == "all":
        return [(region.start, region.end)], None
    if requested not in {"exon", "intron", "cds", "utr"}:
        raise ValueError(f"Unsupported gene feature: {feature}")

    gene = choose_primary_gene(features, region)
    if gene is None:
        raise ValueError(
            f"No gene/transcript annotation overlaps {region.vcf_label}; cannot apply --gene-feature {requested}"
        )

    # Restrict child features to the selected gene span. This is intentionally permissive
    # across GFF3/GTF attribute conventions while remaining consistent with the plotted gene model.
    children = [
        f for f in features
        if f.start <= gene.end and f.end >= gene.start
        and f.ftype in {"exon", "CDS", "UTR", "five_prime_UTR", "three_prime_UTR"}
    ]
    exons = [(f.start, f.end) for f in children if f.ftype == "exon"]
    cds = [(f.start, f.end) for f in children if f.ftype == "CDS"]
    utrs = [
        (f.start, f.end)
        for f in children
        if f.ftype in {"UTR", "five_prime_UTR", "three_prime_UTR"}
    ]

    if requested == "exon":
        intervals = exons or (cds + utrs)
    elif requested == "cds":
        intervals = cds
    elif requested == "utr":
        intervals = utrs
    else:  # intron
        exon_like = exons or (cds + utrs)
        if not exon_like:
            raise ValueError(
                f"No exon/CDS/UTR records are available for the primary gene in {region.vcf_label}; "
                "cannot infer introns"
            )
        intervals = _subtract_intervals((gene.start, gene.end), exon_like)

    intervals = _merge_intervals(
        (max(s, region.start), min(e, region.end))
        for s, e in intervals
        if e >= region.start and s <= region.end
    )
    if not intervals:
        raise ValueError(
            f"No {requested} intervals overlap {region.vcf_label}; cannot apply --gene-feature {requested}"
        )
    return intervals, gene


def filter_calls_by_gene_feature(
    calls: Sequence[VariantCall],
    gff_file: str,
    region: Region,
    feature: str,
) -> Tuple[List[VariantCall], List[Tuple[int, int]], Optional[AnnotationFeature]]:
    requested = str(feature).strip().lower()
    if requested == "all":
        return list(calls), [(region.start, region.end)], None
    features = read_annotation_features(gff_file, region.chrom, region.start, region.end)
    intervals, gene = feature_intervals(features, region, requested)
    kept = [
        call for call in calls
        if any(s <= int(call.pos) <= e for s, e in intervals)
    ]
    return kept, intervals, gene
