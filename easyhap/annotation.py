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


def read_annotation_features(
    path: str, chrom: str, start: int, end: int
) -> List[AnnotationFeature]:
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
            seqid, _, raw_type, s, e, _, strand, _, attrs = parts[:9]
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
            features.append(
                AnnotationFeature(seqid, si, ei, strand, ftype, _parse_attrs(attrs))
            )
    return features


def _overlap_len(start1: int, end1: int, start2: int, end2: int) -> int:
    return max(0, min(end1, end2) - max(start1, start2) + 1)


def choose_primary_gene(
    features: Sequence[AnnotationFeature], region: Region
) -> Optional[AnnotationFeature]:
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


def _subtract_intervals(
    base: Tuple[int, int], covered: Sequence[Tuple[int, int]]
) -> List[Tuple[int, int]]:
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


def gene_children(features, gene):
    """Resolve GFF3 Parent or GTF gene_id/transcript_id, excluding neighboring genes."""
    primary_ids = {
        gene.attrs[k] for k in ("ID", "gene_id", "transcript_id") if gene.attrs.get(k)
    }
    selected_ids = set(primary_ids)
    for _ in range(len(features) + 1):
        old = set(selected_ids)
        for f in features:
            parents = set(f.attrs.get("Parent", "").split(",")) - {""}
            belongs = bool(
                parents & selected_ids
                or f.attrs.get("gene_id") in primary_ids
                or f.attrs.get("transcript_id") in selected_ids
            )
            if belongs:
                selected_ids.update(
                    f.attrs[k] for k in ("ID", "transcript_id") if f.attrs.get(k)
                )
        if old == selected_ids:
            break
    children = []
    for f in features:
        if f.ftype not in {"exon", "CDS", "UTR", "five_prime_UTR", "three_prime_UTR"}:
            continue
        parents = set(f.attrs.get("Parent", "").split(",")) - {""}
        if (
            parents & selected_ids
            or f.attrs.get("gene_id") in primary_ids
            or f.attrs.get("transcript_id") in selected_ids
        ):
            children.append(f)
    # An unannotated parent relationship is ambiguous; do not merge by coordinate overlap.
    return children


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

    children = gene_children(features, gene)
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
        by_transcript = {}
        source = [f for f in children if f.ftype == "exon"] or [
            f
            for f in children
            if f.ftype in {"CDS", "UTR", "five_prime_UTR", "three_prime_UTR"}
        ]
        for f in source:
            parents = f.attrs.get(
                "Parent", f.attrs.get("transcript_id", f.attrs.get("gene_id", ""))
            ).split(",")
            for parent in parents:
                by_transcript.setdefault(parent, []).append((f.start, f.end))
        intervals = []
        for exon_intervals in by_transcript.values():
            ordered = _merge_intervals(exon_intervals)
            intervals.extend(
                (left[1] + 1, right[0] - 1)
                for left, right in zip(ordered, ordered[1:])
                if right[0] > left[1] + 1
            )

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
    features = read_annotation_features(
        gff_file, region.chrom, region.start, region.end
    )
    intervals, gene = feature_intervals(features, region, requested)
    kept = [
        call for call in calls if any(s <= int(call.pos) <= e for s, e in intervals)
    ]
    return kept, intervals, gene
