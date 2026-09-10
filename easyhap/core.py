from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Tuple
import csv
import os
import platform
from datetime import datetime

import pandas as pd

from .io_utils import ensure_dir, read_group_metadata, read_trait_file, sanitize_filename
from .stats import bh_adjust, connected_component_clusters, fisher_exact_2x2
from .vcf_reader import (
    ABSENCE_ALLELE, MISSING_ALLELE, PlinkBedReader, Region, VCFReader, VariantCall,
    deletion_allele_indices, deletion_end, read_regions, token_for_gt,
)
from .research import write_population_haplotype_statistics, write_trait_association, write_ld_outputs
from .annotation import filter_calls_by_gene_feature
from .performance import PerformanceMetrics, PerformanceMonitor

IUPAC = {
    frozenset({"A", "G"}): "R",
    frozenset({"C", "T"}): "Y",
    frozenset({"G", "C"}): "S",
    frozenset({"A", "T"}): "W",
    frozenset({"G", "T"}): "K",
    frozenset({"A", "C"}): "M",
    frozenset({"A", "C", "G"}): "V",
    frozenset({"A", "C", "T"}): "H",
    frozenset({"A", "G", "T"}): "D",
    frozenset({"C", "G", "T"}): "B",
    frozenset({"A", "C", "G", "T"}): "N",
}

STATE_ALPHABET = list("ACGTNRYSWKMBDHV0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ")


class AnalysisResults(list):
    """List-like analysis result carrying run-level performance metrics.

    This preserves backward compatibility with code that expects a normal list
    while allowing the GUI/CLI to report the exact analysis-window benchmark.
    """

    def __init__(self, values=(), performance: Optional[PerformanceMetrics] = None) -> None:
        super().__init__(values)
        self.performance = performance

    @property
    def runtime_s(self) -> Optional[float]:
        return None if self.performance is None else self.performance.runtime_s

    @property
    def peak_ram_gb(self) -> Optional[float]:
        return None if self.performance is None else self.performance.peak_ram_gb


@dataclass
class HapResult:
    region: Region
    variants: List[VariantCall]
    hap_ids: Dict[Tuple[str, ...], str]
    hap_clusters: Dict[str, str]
    hap_accessions: Dict[str, List[str]]
    hap_sequences: Dict[str, Tuple[str, ...]]
    sample_haps: Dict[str, List[str]]
    sample_class: Dict[str, str]
    sample_cluster_class: Dict[str, str]
    hap_summary_path: str
    hap_group_path: str
    output_prefix: str


def _gt_tokens(call: VariantCall, sample: str) -> List[str]:
    gt = call.genotypes.get(sample, tuple())
    return [token_for_gt(call.allele_tokens, a) for a in gt]


def _consensus_token(tokens: Sequence[str], hetero_policy: str = "slash") -> str:
    valid = [t for t in tokens if t != MISSING_ALLELE]
    if not valid:
        return MISSING_ALLELE
    # Phase-independent genotype states must not depend on the textual GT order:
    # 0/1 and 1/0 are the same multilocus genotype state.
    uniq = sorted(set(valid))
    if len(uniq) == 1:
        return uniq[0]
    if hetero_policy == "missing":
        return MISSING_ALLELE
    if hetero_policy == "iupac":
        bases = set(uniq)
        if all(len(x) == 1 and x in "ACGT" for x in bases):
            return IUPAC.get(frozenset(bases), MISSING_ALLELE)
    return "/".join(uniq)


def _get_ploidy(calls: Sequence[VariantCall], sample: str) -> int:
    ploidy = 0
    for call in calls:
        ploidy = max(ploidy, len(call.genotypes.get(sample, tuple())))
    return max(ploidy, 1)


def _format_hap_id(n: int) -> str:
    return f"Hap{n:03d}"


def _sort_hap_sequences(seq_to_samples: Dict[Tuple[str, ...], List[str]]) -> List[Tuple[Tuple[str, ...], List[str]]]:
    return sorted(seq_to_samples.items(), key=lambda kv: (-len(set(kv[1])), kv[0]))


def _sample_pool(reader_samples: Sequence[str], group_map: Dict[str, str]) -> List[str]:
    samples = [s for s in reader_samples if s in group_map]
    if not samples:
        raise ValueError("None of the samples in the group file are present in the VCF")
    return samples


def filter_variants_by_fisher(
    calls: Sequence[VariantCall],
    samples: Sequence[str],
    group_map: Dict[str, str],
    group1: Optional[str],
    group2: Optional[str],
    alpha: Optional[float],
    adjust: str = "none",
) -> Tuple[List[VariantCall], Optional[pd.DataFrame]]:
    if not group1 or not group2 or alpha is None:
        return list(calls), None
    stats_rows = []
    pvals: List[float] = []
    for call in calls:
        a = b = c = d = 0
        for sample in samples:
            grp = group_map.get(sample)
            if grp not in {group1, group2}:
                continue
            gt = call.genotypes.get(sample, tuple())
            called = [x for x in gt if x is not None and x >= 0]
            if not called:
                continue
            alt_present = any(x > 0 for x in called)
            if grp == group1:
                if alt_present:
                    a += 1
                else:
                    b += 1
            elif grp == group2:
                if alt_present:
                    c += 1
                else:
                    d += 1
        p = fisher_exact_2x2(a, b, c, d) if (a + b + c + d) > 0 else 1.0
        pvals.append(p)
        stats_rows.append(
            {
                "CHROM": call.chrom,
                "POS": call.pos,
                "ID": call.variant_id,
                f"{group1}_ALT": a,
                f"{group1}_REF": b,
                f"{group2}_ALT": c,
                f"{group2}_REF": d,
                "pvalue": p,
            }
        )
    if adjust.lower() in {"bh", "fdr", "padj", "adjusted"}:
        padj = bh_adjust(pvals)
    else:
        padj = pvals
    kept: List[VariantCall] = []
    for row, call, q in zip(stats_rows, calls, padj):
        row["padj"] = q
        row["keep"] = q <= alpha
        if q <= alpha:
            kept.append(call)
    return kept, pd.DataFrame(stats_rows)


def _normalize_mode(mode: str) -> str:
    m = str(mode).strip().lower()
    if m in {"inbred", "genotype"}:
        return "genotype"
    if m in {"hybrid", "copy"}:
        return "copy"
    raise ValueError("--mode should be genotype/copy (legacy aliases: inbred/hybrid)")


def _is_heterozygous_called(gt: Sequence[Optional[int]]) -> bool:
    called = [a for a in gt if a is not None and a >= 0]
    return len(set(called)) > 1


def validate_copy_phase(calls: Sequence[VariantCall], samples: Sequence[str]) -> None:
    """Reject unphased heterozygous GTs before copy-resolved reconstruction."""
    examples: List[str] = []
    total = 0
    for call in calls:
        for sample in samples:
            gt = call.genotypes.get(sample, tuple())
            if _is_heterozygous_called(gt) and not call.phased.get(sample, False):
                total += 1
                if len(examples) < 8:
                    gt_text = "/".join("." if a is None else str(a) for a in gt)
                    examples.append(f"{sample}@{call.chrom}:{call.pos}={gt_text}")
    if total:
        detail = ", ".join(examples)
        raise ValueError(
            f"Copy-resolved mode requires phased heterozygous genotypes. Detected {total} unphased "
            f"heterozygous genotype(s); examples: {detail}. Use --mode genotype for phase-independent "
            "multilocus genotype analysis, or phase the VCF with a dedicated phasing tool first."
        )


def build_structural_absence_mask(
    calls: Sequence[VariantCall], samples: Sequence[str], mode: str
) -> Dict[Tuple[str, int], set]:
    """Map (sample, downstream-call-index) to chromosome copies absent by deletion.

    A mask is inferred only from an explicit deletion allele and a known deletion interval.
    Ordinary missing genotypes are never converted to genomic absence. In genotype mode,
    masking is only inferred when all called chromosome copies carry the deletion; an
    unphased heterozygous deletion is intentionally not copy-assigned.
    """
    resolved = _normalize_mode(mode)
    mask: Dict[Tuple[str, int], set] = {}
    for di, dcall in enumerate(calls):
        del_idx = set(deletion_allele_indices(dcall))
        if not del_idx:
            continue
        end = deletion_end(dcall)
        if end <= dcall.pos:
            continue
        for sample in samples:
            gt = list(dcall.genotypes.get(sample, tuple()))
            called = [a for a in gt if a is not None and a >= 0]
            if not called or not any(a in del_idx for a in called):
                continue
            absent_copies: set = set()
            if resolved == "copy":
                # Phase is irrelevant for homozygous/all-copy deletions; otherwise copy
                # assignment is valid only after validate_copy_phase() has succeeded.
                if all(a in del_idx for a in called):
                    absent_copies = {i for i, a in enumerate(gt) if a is not None}
                elif dcall.phased.get(sample, False):
                    absent_copies = {i for i, a in enumerate(gt) if a in del_idx}
            else:
                if all(a in del_idx for a in called):
                    absent_copies = {i for i in range(max(1, len(gt)))}
            if not absent_copies:
                continue
            for ci, call in enumerate(calls):
                if ci == di or call.chrom != dcall.chrom:
                    continue
                if dcall.pos < call.pos <= end:
                    mask.setdefault((sample, ci), set()).update(absent_copies)
    return mask


def write_variant_overlap_report(prefix: str, calls: Sequence[VariantCall]) -> str:
    """Write explicit overlap relationships and the action taken by EasyHap."""
    path = prefix + ".VariantOverlap.tsv"
    rows = []
    block = 0
    for i, a in enumerate(calls):
        a_end = max(a.pos + len(a.ref) - 1, deletion_end(a))
        for j in range(i + 1, len(calls)):
            b = calls[j]
            if b.chrom != a.chrom:
                continue
            b_end = max(b.pos + len(b.ref) - 1, deletion_end(b))
            if b.pos > a_end and a.pos > b_end:
                continue
            if max(a.pos, b.pos) > min(a_end, b_end):
                continue
            block += 1
            adel = bool(deletion_allele_indices(a)); bdel = bool(deletion_allele_indices(b))
            if adel and a.pos < b.pos <= deletion_end(a):
                rel, action = "contained_in_deletion", "structural_absence_mask_when_genotype_supports_deletion"
            elif bdel and b.pos < a.pos <= deletion_end(b):
                rel, action = "contained_in_deletion", "structural_absence_mask_when_genotype_supports_deletion"
            elif a.pos == b.pos:
                rel, action = "shared_start", "retained_separately"
            else:
                rel, action = "overlap_or_nested", "retained_separately_with_report"
            rows.append((f"OB{block:04d}", a.variant_id, b.variant_id, a.chrom, a.pos, a_end, b.pos, b_end, rel, action))
    with open(path, "w", encoding="utf-8") as out:
        out.write("BlockID\tVariant1\tVariant2\tCHROM\tStart1\tEnd1\tStart2\tEnd2\tRelationship\tAction\n")
        for row in rows:
            out.write("\t".join(map(str, row)) + "\n")
    return path


def build_haplotypes(
    calls: Sequence[VariantCall],
    samples: Sequence[str],
    mode: str = "inbred",
    hetero_policy: str = "slash",
) -> Tuple[
    Dict[Tuple[str, ...], str],
    Dict[str, List[str]],
    Dict[str, Tuple[str, ...]],
    Dict[str, List[str]],
]:
    """Return seq->hap_id, hap_id->accessions, hap_id->seq, sample->hap_ids."""
    seq_to_samples: Dict[Tuple[str, ...], List[str]] = {}
    sample_hap_sequences: Dict[str, List[Tuple[str, ...]]] = {}

    resolved_mode = _normalize_mode(mode)
    if resolved_mode == "copy":
        validate_copy_phase(calls, samples)
    absence_mask = build_structural_absence_mask(calls, samples, resolved_mode)

    if resolved_mode == "genotype":
        for sample in samples:
            seq_tokens: List[str] = []
            for ci, call in enumerate(calls):
                gt_tokens = _gt_tokens(call, sample)
                masked = absence_mask.get((sample, ci), set())
                if masked and gt_tokens:
                    gt_tokens = [ABSENCE_ALLELE if i in masked else tok for i, tok in enumerate(gt_tokens)]
                seq_tokens.append(_consensus_token(gt_tokens, hetero_policy))
            seq = tuple(seq_tokens)
            sample_hap_sequences[sample] = [seq]
            seq_to_samples.setdefault(seq, []).append(sample)
    else:
        for sample in samples:
            ploidy = _get_ploidy(calls, sample)
            copy_seqs: List[List[str]] = [[] for _ in range(ploidy)]
            for ci, call in enumerate(calls):
                gt = list(call.genotypes.get(sample, tuple()))
                if len(gt) < ploidy:
                    gt.extend([None] * (ploidy - len(gt)))
                masked = absence_mask.get((sample, ci), set())
                for i in range(ploidy):
                    tok = ABSENCE_ALLELE if i in masked else token_for_gt(call.allele_tokens, gt[i])
                    copy_seqs[i].append(tok)
            seqs = [tuple(x) for x in copy_seqs]
            sample_hap_sequences[sample] = seqs
            for seq in sorted(set(seqs)):
                seq_to_samples.setdefault(seq, []).append(sample)

    seq_to_hap: Dict[Tuple[str, ...], str] = {}
    hap_to_samples: Dict[str, List[str]] = {}
    hap_to_seq: Dict[str, Tuple[str, ...]] = {}
    for idx, (seq, accs) in enumerate(_sort_hap_sequences(seq_to_samples), 1):
        hid = _format_hap_id(idx)
        seq_to_hap[seq] = hid
        hap_to_samples[hid] = sorted(set(accs))
        hap_to_seq[hid] = seq

    sample_haps = {
        sample: [seq_to_hap[seq] for seq in seqs]
        for sample, seqs in sample_hap_sequences.items()
    }
    return seq_to_hap, hap_to_samples, hap_to_seq, sample_haps


def _hap_class(haps: Sequence[str], mode: str) -> str:
    if _normalize_mode(mode) == "genotype":
        return haps[0] if haps else "NA"
    # For hybrids, retain one entry for each haplotype copy, but sort for a stable diplotype/multiplotype label.
    return "_".join(sorted(haps)) if haps else "NA"


def _cluster_class(haps: Sequence[str], hap_clusters: Dict[str, str], mode: str) -> str:
    clusters = [hap_clusters.get(h, "NA") for h in haps]
    if _normalize_mode(mode) == "genotype":
        return clusters[0] if clusters else "NA"
    return "_".join(sorted(dict.fromkeys(clusters))) if clusters else "NA"


def write_hap_summary(
    path: str,
    calls: Sequence[VariantCall],
    hap_to_seq: Dict[str, Tuple[str, ...]],
    hap_clusters: Dict[str, str],
    hap_to_samples: Dict[str, List[str]],
) -> None:
    header = ["Hap", "ClusterID"] + [str(c.pos) for c in calls] + ["Accession", "Number"]
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerow(["CHR", ""] + [c.chrom for c in calls] + ["", "NA"])
        writer.writerow(["POS", ""] + [str(c.pos) for c in calls] + ["", "NA"])
        for hid in sorted(hap_to_seq):
            accs = hap_to_samples[hid]
            writer.writerow([hid, hap_clusters.get(hid, "NA")] + list(hap_to_seq[hid]) + [";".join(accs), len(accs)])


def write_hap_group(
    path: str,
    sample_haps: Dict[str, List[str]],
    hap_clusters: Dict[str, str],
    group_map: Dict[str, str],
    traits: Optional[pd.DataFrame],
    mode: str,
) -> pd.DataFrame:
    trait_map: Dict[str, Dict[str, object]] = {}
    trait_cols: List[str] = []
    if traits is not None:
        trait_cols = [c for c in traits.columns if c != "Accession"]
        trait_map = traits.set_index("Accession").to_dict(orient="index")
    rows = []
    for sample in sorted(sample_haps):
        haps = sample_haps[sample]
        row = {
            "Hap": _hap_class(haps, mode),
            "ClusterID": _cluster_class(haps, hap_clusters, mode),
            "Accession": sample,
            "Type": group_map.get(sample, "NA"),
        }
        for col in trait_cols:
            row[col] = trait_map.get(sample, {}).get(col, pd.NA)
        rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(path, sep="\t", index=False)
    return df


def write_processed_tables(prefix: str, calls: Sequence[VariantCall], samples: Sequence[str]) -> Tuple[str, str]:
    variant_path = prefix + ".ProcessedVariants.tsv"
    genotype_path = prefix + ".SampleGenotypeTokens.tsv"
    with open(variant_path, "w", encoding="utf-8") as out:
        out.write("CHROM\tPOS\tID\tREF\tALT\tAlleleTokens\tSVTYPE\tEND\tSVLEN\tDeletionAlleleIndices\n")
        for c in calls:
            dins = ",".join(map(str, deletion_allele_indices(c)))
            out.write(
                f"{c.chrom}\t{c.pos}\t{c.variant_id}\t{c.ref}\t{','.join(c.alts)}\t{','.join(c.allele_tokens)}\t"
                f"{c.info.get('SVTYPE','')}\t{c.info.get('END','')}\t{c.info.get('SVLEN','')}\t{dins}\n"
            )
    with open(genotype_path, "w", encoding="utf-8") as out:
        out.write("CHROM\tPOS\tID\t" + "\t".join(samples) + "\n")
        for c in calls:
            vals = []
            for s in samples:
                sep = "|" if c.phased.get(s, False) else "/"
                vals.append(sep.join(_gt_tokens(c, s)))
            out.write(f"{c.chrom}\t{c.pos}\t{c.variant_id}\t" + "\t".join(vals) + "\n")
    return variant_path, genotype_path


def _state_encode_sequences(sequences: Dict[str, Sequence[str]]) -> Tuple[Dict[str, str], List[Dict[str, str]]]:
    if not sequences:
        return {}, []
    names = list(sequences)
    n_sites = len(sequences[names[0]])
    per_site_maps: List[Dict[str, str]] = []
    encoded = {name: [] for name in names}
    for i in range(n_sites):
        tokens = sorted({sequences[name][i] for name in names})
        site_map: Dict[str, str] = {}
        used = set()
        # Preserve simple DNA symbols when possible.
        for tok in tokens:
            if len(tok) == 1 and tok in "ACGTN" and tok not in used:
                site_map[tok] = tok
                used.add(tok)
        for tok in tokens:
            if tok in site_map:
                continue
            for ch in STATE_ALPHABET:
                if ch not in used:
                    site_map[tok] = ch
                    used.add(ch)
                    break
            else:
                raise ValueError("Too many allele states at one site to encode as single characters")
        per_site_maps.append(site_map)
        for name in names:
            encoded[name].append(site_map[sequences[name][i]])
    return {name: "".join(chars) for name, chars in encoded.items()}, per_site_maps


def write_alignment_files(
    prefix: str,
    calls: Sequence[VariantCall],
    hap_to_seq: Dict[str, Tuple[str, ...]],
) -> None:
    hap_encoded, maps = _state_encode_sequences(hap_to_seq)
    fasta_path = prefix + ".Haplotype.fa"
    phy_path = prefix + ".Haplotype.phy"
    nex_path = prefix + ".Haplotype.nex"
    map_path = prefix + ".AlleleStateMap.tsv"

    with open(fasta_path, "w", encoding="utf-8") as fh:
        for name, seq in hap_encoded.items():
            fh.write(f">{name}\n{seq}\n")

    n_tax = len(hap_encoded)
    n_char = len(calls)
    with open(phy_path, "w", encoding="utf-8") as fh:
        fh.write(f"{n_tax} {n_char}\n")
        for name, seq in hap_encoded.items():
            fh.write(f"{name} {seq}\n")

    symbols = "".join(sorted(set("".join(hap_encoded.values())))) or "ACGT"
    with open(nex_path, "w", encoding="utf-8") as fh:
        fh.write("#NEXUS\n\nBEGIN DATA;\n")
        fh.write(f"  DIMENSIONS NTAX={n_tax} NCHAR={n_char};\n")
        fh.write(f"  FORMAT DATATYPE=STANDARD SYMBOLS=\"{symbols}\" MISSING=N GAP=-;\n")
        fh.write("  MATRIX\n")
        for name, seq in hap_encoded.items():
            fh.write(f"  {name} {seq}\n")
        fh.write("  ;\nEND;\n")

    with open(map_path, "w", encoding="utf-8") as fh:
        fh.write("CHROM\tPOS\tOriginalToken\tEncodedState\n")
        for call, site_map in zip(calls, maps):
            for tok, ch in sorted(site_map.items(), key=lambda kv: kv[1]):
                fh.write(f"{call.chrom}\t{call.pos}\t{tok}\t{ch}\n")



def run_region_analysis(
    vcf_path: Optional[str],
    group_file: Optional[str],
    region: Region,
    outdir: str,
    mode: str = "inbred",
    hetero_policy: str = "slash",
    trait_file: Optional[str] = None,
    fisher_group1: Optional[str] = None,
    fisher_group2: Optional[str] = None,
    fisher_alpha: Optional[float] = None,
    fisher_adjust: str = "none",
    cluster_threshold: float = 0.15,
    vcf_backend: str = "auto",
    write_processed: bool = True,
    make_plots: bool = False,
    gff_file: Optional[str] = None,
    gene_feature: str = "all",
    plot_formats: Sequence[str] = ("pdf",),
    traits_to_plot: Optional[Sequence[str]] = None,
    plot_hap_level: str = "hap",
    plot_min_count: int = 1,
    min_variants: int = 2,
    ref_color: str = "#70AD47",
    alt_color: str = "#4472C4",
    missing_color: str = "#D9D9D9",
    absence_color: str = "#FFFFFF",
    make_ld: bool = True,
    make_network: bool = True,
    map_style: str = "auto",
    hap_palette: Optional[Sequence[str]] = None,
    ld_cmap: Optional[str] = None,
    bfile: Optional[str] = None,
    reader_override=None,
) -> HapResult:
    ensure_dir(outdir)
    if bool(vcf_path) == bool(bfile):
        raise ValueError("Provide exactly one genotype source: --vcf or --bfile")
    if bfile and _normalize_mode(mode) != "genotype":
        raise ValueError("PLINK BED/BIM/FAM is unphased and is supported only with --mode genotype")
    reader = reader_override if reader_override is not None else (PlinkBedReader(bfile) if bfile else VCFReader(str(vcf_path), prefer=vcf_backend))
    group_metadata = read_group_metadata(group_file)
    if not group_metadata.empty:
        group_map = dict(zip(group_metadata["Accession"].astype(str), group_metadata["Type"].astype(str)))
        samples = _sample_pool(reader.samples, group_map)
    else:
        samples = list(reader.samples)
        group_map = {s: "All" for s in samples}
        group_metadata = pd.DataFrame({
            "Accession": samples, "Type": ["All"] * len(samples),
            "Latitude": [pd.NA] * len(samples), "Longitude": [pd.NA] * len(samples),
            "Location": [pd.NA] * len(samples),
        })
    traits = read_trait_file(trait_file)

    calls = list(reader.iter_region(region))
    original_call_count = len(calls)
    requested_feature = str(gene_feature or "all").strip().lower()
    if requested_feature != "all":
        if not gff_file:
            raise ValueError(f"--gene-feature {requested_feature} requires --gff")
        calls, feature_intervals, primary_gene = filter_calls_by_gene_feature(
            calls, gff_file, region, requested_feature
        )
    else:
        feature_intervals, primary_gene = [(region.start, region.end)], None
    feature_filtered_count = len(calls)

    if len(calls) < max(1, int(min_variants)):
        detail = (
            f" after --gene-feature {requested_feature} filtering (from {original_call_count} regional variants)"
            if requested_feature != "all" else ""
        )
        raise ValueError(
            f"Insufficient variants in region {region.vcf_label}{detail}: found {len(calls)}, "
            f"minimum required is {max(1, int(min_variants))}"
        )

    calls, fisher_df = filter_variants_by_fisher(
        calls, samples, group_map, fisher_group1, fisher_group2, fisher_alpha, fisher_adjust
    )
    if len(calls) < max(1, int(min_variants)):
        raise ValueError(
            f"Insufficient variants after filtering in region {region.vcf_label}: found {len(calls)}, "
            f"minimum required is {max(1, int(min_variants))}"
        )

    safe_label = sanitize_filename(region.label)
    prefix = os.path.join(outdir, safe_label)
    if requested_feature != "all":
        with open(prefix + ".GeneFeatureFilter.tsv", "w", encoding="utf-8") as fh:
            fh.write("Region\tGeneFeature\tOriginalVariants\tFeatureRetainedVariants\tFinalRetainedVariants\tIntervals\n")
            interval_text = ";".join(f"{s}-{e}" for s, e in feature_intervals)
            fh.write(
                f"{region.vcf_label}\t{requested_feature}\t{original_call_count}\t{feature_filtered_count}\t{len(calls)}\t{interval_text}\n"
            )
    if fisher_df is not None:
        fisher_df.to_csv(prefix + ".FisherFilter.tsv", sep="\t", index=False)
    if write_processed:
        write_processed_tables(prefix, calls, samples)
    write_variant_overlap_report(prefix, calls)

    seq_to_hap, hap_to_samples, hap_to_seq, sample_haps = build_haplotypes(calls, samples, mode, hetero_policy)
    cluster_labels = connected_component_clusters([hap_to_seq[h] for h in sorted(hap_to_seq)], threshold=cluster_threshold)
    hap_clusters = {hid: cid for hid, cid in zip(sorted(hap_to_seq), cluster_labels)}

    hap_summary_path = prefix + ".HapSummary.tsv"
    hap_group_path = prefix + ".HapGroup.tsv"
    write_hap_summary(hap_summary_path, calls, hap_to_seq, hap_clusters, hap_to_samples)
    hap_group_df = write_hap_group(hap_group_path, sample_haps, hap_clusters, group_map, traits, mode)
    write_alignment_files(prefix, calls, hap_to_seq)

    result = HapResult(
        region=region,
        variants=list(calls),
        hap_ids=seq_to_hap,
        hap_clusters=hap_clusters,
        hap_accessions=hap_to_samples,
        hap_sequences=hap_to_seq,
        sample_haps=sample_haps,
        sample_class={s: _hap_class(hs, mode) for s, hs in sample_haps.items()},
        sample_cluster_class={s: _cluster_class(hs, hap_clusters, mode) for s, hs in sample_haps.items()},
        hap_summary_path=hap_summary_path,
        hap_group_path=hap_group_path,
        output_prefix=prefix,
    )

    # Research outputs added in 1.1.0.
    write_population_haplotype_statistics(prefix, sample_haps, group_map, mode)
    write_trait_association(prefix, hap_group_df, traits_to_plot, class_col=("Hap" if plot_hap_level == "hap" else "ClusterID"))
    if make_ld:
        write_ld_outputs(prefix, calls, samples)

    if make_plots:
        from .plotting import make_all_plots
        make_all_plots(
            result=result,
            hap_group_df=hap_group_df,
            group_map=group_map,
            group_metadata=group_metadata,
            gff_file=gff_file,
            plot_formats=plot_formats,
            traits_to_plot=traits_to_plot,
            plot_hap_level=plot_hap_level,
            plot_min_count=plot_min_count,
            ref_color=ref_color,
            alt_color=alt_color,
            missing_color=missing_color,
            absence_color=absence_color,
            make_ld_plot=make_ld,
            make_network_plot=make_network,
            map_style=map_style,
            hap_palette=hap_palette,
            ld_cmap=ld_cmap,
        )
    return result


def _write_performance_record(
    outdir: str,
    metrics: PerformanceMetrics,
    status: str,
    results: Sequence[HapResult],
    backend: str,
) -> str:
    """Append one machine-readable benchmark row and return its path."""
    path = os.path.join(outdir, "EasyHap.performance.tsv")
    fields = [
        "Timestamp", "OS", "Backend", "Status", "Samples",
        "Completed_regions", "Retained_variants", "Runtime_s", "Peak_RAM_GB",
    ]
    samples = max((len(r.sample_class) for r in results), default=0)
    retained = ",".join(str(len(r.variants)) for r in results) if results else "NA"
    row = {
        "Timestamp": datetime.now().astimezone().isoformat(timespec="seconds"),
        "OS": f"{platform.system()} {platform.release()}",
        "Backend": backend,
        "Status": status,
        "Samples": samples if samples else "NA",
        "Completed_regions": len(results),
        "Retained_variants": retained,
        "Runtime_s": f"{metrics.runtime_s:.6f}",
        "Peak_RAM_GB": f"{metrics.peak_ram_gb:.6f}",
    }
    write_header = not os.path.exists(path) or os.path.getsize(path) == 0
    with open(path, "a", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        if write_header:
            writer.writeheader()
        writer.writerow(row)
    return path


def run_analysis(
    vcf_path: Optional[str] = None,
    group_file: Optional[str] = None,
    outdir: str = "EasyHap_results",
    region: Optional[str] = None,
    region_file: Optional[str] = None,
    min_variants: int = 2,
    **kwargs,
) -> AnalysisResults:
    """Run one or more regions and record cross-platform performance.

    Runtime is wall-clock time measured strictly inside the analysis call, so
    GUI parameter-entry time is excluded. Peak RAM is the maximum simultaneous resident memory/working set of the
    EasyHap process tree observed during the same analysis window.
    """
    monitor = PerformanceMonitor(sample_interval_s=0.01).start()
    results: List[HapResult] = []
    status = "ERROR"
    log_path = os.path.join(outdir, "EasyHap.log")

    try:
        ensure_dir(outdir)
        regions = read_regions(region, region_file)
        # PLINK BIM/FAM parsing and the coordinate index can be expensive for very
        # large datasets. Build the PLINK reader once per batch and reuse it across
        # regions; each region then performs only indexed BED seeks.
        shared_reader = None
        bfile = kwargs.get("bfile")
        if bool(vcf_path) == bool(bfile):
            raise ValueError("Provide exactly one genotype source: --vcf or --bfile")
        if bfile:
            mode = kwargs.get("mode", "genotype")
            if _normalize_mode(mode) != "genotype":
                raise ValueError("PLINK BED/BIM/FAM is unphased and is supported only with --mode genotype")
            shared_reader = PlinkBedReader(bfile)

        with open(log_path, "a", encoding="utf-8") as log:
            for reg in regions:
                log.write(f"[INFO] Processing {reg.vcf_label}\n")
                log.flush()
                try:
                    res = run_region_analysis(
                        vcf_path, group_file, reg, outdir,
                        min_variants=min_variants, reader_override=shared_reader, **kwargs
                    )
                    results.append(res)
                    log.write(f"[INFO] Completed {reg.vcf_label}; variants={len(res.variants)}; haplotypes={len(res.hap_sequences)}\n")
                except ValueError as exc:
                    msg = str(exc)
                    if "variants" in msg.lower() or "filter" in msg.lower() or "--gene-feature" in msg.lower():
                        log.write(f"[WARNING] Skipped {reg.vcf_label}: {msg}\n")
                        print(f"[EasyHap] SKIP {reg.vcf_label}: {msg}")
                        continue
                    log.write(f"[ERROR] {reg.vcf_label}: {msg}\n")
                    raise
                except Exception as exc:
                    # A malformed single region should be visible in the log; unexpected errors remain fatal.
                    log.write(f"[ERROR] {reg.vcf_label}: {type(exc).__name__}: {exc}\n")
                    raise
                finally:
                    log.flush()
        status = "OK"
    finally:
        metrics = monitor.stop()
        backend = "plink-bed" if kwargs.get("bfile") else str(kwargs.get("vcf_backend", "auto"))
        perf_path = None
        try:
            ensure_dir(outdir)
            perf_path = _write_performance_record(outdir, metrics, status, results, backend)
            with open(log_path, "a", encoding="utf-8") as log:
                log.write(
                    f"[PERFORMANCE] Status={status}; Runtime_s={metrics.runtime_s:.6f}; "
                    f"Peak_RAM_GB={metrics.peak_ram_gb:.6f}\n"
                )
        except Exception as perf_exc:
            print(f"[EasyHap] WARNING: could not write performance record: {perf_exc}")
        print(f"[EasyHap] Runtime_s={metrics.runtime_s:.6f}")
        print(f"[EasyHap] Peak_RAM_GB={metrics.peak_ram_gb:.6f}")
        if perf_path:
            print(f"[EasyHap] Performance_record={perf_path}")

    return AnalysisResults(results, performance=metrics)


def prepare_vcf_tables(
    vcf_path: str,
    outdir: str,
    region: Optional[str] = None,
    region_file: Optional[str] = None,
    vcf_backend: str = "auto",
) -> List[Tuple[str, str]]:
    ensure_dir(outdir)
    reader = VCFReader(vcf_path, prefer=vcf_backend)
    outputs = []
    for reg in read_regions(region, region_file):
        calls = list(reader.iter_region(reg))
        if not calls:
            continue
        prefix = os.path.join(outdir, sanitize_filename(reg.label))
        outputs.append(write_processed_tables(prefix, calls, reader.samples))
    return outputs
