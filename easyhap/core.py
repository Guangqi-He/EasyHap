from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Tuple, Union
import csv
import gzip
import os
import re
import platform
import hashlib
import json
import time
import heapq
from functools import lru_cache
from datetime import datetime

import pandas as pd

from .io_utils import (
    ensure_dir,
    read_group_metadata,
    read_trait_file,
    sanitize_filename,
)
from .stats import bh_adjust, connected_component_clusters, fisher_exact_2x2
from .vcf_reader import (
    ABSENCE_ALLELE,
    MISSING_ALLELE,
    PlinkBedReader,
    Plink2PgenReader,
    HapMapReader,
    open_genotype_reader,
    Region,
    VCFReader,
    VariantCall,
    deletion_allele_indices,
    deletion_end,
    deletion_interval,
    read_regions,
    token_for_gt,
    variant_selected,
)
from .research import (
    write_population_haplotype_statistics,
    write_trait_association,
    write_ld_outputs,
)
from .annotation import filter_calls_by_gene_feature
from .performance import PerformanceMetrics, PerformanceMonitor
from . import __version__

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

STATE_ALPHABET = list(
    dict.fromkeys("ACGTRYSWKMBDHV0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ")
)


class SkipRegionError(ValueError):
    """Expected empty/insufficient marker selection; other failures are errors."""


@lru_cache(maxsize=32)
def _file_hash_cached(path, size, mtime_ns):
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def file_fingerprint(path):
    stat = os.stat(path)
    return _file_hash_cached(os.path.abspath(path), stat.st_size, stat.st_mtime_ns)


class AnalysisResults(list):
    """List-like analysis result carrying run-level performance metrics.

    This preserves backward compatibility with code that expects a normal list
    while allowing the GUI/CLI to report the exact analysis-window benchmark.
    """

    def __init__(
        self, values=(), performance: Optional[PerformanceMetrics] = None
    ) -> None:
        super().__init__(values)
        self.performance = performance
        self.status = "OK"
        self.region_statuses = []

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
    absence_mask: Dict[Tuple[str, int], set] = field(default_factory=dict)
    mode: str = "genotype"


def _gt_tokens(call: VariantCall, sample: str) -> List[str]:
    gt = call.genotypes.get(sample, tuple())
    return [token_for_gt(call.allele_tokens, a) for a in gt]


def _consensus_token(tokens, hetero_policy="slash"):
    if not tokens or all(t == MISSING_ALLELE for t in tokens):
        return MISSING_ALLELE
    # The default state is the complete unordered allele multiset (dosage retained).
    if hetero_policy == "slash":
        return "/".join(sorted(tokens))
    if MISSING_ALLELE in tokens:
        return MISSING_ALLELE
    uniq = set(tokens)
    if len(uniq) == 1:
        return tokens[0]
    if hetero_policy == "missing":
        return MISSING_ALLELE
    if hetero_policy == "iupac" and all(len(x) == 1 and x in "ACGT" for x in uniq):
        return IUPAC.get(frozenset(uniq), MISSING_ALLELE)
    return "/".join(sorted(tokens))


def _get_ploidy(calls: Sequence[VariantCall], sample: str) -> int:
    ploidy = 0
    for call in calls:
        ploidy = max(ploidy, len(call.genotypes.get(sample, tuple())))
    return max(ploidy, 1)


def _format_hap_id(n: int) -> str:
    return f"Hap{n:03d}"


def _sort_hap_sequences(
    seq_to_samples: Dict[Tuple[str, ...], List[str]],
) -> List[Tuple[Tuple[str, ...], List[str]]]:
    return sorted(seq_to_samples.items(), key=lambda kv: (-len(set(kv[1])), kv[0]))


def _sample_pool(reader_samples: Sequence[str], group_map: Dict[str, str]) -> List[str]:
    samples = [s for s in reader_samples if s in group_map]
    if not samples:
        raise ValueError("None of the samples in the group file are present in the VCF")
    return samples


def load_variant_list(path_or_str: Union[str, Sequence[str]]) -> set[str]:
    """Parse a file or comma/space-separated string of variant IDs/positions."""
    if not path_or_str:
        return set()
    variants = set()
    if isinstance(path_or_str, (list, tuple, set)):
        items = list(path_or_str)
    elif os.path.exists(str(path_or_str)):
        opener = gzip.open if str(path_or_str).endswith(".gz") else open
        with opener(str(path_or_str), "rt", encoding="utf-8", errors="replace") as fh:
            items = []
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                if line.upper() in {
                    "ID",
                    "POS",
                    "VARIANT",
                    "SNP",
                } or line.upper().startswith("CHROM"):
                    continue
                parts = re.split(r"[\s,\t]+", line)
                if len(parts) >= 2 and parts[1].isdigit():
                    items.append(f"{parts[0]}:{parts[1]}")
                    items.append(f"{parts[0]}_{parts[1]}")
                    if len(parts) >= 3 and not parts[2].isdigit():
                        items.append(parts[2])
                else:
                    for part in parts:
                        if part:
                            items.append(part)
    else:
        items = re.split(r"[\s,\t]+", str(path_or_str).strip())
    for item in items:
        item = str(item).strip()
        if item:
            variants.add(item)
    return variants


def filter_calls_by_variant_list(
    calls: Sequence[VariantCall], variant_list: Optional[Union[str, Sequence[str]]]
) -> Tuple[List[VariantCall], int]:
    if not variant_list:
        return list(calls), len(calls)
    allowed = load_variant_list(variant_list)
    retained = [c for c in calls if variant_selected(c.chrom, c.pos, c.vid, allowed)]
    return retained, len(retained)


def filter_variants_by_qc(
    calls: Sequence[VariantCall],
    samples: Sequence[str],
    max_missing: Optional[float] = None,
    min_maf: Optional[float] = None,
) -> Tuple[List[VariantCall], Optional[pd.DataFrame]]:
    """Filter variants in region by missing rate and minor allele frequency (MAF)."""
    if max_missing is None and min_maf is None:
        return list(calls), None

    n_samples = max(1, len(samples))
    retained: List[VariantCall] = []
    qc_rows = []

    for c in calls:
        n_missing = 0
        allele_counts: Dict[int, int] = {}
        for s in samples:
            gt = c.genotypes.get(s, ())
            if not gt or any(a is None or a < 0 for a in gt):
                n_missing += 1
            else:
                for a in gt:
                    allele_counts[a] = allele_counts.get(a, 0) + 1

        missing_rate = n_missing / float(n_samples)
        total_alleles = sum(allele_counts.values())
        if total_alleles == 0:
            maf = 0.0
        else:
            freqs = sorted(
                [cnt / float(total_alleles) for cnt in allele_counts.values()],
                reverse=True,
            )
            maf = freqs[1] if len(freqs) > 1 else 0.0

        pass_miss = (max_missing is None) or (missing_rate <= float(max_missing))
        pass_maf = (min_maf is None) or (maf >= float(min_maf))

        qc_status = "PASS"
        if not pass_miss and not pass_maf:
            qc_status = "FAIL_MISSING_AND_MAF"
        elif not pass_miss:
            qc_status = "FAIL_MISSING"
        elif not pass_maf:
            qc_status = "FAIL_MAF"

        qc_rows.append(
            {
                "CHROM": c.chrom,
                "POS": c.pos,
                "ID": c.variant_id,
                "MissingRate": missing_rate,
                "MAF": maf,
                "Filter": qc_status,
            }
        )
        if pass_miss and pass_maf:
            retained.append(c)

    qc_df = pd.DataFrame(qc_rows)
    return retained, qc_df


def filter_variants_by_fisher(
    calls, samples, group_map, group1, group2, alpha, adjust="none", test="allele"
):
    if alpha is None:
        return list(calls), None
    if not group1 or not group2 or group1 == group2:
        raise ValueError("Fisher filtering requires two distinct groups")
    if not 0 < alpha <= 1 or test not in {"allele", "carrier"}:
        raise ValueError("Invalid Fisher alpha/test")
    rows, pvals = [], []
    for call in calls:
        # Each sample contributes one observation, avoiding pseudoreplication of copies.
        alleles = (
            [None]
            if test == "carrier"
            else sorted(
                {
                    a
                    for s in samples
                    for a in call.genotypes.get(s, ())
                    if a is not None and a > 0
                }
            )
        )
        results = []
        for allele in alleles:
            table = [[0, 0], [0, 0]]
            for sample in samples:
                group = group_map.get(sample)
                gt = call.genotypes.get(sample, ())
                if (
                    group not in {group1, group2}
                    or not gt
                    or any(a is None for a in gt)
                ):
                    continue
                present = any(a > 0 for a in gt) if allele is None else allele in gt
                table[int(group == group2)][int(not present)] += 1
            a, b = table[0]
            c, d = table[1]
            p = fisher_exact_2x2(a, b, c, d) if a + b and c + d else float("nan")
            results.append((p, allele, a, b, c, d))
        valid = [r for r in results if pd.notna(r[0])]
        best = (
            min(valid, key=lambda x: x[0])
            if valid
            else (float("nan"), None, 0, 0, 0, 0)
        )
        # Bonferroni inside a multiallelic site; BH, if requested, is across sites.
        p = min(1.0, best[0] * len(valid)) if valid else float("nan")
        pvals.append(p)
        rows.append(
            {
                "CHROM": call.chrom,
                "POS": call.pos,
                "ID": call.variant_id,
                "Test": test,
                "BestAlleleIndex": best[1],
                "N_allele_tests": len(valid),
                f"{group1}_Present": best[2],
                f"{group1}_Absent": best[3],
                f"{group2}_Present": best[4],
                f"{group2}_Absent": best[5],
                "pvalue": p,
            }
        )
    qs = (
        bh_adjust(pvals)
        if adjust.lower() in {"bh", "fdr", "padj", "adjusted"}
        else pvals
    )
    kept = []
    for row, call, q in zip(rows, calls, qs):
        row["padj"] = q
        row["keep"] = q <= alpha
        if row["keep"]:
            kept.append(call)
    return kept, pd.DataFrame(rows)


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


def phase_anchor(calls, sample):
    anchors = set()
    for call in calls:
        gt = call.genotypes.get(sample, ())
        informative = len(set(a for a in gt if a is not None)) > 1 or (
            any(a is None for a in gt) and any(a is not None for a in gt)
        )
        if informative:
            anchors.add(call.phase_sets.get(sample) or "IMPLICIT")
    return (
        next(iter(anchors))
        if len(anchors) == 1
        else ("UNINFORMATIVE" if not anchors else "MULTIPLE")
    )


def validate_copy_phase(calls, samples):
    if len({c.chrom for c in calls}) > 1:
        raise ValueError(
            "Copy reconstruction cannot concatenate different chromosomes; analyze separate regions"
        )
    for sample in samples:
        ploidies = {
            len(c.genotypes[sample])
            for c in calls
            if any(a is not None for a in c.genotypes.get(sample, ()))
        }
        if len(ploidies) > 1:
            raise ValueError(
                f"Inconsistent called ploidy for {sample}; split regions with different ploidies"
            )
        for call in calls:
            gt = call.genotypes.get(sample, ())
            informative = len(set(a for a in gt if a is not None)) > 1 or (
                any(a is None for a in gt) and any(a is not None for a in gt)
            )
            if informative and not call.phased.get(sample, False):
                raise ValueError(
                    f"Copy mode requires fully phased heterozygous/partial GTs: {sample}@{call.chrom}:{call.pos}; use --mode genotype"
                )
        if phase_anchor(calls, sample) == "MULTIPLE":
            raise ValueError(
                f"Multiple phase sets for {sample}; split this region by PS or use --mode genotype"
            )


def build_structural_absence_mask(calls, samples, mode, context_calls=()):
    """Infer only removed reference bases, preserving uncertain copy assignment."""
    import bisect

    resolved = _normalize_mode(mode)
    by_chrom = {}
    for ci, call in enumerate(calls):
        by_chrom.setdefault(call.chrom, []).append((call.pos, ci))
    for entries in by_chrom.values():
        entries.sort()
    mask = {}
    for event in list(context_calls) + list(calls):
        intervals = {
            a: deletion_interval(event, a) for a in deletion_allele_indices(event)
        }
        if not intervals:
            continue
        entries = by_chrom.get(event.chrom, [])
        for sample in samples:
            gt = event.genotypes.get(sample, ())
            if not gt:
                continue
            complete = all(a is not None for a in gt)
            event_masks = {}
            for allele, (start, end) in intervals.items():
                copies = {i for i, a in enumerate(gt) if a == allele}
                if not copies:
                    continue
                for pos, ci in entries[
                    bisect.bisect_left(entries, (start, -1)) : bisect.bisect_right(
                        entries, (end, len(calls))
                    )
                ]:
                    call = calls[ci]
                    if call is event:
                        continue
                    target_gt = call.genotypes.get(sample, ())
                    if any(a is not None for a in target_gt) and len(target_gt) != len(
                        gt
                    ):
                        continue
                    if resolved == "copy" and len(copies) < len(gt):
                        if (
                            not event.phased.get(sample, False)
                            or phase_anchor([event, call], sample) == "MULTIPLE"
                        ):
                            continue
                    event_masks.setdefault(ci, set()).update(copies)
            # In genotype mode, multiple DEL alleles at the SAME record may jointly
            # remove all copies. Partial events at different records cannot be aligned.
            for ci, copies in event_masks.items():
                if resolved == "genotype" and not (complete and len(copies) == len(gt)):
                    continue
                mask.setdefault((sample, ci), set()).update(copies)
    return mask


def effective_tokens(call, sample, index, absence_mask, ploidy=None):
    tokens = _gt_tokens(call, sample)
    ploidy = max(ploidy or 0, max(absence_mask.get((sample, index), {-1})) + 1)
    if len(tokens) < ploidy:
        tokens.extend([MISSING_ALLELE] * (ploidy - len(tokens)))
    return [
        ABSENCE_ALLELE if i in absence_mask.get((sample, index), set()) else t
        for i, t in enumerate(tokens)
    ]


def write_variant_overlap_report(
    prefix: str, calls: Sequence[VariantCall], context_calls: Sequence[VariantCall] = ()
) -> str:
    """Write explicit overlap relationships and the action taken by EasyHap."""
    path = prefix + ".VariantOverlap.tsv"
    n_targets = len(calls)
    calls = list(calls) + list(context_calls)
    # Check actual overlapping pairs involving at least one target marker.
    ends = [max(c.pos + len(c.ref) - 1, deletion_end(c)) for c in calls]
    pairs = []
    active, expiry, last_chrom = {}, [], None
    for j in sorted(range(len(calls)), key=lambda i: (calls[i].chrom, calls[i].pos, i)):
        b = calls[j]
        if b.chrom != last_chrom:
            active, expiry, last_chrom = {}, [], b.chrom
        while expiry and expiry[0][0] < b.pos:
            _, expired = heapq.heappop(expiry)
            active.pop(expired, None)
        for i in active:
            if i < n_targets or j < n_targets:
                pairs.append((min(i, j), max(i, j)))
        active[j] = ends[j]
        heapq.heappush(expiry, (ends[j], j))
    rows = []
    for block, (i, j) in enumerate(sorted(pairs), 1):
        a, b = calls[i], calls[j]
        a_end, b_end = ends[i], ends[j]
        adel = bool(deletion_allele_indices(a))
        bdel = bool(deletion_allele_indices(b))
        if adel and a.pos < b.pos <= deletion_end(a):
            rel, action = (
                "contained_in_deletion",
                "structural_absence_mask_when_genotype_supports_deletion",
            )
        elif bdel and b.pos < a.pos <= deletion_end(b):
            rel, action = (
                "contained_in_deletion",
                "structural_absence_mask_when_genotype_supports_deletion",
            )
        elif a.pos == b.pos:
            rel, action = "shared_start", "retained_separately"
        else:
            rel, action = "overlap_or_nested", "retained_separately_with_report"
        rows.append(
            (f"OB{block:04d}", a.variant_id, b.variant_id, a.chrom,
             a.pos, a_end, b.pos, b_end, rel, action)
        )
    with open(path, "w", encoding="utf-8") as out:
        out.write(
            "BlockID\tVariant1\tVariant2\tCHROM\tStart1\tEnd1\tStart2\tEnd2\tRelationship\tAction\n"
        )
        for row in rows:
            out.write("\t".join(map(str, row)) + "\n")
    return path


def build_haplotypes(
    calls: Sequence[VariantCall],
    samples: Sequence[str],
    mode: str = "inbred",
    hetero_policy: str = "slash",
    context_calls: Sequence[VariantCall] = (),
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
    absence_mask = build_structural_absence_mask(
        calls, samples, resolved_mode, context_calls
    )

    if resolved_mode == "genotype":
        for sample in samples:
            seq_tokens: List[str] = []
            for ci, call in enumerate(calls):
                gt_tokens = effective_tokens(call, sample, ci, absence_mask)
                seq_tokens.append(_consensus_token(gt_tokens, hetero_policy))
            seq = tuple(seq_tokens)
            sample_hap_sequences[sample] = [seq]
            seq_to_samples.setdefault(seq, []).append(sample)
    else:
        for sample in samples:
            ploidy = _get_ploidy(list(calls) + list(context_calls), sample)
            copy_seqs: List[List[str]] = [[] for _ in range(ploidy)]
            for ci, call in enumerate(calls):
                gt = list(call.genotypes.get(sample, tuple()))
                if len(gt) < ploidy:
                    gt.extend([None] * (ploidy - len(gt)))
                masked = absence_mask.get((sample, ci), set())
                for i in range(ploidy):
                    tok = (
                        ABSENCE_ALLELE
                        if i in masked
                        else token_for_gt(call.allele_tokens, gt[i])
                    )
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
    header = (
        ["Hap", "ClusterID"] + [str(c.pos) for c in calls] + ["Accession", "Number"]
    )
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerow(["CHR", ""] + [c.chrom for c in calls] + ["", "NA"])
        writer.writerow(["POS", ""] + [str(c.pos) for c in calls] + ["", "NA"])
        for hid in sorted(hap_to_seq):
            accs = hap_to_samples[hid]
            writer.writerow(
                [hid, hap_clusters.get(hid, "NA")]
                + list(hap_to_seq[hid])
                + [";".join(accs), len(accs)]
            )


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
        traits_copy = traits.copy()
        traits_copy["Accession"] = traits_copy["Accession"].astype(str).str.strip()
        trait_cols = [c for c in traits_copy.columns if c != "Accession"]
        trait_map = traits_copy.set_index("Accession").to_dict(orient="index")
    rows = []
    for sample in sorted(sample_haps):
        sample_str = str(sample).strip()
        haps = sample_haps[sample]
        row = {
            "Hap": _hap_class(haps, mode),
            "ClusterID": _cluster_class(haps, hap_clusters, mode),
            "Accession": sample_str,
            "Type": group_map.get(sample, "NA"),
        }
        for col in trait_cols:
            row[col] = trait_map.get(sample_str, {}).get(col, pd.NA)
        rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(path, sep="	", index=False)
    return df


def write_processed_tables(prefix, calls, samples, mode="genotype", context_calls=()):
    variant_path = prefix + ".ProcessedVariants.tsv"
    genotype_path = prefix + ".SampleGenotypeTokens.tsv"
    mask = build_structural_absence_mask(calls, samples, mode, context_calls)
    pd.DataFrame(
        [
            {
                "CHROM": c.chrom,
                "POS": c.pos,
                "ID": c.variant_id,
                "REF": c.ref,
                "ALT": ",".join(c.alts),
                "AlleleTokens": ",".join(c.allele_tokens),
                "SVTYPE": c.info.get("SVTYPE", ""),
                "END": c.info.get("END", ""),
                "SVLEN": c.info.get("SVLEN", ""),
                "DeletionAlleleIndices": ",".join(map(str, deletion_allele_indices(c))),
            }
            for c in calls
        ]
    ).to_csv(variant_path, sep="\t", index=False)
    rows = []
    for ci, c in enumerate(calls):
        row = {"CHROM": c.chrom, "POS": c.pos, "ID": c.variant_id}
        for sample in samples:
            ts = effective_tokens(c, sample, ci, mask)
            row[sample] = ("|" if c.phased.get(sample, False) else "/").join(ts)
        rows.append(row)
    pd.DataFrame(rows).to_csv(genotype_path, sep="\t", index=False)
    # Remove the redundant output from an earlier run of this same region.
    obsolete_path = prefix + ".RawSampleGenotypeTokens.tsv"
    if os.path.isfile(obsolete_path):
        os.remove(obsolete_path)
    return variant_path, genotype_path


def write_analysis_metadata(
    prefix, calls, samples, sample_haps, mode, reader, context_calls=(),
    overlap_flank_bp=0, overlap_window=None, plot_min_count=1, plot_hap_level="hap"
):
    resolved = _normalize_mode(mode)
    sources = [
        getattr(reader, k, None)
        for k in ("path", "bed", "bim", "fam", "pgen", "pvar", "psam")
    ]
    sources = list(dict.fromkeys(p for p in sources if p and os.path.isfile(p)))
    hashes = {os.path.abspath(p): file_fingerprint(p) for p in sources}
    input_hash = hashlib.sha256("".join(sorted(hashes.values())).encode()).hexdigest()
    rows = []
    qc = []
    mask = build_structural_absence_mask(calls, samples, mode, context_calls)
    for sample in samples:
        n_effective = sum(
            any(t != MISSING_ALLELE for t in effective_tokens(c, sample, i, mask))
            for i, c in enumerate(calls)
        )
        n_called = sum(
            any(a is not None for a in c.genotypes.get(sample, ())) for c in calls
        )
        qc.append(
            {
                "Accession": sample,
                "N_markers": len(calls),
                "N_called_markers": n_called,
                "N_effective_markers": n_effective,
                "AllMissing": n_effective == 0,
            }
        )
        if resolved == "copy":
            for copy, hap in enumerate(sample_haps[sample], 1):
                rows.append(
                    {
                        "Accession": sample,
                        "Copy": copy,
                        "Hap": hap,
                        "Chrom": calls[0].chrom,
                        "PhaseSet": phase_anchor(
                            list(calls) + list(context_calls), sample
                        ),
                        "InputSHA256": input_hash,
                        "Mode": resolved,
                    }
                )
    pd.DataFrame(
        rows,
        columns=[
            "Accession",
            "Copy",
            "Hap",
            "Chrom",
            "PhaseSet",
            "InputSHA256",
            "Mode",
        ],
    ).to_csv(prefix + ".HaplotypeCopies.tsv", sep="\t", index=False)
    pd.DataFrame(qc).to_csv(prefix + ".SampleQC.tsv", sep="\t", index=False)
    mask_rows = []
    for (sample, index), copies in sorted(mask.items()):
        call = calls[index]
        raw = _gt_tokens(call, sample)
        mask_rows.append(
            {
                "Accession": sample,
                "CHROM": call.chrom,
                "POS": call.pos,
                "ID": call.variant_id,
                "MaskedCopies": ",".join(str(i + 1) for i in sorted(copies)),
                "RawTokens": "/".join(raw),
                "EffectiveTokens": "/".join(
                    effective_tokens(call, sample, index, mask)
                ),
                "Status": "CALLED_GT_OVERLAPS_DELETION"
                if any(i < len(raw) and raw[i] != MISSING_ALLELE for i in copies)
                else "INFERRED_FROM_DELETION",
            }
        )
    pd.DataFrame(
        mask_rows,
        columns=[
            "Accession",
            "CHROM",
            "POS",
            "ID",
            "MaskedCopies",
            "RawTokens",
            "EffectiveTokens",
            "Status",
        ],
    ).to_csv(prefix + ".StructuralAbsenceAudit.tsv", sep="\t", index=False)
    pd.DataFrame(
        [
            {
                "CHROM": c.chrom,
                "POS": c.pos,
                "ID": c.variant_id,
                "REF": c.ref,
                "ALT": ",".join(c.alts),
                "AlleleTokens": ",".join(c.allele_tokens),
            }
            for c in calls
        ]
    ).to_csv(prefix + ".AnalysisVariants.tsv", sep="\t", index=False)
    metadata = {
        "version": __version__,
        "mode": resolved,
        "backend": reader.backend,
        "input_sha256": input_hash,
        "sources": hashes,
        "n_samples": len(samples),
        "n_variants": len(calls),
        "plot_min_count": plot_min_count,
        "plot_hap_level": plot_hap_level,
        "trait_min_class_n": max(2, plot_min_count),
        "trait_count_basis": "finite_phenotype_per_selected_class",
        "overlap_flank_bp": overlap_flank_bp,
        "overlap_window": overlap_window,
        "context_deletions": [
            {"chrom": c.chrom, "pos": c.pos, "id": c.variant_id} for c in context_calls
        ],
    }
    with open(prefix + ".AnalysisManifest.json", "w", encoding="utf-8") as fh:
        json.dump(metadata, fh, ensure_ascii=False, indent=2)


def _state_encode_sequences(
    sequences: Dict[str, Sequence[str]],
) -> Tuple[Dict[str, str], List[Dict[str, str]]]:
    if not sequences:
        return {}, []
    names = list(sequences)
    n_sites = len(sequences[names[0]])
    per_site_maps: List[Dict[str, str]] = []
    encoded = {name: [] for name in names}
    for i in range(n_sites):
        tokens = sorted({sequences[name][i] for name in names})
        site_map: Dict[str, str] = {}
        used = {"N"}
        # Preserve simple DNA symbols when possible.
        for tok in tokens:
            if len(tok) == 1 and tok in "ACGTN" and (tok == "N" or tok not in used):
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
                raise ValueError(
                    "Too many allele states at one site to encode as single characters"
                )
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

    symbols = "".join(sorted(set("".join(hap_encoded.values())) - {"N"})) or "ACGT"
    with open(nex_path, "w", encoding="utf-8") as fh:
        fh.write("#NEXUS\n\nBEGIN DATA;\n")
        fh.write(f"  DIMENSIONS NTAX={n_tax} NCHAR={n_char};\n")
        fh.write(f'  FORMAT DATATYPE=STANDARD SYMBOLS="{symbols}" MISSING=N GAP=-;\n')
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
    pfile: Optional[str] = None,
    hmp: Optional[str] = None,
    variant_list: Optional[Union[str, Sequence[str]]] = None,
    max_missing: Optional[float] = None,
    min_maf: Optional[float] = None,
    reader_override=None,
    overlap_flank_bp: int = 0,
    **kwargs,
) -> HapResult:
    if not isinstance(overlap_flank_bp, int) or isinstance(overlap_flank_bp, bool) or overlap_flank_bp < 0:
        raise ValueError("overlap_flank_bp must be a nonnegative integer (bp)")
    explicit_region = region.chrom != "CustomVariants"
    if overlap_flank_bp and not explicit_region:
        raise ValueError("Positive overlap_flank_bp requires --region or --region-file")
    if max_missing is not None and not 0 <= max_missing <= 1:
        raise ValueError("max_missing must be between 0 and 1")
    if min_maf is not None and not 0 <= min_maf <= 0.5:
        raise ValueError("min_maf must be between 0 and 0.5")
    if not 0 <= cluster_threshold <= 1 or min_variants < 1:
        raise ValueError("Invalid clustering threshold/minimum marker count")
    if not isinstance(plot_min_count, int) or isinstance(plot_min_count, bool) or plot_min_count < 1:
        raise ValueError("plot_min_count must be a positive integer")
    ensure_dir(outdir)
    trait_file = trait_file or kwargs.get("trait_file") or kwargs.get("traits")
    pfile = pfile or kwargs.get("pfile")
    hmp = hmp or kwargs.get("hmp")
    variant_list = variant_list or kwargs.get("variant_list")
    max_missing = max_missing if max_missing is not None else kwargs.get("max_missing")
    min_maf = min_maf if min_maf is not None else kwargs.get("min_maf")
    sources = [s for s in (vcf_path, bfile, pfile, hmp) if s is not None]
    if len(sources) != 1 and reader_override is None:
        raise ValueError(
            "Provide exactly one genotype source: --vcf, --bfile, --pfile, or --hmp"
        )
    if (bfile or pfile) and _normalize_mode(mode) != "genotype":
        raise ValueError(
            "PLINK binary filesets are unphased and are supported only with --mode genotype"
        )
    reader = (
        reader_override
        if reader_override is not None
        else open_genotype_reader(
            vcf=vcf_path, bfile=bfile, pfile=pfile, hmp=hmp, prefer=vcf_backend
        )
    )
    group_metadata = read_group_metadata(group_file)
    if not group_metadata.empty:
        group_map = dict(
            zip(
                group_metadata["Accession"].astype(str),
                group_metadata["Type"].astype(str),
            )
        )
        samples = _sample_pool(reader.samples, group_map)
    else:
        samples = list(reader.samples)
        group_map = {s: "All" for s in samples}
        group_metadata = pd.DataFrame(
            {
                "Accession": samples,
                "Type": ["All"] * len(samples),
                "Latitude": [pd.NA] * len(samples),
                "Longitude": [pd.NA] * len(samples),
                "Location": [pd.NA] * len(samples),
            }
        )
    traits = read_trait_file(trait_file)

    safe_label = sanitize_filename(region.label)
    prefix = os.path.join(outdir, safe_label)

    started = time.monotonic()

    def progress(stage, detail=""):
        message = (
            f"[EasyHap] {region.vcf_label} Stage={stage}; "
            f"Elapsed_s={time.monotonic() - started:.3f}"
            + (f"; {detail}" if detail else "")
        )
        print(message, flush=True)
        with open(prefix + ".Progress.log", "a", encoding="utf-8") as log:
            log.write(message + "\n")

    progress("read_variants", f"Backend={reader.backend}; Samples={len(samples)}")

    variant_list = variant_list or kwargs.get("variant_list")
    max_missing = max_missing if max_missing is not None else kwargs.get("max_missing")
    min_maf = min_maf if min_maf is not None else kwargs.get("min_maf")

    allowed = load_variant_list(variant_list) if variant_list else None
    if allowed and hasattr(reader, "iter_variants"):
        calls = list(
            reader.iter_variants(
                allowed, region=None if region.chrom == "CustomVariants" else region
            )
        )
    elif isinstance(reader, VCFReader):
        calls = list(
            reader.iter_region(
                None if region.chrom == "CustomVariants" else region,
                variant_list=None if explicit_region else allowed,
                samples=samples,
            )
        )
    elif region.chrom == "CustomVariants":
        calls = list(reader.iter_all_variants())
    else:
        calls = list(reader.iter_region(region))
    overlap_candidates = list(calls)
    calls, _ = filter_calls_by_variant_list(calls, variant_list)
    if region.chrom == "CustomVariants" and calls:
        chroms = {c.chrom for c in calls}
        if len(chroms) > 1:
            if _normalize_mode(mode) == "copy" or str(gene_feature).lower() != "all":
                raise ValueError(
                    "Copy/gene-feature analysis requires variants on one chromosome"
                )
            region = Region("CustomVariants", 1, len(calls))
        else:
            region = Region(
                calls[0].chrom, min(c.pos for c in calls), max(c.pos for c in calls)
            )
        prefix = os.path.join(outdir, sanitize_filename(region.label))
    original_call_count = len(calls)
    requested_feature = str(gene_feature or "all").strip().lower()
    feature_intervals, primary_gene = [(region.start, region.end)], None
    if requested_feature != "all":
        if not gff_file:
            raise ValueError(f"--gene-feature {requested_feature} requires --gff")
        calls, feature_intervals, primary_gene = filter_calls_by_gene_feature(
            calls, gff_file, region, requested_feature
        )
    feature_filtered_count = len(calls)
    calls, qc_df = filter_variants_by_qc(calls, samples, max_missing, min_maf)
    if qc_df is not None:
        qc_df.to_csv(prefix + ".VariantQCFilter.tsv", sep="\t", index=False)
    progress("qc_complete", f"RetainedVariants={len(calls)}")
    overlap_window = None
    if explicit_region:
        window_start = max(1, region.start - overlap_flank_bp)
        window_end = region.end + overlap_flank_bp
        overlap_window = {"chrom": region.chrom, "start": window_start, "end": window_end}
        progress("overlap_context", f"Window={region.chrom}:{window_start}-{window_end}; FlankBP={overlap_flank_bp}")
        # Default 0 reuses the already-read regional variants: no second query.
        flank_regions = []
        if window_start < region.start:
            flank_regions.append(Region(region.chrom, window_start, region.start - 1))
        if window_end > region.end:
            flank_regions.append(Region(region.chrom, region.end + 1, window_end))
        for flank_region in flank_regions:
            progress("read_overlap_flank", f"Query={flank_region.vcf_label}")
            iterator = (reader.iter_region(flank_region, samples=samples)
                        if isinstance(reader, VCFReader) else reader.iter_region(flank_region))
            overlap_candidates.extend(iterator)
    else:
        progress("overlap_context", "SpecifiedVariantsOnly; FlankBP=0")
    progress("overlap_context_complete", f"CandidateVariants={len(overlap_candidates)}")

    # 4. Fisher exact test filtering
    calls, fisher_df = filter_variants_by_fisher(
        calls,
        samples,
        group_map,
        fisher_group1,
        fisher_group2,
        fisher_alpha,
        fisher_adjust,
        kwargs.get("fisher_test", "allele"),
    )
    if len(calls) < max(1, int(min_variants)):
        raise SkipRegionError(
            f"Insufficient variants after filtering in region {region.vcf_label}: found {len(calls)}, "
            f"minimum required is {max(1, int(min_variants))}"
        )
    # Deletions omitted as markers can still establish physical absence.
    # Retained events already contribute through calls and must not be repeated.
    selected_keys = {(c.chrom, c.pos, c.vid, c.ref, c.alts) for c in calls}
    overlap_context_calls = [
        c
        for c in overlap_candidates
        if (c.chrom, c.pos, c.vid, c.ref, c.alts) not in selected_keys
    ]
    context_calls = [c for c in overlap_context_calls if deletion_allele_indices(c)]
    pd.DataFrame([{
        "CHROM": region.chrom, "RegionStart": region.start, "RegionEnd": region.end,
        "FlankBP": overlap_flank_bp,
        "WindowStart": overlap_window["start"] if overlap_window else None,
        "WindowEnd": overlap_window["end"] if overlap_window else None,
        "CandidateVariants": len(overlap_candidates),
        "ContextDeletions": len(context_calls),
    }]).to_csv(prefix + ".OverlapWindow.tsv", sep="\t", index=False)
    if requested_feature != "all":
        with open(prefix + ".GeneFeatureFilter.tsv", "w", encoding="utf-8") as fh:
            fh.write(
                "Region\tGeneFeature\tOriginalVariants\tFeatureRetainedVariants\tFinalRetainedVariants\tIntervals\n"
            )
            interval_text = ";".join(f"{s}-{e}" for s, e in feature_intervals)
            fh.write(
                f"{region.vcf_label}\t{requested_feature}\t{original_call_count}\t{feature_filtered_count}\t{len(calls)}\t{interval_text}\n"
            )
    if fisher_df is not None:
        fisher_df.to_csv(prefix + ".FisherFilter.tsv", sep="\t", index=False)
    progress("write_processed", f"RetainedVariants={len(calls)}")
    if write_processed:
        write_processed_tables(prefix, calls, samples, mode, context_calls)
    write_variant_overlap_report(prefix, calls, overlap_context_calls)

    progress("build_haplotypes")
    seq_to_hap, hap_to_samples, hap_to_seq, sample_haps = build_haplotypes(
        calls, samples, mode, hetero_policy, context_calls
    )
    progress("cluster_haplotypes", f"Classes={len(hap_to_seq)}")
    cluster_labels = connected_component_clusters(
        [hap_to_seq[h] for h in sorted(hap_to_seq)], threshold=cluster_threshold
    )
    hap_clusters = {hid: cid for hid, cid in zip(sorted(hap_to_seq), cluster_labels)}

    hap_summary_path = prefix + ".HapSummary.tsv"
    hap_group_path = prefix + ".HapGroup.tsv"
    write_hap_summary(hap_summary_path, calls, hap_to_seq, hap_clusters, hap_to_samples)
    hap_group_df = write_hap_group(
        hap_group_path, sample_haps, hap_clusters, group_map, traits, mode
    )
    write_alignment_files(prefix, calls, hap_to_seq)
    write_analysis_metadata(
        prefix, calls, samples, sample_haps, mode, reader, context_calls,
        overlap_flank_bp=overlap_flank_bp, overlap_window=overlap_window,
        plot_min_count=plot_min_count, plot_hap_level=plot_hap_level
    )

    result = HapResult(
        region=region,
        variants=list(calls),
        hap_ids=seq_to_hap,
        hap_clusters=hap_clusters,
        hap_accessions=hap_to_samples,
        hap_sequences=hap_to_seq,
        sample_haps=sample_haps,
        sample_class={s: _hap_class(hs, mode) for s, hs in sample_haps.items()},
        sample_cluster_class={
            s: _cluster_class(hs, hap_clusters, mode) for s, hs in sample_haps.items()
        },
        hap_summary_path=hap_summary_path,
        hap_group_path=hap_group_path,
        output_prefix=prefix,
        absence_mask=build_structural_absence_mask(calls, samples, mode, context_calls),
        mode=_normalize_mode(mode),
    )

    progress("statistics")
    # Research outputs added in 1.1.0.
    write_population_haplotype_statistics(prefix, sample_haps, group_map, mode)
    write_trait_association(
        prefix,
        hap_group_df,
        traits_to_plot,
        class_col=("Hap" if plot_hap_level == "hap" else "ClusterID"),
        min_class_n=plot_min_count,
    )
    if make_ld:
        progress("ld")
        write_ld_outputs(
            prefix,
            calls,
            samples,
            build_structural_absence_mask(calls, samples, mode, context_calls),
        )

    if make_plots:
        progress("plots")
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
    progress("complete")
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
        "Timestamp",
        "OS",
        "Backend",
        "Status",
        "Samples",
        "Completed_regions",
        "Retained_variants",
        "Runtime_s",
        "Peak_RAM_GB",
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


_WORKER_READERS = {}


def _run_single_region_worker(task):
    reg = task["region"]
    try:
        if "reader_override" not in task:
            key = tuple(
                task.get(k)
                for k in ("vcf_path", "bfile", "pfile", "hmp", "vcf_backend")
            )
            if key not in _WORKER_READERS:
                _WORKER_READERS[key] = open_genotype_reader(
                    *key[:4], prefer=key[4] or "auto"
                )
            task = dict(task, reader_override=_WORKER_READERS[key])
        return "OK", reg, run_region_analysis(**task), ""
    except SkipRegionError as exc:
        return "SKIP", reg, None, str(exc)
    except Exception as exc:
        return "ERROR", reg, None, f"{type(exc).__name__}: {exc}"


def run_analysis(
    vcf_path=None,
    group_file=None,
    outdir="EasyHap_results",
    region=None,
    region_file=None,
    min_variants=2,
    bfile=None,
    pfile=None,
    hmp=None,
    variant_list=None,
    max_missing=None,
    min_maf=None,
    threads=1,
    **kwargs,
):
    monitor = PerformanceMonitor(sample_interval_s=0.01).start()
    results, records = [], []
    status, backend = "ERROR", "unopened"
    ensure_dir(outdir)
    log_path = os.path.join(outdir, "EasyHap.log")
    try:
        if sum(s is not None for s in (vcf_path, bfile, pfile, hmp)) != 1:
            raise ValueError(
                "Provide exactly one genotype source: --vcf, --bfile, --pfile, or --hmp"
            )
        if int(threads) < 1:
            raise ValueError("threads must be at least 1")
        regions = read_regions(region, region_file, variant_list)
        reader = open_genotype_reader(
            vcf_path, bfile, pfile, hmp, kwargs.get("vcf_backend", "auto")
        )
        backend = reader.backend
        common = dict(
            vcf_path=vcf_path,
            group_file=group_file,
            outdir=outdir,
            min_variants=min_variants,
            bfile=bfile,
            pfile=pfile,
            hmp=hmp,
            variant_list=variant_list,
            max_missing=max_missing,
            min_maf=min_maf,
            **kwargs,
        )
        tasks = [dict(common, region=reg) for reg in regions]
        if threads > 1 and len(tasks) > 1:
            from concurrent.futures import ProcessPoolExecutor
            import multiprocessing

            with ProcessPoolExecutor(
                max_workers=min(int(threads), len(tasks)),
                mp_context=multiprocessing.get_context("spawn"),
            ) as executor:
                outcomes = list(executor.map(_run_single_region_worker, tasks))
        else:
            outcomes = [
                _run_single_region_worker(dict(t, reader_override=reader))
                for t in tasks
            ]
        with open(log_path, "a", encoding="utf-8") as log:
            for flag, reg, result, reason in outcomes:
                records.append(
                    {"Region": reg.vcf_label, "Status": flag, "Reason": reason}
                )
                if result is not None:
                    results.append(result)
                message = f"[EasyHap] {flag} {reg.vcf_label}" + (
                    f": {reason}" if reason else ""
                )
                print(message)
                log.write(message + "\n")
        errors = [r for r in records if r["Status"] == "ERROR"]
        status = (
            ("PARTIAL" if results else "ERROR")
            if errors
            else (
                "OK_WITH_SKIPS"
                if len(results) < len(tasks) and results
                else ("OK" if results else "SKIPPED")
            )
        )
        if errors:
            raise RuntimeError(
                f"{len(errors)} region(s) failed; completed results are retained. "
                + "; ".join(r["Region"] + ": " + r["Reason"] for r in errors)
            )
        if not results:
            raise RuntimeError(
                "No region completed; see EasyHap.RegionStatus.tsv for skipped selections"
            )
    finally:
        metrics = monitor.stop()
        pd.DataFrame(records, columns=["Region", "Status", "Reason"]).to_csv(
            os.path.join(outdir, "EasyHap.RegionStatus.tsv"), sep="\t", index=False
        )
        _write_performance_record(outdir, metrics, status, results, backend)
        with open(log_path, "a", encoding="utf-8") as log:
            log.write(
                f"[PERFORMANCE] Status={status}; Runtime_s={metrics.runtime_s:.6f}; Peak_RAM_GB={metrics.peak_ram_gb:.6f}\n"
            )
        print(
            f"[EasyHap] Status={status}; Runtime_s={metrics.runtime_s:.6f}; Peak_RAM_GB={metrics.peak_ram_gb:.6f}"
        )
    output = AnalysisResults(results, performance=metrics)
    output.status, output.region_statuses = status, records
    return output


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
