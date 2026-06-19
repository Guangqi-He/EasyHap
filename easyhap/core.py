from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Tuple
import csv
import os

import pandas as pd

from .io_utils import ensure_dir, read_group_file, read_trait_file, sanitize_filename
from .stats import bh_adjust, connected_component_clusters, fisher_exact_2x2
from .vcf_reader import MISSING_ALLELE, Region, VCFReader, VariantCall, read_regions, token_for_gt

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
    uniq = list(dict.fromkeys(valid))
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

    if mode == "inbred":
        for sample in samples:
            seq = tuple(_consensus_token(_gt_tokens(call, sample), hetero_policy) for call in calls)
            sample_hap_sequences[sample] = [seq]
            seq_to_samples.setdefault(seq, []).append(sample)
    elif mode == "hybrid":
        for sample in samples:
            ploidy = _get_ploidy(calls, sample)
            copy_seqs: List[List[str]] = [[] for _ in range(ploidy)]
            for call in calls:
                gt = list(call.genotypes.get(sample, tuple()))
                if len(gt) < ploidy:
                    gt.extend([None] * (ploidy - len(gt)))
                for i in range(ploidy):
                    copy_seqs[i].append(token_for_gt(call.allele_tokens, gt[i]))
            seqs = [tuple(x) for x in copy_seqs]
            sample_hap_sequences[sample] = seqs
            # Count accession membership once per unique haplotype sequence.
            for seq in sorted(set(seqs)):
                seq_to_samples.setdefault(seq, []).append(sample)
    else:
        raise ValueError("--mode should be inbred or hybrid")

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
    if mode == "inbred":
        return haps[0] if haps else "NA"
    # For hybrids, retain one entry for each haplotype copy, but sort for a stable diplotype/multiplotype label.
    return "_".join(sorted(haps)) if haps else "NA"


def _cluster_class(haps: Sequence[str], hap_clusters: Dict[str, str], mode: str) -> str:
    clusters = [hap_clusters.get(h, "NA") for h in haps]
    if mode == "inbred":
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
        out.write("CHROM\tPOS\tID\tREF\tALT\tAlleleTokens\n")
        for c in calls:
            out.write(f"{c.chrom}\t{c.pos}\t{c.variant_id}\t{c.ref}\t{','.join(c.alts)}\t{','.join(c.allele_tokens)}\n")
    with open(genotype_path, "w", encoding="utf-8") as out:
        out.write("CHROM\tPOS\tID\t" + "\t".join(samples) + "\n")
        for c in calls:
            vals = []
            for s in samples:
                vals.append("/".join(_gt_tokens(c, s)))
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
    sample_haps: Dict[str, List[str]],
    mode: str,
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

    # Sample-level copy alignment.
    sample_sequences: Dict[str, Tuple[str, ...]] = {}
    hap_id_to_seq = hap_to_seq
    for sample, haps in sample_haps.items():
        if mode == "inbred":
            sample_sequences[f"{sample}|{haps[0]}"] = hap_id_to_seq[haps[0]]
        else:
            for i, hid in enumerate(haps, 1):
                sample_sequences[f"{sample}|copy{i}|{hid}"] = hap_id_to_seq[hid]
    sample_encoded, _ = _state_encode_sequences(sample_sequences)
    with open(prefix + ".Haplotype_sample.fa", "w", encoding="utf-8") as fh:
        for name, seq in sample_encoded.items():
            fh.write(f">{name}\n{seq}\n")


def run_region_analysis(
    vcf_path: str,
    group_file: str,
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
    plot_formats: Sequence[str] = ("pdf",),
    traits_to_plot: Optional[Sequence[str]] = None,
    plot_hap_level: str = "hap",
    plot_min_count: int = 1,
) -> HapResult:
    ensure_dir(outdir)
    reader = VCFReader(vcf_path, prefer=vcf_backend)
    group_map = read_group_file(group_file)
    samples = _sample_pool(reader.samples, group_map)
    traits = read_trait_file(trait_file)

    calls = list(reader.iter_region(region))
    if not calls:
        raise ValueError(f"No variants found in region {region.vcf_label}")

    calls, fisher_df = filter_variants_by_fisher(
        calls, samples, group_map, fisher_group1, fisher_group2, fisher_alpha, fisher_adjust
    )
    if not calls:
        raise ValueError("No variants remained after Fisher exact filtering")

    safe_label = sanitize_filename(region.label)
    prefix = os.path.join(outdir, safe_label)
    if fisher_df is not None:
        fisher_df.to_csv(prefix + ".FisherFilter.tsv", sep="\t", index=False)
    if write_processed:
        write_processed_tables(prefix, calls, samples)

    seq_to_hap, hap_to_samples, hap_to_seq, sample_haps = build_haplotypes(calls, samples, mode, hetero_policy)
    cluster_labels = connected_component_clusters([hap_to_seq[h] for h in sorted(hap_to_seq)], threshold=cluster_threshold)
    hap_clusters = {hid: cid for hid, cid in zip(sorted(hap_to_seq), cluster_labels)}

    hap_summary_path = prefix + ".HapSummary.tsv"
    hap_group_path = prefix + ".HapGroup.tsv"
    write_hap_summary(hap_summary_path, calls, hap_to_seq, hap_clusters, hap_to_samples)
    hap_group_df = write_hap_group(hap_group_path, sample_haps, hap_clusters, group_map, traits, mode)
    write_alignment_files(prefix, calls, hap_to_seq, sample_haps, mode)

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

    if make_plots:
        from .plotting import make_all_plots

        make_all_plots(
            result=result,
            hap_group_df=hap_group_df,
            group_map=group_map,
            gff_file=gff_file,
            plot_formats=plot_formats,
            traits_to_plot=traits_to_plot,
            plot_hap_level=plot_hap_level,
            plot_min_count=plot_min_count,
        )
    return result


def run_analysis(
    vcf_path: str,
    group_file: str,
    outdir: str,
    region: Optional[str] = None,
    region_file: Optional[str] = None,
    **kwargs,
) -> List[HapResult]:
    regions = read_regions(region, region_file)
    results = []
    for reg in regions:
        results.append(run_region_analysis(vcf_path, group_file, reg, outdir, **kwargs))
    return results


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
