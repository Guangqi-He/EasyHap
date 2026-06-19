from __future__ import annotations

import argparse
import sys
from typing import List, Optional

from . import __version__
from .core import prepare_vcf_tables, run_analysis


def _csv_list(text: Optional[str]) -> Optional[List[str]]:
    if not text:
        return None
    return [x.strip() for x in text.split(",") if x.strip()]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="easyhap",
        description="EasyHap 1.0: haplotype analysis for phased VCF regions from fungi, plants, and animals.",
    )
    parser.add_argument("--version", action="version", version=f"EasyHap {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    p_prepare = sub.add_parser("prepare", help="Convert REF/ALT alleles into compact downstream tokens")
    p_prepare.add_argument("--vcf", required=True, help="Phased VCF/VCF.gz/BCF file")
    p_prepare.add_argument("--region", help="Single region, e.g. Chr10:1-500")
    p_prepare.add_argument("--region-file", help="TAB-delimited file: chr start end")
    p_prepare.add_argument("--outdir", required=True, help="Output directory")
    p_prepare.add_argument("--vcf-backend", default="auto", choices=["auto", "cyvcf2", "pysam", "plain"], help="VCF reader backend")

    p = sub.add_parser("analyze", help="Run haplotype summarization, filtering, alignments, and optional plots")
    p.add_argument("--vcf", required=True, help="Phased VCF/VCF.gz/BCF file. Indexed files allow fast region access.")
    p.add_argument("--group", required=True, help="TAB-delimited sample group file: sample group")
    p.add_argument("--region", help="Single region, e.g. Chr10:1-500")
    p.add_argument("--region-file", help="TAB-delimited file with no header: chr start end")
    p.add_argument("--outdir", required=True, help="Output directory")
    p.add_argument("--mode", default="inbred", choices=["inbred", "hybrid"], help="inbred: genotype-level haplotypes; hybrid: phased copy-level haplotypes")
    p.add_argument("--hetero-policy", default="slash", choices=["slash", "iupac", "missing"], help="How to encode heterozygous sites in inbred mode")
    p.add_argument("--traits", help="Optional TAB-delimited trait table with header. First column should be sample/accession.")
    p.add_argument("--trait-cols", help="Comma-separated trait columns to plot")
    p.add_argument("--fisher-groups", help="Two group names for Fisher filtering, e.g. Cultivar,Landrace")
    p.add_argument("--fisher-alpha", type=float, help="P-value or adjusted P-value cutoff for Fisher filtering")
    p.add_argument("--fisher-adjust", default="none", choices=["none", "bh"], help="Multiple-testing adjustment for Fisher filtering")
    p.add_argument("--cluster-threshold", type=float, default=0.15, help="Maximum pairwise normalized Hamming distance for connected-component haplotype clustering")
    p.add_argument("--vcf-backend", default="auto", choices=["auto", "cyvcf2", "pysam", "plain"], help="VCF reader backend")
    p.add_argument("--no-processed", action="store_true", help="Do not write processed variant/genotype token tables")
    p.add_argument("--plot", action="store_true", help="Generate haplotype heatmap, group proportions, and trait plots")
    p.add_argument("--gff", help="Optional GFF3/GTF gene annotation for gene-haplotype plot")
    p.add_argument("--plot-format", default="pdf", help="Comma-separated output image formats: pdf,svg,png")
    p.add_argument("--plot-hap-level", default="hap", choices=["hap", "cluster"], help="Plot by raw haplotype/diplotype labels (hap) or clustered haplotype/diplotype labels (cluster)")
    p.add_argument("--plot-min-count", type=int, default=1, help="Minimum total sample/accession count required for a haplotype/cluster class to be displayed in plots only")
    return parser


def main(argv: Optional[List[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command == "prepare":
        outputs = prepare_vcf_tables(
            vcf_path=args.vcf,
            outdir=args.outdir,
            region=args.region,
            region_file=args.region_file,
            vcf_backend=args.vcf_backend,
        )
        print(f"Prepared {len(outputs)} region(s).")
        for vp, gp in outputs:
            print(vp)
            print(gp)
        return 0

    fisher_group1 = fisher_group2 = None
    if args.fisher_groups:
        parts = [x.strip() for x in args.fisher_groups.split(",") if x.strip()]
        if len(parts) != 2:
            parser.error("--fisher-groups should contain exactly two comma-separated group names")
        fisher_group1, fisher_group2 = parts
    plot_formats = _csv_list(args.plot_format) or ["pdf"]
    results = run_analysis(
        vcf_path=args.vcf,
        group_file=args.group,
        region=args.region,
        region_file=args.region_file,
        outdir=args.outdir,
        mode=args.mode,
        hetero_policy=args.hetero_policy,
        trait_file=args.traits,
        fisher_group1=fisher_group1,
        fisher_group2=fisher_group2,
        fisher_alpha=args.fisher_alpha,
        fisher_adjust=args.fisher_adjust,
        cluster_threshold=args.cluster_threshold,
        vcf_backend=args.vcf_backend,
        write_processed=not args.no_processed,
        make_plots=args.plot,
        gff_file=args.gff,
        plot_formats=plot_formats,
        traits_to_plot=_csv_list(args.trait_cols),
        plot_hap_level=args.plot_hap_level,
        plot_min_count=args.plot_min_count,
    )
    print(f"Finished {len(results)} region(s).")
    for r in results:
        print(f"[{r.region.vcf_label}]")
        print(f"  HapSummary: {r.hap_summary_path}")
        print(f"  HapGroup:   {r.hap_group_path}")
        print(f"  Prefix:     {r.output_prefix}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
