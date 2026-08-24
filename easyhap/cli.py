from __future__ import annotations

import argparse
from typing import List, Optional

from . import __version__
from .core import prepare_vcf_tables, run_analysis


def _csv_list(text: Optional[str]) -> Optional[List[str]]:
    if not text:
        return None
    return [x.strip() for x in text.split(",") if x.strip()]


def _resolve_palette(args):
    ref, alt, missing = args.ref_color, args.alt_color, args.missing_color
    if args.palette:
        parts = [x.strip() for x in args.palette.split(",") if x.strip()]
        if len(parts) not in {2, 3}:
            raise ValueError("--palette should contain REF,ALT or REF,ALT,MISSING colors")
        ref, alt = parts[0], parts[1]
        if len(parts) == 3:
            missing = parts[2]
    return ref, alt, missing


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="easyhap",
        description="EasyHap 1.1.0: ploidy-aware regional haplotype analysis, population genetics, trait association and visualization.",
    )
    parser.add_argument("--version", action="version", version=f"EasyHap {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    p_prepare = sub.add_parser("prepare", help="Convert REF/ALT alleles into compact downstream tokens")
    p_prepare.add_argument("--vcf", required=True, help="Phased VCF/VCF.gz/BCF file")
    p_prepare.add_argument("--region", help="Single region, e.g. Chr10:1-500")
    p_prepare.add_argument("--region-file", help="TAB-delimited file: chr start end")
    p_prepare.add_argument("--outdir", default="EasyHap_prepare", help="Output directory [EasyHap_prepare]")
    p_prepare.add_argument("--vcf-backend", default="auto", choices=["auto", "cyvcf2", "pysam", "plain"])

    p = sub.add_parser("analyze", help="Run regional haplotype, population, trait, LD and visualization analyses")
    p.add_argument("--vcf", required=True, help="Phased VCF/VCF.gz/BCF file")
    rg = p.add_mutually_exclusive_group(required=True)
    rg.add_argument("--region", help="Single region, e.g. Chr10:1-500")
    rg.add_argument("--region-file", help="TAB-delimited batch region file: chr start end")
    p.add_argument("--outdir", default="EasyHap_results", help="Output directory [EasyHap_results]")

    optional = p.add_argument_group("optional biological inputs")
    optional.add_argument("--group", help="TAB-delimited sample-group file without a header. If omitted, all VCF samples are assigned to group 'All'.")
    optional.add_argument("--traits", help="TAB-delimited trait table with a header; the first column is sample/accession")
    optional.add_argument("--trait-cols", help="Comma-separated trait columns to analyze/plot; blank = all trait columns")
    optional.add_argument("--gff", help="GFF3/GTF annotation for strand-aware gene structure visualization")

    analysis = p.add_argument_group("analysis options")
    analysis.add_argument("--mode", default="inbred", choices=["inbred", "hybrid"], help="Haplotype reconstruction mode [inbred]")
    analysis.add_argument("--hetero-policy", default="slash", choices=["slash", "iupac", "missing"], help="Heterozygous-site encoding used in inbred mode [slash]")
    analysis.add_argument("--min-variants", type=int, default=2, help="Skip regions with fewer variants [2]")
    analysis.add_argument("--cluster-threshold", type=float, default=0.15, help="Sequence-distance threshold for haplotype clustering [0.15]")
    analysis.add_argument("--fisher-groups", help="Two group names, e.g. Cultivar,Landrace")
    analysis.add_argument("--fisher-alpha", type=float, help="Significance threshold for optional Fisher variant filtering")
    analysis.add_argument("--fisher-adjust", default="none", choices=["none", "bh"], help="Multiple-testing adjustment for Fisher filtering [none]")
    analysis.add_argument("--vcf-backend", default="auto", choices=["auto", "cyvcf2", "pysam", "plain"], help="VCF reader backend [auto]")
    analysis.add_argument("--no-ld", action="store_true", help="Disable LD r² calculation and LD heatmap output")
    analysis.add_argument("--no-processed", action="store_true", help="Do not write processed allele/genotype token tables")

    plots = p.add_argument_group("visualization options")
    plots.add_argument("--plot", action="store_true", help="Generate standalone figures")
    plots.add_argument("--plot-format", default="pdf", help="Comma-separated figure formats: pdf,svg,png [pdf]")
    plots.add_argument("--plot-hap-level", default="hap", choices=["hap", "cluster"], help="Plot individual haplotypes or sequence clusters [hap]")
    plots.add_argument("--plot-min-count", type=int, default=1, help="Minimum displayed haplotype/cluster count [1]")
    plots.add_argument("--palette", help="Custom allele palette: REF,ALT[,MISSING], e.g. '#70AD47,#4472C4,#D9D9D9'")
    plots.add_argument("--hap-palette", help="Comma-separated haplotype colors shared by pie, stacked-bar and trait boxplots")
    plots.add_argument("--ld-cmap", default="viridis", help="Matplotlib colormap for the LD heatmap, e.g. viridis, magma, coolwarm [viridis]")
    plots.add_argument("--ref-color", default="#70AD47", help="REF cell color in haplotype heatmaps [#70AD47]")
    plots.add_argument("--alt-color", default="#4472C4", help="ALT cell color in haplotype heatmaps [#4472C4]")
    plots.add_argument("--missing-color", default="#D9D9D9", help="Missing-data cell color in haplotype heatmaps [#D9D9D9]")
    return parser


def main(argv: Optional[List[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command == "prepare":
        outputs = prepare_vcf_tables(args.vcf, args.outdir, args.region, args.region_file, args.vcf_backend)
        print(f"Prepared {len(outputs)} region(s).")
        return 0

    if args.min_variants < 1:
        parser.error("--min-variants must be >=1")
    fisher_group1 = fisher_group2 = None
    if args.fisher_groups:
        parts = [x.strip() for x in args.fisher_groups.split(",") if x.strip()]
        if len(parts) != 2:
            parser.error("--fisher-groups requires exactly two comma-separated group names")
        if not args.group:
            parser.error("--fisher-groups requires --group")
        fisher_group1, fisher_group2 = parts
    if args.fisher_alpha is not None and not args.fisher_groups:
        parser.error("--fisher-alpha requires --fisher-groups")
    try:
        ref_color, alt_color, missing_color = _resolve_palette(args)
    except ValueError as exc:
        parser.error(str(exc))

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
        plot_formats=_csv_list(args.plot_format) or ["pdf"],
        traits_to_plot=_csv_list(args.trait_cols),
        plot_hap_level=args.plot_hap_level,
        plot_min_count=args.plot_min_count,
        min_variants=args.min_variants,
        ref_color=ref_color,
        alt_color=alt_color,
        missing_color=missing_color,
        make_ld=not args.no_ld,
        hap_palette=_csv_list(args.hap_palette),
        ld_cmap=args.ld_cmap,
    )
    print(f"Finished {len(results)} region(s). See {args.outdir}/EasyHap.log for processed/skipped regions.")
    for r in results:
        print(f"[{r.region.vcf_label}] {r.output_prefix}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
