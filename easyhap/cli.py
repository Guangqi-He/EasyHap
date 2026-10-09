from __future__ import annotations

import argparse
from typing import List, Optional

from . import __version__
from .core import prepare_vcf_tables, run_analysis
from .plotting import (
    DEFAULT_ALT_COLOR,
    DEFAULT_MISSING_COLOR,
    DEFAULT_REF_COLOR,
    DEFAULT_ABSENCE_COLOR,
)


def _csv_list(text: Optional[str]) -> Optional[List[str]]:
    if not text:
        return None
    return [x.strip() for x in text.split(",") if x.strip()]


def _resolve_palette(args):
    ref, alt, missing = args.ref_color, args.alt_color, args.missing_color
    if getattr(args, "palette", None):
        parts = [x.strip() for x in args.palette.split(",") if x.strip()]
        if len(parts) not in {2, 3}:
            raise ValueError(
                "--palette should contain REF,ALT or REF,ALT,MISSING colors"
            )
        ref, alt = parts[0], parts[1]
        if len(parts) == 3:
            missing = parts[2]
    return ref, alt, missing


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="easyhap",
        description=f"EasyHap {__version__}: ploidy-aware regional haplotype analysis, population genetics, trait association and visualization.",
    )
    parser.add_argument("--version", action="version", version=f"EasyHap {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    # 1. prepare
    p_prepare = sub.add_parser(
        "prepare", help="Convert REF/ALT alleles into compact downstream tokens"
    )
    p_prepare.add_argument("--vcf", required=True, help="Phased VCF/VCF.gz/BCF file")
    p_prepare.add_argument("--region", help="Single region, e.g. Chr10:1-500")
    p_prepare.add_argument("--region-file", help="TAB-delimited file: chr start end")
    p_prepare.add_argument(
        "--outdir", default="EasyHap_prepare", help="Output directory [EasyHap_prepare]"
    )
    p_prepare.add_argument(
        "--vcf-backend",
        default="auto",
        choices=["auto", "bcftools", "cyvcf2", "pysam", "plain"],
    )

    # 2. analyze
    p = sub.add_parser(
        "analyze",
        help="Run regional haplotype, population, trait, LD and visualization analyses",
    )
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument(
        "--vcf",
        help="VCF/VCF.gz/BCF genotype file; phased heterozygotes are required only in copy mode",
    )
    src.add_argument(
        "--bfile",
        help="PLINK 1 binary prefix for PREFIX.bed/.bim/.fam; genotype mode only",
    )
    src.add_argument(
        "--pfile",
        help="PLINK 2 binary prefix for PREFIX.pgen/.pvar/.psam; genotype mode only",
    )
    src.add_argument("--hmp", help="HapMap file path (.hmp, .hmp.txt, .hmp.gz)")
    rg = p.add_mutually_exclusive_group(required=False)
    rg.add_argument("--region", help="Single region, e.g. Chr10:1-500")
    rg.add_argument(
        "--region-file", help="TAB-delimited batch region file: chr start end"
    )
    p.add_argument(
        "--outdir", default="EasyHap_results", help="Output directory [EasyHap_results]"
    )

    optional = p.add_argument_group("optional biological inputs")
    optional.add_argument(
        "--group",
        help="TAB-delimited sample-group file: sample, group, and optional latitude/longitude/location columns. If omitted, all VCF samples are assigned to group 'All'. Supports .gz.",
    )
    optional.add_argument(
        "--traits",
        help="TAB-delimited trait table with a header; the first column is sample/accession. Supports .gz.",
    )
    optional.add_argument(
        "--trait-cols",
        help="Comma-separated trait columns to analyze/plot; blank = all trait columns",
    )
    optional.add_argument(
        "--gff",
        help="GFF3/GTF annotation for gene structure visualization and optional variant filtering. Supports .gz.",
    )

    analysis = p.add_argument_group("analysis options")
    analysis.add_argument(
        "--mode",
        default="genotype",
        choices=["genotype", "copy", "inbred", "hybrid"],
        help="Analysis mode: genotype=phase-independent multilocus genotype; copy=copy-resolved haplotypes requiring phased heterozygotes. Legacy aliases: inbred, hybrid [genotype]",
    )
    analysis.add_argument(
        "--hetero-policy",
        default="slash",
        choices=["slash", "iupac", "missing"],
        help="Heterozygous-site encoding used in inbred mode [slash]",
    )
    analysis.add_argument(
        "--min-variants",
        type=int,
        default=2,
        help="Skip regions with fewer variants [2]",
    )
    analysis.add_argument(
        "--threads",
        "-t",
        type=int,
        default=1,
        help="Number of parallel processes/threads for multi-region analysis [1]",
    )
    analysis.add_argument(
        "--variant-list",
        "-l",
        "--var-list",
        dest="variant_list",
        help="Path to variant list file or comma-separated list of variant IDs/positions to restrict the analysis to",
    )
    analysis.add_argument(
        "--max-missing",
        "--missing",
        dest="max_missing",
        type=float,
        default=None,
        help="Maximum allowed missing genotype rate per variant in region (0.0 - 1.0)",
    )
    analysis.add_argument(
        "--min-maf",
        "--maf",
        dest="min_maf",
        type=float,
        default=None,
        help="Minimum minor allele frequency (MAF) per variant in region (0.0 - 0.5)",
    )
    analysis.add_argument(
        "--overlap-flank-bp",
        type=int,
        default=0,
        help="Include variants with POS up to this many bp upstream/downstream of the region in overlap checks; 0 checks only the supplied interval [0]",
    )
    analysis.add_argument(
        "--cluster-threshold",
        type=float,
        default=0.15,
        help="Sequence-distance threshold for haplotype clustering [0.15]",
    )
    analysis.add_argument(
        "--fisher-groups", help="Two group names, e.g. Cultivar,Landrace"
    )
    analysis.add_argument(
        "--fisher-alpha",
        type=float,
        help="Significance threshold for optional Fisher variant filtering",
    )
    analysis.add_argument(
        "--fisher-test",
        default="allele",
        choices=["allele", "carrier"],
        help="Per-allele carrier tests with within-site correction, or legacy any-ALT carrier test [allele]",
    )
    analysis.add_argument(
        "--fisher-adjust",
        default="none",
        choices=["none", "bh"],
        help="Multiple-testing adjustment for Fisher filtering [none]",
    )
    analysis.add_argument(
        "--vcf-backend",
        default="auto",
        choices=["auto", "bcftools", "cyvcf2", "pysam", "plain"],
        help="VCF reader backend [auto]",
    )
    analysis.add_argument(
        "--no-ld",
        action="store_true",
        help="Disable LD r² calculation and LD heatmap output",
    )
    analysis.add_argument(
        "--no-processed",
        action="store_true",
        help="Do not write processed allele/genotype token tables",
    )
    analysis.add_argument(
        "--gene-feature",
        default="all",
        choices=["all", "exon", "intron", "cds", "utr"],
        help="Restrict variants to the selected feature of the primary overlapping gene; non-all choices require --gff [all]",
    )

    plots = p.add_argument_group("visualization options")
    plots.add_argument(
        "--plot", action="store_true", help="Generate standalone figures"
    )
    plots.add_argument(
        "--plot-format",
        default="pdf",
        help="Comma-separated figure formats: pdf,svg,png [pdf]",
    )
    plots.add_argument(
        "--plot-hap-level",
        default="hap",
        choices=["hap", "cluster"],
        help="Plot individual haplotypes or sequence clusters [hap]",
    )
    plots.add_argument(
        "--plot-min-count",
        type=int,
        default=1,
        help="Minimum class count for plots and finite phenotype N for trait tests/candidates (at least 2 for tests/candidates) [1]",
    )
    plots.add_argument(
        "--cell-text",
        default="auto",
        choices=["auto", "always", "never"],
        help=(
            "Allele/base labels inside haplotype heatmap cells: "
            "auto=dynamically show when readable, always=force labels, never=hide labels [auto]"
        ),
    )
    plots.add_argument(
        "--palette",
        help="Custom allele palette: REF,ALT[,MISSING], e.g. '#70AD47,#4472C4,#D9D9D9'",
    )
    plots.add_argument(
        "--hap-palette",
        help="Comma-separated haplotype colors shared by pie, stacked-bar and trait boxplots",
    )
    plots.add_argument(
        "--ld-cmap",
        default="viridis",
        help="Matplotlib colormap for the LD heatmap, e.g. viridis, magma, coolwarm [viridis]",
    )
    plots.add_argument(
        "--ref-color",
        default="#70AD47",
        help="REF cell color in haplotype heatmaps [#70AD47]",
    )
    plots.add_argument(
        "--alt-color",
        default="#4472C4",
        help="ALT cell color in haplotype heatmaps [#4472C4]",
    )
    plots.add_argument(
        "--missing-color",
        default="#D9D9D9",
        help="Missing-genotype cell color in heatmaps [#D9D9D9]",
    )
    plots.add_argument(
        "--absence-color",
        default="#FFFFFF",
        help="Structural/genomic-absence cell color in heatmaps [#FFFFFF]",
    )
    plots.add_argument(
        "--map-style",
        default="auto",
        choices=["auto", "pie", "bar", "both", "none"],
        help="Geographic haplotype overlay when sample-group columns 3-4 contain latitude/longitude: auto=pie [auto]",
    )
    plots.add_argument(
        "--no-network",
        action="store_true",
        help="Disable minimum-spanning haplotype network output",
    )

    # 3. aggregate
    p_agg = sub.add_parser(
        "aggregate",
        help="Multi-locus haplotype combination aggregation and trait association",
    )
    p_agg.add_argument(
        "--indir",
        "--dir",
        "--results-dir",
        dest="indir",
        default="EasyHap_results",
        help="Directory containing easyhap analyze results [EasyHap_results]",
    )
    p_agg.add_argument(
        "--regions",
        default=None,
        help="Comma-separated or space-separated list of 2 or more regions to aggregate, e.g. 'Chr1:100-500,Chr2:200-600'",
    )
    p_agg.add_argument(
        "--region-file",
        default=None,
        help="Path to file containing combinations of regions to aggregate (one combination per line, regions tab-delimited)",
    )
    p_agg.add_argument(
        "--outdir",
        default="EasyHap_aggregate",
        help="Output directory [EasyHap_aggregate]",
    )
    p_agg.add_argument(
        "--traits",
        help="Optional TAB-delimited trait file (if omitted, traits are read from analyze HapGroup.tsv results)",
    )
    p_agg.add_argument(
        "--group",
        help="Optional sample-group TSV; latitude/longitude enable combination maps. Supplied groups override embedded labels. Supports .gz.",
    )
    p_agg.add_argument(
        "--map-style", choices=["auto", "pie", "bar", "both", "none"], default="auto",
        help="Combination maps when --group supplies coordinates: auto=pie [auto]",
    )
    p_agg.add_argument(
        "--trait-cols",
        help="Comma-separated trait columns to analyze; blank = all traits",
    )
    p_agg.add_argument(
        "--min-count",
        type=int,
        default=1,
        help="Minimum sample count per combination to display in plots [1]",
    )
    p_agg.add_argument(
        "--test-min-count",
        type=int,
        default=2,
        help="Minimum nonmissing phenotype N per class for tests [2]",
    )
    p_agg.add_argument(
        "--candidate-min-count",
        type=int,
        default=5,
        help="Minimum phenotype N for candidate ranking [5]",
    )
    p_agg.add_argument(
        "--trait-direction",
        choices=["higher", "lower"],
        default="higher",
        help="Descriptive ranking direction [higher]",
    )
    p_agg.add_argument(
        "--trait-directions",
        help="Per-trait overrides, e.g. Height=higher,Disease=lower",
    )
    p_agg.add_argument(
        "--fdr-scope",
        choices=["run", "trait"],
        default="run",
        help="BH correction across combinations/traits, with separate test families [run]",
    )
    p_agg.add_argument(
        "--phase-link",
        action="store_true",
        help="Connect known phased copies on the same chromosome/PS; unresolved samples are reported",
    )
    p_agg.add_argument(
        "--interaction",
        action="store_true",
        help="Nested categorical OLS two-region interaction test without covariates",
    )
    p_agg.add_argument(
        "--plot",
        action="store_true",
        default=True,
        help="Generate combination trait/group figures and coordinate maps [default: True]",
    )
    p_agg.add_argument(
        "--no-plot",
        action="store_false",
        dest="plot",
        help="Disable all aggregate figures, including maps",
    )
    p_agg.add_argument(
        "--plot-format",
        default="pdf",
        help="Comma-separated figure formats: pdf,svg,png [pdf]",
    )
    p_agg.add_argument(
        "--hap-palette", help="Comma-separated color palette for combinations"
    )

    # 4. heatmap
    p_hm = sub.add_parser(
        "heatmap",
        help="All-sample haplotype color heatmap with population group annotation and clustering",
    )
    p_hm.add_argument(
        "--indir",
        "--dir",
        "--results-dir",
        dest="indir",
        default="EasyHap_results",
        help="Directory containing easyhap analyze results [EasyHap_results]",
    )
    p_hm.add_argument(
        "--region",
        help="Region to plot (e.g. Chr10:1-500); if omitted and only one region exists, it is selected automatically",
    )
    p_hm.add_argument(
        "--region-file",
        help="Path to file containing regions to plot, one region per line",
    )
    p_hm.add_argument(
        "--group",
        help="Optional TAB-delimited sample-group file (if omitted, read from HapGroup.tsv)",
    )
    p_hm.add_argument("--outdir", help="Output directory [default: same as indir]")
    p_hm.add_argument(
        "--plot-format",
        default="pdf",
        help="Comma-separated figure formats: pdf,svg,png [pdf]",
    )
    p_hm.add_argument(
        "--cluster-method",
        default="average",
        choices=["average", "complete", "single", "ward"],
        help="Clustering method for ordering samples within groups [average]",
    )
    p_hm.add_argument("--palette", help="Custom allele palette: REF,ALT[,MISSING]")
    p_hm.add_argument(
        "--ref-color", default=DEFAULT_REF_COLOR, help="REF cell color [#70AD47]"
    )
    p_hm.add_argument(
        "--alt-color", default=DEFAULT_ALT_COLOR, help="ALT cell color [#4472C4]"
    )
    p_hm.add_argument(
        "--het-color", default="#F4B183", help="Heterozygous cell color [#F4B183]"
    )
    p_hm.add_argument(
        "--missing-color",
        default=DEFAULT_MISSING_COLOR,
        help="Missing genotype cell color [#D9D9D9]",
    )
    p_hm.add_argument(
        "--absence-color",
        default=DEFAULT_ABSENCE_COLOR,
        help="Absence cell color [#FFFFFF]",
    )
    p_hm.add_argument("--hap-palette", help="Comma-separated group annotation colors")

    return parser


def main(argv: Optional[List[str]] = None) -> int:
    import multiprocessing

    multiprocessing.freeze_support()
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command == "prepare":
        outputs = prepare_vcf_tables(
            args.vcf, args.outdir, args.region, args.region_file, args.vcf_backend
        )
        print(f"Prepared {len(outputs)} region(s).")
        return 0

    if args.command == "aggregate":
        from .aggregate import run_aggregate

        regions_list = None
        if args.regions:
            regions_list = [
                r.strip()
                for r in args.regions.replace(";", ",").split(",")
                if r.strip()
            ]
            if len(regions_list) < 2:
                regions_list = [r.strip() for r in args.regions.split() if r.strip()]
        if not regions_list and not args.region_file:
            parser.error(
                "Must provide either --regions (at least two regions) or --region-file containing combinations."
            )
        try:
            direction = args.trait_direction
            if args.trait_directions:
                direction = {"*": args.trait_direction}
                for item in args.trait_directions.split(","):
                    name, value = item.split("=", 1)
                    if value.strip() not in {"higher", "lower"}:
                        raise ValueError(
                            "Trait directions require Trait=higher or Trait=lower"
                        )
                    direction[name.strip()] = value.strip()
            run_aggregate(
                indir=args.indir,
                regions=regions_list,
                region_file=args.region_file,
                outdir=args.outdir,
                trait_file=args.traits,
                group_file=args.group,
                map_style=args.map_style,
                trait_cols=_csv_list(args.trait_cols),
                make_plots=args.plot,
                plot_formats=_csv_list(args.plot_format) or ["pdf"],
                min_count=args.min_count,
                test_min_count=args.test_min_count,
                candidate_min_count=args.candidate_min_count,
                trait_direction=direction,
                fdr_scope=args.fdr_scope,
                phase_link=args.phase_link,
                interaction=args.interaction,
                hap_palette=_csv_list(args.hap_palette),
            )
        except Exception as exc:
            parser.error(str(exc))
        return 0

    if args.command == "heatmap":
        from .all_sample_heatmap import run_all_sample_heatmap

        try:
            ref_color, alt_color, missing_color = _resolve_palette(args)
        except ValueError as exc:
            parser.error(str(exc))
        try:
            run_all_sample_heatmap(
                indir=args.indir,
                region=args.region,
                region_file=args.region_file,
                group_file=args.group,
                outdir=args.outdir,
                plot_formats=_csv_list(args.plot_format) or ["pdf"],
                cluster_method=args.cluster_method,
                ref_color=ref_color,
                alt_color=alt_color,
                het_color=args.het_color,
                missing_color=missing_color,
                absence_color=args.absence_color,
                hap_palette=_csv_list(args.hap_palette),
            )
        except Exception as exc:
            parser.error(str(exc))
        return 0

    # analyze command:
    if (args.bfile or args.pfile) and args.mode in {"copy", "hybrid"}:
        parser.error(
            "--bfile and --pfile contain unphased genotypes and are supported only with --mode genotype/inbred"
        )

    if not args.region and not args.region_file and not args.variant_list:
        parser.error("Please provide --region, --region-file, or --variant-list")
    if args.min_variants < 1:
        parser.error("--min-variants must be >=1")
    if args.threads < 1:
        parser.error("--threads must be >= 1")
    if args.max_missing is not None and not (0.0 <= args.max_missing <= 1.0):
        parser.error("--max-missing must be between 0.0 and 1.0")
    if args.min_maf is not None and not (0.0 <= args.min_maf <= 0.5):
        parser.error("--min-maf must be between 0.0 and 0.5")
    if args.overlap_flank_bp < 0:
        parser.error("--overlap-flank-bp must be a nonnegative integer (bp)")
    if args.overlap_flank_bp and not (args.region or args.region_file):
        parser.error("Positive --overlap-flank-bp requires --region or --region-file")

    if args.gene_feature != "all" and not args.gff:
        parser.error(f"--gene-feature {args.gene_feature} requires --gff")
    fisher_group1 = fisher_group2 = None
    if args.fisher_groups:
        parts = [x.strip() for x in args.fisher_groups.split(",") if x.strip()]
        if len(parts) != 2:
            parser.error(
                "--fisher-groups requires exactly two comma-separated group names"
            )
        if not args.group:
            parser.error("--fisher-groups requires --group")
        fisher_group1, fisher_group2 = parts
    if args.fisher_alpha is not None and not args.fisher_groups:
        parser.error("--fisher-alpha requires --fisher-groups")
    try:
        ref_color, alt_color, missing_color = _resolve_palette(args)
    except ValueError as exc:
        parser.error(str(exc))

    if args.plot:
        from .plotting import set_cell_text_mode

        set_cell_text_mode(args.cell_text)

    try:
        results = run_analysis(
            vcf_path=args.vcf,
            bfile=args.bfile,
            pfile=args.pfile,
            hmp=args.hmp,
            group_file=args.group,
            region=args.region,
            region_file=args.region_file,
            outdir=args.outdir,
            mode=args.mode,
            hetero_policy=args.hetero_policy,
            variant_list=args.variant_list,
            max_missing=args.max_missing,
            min_maf=args.min_maf,
            overlap_flank_bp=args.overlap_flank_bp,
            threads=args.threads,
            trait_file=args.traits,
            fisher_group1=fisher_group1,
            fisher_group2=fisher_group2,
            fisher_alpha=args.fisher_alpha,
            fisher_adjust=args.fisher_adjust,
            fisher_test=args.fisher_test,
            cluster_threshold=args.cluster_threshold,
            vcf_backend=args.vcf_backend,
            write_processed=not args.no_processed,
            make_plots=args.plot,
            gff_file=args.gff,
            gene_feature=args.gene_feature,
            plot_formats=_csv_list(args.plot_format) or ["pdf"],
            traits_to_plot=_csv_list(args.trait_cols),
            plot_hap_level=args.plot_hap_level,
            plot_min_count=args.plot_min_count,
            min_variants=args.min_variants,
            ref_color=ref_color,
            alt_color=alt_color,
            missing_color=missing_color,
            absence_color=args.absence_color,
            make_ld=not args.no_ld,
            make_network=not args.no_network,
            map_style=args.map_style,
            hap_palette=_csv_list(args.hap_palette),
            ld_cmap=args.ld_cmap,
        )
    except (ValueError, RuntimeError, OSError) as exc:
        parser.error(str(exc))

    print(
        f"Finished {len(results)} region(s). See {args.outdir}/EasyHap.log for processed/skipped regions."
    )
    for r in results:
        print(f"[{r.region.vcf_label}] {r.output_prefix}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
