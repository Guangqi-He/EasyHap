# EasyHap changelog

## 1.1.0

### Workflow and stability
- Batch analysis now skips regions with fewer than `--min-variants` retained variants instead of terminating the complete run.
- Added `EasyHap.log` for completed, skipped, and failed regions.
- Restored direct handling of unindexed plain-text VCF files with `--vcf-backend auto`.
- Added Matplotlib-version-aware boxplot arguments for compatibility with Matplotlib 3.9 and later.
- The minimum analysis input is now a phased VCF/BCF plus either `--region` or `--region-file`; group, trait, and GFF/GTF files are optional.

### Haplotype and population analysis
- Added group-wise haplotype frequency and unbiased haplotype diversity (Hd).
- Added private-haplotype detection when multiple biological groups are supplied.
- Added FASTA, PHYLIP, NEXUS, and allele-state-map exports for representative haplotypes.

### Trait analysis
- Added haplotype-wise trait summaries and descriptive superior-haplotype candidate ranking.
- Added overall Kruskal-Wallis tests.
- Added pairwise two-sided Mann-Whitney U tests with Benjamini-Hochberg correction.
- Trait figures now use colored boxplots without overlaid raw-data scatter points.

### LD analysis
- Added generalized allele-dosage pairwise LD `r²` calculation.
- Added `*.LD_r2_matrix.tsv` output.
- Added inverted-triangle LD heatmaps with a gene/region track and variant connectors.
- Added configurable LD heatmap colormap through `--ld-cmap` and the GUI.

### Visualization
- Haplotype heatmaps use REF/ALT/missing states with square black-bordered cells and allele labels.
- Added gene-to-variant connectors between the gene structure and haplotype heatmap.
- Gene and haplotype panels are aligned to the same plotting width.
- Gene models are strand-aware and displayed 5'→3' from left to right.
- UTRs are gray, coding/exon blocks are light blue, and introns are black; UTR and coding blocks are non-overlapping.
- Added group haplotype pie charts and stacked barplots with displayed sample counts (`n`).
- Added a shared haplotype color palette for pie charts, stacked barplots, and trait boxplots through `--hap-palette` and the GUI.
- Added configurable REF, ALT, and missing-state colors.
- PDF, SVG, and 600-dpi PNG output are supported.

### Output cleanup
- Removed `*.LD_r2_pairs.tsv`.
- Removed `*.Haplotype_sample.fa`.
- Removed integrated plot/report output; figures are produced as standalone files.

### GUI and documentation
- Reorganized the GUI into Inputs, Analysis, and Visualization tabs.
- Updated GUI labels to match CLI terminology.
- Added GUI controls for haplotype palette and LD heatmap palette.
- Refocused the README on project overview, quick start, core parameters, and output overview; detailed documentation is intended for the GitHub Wiki.

### Packaging
- Corrected PEP 621 dependency metadata in `pyproject.toml`.
- `pysam>=0.21` and `cyvcf2>=0.30` are standard runtime dependencies.
