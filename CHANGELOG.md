# EasyHap changelog

## 1.2.0

### Geographic and network visualization
- Extended the sample-group file to optionally accept latitude, longitude, and location columns while retaining compatibility with the historical two-column format.
- Added geographic haplotype maps with pie or stacked-bar overlays through `--map-style auto|pie|bar|both|none`.
- Added minimum-spanning haplotype networks with node size proportional to haplotype copy count, group composition shown as node pies, and mutational distance on edges.
- Added `*.HaplotypeNetworkEdges.tsv` for the network topology and edge distances.

### Gene-structure variant filtering
- Added `--gene-feature all|exon|intron|cds|utr` to determine which variants participate in haplotype reconstruction and downstream analyses.
- Non-`all` feature filtering requires `--gff`; introns are inferred from the primary overlapping gene as gene sequence not covered by exons.
- Added `*.GeneFeatureFilter.tsv` to report original and retained variant counts and the intervals used for filtering.

### GUI and examples
- Added GUI controls for geographic map style, haplotype network output, and gene-feature filtering.
- Updated the demonstration `sample_group.tsv` with reproducible geographic coordinates and location labels.

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
