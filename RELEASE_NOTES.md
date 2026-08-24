# EasyHap 1.1.0 release notes

EasyHap 1.1.0 extends the regional haplotype workflow with more robust batch processing, population statistics, trait association, LD analysis, and publication-ready standalone figures.

## Highlights

- Batch regions with insufficient retained variants are skipped and recorded in `EasyHap.log`.
- Group-wise haplotype frequency, haplotype diversity, and private haplotypes are reported when appropriate.
- Trait analysis includes Kruskal-Wallis testing and BH-adjusted pairwise Mann-Whitney U tests.
- LD is reported as a pairwise `r²` matrix and can be visualized as an inverted-triangle heatmap with gene-to-variant connectors.
- Gene-haplotype plots use aligned gene/heatmap widths, compact square cells, black cell borders, and strand-aware gene structure.
- Group composition is available as both stacked barplots and pie charts with displayed sample counts.
- Trait plots use haplotype-colored boxplots without raw-data scatter overlays.
- Haplotype colors used across pie, stacked-bar, and trait plots can be changed with `--hap-palette` or in the GUI.
- LD heatmap colors can be changed with `--ld-cmap` or in the GUI.
- Integrated plot/report output, `*.LD_r2_pairs.tsv`, and `*.Haplotype_sample.fa` have been removed to keep the output set focused.

For usage details, see the GitHub Wiki: https://github.com/Guangqi-He/EasyHap/wiki
