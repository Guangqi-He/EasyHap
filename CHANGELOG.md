# EasyHap 1.4.0 — 2026-10-08

- Restore a compact gene track and small fixed gap above the haplotype heatmap, shortening guides in both haplotype and cluster views.
- Export only effective SampleGenotypeTokens.tsv and remove the redundant raw table on reruns.
- Apply analyze --plot-min-count to finite phenotype N for summaries, candidates, tests and boxplot annotations; omit singleton candidates and untestable empty comparisons.
- Add aggregate --group and --map-style for coordinate-based combination maps, retaining rare classes in Others and sharing colors with group composition figures.
- Support named/reversed coordinate headers, missing coordinates, and group overrides; expose map inputs in the GUI.
- Use buffered PNG output; validate the generated demonstration images for completeness.

- Compact diploid homozygote labels in genotype figures; preserve heterozygotes, partial calls and polyploid dosage.
- Size each variant column to its longest visible label and use that width for every row in the column; columns containing only single characters remain square, with fixed row height.
- Use actual label widths for figure sizing, including gene/cluster views, without changing exported genotype states.


- Version number set to 1.4.0 as requested, retaining prior format/processing/aggregate fixes.
- Overlap checks default to variants within the supplied interval (POS selection).
- Add `--overlap-flank-bp` and matching GUI input; default 0, symmetric upstream/downstream bp window, no END/SVLEN-based automatic extension.
- Reuse in-region variants at 0 bp; positive values query only the two explicit flanks.
- Flank candidates participate in overlap checks and deletion evidence, while haplotype markers remain in the requested region.
- Export explicit windows in OverlapWindow.tsv, progress logs and manifests.
- Retain interval sweep overlap reporting and all previous genotype/aggregate corrections.
