# EasyHap 1.3.0 reviewer-response source revision

EasyHap 1.3.0 incorporates the P0-P2 reviewer-response changes and the PLINK/Windows usability and performance improvements listed below.

## P0: methodological correctness

1. Phase-independent `genotype` mode and copy-resolved `copy` mode were made explicit; legacy `inbred`/`hybrid` aliases remain accepted.
2. Copy mode now validates phase and rejects unphased heterozygous genotypes instead of treating the textual allele order of `0/1` as chromosome-copy phase.
3. Genotype mode canonicalizes unphased allele order (`0/1` and `1/0` become the same multilocus state).
4. VCF spanning-deletion `*` is represented as `ABS`, distinct from ordinary missing genotype `N`.
5. Explicit deletion intervals can mask downstream markers as `ABS` when the genotype supports genomic absence. Ordinary `./.` is never converted to absence without structural evidence.
6. Multiallelic indel/SV tokens preserve VCF ALT identity, preventing equal-length but sequence-distinct ALT alleles from collapsing.

## P1: SV/overlap handling and validation

1. Haplotype distance treats deletion-derived downstream `ABS` states as non-comparable, preventing a single large deletion from being counted repeatedly through all covered SNPs; the deletion event itself remains informative.
2. `*.VariantOverlap.tsv` reports same-start, nested/overlapping, and deletion-contained records and the action taken.
3. Heatmaps distinguish REF, ALT, missing GT, and genomic absence; `--absence-color` was added.
4. Synthetic reviewer-validation examples were added under `examples/reviewer_validation/`.

## P2: usability/input compatibility

1. Added PLINK 1 binary `.bed/.bim/.fam` input through `--bfile PREFIX` in genotype mode.
2. README now explains when phased input is and is not required.
3. README includes current example commands for Beagle 5.5, WhatsHap (diploid and polyploid), and SHAPEIT5, and states that EasyHap itself does not perform phasing.

## EasyHap 1.3 PLINK/Windows additions

1. Added a GUI genotype-source selector for VCF/BCF versus PLINK 1 BED/BIM/FAM, enabling the Windows GUI executable to expose the PLINK reader already available in the analysis core.
2. Browsing a `.bed`, `.bim`, or `.fam` file in the GUI automatically derives the common PLINK prefix and checks that all three files are present.
3. PLINK input is guarded against `copy`/`hybrid` mode because BED contains unphased genotypes.
4. Replaced whole-BED scanning for every region with a chromosome/position index built from BIM plus fixed-width SNP-major BED `seek()` access.
5. Batch analyses reuse one PLINK reader/index across regions instead of repeatedly reparsing BIM/FAM.
6. PLINK reader input now accepts a bare prefix or a `.bed`/`.bim`/`.fam` path and normalizes it to the common prefix.
7. Package, CLI, GUI, README, VERSION, and project metadata were advanced to 1.3.0.

## Verification performed for 1.3.0

- Python compilation/import checks passed for the package, CLI, and GUI modules.
- A synthetic PLINK BED/BIM/FAM dataset with deliberately unsorted BIM coordinates confirmed that indexed regional access returns the correct variants in genomic order.
- End-to-end batch PLINK analysis passed for multiple regions, including use of a `.bed` path instead of a bare prefix.
- Batch instrumentation confirmed that the PLINK reader/BIM index is constructed once and reused across regions.
- The CLI correctly rejects PLINK input with `copy` mode.
- A plain-VCF genotype-mode smoke test passed after the PLINK/core changes, providing a regression check for the existing VCF path.
- The project successfully built an EasyHap 1.3.0 wheel using the local build environment (`--no-build-isolation`).

## Reviewer-validation GFF3 and plotting examples

- Added one matching GFF3 annotation for every reviewer-validation dataset.
- Each GFF3 contains gene, mRNA, exon, CDS, 5-prime UTR, and 3-prime UTR records for gene-structure visualization; the allotetraploid annotation contains separate models on ChrA01 and ChrB01.
- Updated reviewer-validation commands to include `--group`, `--gff`, `--plot`, and `--plot-format pdf,png`.
- Updated the root README plotting examples to use the same PDF+PNG convention.
- Revalidated all synthetic targets: each still resolves to exactly six abundant haplotypes in its intended mode, with each haplotype present in >5 accessions and at >5% frequency.
- Executed the documented plotting commands and confirmed PDF and PNG output, including the gene-structure/haplotype figure.
