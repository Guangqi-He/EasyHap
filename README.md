# HaploFlow

HaploFlow is a first working prototype for **regional haplotype analysis from phased VCF files** in animal and plant population resequencing projects.

It focuses on candidate regions or candidate genes. It can extract phased genotypes, filter variants, optionally select variants by Fisher exact test between two groups, build haplotypes, count haplotypes by group, export network-compatible sequence files, and generate basic plots.

## Main functions

1. Read phased VCF/VCF.GZ plus sample group information.
2. Accept a single region such as `Chr01:100000-150000` or a BED-like file of multiple regions.
3. Filter variants by MAF, missing rate, QUAL, PASS status, SNP/INDEL, biallelic/multiallelic status.
4. Optionally compare two groups and compute per-site Fisher exact test, BH-FDR, odds ratio and allele-frequency difference.
5. Construct haplotypes from phased genotypes.
6. Count each haplotype in each sample group and list the sample names.
7. Name haplotypes as `Hap01`, `Hap02`, etc. by descending frequency.
8. Optionally cluster haplotypes by Hamming distance.
9. Export FASTA, NEXUS and PHYLIP files for downstream haplotype-network analysis.
10. Plot gene structure + haplotype allele matrix, group haplotype composition, and trait boxplots.
11. Provide a simple Windows-compatible GUI through `tkinter`.

## Important biological assumptions

For a diploid phased genotype such as:

```text
0|1
```

HaploFlow treats the sample as two haplotype copies:

```text
sample_H1 = allele 0
sample_H2 = allele 1
```

For phenotype/trait plotting, the safest default category is the **sample-level diplotype**, for example:

```text
Hap01/Hap03
```

This avoids incorrectly assigning one phenotype value to only one chromosome copy. For highly homozygous inbred materials, `Hap01/Hap01` can usually be interpreted as the corresponding sample haplotype.

## Installation

From the project directory:

```bash
pip install -e .
```

For a clean conda environment:

```bash
conda create -n haploflow python=3.10 -y
conda activate haploflow
pip install -e .
```

## Required input files

### 1. Phased VCF

The VCF should contain phased genotypes such as:

```text
0|0
0|1
1|0
1|1
```

Unphased genotypes such as `0/1` are treated as missing by default. To keep them, use:

```bash
--keep-unphased
```

### 2. Sample group file

Tab- or comma-delimited file:

```text
sample  group
S1      wild
S2      wild
S3      cultivar
S4      cultivar
```

### 3. Region input

Single region:

```bash
--region Chr01:100000-150000
```

or BED-like file:

```text
Chr01   99999   150000   GeneA
Chr05   200000  250000   GeneB
```

The BED file uses 0-based start and 1-based end, as standard BED format.

### 4. Optional trait file

```text
sample  plant_height  seed_weight
S1      82.1          21.5
S2      76.4          19.8
S3      95.2          25.1
S4      88.9          23.7
```

## Basic usage

```bash
haploflow run \
  --vcf examples/demo.vcf \
  --group examples/sample_group.tsv \
  --region Chr01:100-250 \
  --out demo_out
```

## With MAF and missing filtering

```bash
haploflow run \
  --vcf input.vcf.gz \
  --group sample_group.tsv \
  --region Chr01:100000-150000 \
  --maf 0.05 \
  --max-missing 0.2 \
  --out result
```

## With Fisher exact test filtering

```bash
haploflow run \
  --vcf input.vcf.gz \
  --group sample_group.tsv \
  --region Chr01:100000-150000 \
  --compare wild,cultivar \
  --fisher-p 0.01 \
  --min-delta-af 0.3 \
  --out result
```

If `--compare` is provided but no `--fisher-p`, `--fisher-q` or `--min-delta-af` threshold is given, HaploFlow calculates Fisher statistics but keeps all sites that pass the basic filters.

## With gene annotation and trait plotting

```bash
haploflow run \
  --vcf input.vcf.gz \
  --group sample_group.tsv \
  --region Chr01:100000-150000 \
  --gff genome.gff3 \
  --trait traits.tsv \
  --trait-name plant_height \
  --out result
```

## Windows GUI

After installation, run:

```bash
haploflow-gui
```

or:

```bash
python -m haploflow.gui
```

The GUI is a thin wrapper around the command-line workflow, so the GUI and CLI use the same core code.

## Output files

For each region, HaploFlow creates a subdirectory containing:

```text
variants.tsv
```

Per-site variant statistics, including MAF, missing rate, AF, Fisher p value, FDR and delta AF when applicable.

```text
genotype_matrix.tsv
```

Sample-by-variant genotype matrix.

```text
haplotype_matrix.tsv
```

One row per haplotype copy, such as `S1_H1` and `S1_H2`.

```text
haplotype_summary.tsv
```

Unique haplotype summary, including counts by group and sample names.

```text
sample_haplotype_assignment.tsv
```

Sample-level assignment of `hap1`, `hap2`, `diplotype`, and `primary_haplotype`.

```text
haplotypes.unique.fa
haplotypes.samples.fa
haplotypes.unique.nex
haplotypes.unique.phy
popart_group_counts.tsv
```

Sequence files for downstream network analysis.

```text
haplotype_structure.png
haplotype_group_stacked_bar.png
haplotype_group_pies.png
trait_<trait_name>_boxplot.png
```

Basic plots.

## Notes on SNPs and INDELs

By default, HaploFlow keeps only biallelic SNPs. This is recommended for haplotype-network export because sequence alignments require equal-length character strings.

To include INDELs:

```bash
--include-indels
```

For non-SNP variants, HaploFlow encodes alleles as single-character allele indices such as `0` and `1` to keep the alignment length consistent. For strict DNA-only network software, using only SNPs is safer.

## Limitations of this first prototype

This is a functional prototype, not yet a fully optimized production release.

Current limitations:

1. The built-in VCF reader streams through VCF/VCF.GZ files. This is robust but slow for very large files. A future version should add a `pysam`/`htslib` indexed backend for fast random access.
2. Phase-block checking by `PS` tag is not yet implemented.
3. The gene-structure plot is intentionally simple and may need visual tuning for complex gene models.
4. The trait plot uses sample-level diplotypes by default, which is statistically safer than assigning phenotype values directly to haplotype copies.
5. The GUI is a basic wrapper for Windows users and can be improved later with better parameter validation and result preview.

## Suggested next development steps

1. Add indexed VCF reading through `pysam` or `cyvcf2`.
2. Add phase-block checking and automatic splitting by `PS` block.
3. Add HTML report output.
4. Add LD and local diversity statistics.
5. Add support for SV/TE presence-absence haplotypes.
6. Add a richer Qt-based GUI if the software will be distributed to non-command-line users.

