"""All-sample haplotype color heatmap with population group annotation and clustering."""

from __future__ import annotations

import glob
import os
import re
from typing import Dict, List, Optional, Sequence, Tuple

import matplotlib

matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
from scipy.cluster import hierarchy
from scipy.spatial import distance

from .io_utils import ensure_dir, read_group_metadata, sanitize_filename
from .plotting import (
    DEFAULT_ALT_COLOR,
    DEFAULT_MISSING_COLOR,
    DEFAULT_REF_COLOR,
    DEFAULT_ABSENCE_COLOR,
    _hap_color_map,
    _natural_key,
    _save_formats,
)


def _find_target_prefix(
    indir: str, region_str: Optional[str] = None
) -> Tuple[str, str]:
    """Find the output prefix and label for the requested region in indir."""
    group_files = glob.glob(os.path.join(indir, "**", "*.HapGroup.tsv"), recursive=True)
    if not group_files:
        group_files = glob.glob(os.path.join(indir, "*.HapGroup.tsv"))

    if not group_files:
        raise FileNotFoundError(f"No .HapGroup.tsv files found in directory: {indir}")

    if region_str:
        from .aggregate import find_region_hap_group_file

        path = find_region_hap_group_file(indir, region_str)
        prefix = path[: -len(".HapGroup.tsv")]
        return prefix, os.path.basename(prefix)
    if len(group_files) > 1:
        raise ValueError("Multiple regions exist; specify --region or --region-file")

    first = sorted(group_files)[0]
    prefix = first[: -len(".HapGroup.tsv")]
    return prefix, os.path.basename(prefix)


def _load_sample_variant_matrix(
    prefix: str,
) -> Tuple[List[str], List[int], np.ndarray, pd.DataFrame]:
    """Load sample by variant matrix and return (samples, positions, matrix_codes, sample_info_df)."""
    token_file = prefix + ".SampleGenotypeTokens.tsv"
    hap_group_file = prefix + ".HapGroup.tsv"
    hap_summary_file = prefix + ".HapSummary.tsv"
    var_file = prefix + ".ProcessedVariants.tsv"

    hap_group_df = pd.read_csv(
        hap_group_file, sep="\t", dtype=str, keep_default_na=False
    )
    if hap_group_df.Accession.duplicated().any():
        raise ValueError("Duplicate sample IDs in heatmap HapGroup")
    if "Type" not in hap_group_df:
        hap_group_df["Type"] = "All"
    if not os.path.isfile(var_file):
        var_file = prefix + ".AnalysisVariants.tsv"
    vdf = (
        pd.read_csv(var_file, sep="\t", dtype=str, keep_default_na=False)
        if os.path.isfile(var_file)
        else None
    )
    if os.path.isfile(token_file):
        tokens = pd.read_csv(token_file, sep="\t", dtype=str, keep_default_na=False)
        samples = [c for c in tokens if c not in {"CHROM", "POS", "ID", "REF", "ALT"}]
        positions = tokens.POS.astype(int).tolist()
        values = tokens[samples].to_numpy().T
    else:
        summary = pd.read_csv(
            hap_summary_file, sep="\t", dtype=str, keep_default_na=False
        )
        columns = [
            c for c in summary if c not in {"Hap", "ClusterID", "Accession", "Number"}
        ]
        position_row = summary[summary.Hap == "POS"]
        positions = position_row[columns].iloc[0].astype(int).tolist()
        hap_states = {
            row.Hap: row[columns].to_numpy()
            for _, row in summary.iterrows()
            if row.Hap not in {"CHR", "POS"}
        }
        samples = hap_group_df.Accession.tolist()
        values = []
        for hap in hap_group_df.Hap:
            copies = hap.split("_")
            if any(h not in hap_states for h in copies):
                raise ValueError(f"HapGroup category {hap} is absent from HapSummary")
            values.append(
                [
                    "/".join(str(hap_states[h][j]) for h in copies)
                    for j in range(len(columns))
                ]
            )
        values = np.asarray(values)
    if vdf is not None and len(vdf) != len(positions):
        raise ValueError("Variant metadata and heatmap marker counts disagree")
    mat = np.zeros(values.shape, dtype=int)
    for j in range(len(positions)):
        ref = str(vdf.iloc[j].AlleleTokens).split(",")[0] if vdf is not None else None
        for i in range(len(samples)):
            ts = re.split(r"[|/]", str(values[i, j]))
            if any(t in {"N", ".", "", "NA", "<NA>"} for t in ts):
                mat[i, j] = 3
            elif "ABS" in ts:
                mat[i, j] = 4
            elif len(set(ts)) > 1:
                mat[i, j] = 2
            elif ref is None:
                # Legacy no-processed output cannot establish reference identity.
                mat[i, j] = 1
            else:
                mat[i, j] = 0 if ts[0] == ref else 1
    return samples, positions, mat, hap_group_df


def _generate_single_heatmap(
    indir: str,
    region: Optional[str] = None,
    group_file: Optional[str] = None,
    outdir: Optional[str] = None,
    plot_formats: Sequence[str] = ("pdf",),
    cluster_method: str = "average",
    ref_color: str = DEFAULT_REF_COLOR,
    alt_color: str = DEFAULT_ALT_COLOR,
    het_color: str = "#F4B183",
    missing_color: str = DEFAULT_MISSING_COLOR,
    absence_color: str = DEFAULT_ABSENCE_COLOR,
    hap_palette: Optional[Sequence[str]] = None,
) -> str:
    """Generate all-sample haplotype color heatmap for a single region."""
    if outdir is None:
        outdir = indir
    ensure_dir(outdir)

    prefix, label = _find_target_prefix(indir, region)
    samples, pos_list, mat, hap_group_df = _load_sample_variant_matrix(prefix)

    if group_file:
        group_meta = read_group_metadata(group_file)
        group_map = dict(
            zip(group_meta["Accession"].astype(str), group_meta["Type"].astype(str))
        )
    else:
        group_map = dict(
            zip(hap_group_df["Accession"].astype(str), hap_group_df["Type"].astype(str))
        )

    sample_groups = [group_map.get(s, "NA") for s in samples]
    sample_groups = [
        g if g not in {"", "nan", "None", "<NA>"} else "NA" for g in sample_groups
    ]

    unique_groups = sorted(set(sample_groups), key=_natural_key)
    group_colors = _hap_color_map(unique_groups, hap_palette)

    ordered_indices = []
    group_slices = []
    current_start = 0

    for g in unique_groups:
        idx_in_g = [i for i, grp in enumerate(sample_groups) if grp == g]
        if len(idx_in_g) > 2:
            g_mat = mat[idx_in_g, :]
            pdist = distance.pdist(g_mat, metric="hamming")
            if np.all(np.isfinite(pdist)) and not np.all(pdist == 0):
                link = hierarchy.linkage(pdist, method=cluster_method)
                order = hierarchy.leaves_list(link)
                sorted_g_idx = [idx_in_g[o] for o in order]
            else:
                sorted_g_idx = idx_in_g
        else:
            sorted_g_idx = idx_in_g

        ordered_indices.extend(sorted_g_idx)
        group_slices.append((g, current_start, current_start + len(sorted_g_idx)))
        current_start += len(sorted_g_idx)

    ordered_mat = mat[ordered_indices, :]
    ordered_samples = [samples[i] for i in ordered_indices]
    ordered_groups = [sample_groups[i] for i in ordered_indices]

    n_samples, n_vars = ordered_mat.shape

    # Proportional, compact figure dimensions
    fig_w = max(5.5, min(14.0, 2.2 + 0.16 * n_vars))
    fig_h = max(4.0, min(12.0, 1.8 + 0.015 * n_samples))
    fig = plt.figure(figsize=(fig_w, fig_h))

    gs = fig.add_gridspec(1, 2, width_ratios=[0.05, 0.95], wspace=0.02)
    ax_grp = fig.add_subplot(gs[0])
    ax_heat = fig.add_subplot(gs[1])

    # 1. Group annotation bar
    grp_to_int = {g: i for i, g in enumerate(unique_groups)}
    grp_codes = np.array([[grp_to_int[g]] for g in ordered_groups])
    grp_cmap = ListedColormap([group_colors[g] for g in unique_groups])

    ax_grp.pcolormesh(
        np.array([0, 1]),
        np.arange(n_samples + 1),
        grp_codes,
        cmap=grp_cmap,
        vmin=-0.5,
        vmax=len(unique_groups) - 0.5,
        shading="flat",
        edgecolors="none",
    )
    ax_grp.set_xlim(0, 1)
    ax_grp.set_ylim(n_samples, 0)
    ax_grp.set_xticks([])
    ax_grp.set_yticks([])
    ax_grp.set_ylabel("Samples", fontsize=11)

    for g, s, e in group_slices:
        ax_grp.plot([0, 1], [e, e], color="white", linewidth=1.0, zorder=5)

    # 2. Main Genotype Heatmap (colors only, no text letters)
    gt_cmap = ListedColormap(
        [ref_color, alt_color, het_color, missing_color, absence_color]
    )
    ax_heat.pcolormesh(
        np.arange(n_vars + 1),
        np.arange(n_samples + 1),
        ordered_mat,
        cmap=gt_cmap,
        vmin=-0.5,
        vmax=4.5,
        shading="flat",
        edgecolors="none",
        linewidth=0,
    )
    ax_heat.set_xlim(0, n_vars)
    ax_heat.set_ylim(n_samples, 0)

    ax_heat.set_xticks(np.arange(n_vars) + 0.5)
    ax_heat.set_xticklabels(
        [str(p) for p in pos_list],
        rotation=90,
        fontsize=min(8.5, max(5.0, int(100 / max(1, n_vars)))),
    )
    ax_heat.set_xlabel("Variant position", fontsize=11)

    ax_heat.set_yticks([])
    for g, s, e in group_slices:
        ax_heat.plot([0, n_vars], [e, e], color="white", linewidth=1.0, zorder=5)

    ax_heat.set_title(
        f"All-Sample Haplotype Heatmap: {label} (N={n_samples:,})", fontsize=13, pad=10
    )

    gt_handles = [
        Patch(facecolor=ref_color, edgecolor="black", label="REF"),
        Patch(facecolor=alt_color, edgecolor="black", label="ALT"),
        Patch(facecolor=het_color, edgecolor="black", label="HET"),
        Patch(facecolor=missing_color, edgecolor="black", label="Missing"),
        Patch(facecolor=absence_color, edgecolor="black", label="Absence"),
    ]
    group_counts = {g: e - s for g, s, e in group_slices}
    grp_handles = [
        Patch(
            facecolor=group_colors[g],
            edgecolor="black",
            label=f"{g} (n={group_counts[g]:,})",
        )
        for g in unique_groups
    ]

    leg1 = ax_heat.legend(
        handles=gt_handles,
        title="Genotype",
        frameon=False,
        bbox_to_anchor=(1.01, 1.0),
        loc="upper left",
        fontsize=8.5,
        title_fontsize=9.5,
    )
    ax_heat.add_artist(leg1)
    ax_heat.legend(
        handles=grp_handles,
        title="Group",
        frameon=False,
        bbox_to_anchor=(1.01, 0.62),
        loc="upper left",
        fontsize=8.5,
        title_fontsize=9.5,
    )

    fig.subplots_adjust(left=0.08, right=0.76, bottom=0.14, top=0.92)

    plot_stem = os.path.join(outdir, f"{label}.AllSampleHaplotypeHeatmap")
    _save_formats(fig, plot_stem, plot_formats)
    plt.close(fig)

    out_file = plot_stem + "." + list(plot_formats)[0].lower()
    print(f"[EasyHap Heatmap] All-sample haplotype heatmap generated: {out_file}")
    return out_file


def run_all_sample_heatmap(
    indir: str,
    region: Optional[str] = None,
    region_file: Optional[str] = None,
    group_file: Optional[str] = None,
    outdir: Optional[str] = None,
    plot_formats: Sequence[str] = ("pdf",),
    cluster_method: str = "average",
    ref_color: str = DEFAULT_REF_COLOR,
    alt_color: str = DEFAULT_ALT_COLOR,
    het_color: str = "#F4B183",
    missing_color: str = DEFAULT_MISSING_COLOR,
    absence_color: str = DEFAULT_ABSENCE_COLOR,
    hap_palette: Optional[Sequence[str]] = None,
) -> List[str]:
    """Generate all-sample haplotype color heatmap for one or more regions."""
    regions_to_run: List[Optional[str]] = []
    if region_file:
        if not os.path.isfile(region_file):
            raise FileNotFoundError(f"Region file not found: {region_file}")
        with open(region_file, "r", encoding="utf-8") as rf:
            for line in rf:
                r = line.strip()
                if not r or r.startswith("#"):
                    continue
                parts = r.split()
                if len(parts) == 3:
                    from .vcf_reader import Region

                    regions_to_run.append(
                        Region(parts[0], int(parts[1]), int(parts[2])).vcf_label
                    )
                elif len(parts) == 1:
                    regions_to_run.append(parts[0])
                else:
                    raise ValueError(
                        "Heatmap region rows require one region label or chrom/start/end"
                    )
    if region:
        if region not in regions_to_run:
            regions_to_run.append(region.strip())

    if not regions_to_run:
        # Auto-detect if single region in indir
        regions_to_run = [None]

    generated_files = []
    for reg in regions_to_run:
        out = _generate_single_heatmap(
            indir=indir,
            region=reg,
            group_file=group_file,
            outdir=outdir,
            plot_formats=plot_formats,
            cluster_method=cluster_method,
            ref_color=ref_color,
            alt_color=alt_color,
            het_color=het_color,
            missing_color=missing_color,
            absence_color=absence_color,
            hap_palette=hap_palette,
        )
        generated_files.append(out)

    return generated_files
