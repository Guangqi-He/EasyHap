from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple
import gzip
import itertools
import math
import re

import matplotlib
matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import ConnectionPatch, Rectangle
import numpy as np
import pandas as pd

from .core import HapResult
from .io_utils import sanitize_filename
from .stats import bh_adjust

# A unified soft categorical palette used by heatmaps, pies, stacked bars, and boxplots.
# The colors are intentionally muted/pastel so dense haplotype figures are not visually harsh.
SOFT_PALETTE = [
    "#A6CEE3",  # soft blue
    "#B2DF8A",  # soft green
    "#FDBF6F",  # soft orange
    "#CAB2D6",  # soft purple
    "#FB9A99",  # soft red
    "#CCEBC5",  # pale green
    "#B3CDE3",  # pale blue
    "#FBB4AE",  # pale salmon
    "#DECBE4",  # pale lavender
    "#FED9A6",  # pale peach
    "#E5D8BD",  # pale brown
    "#D9D9D9",  # neutral grey
]
NEUTRAL_EDGE = "black"
GENE_LINE_COLOR = "black"
VARIANT_LINE_COLOR = "black"


@dataclass
class GeneFeature:
    seqid: str
    start: int
    end: int
    strand: str
    ftype: str
    attrs: Dict[str, str]


def _save_formats(fig, prefix: str, formats: Sequence[str]) -> None:
    for fmt in formats:
        fmt = fmt.lower().lstrip(".")
        save_kwargs = {"bbox_inches": "tight"}
        if fmt == "png":
            save_kwargs["dpi"] = 300
        fig.savefig(prefix + "." + fmt, **save_kwargs)
    plt.close(fig)


def _soft_cmap(n_states: int) -> ListedColormap:
    n = max(1, n_states)
    colors = [SOFT_PALETTE[i % len(SOFT_PALETTE)] for i in range(n)]
    return ListedColormap(colors)


def _color_map(labels: Sequence[str]) -> Dict[str, str]:
    return {label: SOFT_PALETTE[i % len(SOFT_PALETTE)] for i, label in enumerate(labels)}


def _natural_key(text: str) -> List[object]:
    return [int(x) if x.isdigit() else x for x in re.split(r"(\d+)", str(text))]


def _plot_level_column(plot_hap_level: str) -> str:
    if plot_hap_level == "hap":
        return "Hap"
    if plot_hap_level == "cluster":
        return "ClusterID"
    raise ValueError("plot_hap_level must be 'hap' or 'cluster'")


def _plot_level_title(plot_hap_level: str) -> str:
    return "Haplotype/diplotype" if plot_hap_level == "hap" else "Clustered haplotype/diplotype"


def _plot_stem(result: HapResult, plot_hap_level: str, stem: str) -> str:
    if plot_hap_level == "hap":
        return result.output_prefix + "." + stem
    return result.output_prefix + ".Cluster" + stem


def _parse_attrs(text: str) -> Dict[str, str]:
    attrs: Dict[str, str] = {}
    for item in text.split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" in item:
            k, v = item.split("=", 1)
        elif " " in item:
            k, v = item.split(" ", 1)
            v = v.strip('"')
        else:
            continue
        attrs[k] = v
    return attrs


def read_gff_features(path: str, chrom: str, start: int, end: int) -> List[GeneFeature]:
    features: List[GeneFeature] = []
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:  # type: ignore[arg-type]
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            seqid, _source, ftype_raw, s, e, _score, strand, _phase, attrs = parts
            s_i, e_i = int(s), int(e)
            ftype_l = ftype_raw.lower()
            type_map = {
                "gene": "gene",
                "mrna": "mRNA",
                "transcript": "transcript",
                "exon": "exon",
                "cds": "CDS",
                "utr": "UTR",
                "five_prime_utr": "five_prime_UTR",
                "5utr": "five_prime_UTR",
                "5_prime_utr": "five_prime_UTR",
                "three_prime_utr": "three_prime_UTR",
                "3utr": "three_prime_UTR",
                "3_prime_utr": "three_prime_UTR",
            }
            ftype = type_map.get(ftype_l, ftype_raw)
            if (
                seqid == chrom
                and e_i >= start
                and s_i <= end
                and ftype in {"gene", "mRNA", "transcript", "exon", "CDS", "UTR", "five_prime_UTR", "three_prime_UTR"}
            ):
                features.append(GeneFeature(seqid, s_i, e_i, strand, ftype, _parse_attrs(attrs)))
    return features


def _feature_name(feat: GeneFeature) -> str:
    for key in ("Name", "gene", "gene_name", "ID", "locus_tag"):
        if key in feat.attrs:
            return feat.attrs[key]
    return feat.ftype


def _split_composite_token(token: str) -> List[str]:
    """Split slash-joined genotype tokens and remove empty fragments.

    Cluster heatmaps may combine haplotypes that already contain heterozygous
    states such as A/G. Splitting first prevents labels such as A/G/G.
    """
    if token in {"", "NA"}:
        return []
    return [x for x in str(token).split("/") if x]


def _cluster_consensus_token(tokens: Sequence[str]) -> str:
    flattened: List[str] = []
    for tok in tokens:
        flattened.extend(_split_composite_token(tok))
    # dict.fromkeys keeps first occurrence; sorting makes output stable.
    uniq = sorted(dict.fromkeys(flattened), key=_natural_key)
    if not uniq:
        return "N"
    if len(uniq) == 1:
        return uniq[0]
    return "/".join(uniq)


def _token_matrix(
    result: HapResult,
    plot_hap_level: str = "hap",
    min_count: int = 1,
) -> Tuple[List[str], Dict[str, int], List[str], np.ndarray, List[List[str]]]:
    """Return row ids, accession counts, allele tokens, numeric matrix, text matrix.

    For plot_hap_level='cluster', each row is a haplotype cluster. A site is shown as a
    consensus state when all member haplotypes agree, or as a slash-joined composite
    state when the cluster contains multiple allele states at that site.
    """
    if plot_hap_level == "hap":
        row_ids = sorted(result.hap_sequences, key=_natural_key)
        row_count = {hid: len(result.hap_accessions.get(hid, [])) for hid in row_ids}
        text = [list(result.hap_sequences[hid]) for hid in row_ids]
    elif plot_hap_level == "cluster":
        cluster_to_haps: Dict[str, List[str]] = {}
        for hid, cid in result.hap_clusters.items():
            cluster_to_haps.setdefault(cid, []).append(hid)
        row_ids = sorted(cluster_to_haps, key=_natural_key)
        row_count: Dict[str, int] = {}
        text = []
        for cid in row_ids:
            hids = sorted(cluster_to_haps[cid], key=_natural_key)
            seqs = [result.hap_sequences[hid] for hid in hids]
            accessions = sorted({acc for hid in hids for acc in result.hap_accessions.get(hid, [])})
            row_count[cid] = len(accessions)
            if not seqs:
                text.append([])
            else:
                n_sites = len(seqs[0])
                row = [_cluster_consensus_token([seq[i] for seq in seqs]) for i in range(n_sites)]
                text.append(row)
    else:
        raise ValueError("plot_hap_level must be 'hap' or 'cluster'")

    min_count = max(1, int(min_count or 1))
    if min_count > 1 and row_ids:
        keep_idx = [i for i, row in enumerate(row_ids) if row_count.get(row, 0) >= min_count]
        row_ids = [row_ids[i] for i in keep_idx]
        text = [text[i] for i in keep_idx]
        row_count = {row: row_count[row] for row in row_ids}

    tokens = sorted({tok for row in text for tok in row}, key=_natural_key)
    tok_to_num = {t: i for i, t in enumerate(tokens)}
    matrix = np.array([[tok_to_num[t] for t in row] for row in text], dtype=float) if text else np.empty((0, 0))
    return row_ids, row_count, tokens, matrix, text


def _heatmap_figsize(n_cols: int, n_rows: int, extra_left: float = 3.8, extra_bottom: float = 2.2) -> Tuple[float, float]:
    # A near-square cell geometry. The exact rendered cell size depends on tick labels,
    # but equal axis aspect keeps the plotting area itself square-celled.
    cell = 0.56
    width = max(5.5, min(extra_left + n_cols * cell, 36.0))
    height = max(3.8, min(extra_bottom + n_rows * cell, 32.0))
    return width, height


def _draw_token_heatmap(
    ax: plt.Axes,
    matrix: np.ndarray,
    text_matrix: List[List[str]],
    hap_labels: Sequence[str],
    variant_labels: Sequence[str],
    title: Optional[str] = None,
    show_xlabel: bool = True,
    show_ylabel: bool = True,
    square_cells: bool = True,
) -> None:
    if matrix.size == 0:
        return
    n_rows, n_cols = matrix.shape
    n_states = int(np.nanmax(matrix)) + 1 if matrix.size else 1
    ax.pcolormesh(
        np.arange(n_cols + 1),
        np.arange(n_rows + 1),
        matrix,
        cmap=_soft_cmap(n_states),
        vmin=-0.5,
        vmax=max(n_states - 0.5, 0.5),
        edgecolors=NEUTRAL_EDGE,
        linewidth=0.75,
        shading="flat",
    )
    ax.set_xlim(0, n_cols)
    ax.set_ylim(n_rows, 0)
    if square_cells:
        ax.set_aspect("equal", adjustable="box")
    else:
        ax.set_aspect("auto")
    ax.set_xticks(np.arange(n_cols) + 0.5)
    ax.set_xticklabels(variant_labels, rotation=90, fontsize=9)
    ax.set_yticks(np.arange(n_rows) + 0.5)
    ax.set_yticklabels(hap_labels, fontsize=9)
    if show_xlabel:
        ax.set_xlabel("Variant position")
    if show_ylabel:
        ax.set_ylabel("Haplotype/cluster (sample count)")
    if title:
        ax.set_title(title)

    # Keep letters readable. For extremely large matrices, reduce slightly but avoid tiny text.
    max_token_len = max((len(tok) for row in text_matrix for tok in row), default=1)
    if n_cols <= 35 and n_rows <= 40 and max_token_len <= 4:
        font_size = 11
    elif max_token_len <= 7:
        font_size = 9.5
    else:
        font_size = 8.5
    for i in range(n_rows):
        for j in range(n_cols):
            ax.text(j + 0.5, i + 0.5, text_matrix[i][j], ha="center", va="center", fontsize=font_size, color="black")

    for spine in ax.spines.values():
        spine.set_linewidth(0.8)
        spine.set_edgecolor(NEUTRAL_EDGE)


def plot_haplotype_heatmap(result: HapResult, plot_formats: Sequence[str], plot_hap_level: str = "hap", plot_min_count: int = 1) -> None:
    row_ids, row_count, _tokens, matrix, text_matrix = _token_matrix(result, plot_hap_level, min_count=plot_min_count)
    variants = result.variants
    if not row_ids or not variants:
        return
    width, height = _heatmap_figsize(len(variants), len(row_ids))
    fig, ax = plt.subplots(figsize=(width, height))
    row_labels = [f"{row} ({row_count.get(row, 0)})" for row in row_ids]
    variant_labels = [str(v.pos) for v in variants]
    _draw_token_heatmap(
        ax,
        matrix,
        text_matrix,
        row_labels,
        variant_labels,
        title=f"{_plot_level_title(plot_hap_level)} heatmap: {result.region.vcf_label}",
    )
    _save_formats(fig, _plot_stem(result, plot_hap_level, "HaplotypeHeatmap"), plot_formats)


def plot_gene_structure_with_haps(
    result: HapResult,
    gff_file: Optional[str],
    plot_formats: Sequence[str],
    plot_hap_level: str = "hap",
    plot_min_count: int = 1,
) -> None:
    if not gff_file:
        return
    features = read_gff_features(gff_file, result.region.chrom, result.region.start, result.region.end)
    if not features:
        return
    row_ids, row_count, _tokens, matrix, text_matrix = _token_matrix(result, plot_hap_level, min_count=plot_min_count)
    variants = result.variants
    if not row_ids or not variants:
        return

    # Keep the heatmap panel and gene-structure panel at the same overall width.
    # In this combined figure, the heatmap fills the horizontal axis so the two
    # panels align visually; standalone heatmaps still use near-square cells.
    width = max(7.0, min(4.2 + len(variants) * 0.72, 36.0))
    heat_height = max(2.2, min(1.4 + len(row_ids) * 0.45, 24.0))
    gene_panel_height = 1.125  # 50% taller gene model panel for better readability
    fig = plt.figure(figsize=(width, min(30.0, heat_height + gene_panel_height + 0.55)))
    gs = fig.add_gridspec(2, 1, height_ratios=[gene_panel_height, heat_height], hspace=0.08)
    ax_gene = fig.add_subplot(gs[0, 0])
    ax_heat = fig.add_subplot(gs[1, 0])

    gene_feats = [f for f in features if f.ftype == "gene"]
    if not gene_feats:
        # GTF files sometimes lack explicit gene records in a small interval.
        span_start = min(f.start for f in features)
        span_end = max(f.end for f in features)
        gene_feats = [GeneFeature(result.region.chrom, span_start, span_end, ".", "gene", {"Name": "gene_model"})]

    for idx, gene in enumerate(gene_feats):
        y = idx
        # Region-level baseline: if the requested region is wider than the gene,
        # the non-gene flanking part is shown as a thin line instead of blank space.
        ax_gene.hlines(y, result.region.start, result.region.end, linewidth=0.65, color="black", zorder=1)
        ax_gene.hlines(
            y,
            max(gene.start, result.region.start),
            min(gene.end, result.region.end),
            linewidth=1.2,
            color=GENE_LINE_COLOR,
            zorder=2,
        )
        ax_gene.text(max(gene.start, result.region.start), y + 0.34, _feature_name(gene), fontsize=8, ha="left", color="black")

        # Draw UTRs in grey. CDS/exon blocks remain soft and low-saturation.
        child_feats = [
            f for f in features
            if f.ftype in {"exon", "CDS", "UTR", "five_prime_UTR", "three_prime_UTR"}
            and f.end >= gene.start and f.start <= gene.end
        ]
        # Draw exons first, then UTR/CDS on top, so UTR information remains visible.
        feature_priority = {"exon": 0, "UTR": 1, "five_prime_UTR": 1, "three_prime_UTR": 1, "CDS": 2}
        for feat in sorted(child_feats, key=lambda f: (feature_priority.get(f.ftype, 0), f.start, f.end)):
            left = max(feat.start, result.region.start)
            right = min(feat.end, result.region.end)
            if right < left:
                continue
            if feat.ftype in {"UTR", "five_prime_UTR", "three_prime_UTR"}:
                height = 0.24
                face = "#D9D9D9"
            elif feat.ftype == "CDS":
                height = 0.45
                face = SOFT_PALETTE[0]
            else:  # exon when CDS/UTR subfeatures are absent
                height = 0.33
                face = "#F2F2F2"
            ax_gene.add_patch(
                Rectangle(
                    (left, y - height / 2),
                    max(1, right - left + 1),
                    height,
                    facecolor=face,
                    edgecolor=NEUTRAL_EDGE,
                    linewidth=0.65,
                    zorder=3,
                )
            )

    for j, v in enumerate(variants):
        ax_gene.axvline(v.pos, linestyle="--", linewidth=0.5, color=VARIANT_LINE_COLOR, alpha=1.0, zorder=0)
        con = ConnectionPatch(
            xyA=(v.pos, -0.68),
            coordsA=ax_gene.transData,
            xyB=(j + 0.5, 0),
            coordsB=ax_heat.transData,
            color=VARIANT_LINE_COLOR,
            linewidth=0.45,
            alpha=1.0,
        )
        fig.add_artist(con)

    ax_gene.set_xlim(result.region.start, result.region.end)
    ax_gene.set_ylim(-0.88, max(0.88, len(gene_feats) - 0.12))
    ax_gene.set_yticks([])
    ax_gene.set_ylabel("Gene", fontsize=9)
    ax_gene.set_title(f"Gene structure and {_plot_level_title(plot_hap_level).lower()}: {result.region.vcf_label}", fontsize=11)
    # The gene panel is an annotation track; remove genomic x-axis ticks/labels
    # and the outer frame to keep the combined figure cleaner.
    ax_gene.tick_params(axis="x", which="both", bottom=False, top=False, labelbottom=False)
    ax_gene.tick_params(axis="y", which="both", left=False, right=False)
    for spine in ax_gene.spines.values():
        spine.set_visible(False)

    row_labels = [f"{row} ({row_count.get(row, 0)})" for row in row_ids]
    variant_labels = [str(v.pos) for v in variants]
    _draw_token_heatmap(
        ax_heat,
        matrix,
        text_matrix,
        row_labels,
        variant_labels,
        title=None,
        show_xlabel=True,
        show_ylabel=True,
        square_cells=False,
    )
    stem = "GeneHaplotype" if plot_hap_level == "hap" else "GeneClusterHaplotype"
    _save_formats(fig, result.output_prefix + "." + stem, plot_formats)


def _filter_classes_by_min_count(df: pd.DataFrame, class_col: str, min_count: int) -> Tuple[pd.DataFrame, List[str]]:
    """Filter sample-level plotting classes by total population count."""
    if df.empty or class_col not in df.columns:
        return df, []
    min_count = max(1, int(min_count or 1))
    counts = df[class_col].value_counts()
    keep = counts[counts >= min_count].index.tolist()
    if min_count <= 1:
        keep = counts.index.tolist()
    return df[df[class_col].isin(keep)].copy(), keep


def plot_group_distribution(
    result: HapResult,
    hap_group_df: pd.DataFrame,
    plot_formats: Sequence[str],
    plot_hap_level: str = "hap",
    plot_min_count: int = 1,
) -> None:
    class_col = _plot_level_column(plot_hap_level)
    if hap_group_df.empty or "Type" not in hap_group_df.columns or class_col not in hap_group_df.columns:
        return

    plot_df, _keep = _filter_classes_by_min_count(hap_group_df, class_col, plot_min_count)
    if plot_df.empty:
        return

    counts = plot_df.groupby(["Type", class_col]).size().unstack(fill_value=0)
    # Order classes by their total abundance so the colors and legend are stable
    # across the stacked bar chart and the combined group pie figure.
    class_order = counts.sum(axis=0).sort_values(ascending=False).index.tolist()
    if not class_order:
        return
    counts = counts.loc[:, class_order]
    props = counts.div(counts.sum(axis=1), axis=0)
    cmap = _color_map(class_order)
    label_title = _plot_level_title(plot_hap_level)
    filter_note = f" (n≥{plot_min_count})" if int(plot_min_count or 1) > 1 else ""

    fig, ax = plt.subplots(figsize=(max(7, 0.65 * len(class_order) + 3), 5))
    bottom = np.zeros(len(props.index))
    x = np.arange(len(props.index))
    for cls in class_order:
        vals = props[cls].values
        ax.bar(x, vals, bottom=bottom, label=cls, color=cmap[cls], edgecolor=NEUTRAL_EDGE, linewidth=0.5)
        bottom += vals
    ax.set_xticks(x)
    group_sample_counts = counts.sum(axis=1).astype(int).to_dict()
    group_labels = [f"{grp}\n(n={group_sample_counts.get(grp, 0)})" for grp in props.index]
    ax.set_xticklabels(group_labels, rotation=0, ha="center")
    ax.set_ylabel("Proportion among displayed classes")
    ax.set_xlabel("Group")
    ax.set_ylim(0, 1)
    ax.set_title(f"{label_title} proportion by group{filter_note}")
    ax.legend(fontsize=7, bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False, title=label_title)
    _save_formats(fig, _plot_stem(result, plot_hap_level, "GroupStackedBar"), plot_formats)

    # Combined pie chart: all groups are shown in one figure instead of generating
    # one independent pie-chart file per group. This makes group-level composition
    # easier to compare and keeps the output directory cleaner.
    group_order = counts.index.tolist()
    n_groups = len(group_order)
    if n_groups == 0:
        return
    n_cols = min(3, n_groups)
    n_rows = int(np.ceil(n_groups / n_cols))
    fig_w = max(5.2 * n_cols, 6.0)
    fig_h = max(4.6 * n_rows, 4.8)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(fig_w, fig_h), squeeze=False)

    def _autopct(pct: float) -> str:
        # Hide tiny labels to avoid clutter in dense haplotype compositions.
        return f"{pct:.1f}%" if pct >= 4 else ""

    for ax, grp in zip(axes.ravel(), group_order):
        vc = counts.loc[grp]
        vc = vc[vc > 0]
        colors = [cmap[c] for c in vc.index]
        _wedges, _texts, autotexts = ax.pie(
            vc.values,
            labels=None,
            autopct=_autopct,
            startangle=90,
            counterclock=False,
            colors=colors,
            wedgeprops={"edgecolor": NEUTRAL_EDGE, "linewidth": 0.6},
            textprops={"fontsize": 9, "color": "black"},
        )
        for t in autotexts:
            t.set_fontsize(8.5)
        ax.set_title(f"{grp} (displayed n={int(vc.sum())})", fontsize=11)
        ax.axis("equal")

    for ax in axes.ravel()[n_groups:]:
        ax.axis("off")

    from matplotlib.patches import Patch

    handles = [Patch(facecolor=cmap[c], edgecolor=NEUTRAL_EDGE, linewidth=0.5, label=c) for c in class_order]
    fig.legend(
        handles=handles,
        title=label_title,
        fontsize=8,
        title_fontsize=9,
        bbox_to_anchor=(1.01, 0.5),
        loc="center left",
        frameon=False,
    )
    fig.suptitle(f"{label_title} proportion by group{filter_note}", fontsize=12, y=0.98)
    fig.subplots_adjust(left=0.04, right=0.82, top=0.88, bottom=0.06, wspace=0.20, hspace=0.28)
    _save_formats(fig, _plot_stem(result, plot_hap_level, "GroupPie"), plot_formats)


def _p_to_symbol(p: object) -> str:
    try:
        pf = float(p)  # type: ignore[arg-type]
    except Exception:
        return "NA"
    if math.isnan(pf):
        return "NA"
    if pf <= 0.001:
        return "***"
    if pf <= 0.01:
        return "**"
    if pf <= 0.05:
        return "*"
    return "ns"


def _overall_trait_test(
    tmp: pd.DataFrame,
    class_col: str,
    trait: str,
    order: Sequence[str],
    plot_hap_level: str,
) -> Dict[str, object]:
    row: Dict[str, object] = {
        "Trait": trait,
        "PlotLevel": plot_hap_level,
        "ComparisonType": "overall",
        "Class1": "ALL",
        "Class2": "",
        "N1": int(tmp[trait].notna().sum()),
        "N2": "",
        "Mean1": "",
        "Mean2": "",
        "Test": "Kruskal-Wallis",
        "pvalue": np.nan,
        "padj_BH": np.nan,
        "Significance": "NA",
        "SkippedReason": "",
    }
    groups = [tmp.loc[tmp[class_col] == cls, trait].dropna().astype(float) for cls in order]
    groups = [g for g in groups if len(g) > 0]
    if len(groups) < 2:
        row["SkippedReason"] = "fewer than two classes with numeric trait values"
        return row
    try:
        from scipy.stats import kruskal  # type: ignore
    except Exception:
        row["SkippedReason"] = "scipy is not installed; install scipy>=1.10 for significance testing"
        return row
    try:
        pvalue = float(kruskal(*groups).pvalue)
        row["pvalue"] = pvalue
        row["padj_BH"] = pvalue
        row["Significance"] = _p_to_symbol(pvalue)
    except Exception as exc:
        row["SkippedReason"] = f"Kruskal-Wallis failed: {exc}"
    return row


def _pairwise_trait_tests(
    tmp: pd.DataFrame,
    class_col: str,
    trait: str,
    order: Sequence[str],
    plot_hap_level: str,
) -> pd.DataFrame:
    rows: List[Dict[str, object]] = []
    valid_row_idx: List[int] = []
    pvals: List[float] = []
    try:
        from scipy.stats import mannwhitneyu  # type: ignore
        scipy_error = ""
    except Exception:
        mannwhitneyu = None  # type: ignore[assignment]
        scipy_error = "scipy is not installed; install scipy>=1.10 for significance testing"

    for c1, c2 in itertools.combinations(order, 2):
        x = tmp.loc[tmp[class_col] == c1, trait].dropna().astype(float)
        y = tmp.loc[tmp[class_col] == c2, trait].dropna().astype(float)
        row: Dict[str, object] = {
            "Trait": trait,
            "PlotLevel": plot_hap_level,
            "ComparisonType": "pairwise",
            "Class1": c1,
            "Class2": c2,
            "N1": int(len(x)),
            "N2": int(len(y)),
            "Mean1": float(np.mean(x)) if len(x) else np.nan,
            "Mean2": float(np.mean(y)) if len(y) else np.nan,
            "Test": "two-sided Mann-Whitney U",
            "pvalue": np.nan,
            "padj_BH": np.nan,
            "Significance": "NA",
            "SkippedReason": "",
        }
        if len(x) == 0 or len(y) == 0:
            row["SkippedReason"] = "one or both classes have no numeric trait values"
            rows.append(row)
            continue
        if scipy_error:
            row["SkippedReason"] = scipy_error
            rows.append(row)
            continue
        try:
            try:
                pvalue = float(mannwhitneyu(x, y, alternative="two-sided", method="auto").pvalue)  # type: ignore[misc]
            except TypeError:
                pvalue = float(mannwhitneyu(x, y, alternative="two-sided").pvalue)  # type: ignore[misc]
            row["pvalue"] = pvalue
            pvals.append(pvalue)
            valid_row_idx.append(len(rows))
        except Exception as exc:
            row["SkippedReason"] = f"Mann-Whitney U failed: {exc}"
        rows.append(row)

    if pvals:
        padj = bh_adjust(pvals)
        for idx, q in zip(valid_row_idx, padj):
            rows[idx]["padj_BH"] = q
            rows[idx]["Significance"] = _p_to_symbol(q)
    return pd.DataFrame(rows)


def _annotate_pairwise_significance(
    ax: plt.Axes,
    pairwise_df: pd.DataFrame,
    order: Sequence[str],
    max_annotations: int = 10,
) -> None:
    if pairwise_df.empty or "padj_BH" not in pairwise_df.columns:
        if len(order) >= 2:
            ax.text(0.98, 0.98, "no pairwise test", transform=ax.transAxes, ha="right", va="top", fontsize=8)
        return
    tested = pairwise_df[pd.to_numeric(pairwise_df["padj_BH"], errors="coerce").notna()].copy()
    if tested.empty:
        if len(order) >= 2:
            ax.text(0.98, 0.98, "tests skipped", transform=ax.transAxes, ha="right", va="top", fontsize=8)
        return
    sig = tested.loc[pd.to_numeric(tested["padj_BH"], errors="coerce") <= 0.05].copy()
    if sig.empty:
        if len(order) >= 2:
            ax.text(0.98, 0.98, "pairwise ns", transform=ax.transAxes, ha="right", va="top", fontsize=8)
        return
    sig["padj_BH_numeric"] = pd.to_numeric(sig["padj_BH"], errors="coerce")
    sig["pvalue_numeric"] = pd.to_numeric(sig["pvalue"], errors="coerce")
    sig = sig.sort_values(["padj_BH_numeric", "pvalue_numeric"]).head(max_annotations)
    x_pos = {cls: i + 1 for i, cls in enumerate(order)}
    y_min, y_max = ax.get_ylim()
    y_range = y_max - y_min
    if y_range <= 0:
        y_range = 1.0
    base = y_max + 0.05 * y_range
    h = 0.025 * y_range
    step = 0.075 * y_range
    current = base
    for _, row in sig.iterrows():
        c1, c2 = row["Class1"], row["Class2"]
        if c1 not in x_pos or c2 not in x_pos:
            continue
        x1, x2 = x_pos[c1], x_pos[c2]
        if x1 > x2:
            x1, x2 = x2, x1
        ax.plot([x1, x1, x2, x2], [current, current + h, current + h, current], color=NEUTRAL_EDGE, linewidth=0.8)
        ax.text((x1 + x2) / 2, current + h, str(row["Significance"]), ha="center", va="bottom", fontsize=8.5)
        current += step
    ax.set_ylim(y_min, current + 0.08 * y_range)


def plot_trait_boxplots(
    result: HapResult,
    hap_group_df: pd.DataFrame,
    traits_to_plot: Optional[Sequence[str]],
    plot_formats: Sequence[str],
    plot_hap_level: str = "hap",
    plot_min_count: int = 1,
) -> None:
    class_col = _plot_level_column(plot_hap_level)
    trait_cols = [c for c in hap_group_df.columns if c not in {"Hap", "ClusterID", "Accession", "Type"}]
    if traits_to_plot:
        trait_cols = [c for c in trait_cols if c in set(traits_to_plot)]
    if not trait_cols or class_col not in hap_group_df.columns:
        return

    plot_df, _keep = _filter_classes_by_min_count(hap_group_df, class_col, plot_min_count)
    if plot_df.empty:
        return

    # Keep the overall order consistent across traits, but skip classes with no numeric
    # observations within a specific trait panel.
    base_order = plot_df[class_col].value_counts().index.tolist()
    cmap = _color_map(base_order)
    rng = np.random.default_rng(42)

    valid_traits: List[str] = []
    for trait in trait_cols:
        tmp = plot_df[[class_col, trait]].copy()
        tmp[trait] = pd.to_numeric(tmp[trait], errors="coerce")
        if tmp[trait].notna().any():
            valid_traits.append(trait)
    if not valid_traits:
        return

    n = len(valid_traits)
    fig, axes = plt.subplots(1, n, figsize=(max(5.2 * n, 6), 5.8), squeeze=False)
    axes_list = axes[0]
    significance_rows: List[pd.DataFrame] = []
    level_title = _plot_level_title(plot_hap_level)
    filter_note = f"; displayed classes n≥{plot_min_count}" if int(plot_min_count or 1) > 1 else ""

    for ax, trait in zip(axes_list, valid_traits):
        tmp = plot_df[[class_col, trait]].copy()
        tmp[trait] = pd.to_numeric(tmp[trait], errors="coerce")
        order = [cls for cls in base_order if tmp.loc[tmp[class_col] == cls, trait].notna().sum() > 0]
        if not order:
            ax.axis("off")
            continue
        data = [tmp.loc[tmp[class_col] == cls, trait].dropna() for cls in order]
        sample_counts = {cls: int(tmp.loc[tmp[class_col] == cls, trait].notna().sum()) for cls in order}
        display_labels = [f"{cls}\n(n={sample_counts.get(cls, 0)})" for cls in order]
        bp = ax.boxplot(
            data,
            labels=display_labels,
            showfliers=False,
            patch_artist=True,
            medianprops={"color": "black", "linewidth": 1.0},
            boxprops={"edgecolor": NEUTRAL_EDGE, "linewidth": 0.9},
            whiskerprops={"color": NEUTRAL_EDGE, "linewidth": 0.8},
            capprops={"color": NEUTRAL_EDGE, "linewidth": 0.8},
        )
        for patch, cls in zip(bp["boxes"], order):
            patch.set_facecolor(cmap[cls])
            patch.set_alpha(0.38)
        for i, (cls, vals) in enumerate(zip(order, data), 1):
            if vals.empty:
                continue
            jitter = rng.normal(i, 0.035, size=len(vals))
            ax.scatter(
                jitter,
                vals.values,
                s=12,
                alpha=0.68,
                color=cmap[cls],
                edgecolors="none",
                linewidth=0,
                zorder=3,
            )

        overall = _overall_trait_test(tmp, class_col, trait, order, plot_hap_level)
        p = overall.get("pvalue", np.nan)
        title = trait if pd.isna(p) else f"{trait}\nKruskal-Wallis p={float(p):.3g}"
        ax.set_title(title)
        ax.set_xlabel(level_title)
        ax.set_ylabel(trait)
        ax.tick_params(axis="x", rotation=60, labelsize=8)

        pairwise = _pairwise_trait_tests(tmp, class_col, trait, order, plot_hap_level)
        sig_parts = [pd.DataFrame([overall])]
        if not pairwise.empty:
            sig_parts.append(pairwise)
            _annotate_pairwise_significance(ax, pairwise, order)
        else:
            _annotate_pairwise_significance(ax, pd.DataFrame(), order)
        significance_rows.append(pd.concat(sig_parts, ignore_index=True))

    fig.subplots_adjust(top=0.76, bottom=0.30, wspace=0.30)
    fig.suptitle(f"Trait distribution by sample-level {level_title.lower()}{filter_note}", y=0.98, fontsize=12)
    _save_formats(fig, _plot_stem(result, plot_hap_level, "TraitBoxplot"), plot_formats)

    if significance_rows:
        sig_df = pd.concat(significance_rows, ignore_index=True)
        sig_path = _plot_stem(result, plot_hap_level, "TraitSignificance") + ".tsv"
        sig_df.to_csv(sig_path, sep="\t", index=False)


def make_all_plots(
    result: HapResult,
    hap_group_df: pd.DataFrame,
    group_map: Dict[str, str],
    gff_file: Optional[str] = None,
    plot_formats: Sequence[str] = ("pdf",),
    traits_to_plot: Optional[Sequence[str]] = None,
    plot_hap_level: str = "hap",
    plot_min_count: int = 1,
) -> None:
    plot_hap_level = plot_hap_level.lower()
    if plot_hap_level not in {"hap", "cluster"}:
        raise ValueError("plot_hap_level must be 'hap' or 'cluster'")
    plot_min_count = max(1, int(plot_min_count or 1))
    plot_haplotype_heatmap(result, plot_formats, plot_hap_level=plot_hap_level, plot_min_count=plot_min_count)
    plot_gene_structure_with_haps(result, gff_file, plot_formats, plot_hap_level=plot_hap_level, plot_min_count=plot_min_count)
    plot_group_distribution(result, hap_group_df, plot_formats, plot_hap_level=plot_hap_level, plot_min_count=plot_min_count)
    plot_trait_boxplots(result, hap_group_df, traits_to_plot, plot_formats, plot_hap_level=plot_hap_level, plot_min_count=plot_min_count)
