from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple
import gzip
import math
import re

import matplotlib
matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import ListedColormap, Normalize, to_rgba
from matplotlib.lines import Line2D
from matplotlib.patches import ConnectionPatch, Patch, Polygon, Rectangle
import numpy as np
import pandas as pd

from .core import HapResult
from .stats import bh_adjust

# Publication-oriented defaults.
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans"],
    "font.size": 11,
    "axes.titlesize": 13,
    "axes.labelsize": 11,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 9,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

DEFAULT_REF_COLOR = "#70AD47"
DEFAULT_ALT_COLOR = "#4472C4"
DEFAULT_MISSING_COLOR = "#D9D9D9"
GENE_CDS_COLOR = "black"
GENE_UTR_COLOR = "#BFBFBF"
DEFAULT_HAP_PALETTE = [
    "#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F",
    "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC",
]
DEFAULT_LD_CMAP = "viridis"


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
        kw = {"bbox_inches": "tight"}
        if fmt == "png":
            kw["dpi"] = 600
        fig.savefig(prefix + "." + fmt, **kw)
    plt.close(fig)


def _parse_attrs(text: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
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
        out[k] = v.strip()
    return out


def read_gff_features(path: str, chrom: str, start: int, end: int) -> List[GeneFeature]:
    features: List[GeneFeature] = []
    opener = gzip.open if path.endswith(".gz") else open
    type_map = {
        "gene": "gene", "mrna": "mRNA", "transcript": "transcript", "exon": "exon", "cds": "CDS",
        "utr": "UTR", "five_prime_utr": "five_prime_UTR", "5utr": "five_prime_UTR", "5_prime_utr": "five_prime_UTR",
        "three_prime_utr": "three_prime_UTR", "3utr": "three_prime_UTR", "3_prime_utr": "three_prime_UTR",
    }
    with opener(path, "rt") as fh:  # type: ignore[arg-type]
        for raw in fh:
            if not raw or raw.startswith("#"):
                continue
            p = raw.rstrip("\n").split("\t")
            if len(p) < 9:
                continue
            seqid, _, ft, s, e, _, strand, _, attrs = p
            try:
                si, ei = int(s), int(e)
            except ValueError:
                continue
            ft = type_map.get(ft.lower(), ft)
            if (
                seqid == chrom and ei >= start and si <= end
                and ft in {"gene", "mRNA", "transcript", "exon", "CDS", "UTR", "five_prime_UTR", "three_prime_UTR"}
            ):
                features.append(GeneFeature(seqid, si, ei, strand, ft, _parse_attrs(attrs)))
    return features


def _feature_name(f: GeneFeature) -> str:
    for k in ("Name", "gene_name", "gene", "ID", "locus_tag"):
        if f.attrs.get(k):
            return f.attrs[k]
    return "gene"


def _natural_key(x: str):
    return [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", str(x))]


def _class_rows(result: HapResult, level: str, min_count: int) -> Tuple[List[str], Dict[str, int], List[List[str]]]:
    if level == "hap":
        ids = sorted(result.hap_sequences, key=_natural_key)
        counts = {h: len(result.hap_accessions.get(h, [])) for h in ids}
        text = [list(result.hap_sequences[h]) for h in ids]
    else:
        c2h: Dict[str, List[str]] = {}
        for h, c in result.hap_clusters.items():
            c2h.setdefault(c, []).append(h)
        ids = sorted(c2h, key=_natural_key)
        counts = {}
        text = []
        for c in ids:
            hs = c2h[c]
            accs = {a for h in hs for a in result.hap_accessions.get(h, [])}
            counts[c] = len(accs)
            seqs = [result.hap_sequences[h] for h in hs]
            row = []
            for i in range(len(result.variants)):
                toks = sorted({str(s[i]) for s in seqs}, key=_natural_key)
                row.append("/".join(toks))
            text.append(row)
    keep = [i for i, h in enumerate(ids) if counts[h] >= max(1, min_count)]
    return [ids[i] for i in keep], {ids[i]: counts[ids[i]] for i in keep}, [text[i] for i in keep]


def _state_for_token(token: str, ref_token: str) -> int:
    if token in {"", "NA", "N"}:
        return 2
    parts = [x for x in str(token).split("/") if x]
    return 0 if parts and all(x == ref_token for x in parts) else 1


def _allele_matrix(result: HapResult, level: str, min_count: int):
    rows, counts, text = _class_rows(result, level, min_count)
    mat = []
    for row in text:
        mat.append([_state_for_token(tok, call.allele_tokens[0]) for tok, call in zip(row, result.variants)])
    return rows, counts, text, np.array(mat, dtype=float) if mat else np.empty((0, 0))


def _ordered_heatmap_data(result: HapResult, level: str, min_count: int, reverse_variants: bool = False):
    rows, counts, text, matrix = _allele_matrix(result, level, min_count)
    calls = list(result.variants)
    if reverse_variants and matrix.size:
        matrix = matrix[:, ::-1]
        text = [list(reversed(row)) for row in text]
        calls = list(reversed(calls))
    return rows, counts, text, matrix, calls


def _draw_heatmap(
    ax,
    result: HapResult,
    level: str,
    min_count: int,
    ref_color: str,
    alt_color: str,
    missing_color: str,
    title: Optional[str] = None,
    reverse_variants: bool = False,
    show_legend: bool = True,
):
    rows, counts, text, matrix, calls = _ordered_heatmap_data(result, level, min_count, reverse_variants)
    if matrix.size == 0:
        ax.axis("off")
        return None
    nrow, ncol = matrix.shape
    cmap = ListedColormap([ref_color, alt_color, missing_color])
    # pcolormesh gives every cell an explicit black border. Equal aspect guarantees square cells.
    ax.pcolormesh(
        np.arange(ncol + 1), np.arange(nrow + 1), matrix,
        cmap=cmap, vmin=-0.5, vmax=2.5, shading="flat",
        edgecolors="black", linewidth=0.55, antialiased=True,
    )
    ax.set_xlim(0, ncol)
    ax.set_ylim(nrow, 0)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xticks(np.arange(ncol) + 0.5, [str(v.pos) for v in calls], rotation=90)
    ax.set_yticks(np.arange(nrow) + 0.5, [f"{r} ({counts[r]})" for r in rows])
    ax.set_xlabel("Variant position")
    ax.set_ylabel("Haplotype / cluster")
    if title:
        ax.set_title(title)
    if nrow <= 40 and ncol <= 50:
        for i, row in enumerate(text):
            for j, tok in enumerate(row):
                ax.text(j + 0.5, i + 0.5, tok, ha="center", va="center", fontsize=11.5)
    if show_legend:
        ax.legend(
            handles=[
                Patch(facecolor=ref_color, edgecolor="black", label="REF"),
                Patch(facecolor=alt_color, edgecolor="black", label="ALT"),
                Patch(facecolor=missing_color, edgecolor="black", label="Missing"),
            ],
            frameon=False, loc="upper left", bbox_to_anchor=(1.01, 1), borderaxespad=0,
        )
    return {"rows": rows, "counts": counts, "calls": calls, "nrow": nrow, "ncol": ncol}


def _primary_gene(features: List[GeneFeature], start: int, end: int) -> Optional[GeneFeature]:
    genes = [f for f in features if f.ftype == "gene"]
    if genes:
        return max(genes, key=lambda g: max(0, min(g.end, end) - max(g.start, start) + 1))
    tx = [f for f in features if f.ftype in {"mRNA", "transcript"}]
    return max(tx, key=lambda g: max(0, min(g.end, end) - max(g.start, start) + 1)) if tx else None


def _orient(pos: float, start: int, end: int, strand: str) -> float:
    return start + end - pos if strand == "-" else pos


def _subtract_intervals(base_start: int, base_end: int, covered: Sequence[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """Return portions of [base_start, base_end] not covered by the supplied intervals."""
    segs = [(base_start, base_end)]
    for cs, ce in sorted(covered):
        new = []
        for s, e in segs:
            if ce < s or cs > e:
                new.append((s, e))
                continue
            if cs > s:
                new.append((s, cs - 1))
            if ce < e:
                new.append((ce + 1, e))
        segs = new
    return [(s, e) for s, e in segs if s <= e]


def _gene_context(result: HapResult, gff_file: Optional[str]):
    if not gff_file:
        return [], None, "+"
    features = read_gff_features(gff_file, result.region.chrom, result.region.start, result.region.end)
    if not features:
        return [], None, "+"
    gene = _primary_gene(features, result.region.start, result.region.end)
    strand = gene.strand if gene is not None and gene.strand in {"+", "-"} else "+"
    return features, gene, strand


def _draw_gene_structure(ax, result: HapResult, features: List[GeneFeature], gene: Optional[GeneFeature], strand: str, show_legend: bool = True):
    if gene is None:
        ax.hlines(0, result.region.start, result.region.end, color="black", linewidth=1.2, zorder=1)
        ax.text((result.region.start + result.region.end) / 2, 0.32, result.region.vcf_label, ha="center", va="bottom")
        ax.set_xlim(result.region.start, result.region.end)
        ax.set_ylim(-0.30, 0.72)
        ax.axis("off")
        return

    a = _orient(gene.start, result.region.start, result.region.end, strand)
    b = _orient(gene.end, result.region.start, result.region.end, strand)
    left, right = sorted((a, b))
    # The intron is the black backbone. Exon/UTR patches are drawn without outlines so no black line crosses them.
    ax.hlines(0, left, right, color="black", linewidth=1.6, zorder=1)

    children = [
        f for f in features
        if f.ftype in {"exon", "CDS", "UTR", "five_prime_UTR", "three_prime_UTR"}
        and f.end >= gene.start and f.start <= gene.end
    ]
    exons = [f for f in children if f.ftype == "exon"]
    cds = [f for f in children if f.ftype == "CDS"]
    explicit_utrs = [f for f in children if f.ftype in {"UTR", "five_prime_UTR", "three_prime_UTR"}]

    # Draw coding exon sequence and UTR as non-overlapping blocks.  When CDS records are
    # available, only CDS is blue; UTR is drawn separately in grey.  This prevents the
    # previous "grey UTR on top of a blue exon" appearance.
    coding_source = cds if cds else exons
    for f in coding_source:
        x1 = _orient(f.start, result.region.start, result.region.end, strand)
        x2 = _orient(f.end, result.region.start, result.region.end, strand)
        l, r = sorted((x1, x2))
        ax.add_patch(Rectangle((l, -0.12), max(1, r - l + 1), 0.24, facecolor=GENE_CDS_COLOR, edgecolor="none", linewidth=0, zorder=3))

    # Prefer explicit UTR records. If they are absent but exon/CDS are available, infer UTR
    # as exon sequence outside CDS. UTR blocks are deliberately independent of blue CDS blocks.
    utr_intervals: List[Tuple[int, int]] = []
    if explicit_utrs:
        utr_intervals = [(f.start, f.end) for f in explicit_utrs]
    elif exons and cds:
        cds_intervals = [(f.start, f.end) for f in cds]
        for ex in exons:
            utr_intervals.extend(_subtract_intervals(ex.start, ex.end, cds_intervals))

    for s, e in utr_intervals:
        x1 = _orient(s, result.region.start, result.region.end, strand)
        x2 = _orient(e, result.region.start, result.region.end, strand)
        l, r = sorted((x1, x2))
        ax.add_patch(Rectangle((l, -0.09), max(1, r - l + 1), 0.18, facecolor=GENE_UTR_COLOR, edgecolor="none", linewidth=0, zorder=4))

    ax.text(left, 0.36, f"5'  {_feature_name(gene)}", ha="left", va="bottom")
    ax.text(right, 0.36, "3'", ha="right", va="bottom")
    ax.set_xlim(result.region.start, result.region.end)
    ax.set_ylim(-0.30, 0.72)
    ax.axis("off")
    if show_legend:
        ax.legend(
            handles=[
                Patch(facecolor=GENE_UTR_COLOR, edgecolor="none", label="UTR"),
                Patch(facecolor=GENE_CDS_COLOR, edgecolor="none", label="Exon"),
                Line2D([0], [0], color="black", lw=1.6, label="Intron"),
            ],
            loc="upper center", bbox_to_anchor=(0.62, 1.10), frameon=False, ncol=3,
            handlelength=1.8, columnspacing=1.0, borderaxespad=0,
        )


def _connect_gene_to_axis(fig, gene_ax, target_ax, result: HapResult, strand: str, target_x: Sequence[float], target_y: float = 0.0):
    calls = list(reversed(result.variants)) if strand == "-" else list(result.variants)
    if len(calls) != len(target_x):
        return
    for call, tx in zip(calls, target_x):
        gx = _orient(call.pos, result.region.start, result.region.end, strand)
        gene_ax.vlines(gx, -0.18, -0.13, color="black", linewidth=0.7, zorder=5)
        con = ConnectionPatch(
            xyA=(gx, -0.18), coordsA=gene_ax.transData,
            xyB=(tx, target_y), coordsB=target_ax.transData,
            color="0.35", linewidth=0.6, alpha=0.8, zorder=0, clip_on=False,
        )
        fig.add_artist(con)


def _heatmap_figure_size(nrow: int, ncol: int, with_gene: bool = False) -> Tuple[float, float]:
    # Physical dimensions are selected so equal-aspect heatmap cells remain square without forcing the axes narrower.
    width = max(5.2, min(26.0, 2.8 + 0.42 * max(ncol, 1)))
    axis_width = max(3.2, width * 0.66)
    heat_h = max(1.8, axis_width * max(nrow, 1) / max(ncol, 1))
    heat_h = min(28.0, heat_h)
    height = heat_h + (2.1 if with_gene else 1.8)
    return width, max(4.0, height)


def plot_haplotype_heatmap(
    result: HapResult,
    plot_formats: Sequence[str],
    plot_hap_level="hap",
    plot_min_count=1,
    ref_color=DEFAULT_REF_COLOR,
    alt_color=DEFAULT_ALT_COLOR,
    missing_color=DEFAULT_MISSING_COLOR,
):
    rows, _, _, mat = _allele_matrix(result, plot_hap_level, plot_min_count)
    if mat.size == 0:
        return
    fig, ax = plt.subplots(figsize=_heatmap_figure_size(len(rows), len(result.variants), with_gene=False))
    _draw_heatmap(
        ax, result, plot_hap_level, plot_min_count, ref_color, alt_color, missing_color,
        f"Haplotype heatmap: {result.region.vcf_label}",
    )
    fig.subplots_adjust(left=0.18, right=0.82, bottom=0.18, top=0.90)
    stem = "HaplotypeHeatmap" if plot_hap_level == "hap" else "ClusterHaplotypeHeatmap"
    _save_formats(fig, result.output_prefix + "." + stem, plot_formats)


def plot_gene_structure_with_haps(
    result: HapResult,
    gff_file: Optional[str],
    plot_formats: Sequence[str],
    plot_hap_level="hap",
    plot_min_count=1,
    ref_color=DEFAULT_REF_COLOR,
    alt_color=DEFAULT_ALT_COLOR,
    missing_color=DEFAULT_MISSING_COLOR,
):
    features, gene, strand = _gene_context(result, gff_file)
    if gene is None:
        return
    rows, _, _, mat = _allele_matrix(result, plot_hap_level, plot_min_count)
    if mat.size == 0:
        return
    width, height = _heatmap_figure_size(len(rows), len(result.variants), with_gene=True)
    fig = plt.figure(figsize=(width, height))
    gs = fig.add_gridspec(2, 1, height_ratios=[0.60, max(2.0, height - 2.65)], hspace=0.015)
    ag = fig.add_subplot(gs[0])
    ah = fig.add_subplot(gs[1])
    _draw_gene_structure(ag, result, features, gene, strand, show_legend=True)
    meta = _draw_heatmap(
        ah, result, plot_hap_level, plot_min_count, ref_color, alt_color, missing_color,
        None, reverse_variants=(strand == "-"), show_legend=True,
    )
    # Apply final margins first, then match the gene axis to the actual rendered heatmap width.
    fig.subplots_adjust(left=0.18, right=0.82, bottom=0.12, top=0.95)
    fig.canvas.draw()
    gp = ag.get_position(); hp = ah.get_position()
    ag.set_position([hp.x0, gp.y0, hp.width, gp.height])
    if meta:
        _connect_gene_to_axis(fig, ag, ah, result, strand, np.arange(meta["ncol"]) + 0.5, target_y=0.0)
    ah.set_title(f"Gene structure and haplotypes: {result.region.vcf_label}")
    stem = "GeneHaplotype" if plot_hap_level == "hap" else "GeneClusterHaplotype"
    _save_formats(fig, result.output_prefix + "." + stem, plot_formats)


def _group_count_table(df: pd.DataFrame, plot_hap_level="hap", plot_min_count=1) -> pd.DataFrame:
    col = "Hap" if plot_hap_level == "hap" else "ClusterID"
    if df.empty or col not in df or "Type" not in df:
        return pd.DataFrame()
    counts = df.groupby(["Type", col]).size().unstack(fill_value=0)
    keep = counts.sum(axis=0)
    keep = keep[keep >= max(1, plot_min_count)].index
    return counts.loc[:, keep]


def _resolve_hap_palette(hap_palette: Optional[Sequence[str]]) -> List[Tuple[float, float, float, float]]:
    palette = [str(c).strip() for c in (hap_palette or []) if str(c).strip()]
    if not palette:
        palette = list(DEFAULT_HAP_PALETTE)
    return [to_rgba(c) for c in palette]


def _hap_color_map(columns: Sequence[str], hap_palette: Optional[Sequence[str]] = None) -> Dict[str, object]:
    cols = list(columns)
    palette = _resolve_hap_palette(hap_palette)
    return {h: palette[i % len(palette)] for i, h in enumerate(cols)} if cols else {}


def _resolve_ld_cmap(ld_cmap: Optional[str]):
    try:
        return plt.get_cmap(ld_cmap or DEFAULT_LD_CMAP)
    except Exception:
        return plt.get_cmap(DEFAULT_LD_CMAP)


def _draw_group_stacked(ax, df: pd.DataFrame, plot_hap_level="hap", plot_min_count=1, show_legend=True, color_map=None, hap_palette: Optional[Sequence[str]] = None):
    counts = _group_count_table(df, plot_hap_level, plot_min_count)
    if counts.empty:
        ax.axis("off")
        return
    prop = counts.div(counts.sum(axis=1), axis=0)
    if color_map is None:
        color_map = _hap_color_map(prop.columns, hap_palette)
    bottom = np.zeros(len(prop))
    x = np.arange(len(prop))
    for h in prop.columns:
        ax.bar(x, prop[h].values, bottom=bottom, label=h, color=color_map[h], edgecolor="black", linewidth=0.3)
        bottom += prop[h].values
    # n reflects only the haplotypes actually displayed after plot_min_count filtering.
    display_n = counts.sum(axis=1).astype(int)
    tick_labels = [f"{g}\nn={int(display_n.loc[g])}" for g in prop.index]
    ax.set_xticks(x, tick_labels, rotation=20 if len(prop) > 3 else 0)
    ax.set_ylim(0, 1)
    ax.set_ylabel("Haplotype frequency")
    ax.set_xlabel("Group")
    ax.set_title("Haplotype composition by group")
    if show_legend:
        ax.legend(frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")


def plot_group_distribution(result: HapResult, df: pd.DataFrame, plot_formats: Sequence[str], plot_hap_level="hap", plot_min_count=1, hap_palette: Optional[Sequence[str]] = None):
    counts = _group_count_table(df, plot_hap_level, plot_min_count)
    if counts.empty:
        return
    fig, ax = plt.subplots(figsize=(max(7, 0.6 * len(counts.columns) + 4), 5.5))
    _draw_group_stacked(ax, df, plot_hap_level, plot_min_count, show_legend=True, hap_palette=hap_palette)
    _save_formats(fig, result.output_prefix + ".GroupStackedBar", plot_formats)


def _draw_group_pies(fig, container, df: pd.DataFrame, plot_hap_level="hap", plot_min_count=1, title: Optional[str] = None, color_map=None, hap_palette: Optional[Sequence[str]] = None):
    counts = _group_count_table(df, plot_hap_level, plot_min_count)
    if counts.empty:
        ax = fig.add_subplot(container)
        ax.axis("off")
        return
    if color_map is None:
        color_map = _hap_color_map(counts.columns, hap_palette)
    groups = list(counts.index)
    # n reflects only the haplotypes actually displayed after plot_min_count filtering.
    display_n = counts.sum(axis=1).astype(int)
    n = len(groups)
    ncols = min(3, max(1, n))
    nrows = int(math.ceil(n / ncols))
    sub = container.subgridspec(nrows, ncols, wspace=0.28, hspace=0.34)
    for idx, group in enumerate(groups):
        ax = fig.add_subplot(sub[idx // ncols, idx % ncols])
        s = counts.loc[group]
        s = s[s > 0]
        ax.pie(
            s.values,
            colors=[color_map[h] for h in s.index],
            autopct=lambda p: f"{p:.1f}%" if p >= 4 else "",
            startangle=90,
            counterclock=False,
            wedgeprops={"edgecolor": "white", "linewidth": 0.8},
            textprops={"fontsize": 8},
        )
        ax.set_title(str(group), fontsize=10)
        ax.text(0.5, -0.08, f"n={int(display_n.loc[group])}", transform=ax.transAxes, ha="center", va="top", fontsize=9)
    for idx in range(n, nrows * ncols):
        ax = fig.add_subplot(sub[idx // ncols, idx % ncols])
        ax.axis("off")
    handles = [Patch(facecolor=color_map[h], edgecolor="none", label=h) for h in counts.columns]
    fig.legend(handles=handles, loc="center right", bbox_to_anchor=(0.99, 0.5), frameon=False, title="Haplotype")
    if title:
        fig.suptitle(title, fontsize=13)

def plot_group_pie_chart(result: HapResult, df: pd.DataFrame, plot_formats: Sequence[str], plot_hap_level="hap", plot_min_count=1, hap_palette: Optional[Sequence[str]] = None):
    counts = _group_count_table(df, plot_hap_level, plot_min_count)
    if counts.empty:
        return
    n = len(counts.index)
    fig = plt.figure(figsize=(max(7.5, 3.2 * min(3, n) + 2.0), max(5.0, 3.2 * math.ceil(n / min(3, n)))))
    gs = fig.add_gridspec(1, 1)
    color_map = _hap_color_map(counts.columns, hap_palette)
    _draw_group_pies(fig, gs[0], df, plot_hap_level, plot_min_count, title="Haplotype frequency by group", color_map=color_map, hap_palette=hap_palette)
    fig.subplots_adjust(left=0.05, right=0.82, top=0.90, bottom=0.06)
    _save_formats(fig, result.output_prefix + ".GroupPieChart", plot_formats)


def _boxplot_compat(ax, data, tick_labels, **kwargs):
    """Use the version-appropriate Matplotlib boxplot label parameter."""
    try:
        version = tuple(int(x) for x in matplotlib.__version__.split(".")[:2])
    except Exception:
        version = (0, 0)
    if version >= (3, 9):
        return ax.boxplot(data, tick_labels=tick_labels, **kwargs)
    return ax.boxplot(data, labels=tick_labels, **kwargs)


def _p_text(p: float) -> str:
    if not np.isfinite(p):
        return "NA"
    if p < 1e-4:
        return f"{p:.2e}"
    return f"{p:.4f}"


def _sig_stars(q: float) -> str:
    if q < 1e-4:
        return "****"
    if q < 1e-3:
        return "***"
    if q < 1e-2:
        return "**"
    if q < 0.05:
        return "*"
    return "ns"


def _trait_statistics(tmp: pd.DataFrame, col: str, trait: str, order: Sequence[str]):
    try:
        from scipy.stats import kruskal, mannwhitneyu
    except Exception as exc:
        raise RuntimeError("Trait plotting requires scipy>=1.10") from exc
    valid = [(h, tmp.loc[tmp[col] == h, trait].astype(float).values) for h in order]
    testable = [(h, vals) for h, vals in valid if len(vals) >= 2]
    overall_p = float("nan")
    if len(testable) >= 2:
        overall_p = float(kruskal(*[vals for _, vals in testable]).pvalue)
    pairs = []
    for i in range(len(testable)):
        for j in range(i + 1, len(testable)):
            h1, x = testable[i]
            h2, y = testable[j]
            p = float(mannwhitneyu(x, y, alternative="two-sided", method="auto").pvalue)
            pairs.append({"Class1": h1, "Class2": h2, "pvalue": p})
    if pairs:
        qs = bh_adjust([r["pvalue"] for r in pairs])
        for r, q in zip(pairs, qs):
            r["padj_BH"] = float(q)
    return overall_p, pairs


def _draw_trait_axis(ax, df: pd.DataFrame, trait: str, plot_hap_level="hap", plot_min_count=1, show_title=True, hap_palette: Optional[Sequence[str]] = None):
    col = "Hap" if plot_hap_level == "hap" else "ClusterID"
    tmp = df[[col, trait]].copy()
    tmp[trait] = pd.to_numeric(tmp[trait], errors="coerce")
    tmp = tmp.dropna()
    order = [h for h, n in tmp[col].value_counts().items() if n >= max(1, plot_min_count)]
    if len(order) < 1:
        ax.axis("off")
        return False
    data = [tmp.loc[tmp[col] == h, trait].astype(float).values for h in order]
    bp = _boxplot_compat(
        ax, data, order, showfliers=False, patch_artist=True, widths=0.58,
        medianprops={"color": "black", "linewidth": 1.5},
        whiskerprops={"color": "0.35"}, capprops={"color": "0.35"},
    )
    trait_colors = _hap_color_map(order, hap_palette)
    for hap, b in zip(order, bp["boxes"]):
        b.set_facecolor(trait_colors[hap])
        b.set_alpha(0.72)
        b.set_edgecolor("0.25")
        b.set_linewidth(0.9)

    for i, vals in enumerate(data, 1):
        ax.text(i, 0.015, f"n={len(vals)}", transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=8)

    overall_p, pairs = _trait_statistics(tmp, col, trait, order)
    significant = [r for r in pairs if np.isfinite(r.get("padj_BH", np.nan)) and r["padj_BH"] < 0.05]
    significant = sorted(significant, key=lambda r: r["padj_BH"])[:8]

    finite_vals = np.concatenate([v for v in data if len(v)])
    ymin = float(np.nanmin(finite_vals))
    ymax = float(np.nanmax(finite_vals))
    span = ymax - ymin
    if not np.isfinite(span) or span == 0:
        span = max(abs(ymax), 1.0) * 0.2
    base = ymax + 0.08 * span
    step = 0.10 * span
    pos = {h: i + 1 for i, h in enumerate(order)}
    for k, r in enumerate(significant):
        x1, x2 = pos[r["Class1"]], pos[r["Class2"]]
        y = base + k * step
        h = 0.025 * span
        ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=0.8, color="black", clip_on=False)
        ax.text((x1 + x2) / 2, y + h, _sig_stars(r["padj_BH"]), ha="center", va="bottom", fontsize=9)
    if significant:
        ax.set_ylim(ymin - 0.08 * span, base + (len(significant) + 1.2) * step)
    else:
        ax.set_ylim(ymin - 0.08 * span, ymax + 0.18 * span)

    if show_title:
        ax.set_title(f"{trait} by {col}\nKruskal-Wallis P = {_p_text(overall_p)}")
    ax.set_xlabel(col)
    ax.set_ylabel(trait)
    ax.tick_params(axis="x", rotation=35)
    ax.grid(axis="y", linewidth=0.4, alpha=0.25)
    if significant:
        ax.text(0.99, 0.985, "Pairwise Mann-Whitney U, BH-adjusted",
                transform=ax.transAxes, ha="right", va="top", fontsize=7.5)
    return True


def plot_trait_boxplots(
    result: HapResult,
    df: pd.DataFrame,
    traits_to_plot: Optional[Sequence[str]],
    plot_formats: Sequence[str],
    plot_hap_level="hap",
    plot_min_count=1,
    hap_palette: Optional[Sequence[str]] = None,
):
    traits = [c for c in df.columns if c not in {"Hap", "ClusterID", "Accession", "Type"}]
    if traits_to_plot:
        traits = [c for c in traits if c in set(traits_to_plot)]
    for trait in traits:
        col = "Hap" if plot_hap_level == "hap" else "ClusterID"
        n_groups = 1
        if col in df:
            counts = df.groupby(col).size()
            n_groups = max(1, int((counts >= max(1, plot_min_count)).sum()))
        # Keep trait panels compact; very wide boxplots look disproportionate beside the other outputs.
        fig_w = max(5.4, min(7.8, 3.2 + 0.55 * n_groups))
        fig, ax = plt.subplots(figsize=(fig_w, 5.6))
        if not _draw_trait_axis(ax, df, trait, plot_hap_level, plot_min_count, show_title=True, hap_palette=hap_palette):
            plt.close(fig)
            continue
        fig.subplots_adjust(left=0.12, right=0.97, bottom=0.20, top=0.84)
        _save_formats(fig, result.output_prefix + f".{trait}.TraitBoxplot", plot_formats)


def _draw_ld_triangle(ax, mat: pd.DataFrame, labels: Sequence[str], title: Optional[str] = None, show_colorbar=False, fig=None, ld_cmap: Optional[str] = None):
    n = len(labels)
    cmap = _resolve_ld_cmap(ld_cmap)
    norm = Normalize(vmin=0, vmax=1)
    for i in range(n):
        for j in range(i + 1, n):
            val = mat.iat[i, j]
            if not np.isfinite(val):
                face = (0.92, 0.92, 0.92, 1.0)
            else:
                face = cmap(norm(float(val)))
            cx = (i + j) / 2.0
            cy = -(j - i) / 2.0
            poly = Polygon(
                [(cx, cy + 0.5), (cx + 0.5, cy), (cx, cy - 0.5), (cx - 0.5, cy)],
                closed=True, facecolor=face, edgecolor="white", linewidth=0.55,
            )
            ax.add_patch(poly)
    ax.set_xlim(0, max(1, n - 1))
    ax.set_ylim(-max(1.0, n / 2.0), 0.16)
    ax.set_aspect("equal", adjustable="box")
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position("top")
    ax.set_xticks(np.arange(n), labels, rotation=90)
    ax.tick_params(axis="x", pad=3, length=0)
    ax.set_yticks([])
    for side in ("left", "right", "bottom"):
        ax.spines[side].set_visible(False)
    ax.spines["top"].set_visible(False)
    if title:
        ax.text(0.01, 1.03, title, transform=ax.transAxes, ha="left", va="bottom", fontsize=12)
    sm = ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    if show_colorbar and fig is not None:
        fig.colorbar(sm, ax=ax, label="$r^2$", fraction=0.035, pad=0.025)
    return sm


def _ld_ordered(result: HapResult, strand: str):
    from .research import calculate_ld_r2
    mat, _ = calculate_ld_r2(result.variants, sorted(result.sample_haps))
    if strand == "-":
        order = list(reversed(range(len(result.variants))))
        mat = mat.iloc[order, order]
        calls = list(reversed(result.variants))
    else:
        calls = list(result.variants)
    labels = [str(c.pos) for c in calls]
    return mat, calls, labels


def plot_ld_heatmap(result: HapResult, plot_formats: Sequence[str], gff_file: Optional[str] = None, ld_cmap: Optional[str] = None):
    if len(result.variants) < 2:
        return
    features, gene, strand = _gene_context(result, gff_file)
    mat, calls, labels = _ld_ordered(result, strand)
    n = len(calls)
    width = max(6.5, min(28.0, 3.0 + 0.55 * n))
    # A triangular LD panel has an intrinsic width:height ratio close to 2:1. Reserve a
    # separate colorbar column so the LD axis itself is not narrowed relative to the gene axis.
    main_axis_width = width * 0.76
    ld_h = max(2.2, min(14.0, main_axis_width / 2.0))
    fig = plt.figure(figsize=(width, ld_h + 1.35))
    gs = fig.add_gridspec(2, 2, height_ratios=[0.52, ld_h], width_ratios=[1.0, 0.035], hspace=0.08, wspace=0.08)
    ag = fig.add_subplot(gs[0, 0])
    ald = fig.add_subplot(gs[1, 0])
    cax = fig.add_subplot(gs[1, 1])
    fig.add_subplot(gs[0, 1]).axis("off")
    _draw_gene_structure(ag, result, features, gene, strand, show_legend=(gene is not None))
    sm = _draw_ld_triangle(ald, mat, labels, title="LD ($r^2$)", ld_cmap=ld_cmap)
    # Equal-aspect LD drawing can shrink the axis inside its GridSpec cell. Match the gene
    # axis to the final LD plotting box so their left/right edges are exactly aligned.
    fig.subplots_adjust(left=0.12, right=0.90, bottom=0.08, top=0.95)
    fig.canvas.draw()
    gp = ag.get_position(); lp = ald.get_position()
    ag.set_position([lp.x0, gp.y0, lp.width, gp.height])
    _connect_gene_to_axis(fig, ag, ald, result, strand, np.arange(n), target_y=0.0)
    fig.colorbar(sm, cax=cax, label="$r^2$")
    _save_formats(fig, result.output_prefix + ".LD_r2_Heatmap", plot_formats)



def make_all_plots(
    result: HapResult,
    hap_group_df: pd.DataFrame,
    group_map: Dict[str, str],
    gff_file: Optional[str] = None,
    plot_formats: Sequence[str] = ("pdf",),
    traits_to_plot: Optional[Sequence[str]] = None,
    plot_hap_level="hap",
    plot_min_count=1,
    ref_color=DEFAULT_REF_COLOR,
    alt_color=DEFAULT_ALT_COLOR,
    missing_color=DEFAULT_MISSING_COLOR,
    make_ld_plot: bool = True,
    hap_palette: Optional[Sequence[str]] = None,
    ld_cmap: Optional[str] = None,
):
    # Standalone figures.
    plot_haplotype_heatmap(result, plot_formats, plot_hap_level, plot_min_count, ref_color, alt_color, missing_color)
    plot_gene_structure_with_haps(result, gff_file, plot_formats, plot_hap_level, plot_min_count, ref_color, alt_color, missing_color)
    if group_map:
        plot_group_distribution(result, hap_group_df, plot_formats, plot_hap_level, plot_min_count, hap_palette=hap_palette)
        plot_group_pie_chart(result, hap_group_df, plot_formats, plot_hap_level, plot_min_count, hap_palette=hap_palette)
    plot_trait_boxplots(result, hap_group_df, traits_to_plot, plot_formats, plot_hap_level, plot_min_count, hap_palette=hap_palette)
    if make_ld_plot:
        plot_ld_heatmap(result, plot_formats, gff_file=gff_file, ld_cmap=ld_cmap)

