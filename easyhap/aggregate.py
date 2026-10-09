"""Exploratory multi-region category/phenotype comparisons; no covariates."""

from __future__ import annotations

from . import __version__

import glob
import hashlib
import json
import math
import os
import re
import warnings
from typing import Dict, List, Optional, Sequence, Tuple
import matplotlib

matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch
import numpy as np
import pandas as pd
from scipy import stats

from .io_utils import ensure_dir, read_trait_file, read_group_metadata, sanitize_filename
from .plotting import (
    _save_formats, _resolve_hap_palette, _draw_world_background,
    _draw_geo_pies, _draw_geo_bars,
)
from .stats import bh_adjust
from .core import file_fingerprint
from .multilocus import link_known_phase, fit_categorical_interaction

RESERVED_NON_TRAIT = {
    "accession",
    "sample",
    "type",
    "group",
    "hap",
    "clusterid",
    "hap_combination",
    "copy",
    "chrom",
    "phaseset",
    "inputsha256",
    "mode",
}
MISSING_HAPS = {"", ".", "na", "nan", "none", "<na>", "n"}


def _is_valid_numeric_trait(series, col_name):
    if str(col_name).lower() in RESERVED_NON_TRAIT or str(col_name).startswith("Hap_"):
        return False
    values = pd.to_numeric(series, errors="coerce")
    return int(np.isfinite(values).sum()) >= 2


def _normalize_region_token(r):
    r = r.strip()
    match = re.fullmatch(r"(.+):(\d+)-(\d+)", r) or re.fullmatch(r"(.+)_(\d+)_(\d+)", r)
    if match:
        chrom, start, end = match.groups()
        if int(start) < 1 or int(end) < int(start):
            raise ValueError(f"Invalid region: {r}")
        return sanitize_filename(f"{chrom}_{int(start)}_{int(end)}")
    return sanitize_filename(r)


def find_region_hap_group_file(indir, region_str):
    norm = _normalize_region_token(region_str)
    matches = sorted(
        f
        for f in glob.glob(os.path.join(indir, "**", "*.HapGroup.tsv"), recursive=True)
        if _normalize_region_token(os.path.basename(f)[: -len(".HapGroup.tsv")]) == norm
    )
    if not matches:
        raise FileNotFoundError(
            f"No exact .HapGroup.tsv match for {region_str!r} in {indir}"
        )
    if len(matches) != 1:
        raise ValueError(f"Ambiguous region {region_str!r}: {matches}")
    return matches[0]


def _read_hapgroup(path):
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    if not {"Accession", "Hap"}.issubset(df):
        raise ValueError(f"{path} requires Accession and Hap columns")
    for col in ("Accession", "Hap", "Type"):
        if col in df:
            df[col] = df[col].str.strip()
    if df.Accession.eq("").any() or df.Accession.duplicated().any():
        raise ValueError(f"Empty or duplicate sample IDs in {path}")
    if "Type" not in df:
        df["Type"] = "All"
    df["Type"] = df["Type"].replace("", "All")
    return df.set_index("Accession", drop=False)


def combination_colors(df, palette=None):
    keys = sorted(df.Hap_Combination.unique())
    colors = _resolve_hap_palette(palette)
    result = {key: colors[i % len(colors)] for i, key in enumerate(keys)}
    result["Others"] = (0.78, 0.81, 0.85, 1.0)
    return result


def composition_counts(df, min_count=1, max_display=15):
    """Rare/undisplayed classes stay in Others; group denominators stay intact."""
    counts = df.groupby(["Type", "Hap_Combination"]).size().unstack(fill_value=0)
    totals = counts.sum(axis=0).sort_values(ascending=False, kind="stable")
    keep = list(totals[totals >= min_count].index)
    if len(keep) > max_display:
        keep = keep[: max_display - 1]
    other = [c for c in counts if c not in keep]
    result = counts[keep].copy()
    if other:
        result["Others"] = counts[other].sum(axis=1)
    return result


def _aggregate_geography_table(df, group_metadata, min_count=1):
    """Count retained samples at each coordinate, keeping rare classes in Others."""
    columns = ["Latitude", "Longitude", "Location", "Hap_Combination", "Count", "Frequency"]
    if group_metadata is None or group_metadata.empty:
        return pd.DataFrame(columns=columns)
    metadata = group_metadata.copy()
    for col in ["Latitude", "Longitude"]:
        metadata[col] = pd.to_numeric(metadata[col], errors="coerce")
    metadata = metadata.dropna(subset=["Latitude", "Longitude"])
    if "Location" not in metadata:
        metadata["Location"] = ""
    metadata["Location"] = metadata["Location"].fillna("").astype(str)
    mapped = df[["Accession", "Hap_Combination"]].merge(
        metadata[["Accession", "Latitude", "Longitude", "Location"]],
        on="Accession", how="inner", validate="one_to_one",
    )
    if mapped.empty:
        return pd.DataFrame(columns=columns)
    # Select classes using all retained samples, as in the group composition plots.
    display = [c for c in composition_counts(df, min_count).columns if c != "Others"]
    mapped["Hap_Combination"] = mapped["Hap_Combination"].where(
        mapped["Hap_Combination"].isin(display), "Others"
    )
    coordinates = ["Latitude", "Longitude"]
    labels = mapped.groupby(coordinates)["Location"].agg(
        lambda values: "; ".join(sorted({v.strip() for v in values if v.strip()}))
    ).rename("Location").reset_index()
    table = mapped.groupby(coordinates + ["Hap_Combination"]).size().rename("Count").reset_index()
    table = table.merge(labels, on=coordinates, validate="many_to_one")
    table["Frequency"] = table["Count"] / table.groupby(coordinates)["Count"].transform("sum")
    return table[columns]


def plot_aggregate_geography(table, outdir, plot_formats=("pdf",), colors=None, map_style="auto"):
    if map_style not in {"auto", "pie", "bar", "both", "none"}:
        raise ValueError("map_style must be one of auto, pie, bar, both, none")
    outputs = {"map_pie": None, "map_stacked_bar": None}
    if table.empty or map_style == "none":
        return outputs
    styles = ["pie"] if map_style == "auto" else (["pie", "bar"] if map_style == "both" else [map_style])
    colors = colors or combination_colors(table)
    present = set(table.Hap_Combination)
    classes = [c for c in colors if c in present]
    # Drawing helpers expect the class immediately before Count.
    drawing = table[["Latitude", "Longitude", "Location", "Hap_Combination", "Count"]]
    for style in styles:
        fig, ax = plt.subplots(figsize=(10.5, 6.2))
        _draw_world_background(ax, table.Longitude.to_numpy(dtype=float), table.Latitude.to_numpy(dtype=float))
        if style == "pie":
            _draw_geo_pies(ax, drawing, classes, colors)
            key, stem, label = "map_pie", "Aggregate_CombinationMapPie", "pie"
        else:
            _draw_geo_bars(ax, drawing, classes, colors)
            key, stem, label = "map_stacked_bar", "Aggregate_CombinationMapStackedBar", "stacked bar"
        ax.legend(
            handles=[Patch(facecolor=colors[c], edgecolor="none", label=c) for c in classes],
            title="Haplotype combination", frameon=False, bbox_to_anchor=(1.01, 1),
            loc="upper left", fontsize=9, title_fontsize=11,
        )
        ax.set_title(f"Combined haplotype distribution ({label})\n{int(table.Count.sum())} mapped samples", fontsize=15)
        fig.subplots_adjust(left=0.08, right=0.76, bottom=0.10, top=0.88)
        outputs[key] = os.path.join(outdir, stem)
        _save_formats(fig, outputs[key], plot_formats)
    return outputs


def _direction(trait, direction):
    if isinstance(direction, dict):
        value = direction.get(trait, direction.get("*", "higher"))
    else:
        value = direction
    if value not in {"higher", "lower"}:
        raise ValueError(
            "trait_direction must be higher/lower or a trait-to-direction mapping"
        )
    return value


def _safe_test(fn, *args, **kwargs):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        try:
            result = fn(*args, **kwargs)
            return float(result.statistic), float(result.pvalue)
        except ValueError:
            if all(np.all(np.asarray(x) == np.asarray(args[0])[0]) for x in args):
                return 0.0, 1.0
            return np.nan, np.nan


def _aggregate_single_combination(
    indir,
    regions,
    outdir="EasyHap_aggregate",
    trait_file=None,
    trait_cols=None,
    make_plots=True,
    plot_formats=("pdf",),
    min_count=1,
    hap_palette=None,
    test_min_count=2,
    candidate_min_count=5,
    trait_direction="higher",
    phase_link=False,
    interaction=False,
    fdr_scope="run",
    group_file=None,
    map_style="auto",
    group_metadata=None,
):
    ensure_dir(outdir)
    names = [r.strip() for r in regions if r.strip()]
    if len(names) < 2 or len({_normalize_region_token(r) for r in names}) != len(names):
        raise ValueError("Aggregation requires at least two distinct regions")
    if min_count < 1 or test_min_count < 2 or candidate_min_count < 2:
        raise ValueError(
            "plot min-count >=1; test/candidate min-count >=2 are required"
        )
    if fdr_scope not in {"run", "trait"}:
        raise ValueError("fdr_scope must be run or trait")
    paths = [find_region_hap_group_file(indir, r) for r in names]
    dfs = [_read_hapgroup(p) for p in paths]
    all_ids = sorted(set().union(*(set(df.index) for df in dfs)))
    exclusions, retained = [], []
    all_missing = []
    for path in paths:
        qc_path = path[: -len(".HapGroup.tsv")] + ".SampleQC.tsv"
        if os.path.isfile(qc_path):
            qc = pd.read_csv(qc_path, sep="\t", dtype=str, keep_default_na=False)
            if (
                not {"Accession", "AllMissing"}.issubset(qc)
                or qc.Accession.duplicated().any()
            ):
                raise ValueError(f"Invalid sample QC metadata: {qc_path}")
            all_missing.append(
                set(qc.loc[qc.AllMissing.str.lower().isin(["true", "1"]), "Accession"])
            )
        else:
            all_missing.append(set())
    for sid in all_ids:
        reasons = []
        for name, df, missing in zip(names, dfs, all_missing):
            if sid not in df.index:
                reasons.append(f"{name}:absent_sample")
            elif df.at[sid, "Hap"].lower() in MISSING_HAPS:
                reasons.append(f"{name}:missing_category")
            elif sid in missing:
                reasons.append(f"{name}:all_markers_missing")
        if reasons:
            exclusions.append({"Accession": sid, "Reason": ";".join(reasons)})
        else:
            retained.append(sid)
    pd.DataFrame(exclusions, columns=["Accession", "Reason"]).to_csv(
        os.path.join(outdir, "Aggregate_ExcludedSamples.tsv"), sep="\t", index=False
    )
    if not retained:
        raise ValueError("No complete, valid overlapping samples across regions")
    merged = pd.DataFrame({"Accession": retained})
    hap_cols = [f"Hap_{r}" for r in names]
    metadata = group_metadata if group_metadata is not None else read_group_metadata(group_file)
    supplied_groups = metadata.set_index("Accession")["Type"].to_dict() if not metadata.empty else {}
    if group_file and not set(retained).intersection(supplied_groups):
        raise ValueError("Group file has no matching retained sample IDs")
    groups = []
    for sid in retained:
        if sid in supplied_groups:
            groups.append(supplied_groups[sid])
            continue
        labels = {df.at[sid, "Type"] for df in dfs} - {"All"}
        if len(labels) > 1:
            raise ValueError(f"Conflicting group labels for {sid}: {sorted(labels)}")
        groups.append(next(iter(labels)) if labels else "All")
    merged["Type"] = groups
    for df, col in zip(dfs, hap_cols):
        merged[col] = [df.at[sid, "Hap"] for sid in retained]
    if any(merged[col].str.contains(r"[+|]", regex=True).any() for col in hap_cols):
        raise ValueError(
            "Hap categories cannot contain reserved combination separators + or |"
        )
    merged["Hap_Combination"] = merged[hap_cols].agg("+".join, axis=1)
    # An explicit trait file is authoritative; otherwise coalesce consistent embedded values.
    available = []
    if trait_file:
        traits_df = read_trait_file(trait_file)
        if not set(retained) & set(traits_df.Accession):
            raise ValueError("Trait file has no matching sample IDs")
        traits_df = traits_df.set_index("Accession")
        candidates = [c for c in traits_df if _is_valid_numeric_trait(traits_df[c], c)]
        for col in candidates:
            merged[col] = pd.to_numeric(
                traits_df[col].reindex(retained), errors="coerce"
            ).to_numpy()
        available = candidates
    else:
        candidates = sorted(
            {c for df in dfs for c in df if _is_valid_numeric_trait(df[c], c)}
        )
        for col in candidates:
            values = pd.Series(np.nan, index=retained)
            for df in dfs:
                if col not in df:
                    continue
                new = pd.to_numeric(df[col].reindex(retained), errors="coerce")
                new = new.where(np.isfinite(new))
                common = values.notna() & new.notna()
                if not np.allclose(values[common], new[common], rtol=1e-9, atol=1e-12):
                    raise ValueError(
                        f"Conflicting embedded phenotype values for {col}; use an explicit --traits file"
                    )
                values = values.where(values.notna(), new)
            merged[col] = values.to_numpy()
        available = candidates
    if trait_cols:
        unknown = set(trait_cols) - set(available)
        if unknown:
            raise ValueError(
                f"Requested traits absent or not numeric: {sorted(unknown)}"
            )
        available = list(trait_cols)
    if isinstance(trait_direction, dict) and set(trait_direction) - {"*"} - set(
        available
    ):
        raise ValueError("Trait direction mapping contains unavailable traits")
    for trait in available:
        merged[trait] = pd.to_numeric(merged[trait], errors="coerce").where(
            lambda x: np.isfinite(x)
        )
    combination_path = os.path.join(outdir, "Aggregate_Haplotype_Combinations.tsv")
    merged[["Accession", "Type"] + hap_cols + ["Hap_Combination"] + available].to_csv(
        combination_path, sep="\t", index=False
    )
    freq = merged.groupby(["Type", "Hap_Combination"]).size().rename("N").reset_index()
    freq["Frequency"] = freq.N / freq.groupby("Type").N.transform("sum")
    freq.to_csv(
        os.path.join(outdir, "Aggregate_Combination_Frequency.tsv"),
        sep="\t",
        index=False,
    )
    retention = [
        {"Stage": "union_of_region_samples", "N": len(all_ids)},
        {"Stage": "complete_valid_categories", "N": len(merged)},
        {"Stage": "excluded_categories", "N": len(exclusions)},
    ]
    colors = combination_colors(merged, hap_palette)
    summary, tests = [], []
    for trait in available:
        sub = merged.loc[np.isfinite(merged[trait]), ["Hap_Combination", trait]]
        retention.append({"Stage": f"{trait}:finite_phenotype", "N": len(sub)})
        if sub.empty:
            continue
        overall = float(sub[trait].mean())
        groups = {
            c: g[trait].to_numpy(dtype=float) for c, g in sub.groupby("Hap_Combination")
        }
        eligible = {c: v for c, v in groups.items() if len(v) >= test_min_count}
        test_n = sum(map(len, eligible.values()))
        retention.append({"Stage": f"{trait}:eligible_for_tests", "N": test_n})
        overall_p = {}
        for label, fn in [("Kruskal-Wallis", stats.kruskal), ("ANOVA", stats.f_oneway)]:
            stat, p = (
                _safe_test(fn, *eligible.values())
                if len(eligible) >= 2
                else (np.nan, np.nan)
            )
            tests.append(
                {
                    "Trait": trait,
                    "Comparison": "Overall",
                    "Combination1": "ALL",
                    "Combination2": "",
                    "Test": label,
                    "Primary": label == "Kruskal-Wallis",
                    "N1": test_n,
                    "N2": 0,
                    "N_classes": len(eligible),
                    "Statistic": stat,
                    "pvalue": p,
                    "Status": "OK" if np.isfinite(p) else "NOT_ESTIMABLE",
                }
            )
            overall_p[label] = p
        keys = sorted(eligible)
        for i, c1 in enumerate(keys):
            for c2 in keys[i + 1 :]:
                for label, fn, kwargs in [
                    (
                        "Mann-Whitney U",
                        stats.mannwhitneyu,
                        {"alternative": "two-sided"},
                    ),
                    ("Welch t-test", stats.ttest_ind, {"equal_var": False}),
                ]:
                    stat, p = _safe_test(fn, eligible[c1], eligible[c2], **kwargs)
                    tests.append(
                        {
                            "Trait": trait,
                            "Comparison": "Pairwise",
                            "Combination1": c1,
                            "Combination2": c2,
                            "Test": label,
                            "Primary": label == "Mann-Whitney U",
                            "N1": len(eligible[c1]),
                            "N2": len(eligible[c2]),
                            "N_classes": 2,
                            "Statistic": stat,
                            "pvalue": p,
                            "Status": "OK" if np.isfinite(p) else "NOT_ESTIMABLE",
                        }
                    )
        direction = _direction(trait, trait_direction)
        ranked = sorted(
            groups,
            key=lambda c: (
                (-1 if direction == "higher" else 1) * np.mean(groups[c]),
                -len(groups[c]),
                c,
            ),
        )
        candidate_rank = 0
        for rank, c in enumerate(ranked, 1):
            vals = groups[c]
            n = len(vals)
            mean = float(vals.mean())
            sd = float(vals.std(ddof=1)) if n > 1 else np.nan
            se = sd / math.sqrt(n)
            width = float(stats.t.ppf(0.975, n - 1)) * se if n > 1 else np.nan
            candidate = n >= candidate_min_count
            if candidate:
                candidate_rank += 1
            summary.append(
                {
                    "Trait": trait,
                    "Rank": candidate_rank if candidate else np.nan,
                    "DescriptiveRank": rank,
                    "CandidateRank": candidate_rank if candidate else np.nan,
                    "Candidate": candidate,
                    "Direction": direction,
                    "Hap_Combination": c,
                    "N": n,
                    "Mean": mean,
                    "SD": sd,
                    "SE": se,
                    "CI95_Lower": mean - width,
                    "CI95_Upper": mean + width,
                    "Median": float(np.median(vals)),
                    "IQR": float(stats.iqr(vals)),
                    "Min": float(vals.min()),
                    "Max": float(vals.max()),
                    "Overall_Mean": overall,
                    "Effect": mean - overall,
                    "Effect_reference": "all_retained_samples_with_finite_trait",
                    "Relative_Effect_Percent": 100 * (mean - overall) / abs(overall)
                    if overall
                    else np.nan,
                    "TestEligible": c in eligible,
                    "Test_N": test_n,
                    "Test_classes": len(eligible),
                    "ANOVA_p": overall_p["ANOVA"],
                    "Kruskal_p": overall_p["Kruskal-Wallis"],
                }
            )
        if make_plots:
            shown = sorted(
                (c for c in groups if len(groups[c]) >= min_count),
                key=lambda c: (-len(groups[c]), c),
            )[:20]
            if shown:
                fig, ax = plt.subplots(
                    figsize=(max(5.6, min(16.0, 3 + 0.95 * len(shown))), 5.6)
                )
                bp = ax.boxplot(
                    [groups[c] for c in shown], patch_artist=True, showmeans=True
                )
                rng = np.random.default_rng(2026)
                for i, (patch, c) in enumerate(zip(bp["boxes"], shown), 1):
                    patch.set_facecolor(colors[c])
                    patch.set_alpha(0.85)
                    ax.scatter(
                        i + rng.normal(0, 0.04, len(groups[c])),
                        groups[c],
                        s=15,
                        color=".25",
                        alpha=0.45,
                    )
                ax.set_xticks(
                    range(1, len(shown) + 1),
                    [f"{c}\n(n={len(groups[c])})" for c in shown],
                    rotation=25,
                    ha="right",
                )
                ax.set_ylabel(trait)
                ax.set_xlabel("Multi-region category combination")
                ax.set_title(
                    f"Exploratory comparison: {trait}\nKruskal-Wallis p={overall_p['Kruskal-Wallis']:.3g}; test N={test_n}"
                )
                ax.grid(axis="y", linestyle=":", alpha=0.4)
                fig.subplots_adjust(left=0.12, right=0.96, bottom=0.28, top=0.87)
                _save_formats(
                    fig,
                    os.path.join(
                        outdir, f"Aggregate_{sanitize_filename(trait)}_Boxplot"
                    ),
                    plot_formats,
                )
    summary_path = os.path.join(outdir, "Aggregate_Trait_Summary.tsv")
    tests_path = os.path.join(outdir, "Aggregate_Statistical_Tests.tsv")
    summary_df = pd.DataFrame(
        summary,
        columns=list(summary[0])
        if summary
        else ["Trait", "Rank", "CandidateRank", "Hap_Combination", "N", "Mean"],
    )
    if summary:
        summary_df = summary_df.sort_values(
            ["Trait", "Candidate", "CandidateRank", "DescriptiveRank"],
            ascending=[True, False, True, True],
        )
    summary_df.to_csv(summary_path, sep="\t", index=False)
    pd.DataFrame(
        tests,
        columns=list(tests[0])
        if tests
        else ["Trait", "Comparison", "Test", "pvalue", "Status"],
    ).to_csv(tests_path, sep="\t", index=False)
    pd.DataFrame(retention).to_csv(
        os.path.join(outdir, "Aggregate_SampleRetention.tsv"), sep="\t", index=False
    )
    res = {
        "combinations": combination_path,
        "summary": summary_path,
        "tests": tests_path,
        "outdir": outdir,
        "regions": names,
        "n_samples": len(merged),
        "stacked_bar": None,
        "pie_chart": None,
        "geography": None,
        "map_pie": None,
        "map_stacked_bar": None,
        "n_geographic_samples": 0,
    }
    geography = _aggregate_geography_table(merged, metadata, min_count)
    if not geography.empty:
        res["geography"] = os.path.join(outdir, "Aggregate_Combination_Geography.tsv")
        res["n_geographic_samples"] = int(geography.Count.sum())
        geography.to_csv(res["geography"], sep="\t", index=False)
        print(f"[EasyHap Aggregate] Geography: {res['n_geographic_samples']}/{len(merged)} retained samples have coordinates")
    elif group_file and make_plots and map_style != "none":
        print("[EasyHap Aggregate] Geography skipped: no retained samples with complete coordinates")
    if phase_link:
        res["phase_link"] = link_known_phase(merged, paths, outdir)
    if interaction:
        if len(hap_cols) == 2:
            model, contrasts = fit_categorical_interaction(
                merged, *hap_cols, available, min_cell_n=test_min_count
            )
        else:
            model = pd.DataFrame(
                [
                    {
                        "Trait": t,
                        "Status": "NOT_ESTIMABLE",
                        "Reason": "requires_exactly_two_regions",
                        "pvalue": np.nan,
                    }
                    for t in available
                ],
                columns=["Trait", "Status", "Reason", "pvalue"],
            )
            contrasts = pd.DataFrame(columns=["Trait", "pvalue"])
        res["interaction"] = os.path.join(outdir, "Aggregate_InteractionTests.tsv")
        res["interaction_contrasts"] = os.path.join(
            outdir, "Aggregate_InteractionContrasts.tsv"
        )
        model.to_csv(res["interaction"], sep="\t", index=False)
        contrasts.to_csv(res["interaction_contrasts"], sep="\t", index=False)
    if make_plots:
        res["stacked_bar"] = os.path.join(outdir, "Aggregate_GroupStackedBar")
        res["pie_chart"] = os.path.join(outdir, "Aggregate_GroupPieChart")
        plot_aggregate_group_distribution(
            merged, res["stacked_bar"], plot_formats, min_count, hap_palette, colors
        )
        plot_aggregate_group_pie_chart(
            merged, res["pie_chart"], plot_formats, min_count, hap_palette, colors
        )
        res.update(plot_aggregate_geography(geography, outdir, plot_formats, colors, map_style))
    manifest = {
        "version": __version__,
        "purpose": "exploratory_category_phenotype_comparison",
        "regions": names,
        "source_sha256": {os.path.abspath(p): file_fingerprint(p) for p in paths},
        "trait_source": os.path.abspath(trait_file)
        if trait_file
        else "consistent_embedded_traits",
        "trait_source_sha256": file_fingerprint(trait_file) if trait_file else None,
        "group_source": os.path.abspath(group_file) if group_file else "consistent_embedded_groups",
        "group_source_sha256": file_fingerprint(group_file) if group_file else None,
        "map_style": map_style,
        "n_geographic_samples": res["n_geographic_samples"],
        "map_rare_classes": "Others; same global selection and colors as group composition plots",
        "n_union_samples": len(all_ids),
        "n_retained": len(merged),
        "n_excluded": len(exclusions),
        "plot_min_count": min_count,
        "test_min_count": test_min_count,
        "candidate_min_count": candidate_min_count,
        "trait_direction": trait_direction,
        "fdr_scope": fdr_scope,
        "fdr_families": "separate_test_types; overall_and_pairwise_separate",
        "primary_tests": ["Kruskal-Wallis", "Mann-Whitney U"],
        "secondary_tests": ["ANOVA", "Welch t-test"],
        "phase_link": phase_link,
        "interaction": interaction,
        "covariates": [],
        "all_missing_detection": [
            os.path.isfile(p[: -len(".HapGroup.tsv")] + ".SampleQC.tsv") for p in paths
        ],
    }
    with open(
        os.path.join(outdir, "Aggregate_Manifest.json"), "w", encoding="utf-8"
    ) as fh:
        json.dump(manifest, fh, ensure_ascii=False, indent=2)
    print(
        f"[EasyHap Aggregate] {len(names)} regions; {len(merged)} retained; {len(exclusions)} excluded -> {outdir}"
    )
    return res


def plot_aggregate_group_distribution(
    df: pd.DataFrame,
    out_stem: str,
    plot_formats: Sequence[str] = ("pdf",),
    min_count: int = 1,
    hap_palette: Optional[Sequence[str]] = None,
    color_map=None,
) -> None:
    """Draw 100% stacked bar chart for haplotype combinations across population groups."""
    if df.empty or "Hap_Combination" not in df.columns or "Type" not in df.columns:
        return
    df_clean = df.copy()
    df_clean["Type"] = (
        df_clean["Type"]
        .fillna("NA")
        .astype(str)
        .replace({"": "NA", "nan": "NA", "None": "NA", "<NA>": "NA"})
    )
    df_clean["Hap_Combination"] = (
        df_clean["Hap_Combination"]
        .fillna("NA")
        .astype(str)
        .replace({"": "NA", "nan": "NA", "None": "NA", "<NA>": "NA"})
    )

    counts = (
        df_clean.groupby(["Type", "Hap_Combination"], dropna=False)
        .size()
        .unstack(fill_value=0)
    )
    counts_plot = composition_counts(df_clean, min_count)
    prop = counts_plot.div(counts_plot.sum(axis=1), axis=0)
    color_map = color_map or combination_colors(df_clean, hap_palette)

    fig_width = max(5.5, min(12.0, 3.8 + 1.2 * len(prop.index)))
    fig, ax = plt.subplots(figsize=(fig_width, 5.5))
    bottom = np.zeros(len(prop))
    x = np.arange(len(prop))
    bar_width = 0.55
    for h in prop.columns:
        ax.bar(
            x,
            prop[h].values,
            bottom=bottom,
            width=bar_width,
            label=h,
            color=color_map[h],
            edgecolor="none",
        )
        bottom += prop[h].values
    for xi in x:
        ax.add_patch(
            Rectangle(
                (float(xi) - bar_width / 2.0, 0.0),
                bar_width,
                1.0,
                fill=False,
                edgecolor="black",
                linewidth=0.35,
                zorder=4,
            )
        )
    display_n = counts_plot.sum(axis=1).astype(int)
    ax.set_xticks(x)
    ax.set_xticklabels(
        [f"{g}\n(n={int(display_n.loc[g]):,})" for g in prop.index], fontsize=12
    )
    ax.set_ylim(0, 1)
    ax.set_ylabel("Haplotype combination frequency", fontsize=14)
    ax.set_xlabel("Group", fontsize=14)
    ax.set_title("Haplotype combination composition by group", fontsize=15)
    ax.legend(
        frameon=False,
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        fontsize=10,
        title="Combination",
        title_fontsize=11,
    )
    fig.subplots_adjust(left=0.12, right=0.74, bottom=0.18, top=0.90)
    _save_formats(fig, out_stem, plot_formats)


def plot_aggregate_group_pie_chart(
    df: pd.DataFrame,
    out_stem: str,
    plot_formats: Sequence[str] = ("pdf",),
    min_count: int = 1,
    hap_palette: Optional[Sequence[str]] = None,
    color_map=None,
) -> None:
    """Draw multi-panel pie charts for haplotype combinations across population groups."""
    if df.empty or "Hap_Combination" not in df.columns or "Type" not in df.columns:
        return
    df_clean = df.copy()
    df_clean["Type"] = (
        df_clean["Type"]
        .fillna("NA")
        .astype(str)
        .replace({"": "NA", "nan": "NA", "None": "NA", "<NA>": "NA"})
    )
    df_clean["Hap_Combination"] = (
        df_clean["Hap_Combination"]
        .fillna("NA")
        .astype(str)
        .replace({"": "NA", "nan": "NA", "None": "NA", "<NA>": "NA"})
    )

    counts = (
        df_clean.groupby(["Type", "Hap_Combination"], dropna=False)
        .size()
        .unstack(fill_value=0)
    )
    counts_plot = composition_counts(df_clean, min_count)
    groups = list(counts_plot.index)
    display_n = counts_plot.sum(axis=1).astype(int)
    n_groups = len(groups)
    ncols = min(3, max(1, n_groups))
    nrows = int(math.ceil(n_groups / ncols))

    color_map = color_map or combination_colors(df_clean, hap_palette)

    fig = plt.figure(figsize=(max(7.5, 3.4 * ncols + 2.5), max(5.0, 3.4 * nrows)))
    gs = fig.add_gridspec(1, 1)
    sub = gs[0].subgridspec(nrows, ncols, wspace=0.28, hspace=0.48)
    for idx, group in enumerate(groups):
        ax = fig.add_subplot(sub[idx // ncols, idx % ncols])
        s = counts_plot.loc[group]
        s = s[s > 0]
        ax.pie(
            s.values,
            colors=[color_map[h] for h in s.index],
            autopct=lambda p: f"{p:.1f}%" if p >= 4 else "",
            startangle=90,
            counterclock=False,
            wedgeprops={"edgecolor": "none", "linewidth": 0},
            textprops={"fontsize": 11},
        )
        ax.add_patch(
            plt.Circle(
                (0.0, 0.0), 1.0, fill=False, edgecolor="black", linewidth=0.30, zorder=6
            )
        )
        group_label = str(group)
        n_count = (
            int(display_n.loc[group]) if group in display_n.index else int(s.sum())
        )
        ax.text(
            0.5,
            -0.06,
            f"{group_label}\n(n={n_count:,})",
            transform=ax.transAxes,
            ha="center",
            va="top",
            fontsize=13,
        )

    for idx in range(n_groups, nrows * ncols):
        fig.add_subplot(sub[idx // ncols, idx % ncols]).axis("off")

    handles = [
        Patch(facecolor=color_map[h], edgecolor="none", label=h)
        for h in counts_plot.columns
    ]
    fig.legend(
        handles=handles,
        loc="center right",
        bbox_to_anchor=(0.99, 0.5),
        frameon=False,
        title="Combination",
        fontsize=10,
        title_fontsize=11,
    )
    fig.suptitle("Haplotype combination frequency by group", fontsize=16)
    fig.subplots_adjust(left=0.05, right=0.76, top=0.90, bottom=0.08)
    _save_formats(fig, out_stem, plot_formats)


def _adjust_result_families(results, scope):
    for key in ("tests", "interaction", "interaction_contrasts"):
        pieces = []
        for i, res in enumerate(results):
            if key not in res:
                continue
            df = pd.read_csv(
                res[key], sep="\t", dtype={"Trait": str}, keep_default_na=False
            )
            df["pvalue"] = pd.to_numeric(df["pvalue"], errors="coerce")
            df["_result"] = i
            pieces.append(df)
        if not pieces:
            continue
        combined = pd.concat(pieces, ignore_index=True)
        combined["padj_BH"] = np.nan
        combined["FDR_scope"] = scope
        family_cols = (["Test", "Comparison"] if key == "tests" else []) + (
            ["Trait"] if scope == "trait" else []
        )
        if family_cols:
            for _, indices in combined.groupby(
                family_cols, dropna=False
            ).groups.items():
                combined.loc[indices, "padj_BH"] = bh_adjust(
                    combined.loc[indices, "pvalue"].tolist()
                )
        else:
            combined["padj_BH"] = bh_adjust(combined.pvalue.tolist())
        for i, res in enumerate(results):
            if key in res:
                combined[combined._result == i].drop(columns="_result").to_csv(
                    res[key], sep="\t", index=False
                )


def run_aggregate(
    indir,
    regions=None,
    region_file=None,
    outdir="EasyHap_aggregate",
    trait_file=None,
    trait_cols=None,
    make_plots=True,
    plot_formats=("pdf",),
    min_count=1,
    hap_palette=None,
    test_min_count=2,
    candidate_min_count=5,
    trait_direction="higher",
    phase_link=False,
    interaction=False,
    fdr_scope="run",
    group_file=None,
    map_style="auto",
):
    if map_style not in {"auto", "pie", "bar", "both", "none"}:
        raise ValueError("map_style must be one of auto, pie, bar, both, none")
    group_metadata = read_group_metadata(group_file)
    combinations = []
    if region_file:
        with open(region_file, encoding="utf-8") as fh:
            for line_no, line in enumerate(fh, 1):
                if not line.strip() or line.startswith("#"):
                    continue
                parts = [
                    r.strip() for r in re.split(r"[\t,]+", line.strip()) if r.strip()
                ]
                if len(parts) < 2:
                    raise ValueError(
                        f"{region_file}:{line_no} requires two or more region columns"
                    )
                combinations.append(parts)
    if regions:
        if isinstance(regions, str):
            regions = [x for x in re.split(r"[,\s]+", regions) if x]
        combinations.append(list(regions))
    if not combinations:
        raise ValueError("Provide --regions or --region-file with at least two regions")
    canonical = [
        tuple(sorted(_normalize_region_token(r) for r in c)) for c in combinations
    ]
    if len(set(canonical)) != len(canonical):
        raise ValueError("Duplicate region combinations in aggregate batch")
    ensure_dir(outdir)
    results = []
    for i, combination in enumerate(combinations, 1):
        destination = (
            outdir
            if len(combinations) == 1
            else os.path.join(
                outdir, f"combination_{i}_" + sanitize_filename("_".join(combination))
            )
        )
        results.append(
            _aggregate_single_combination(
                indir,
                combination,
                destination,
                trait_file,
                trait_cols,
                make_plots,
                plot_formats,
                min_count,
                hap_palette,
                test_min_count,
                candidate_min_count,
                trait_direction,
                phase_link,
                interaction,
                fdr_scope,
                group_file=group_file,
                map_style=map_style,
                group_metadata=group_metadata,
            )
        )
    _adjust_result_families(results, fdr_scope)
    index_path = os.path.join(outdir, "Aggregate_BatchIndex.tsv")
    pd.DataFrame(
        [
            {
                "Combination": i,
                "Regions": ",".join(r["regions"]),
                "N_samples": r["n_samples"],
                "Outdir": r["outdir"],
                "Summary": r["summary"],
                "Tests": r["tests"],
                "Geography": r["geography"],
                "MapPie": r["map_pie"],
                "MapStackedBar": r["map_stacked_bar"],
            }
            for i, r in enumerate(results, 1)
        ]
    ).to_csv(index_path, sep="\t", index=False)
    if len(results) == 1:
        return dict(results[0], results=results, index=index_path, status="OK")
    return {"results": results, "index": index_path, "status": "OK"}
