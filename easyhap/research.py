from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple
import math
import os

import numpy as np
import pandas as pd

from .stats import bh_adjust
from .vcf_reader import VariantCall


@dataclass
class ResearchOutputs:
    frequency: str
    diversity: str
    private_haplotypes: str
    trait_summary: Optional[str] = None
    trait_tests: Optional[str] = None
    superior_haplotypes: Optional[str] = None
    ld_matrix: Optional[str] = None
    ld_pairs: Optional[str] = None


def _sample_hap_records(
    sample_haps: Dict[str, List[str]], group_map: Dict[str, str], mode: str
) -> pd.DataFrame:
    rows = []
    for sample, haps in sample_haps.items():
        group = group_map.get(sample, "All")
        if str(mode).lower() in {"inbred", "genotype"}:
            hap_list = haps[:1]
        else:
            hap_list = list(haps)
        for copy_i, hap in enumerate(hap_list, 1):
            rows.append(
                {"Accession": sample, "Group": group, "Copy": copy_i, "Haplotype": hap}
            )
    return pd.DataFrame(rows)


def write_population_haplotype_statistics(
    prefix: str,
    sample_haps: Dict[str, List[str]],
    group_map: Dict[str, str],
    mode: str,
) -> Tuple[str, str, str]:
    rec = _sample_hap_records(sample_haps, group_map, mode)
    freq_path = prefix + ".HaplotypeFrequency.tsv"
    div_path = prefix + ".HaplotypeDiversity.tsv"
    private_path = prefix + ".PrivateHaplotypes.tsv"

    if rec.empty:
        pd.DataFrame(columns=["Group", "Haplotype", "Count", "Frequency"]).to_csv(
            freq_path, sep="\t", index=False
        )
        pd.DataFrame(
            columns=[
                "Group",
                "N_haplotype_copies",
                "N_haplotypes",
                "Haplotype_diversity_Hd",
            ]
        ).to_csv(div_path, sep="\t", index=False)
        pd.DataFrame(columns=["Group", "Haplotype", "Count", "Frequency"]).to_csv(
            private_path, sep="\t", index=False
        )
        return freq_path, div_path, private_path

    counts = rec.groupby(["Group", "Haplotype"]).size().rename("Count").reset_index()
    totals = counts.groupby("Group")["Count"].transform("sum")
    counts["Frequency"] = counts["Count"] / totals
    counts.sort_values(
        ["Group", "Count", "Haplotype"], ascending=[True, False, True]
    ).to_csv(freq_path, sep="\t", index=False)

    div_rows = []
    for group, sub in counts.groupby("Group"):
        n = int(sub["Count"].sum())
        p = sub["Count"].astype(float) / max(n, 1)
        hd = float("nan") if n <= 1 else (n / (n - 1.0)) * (1.0 - float(np.sum(p * p)))
        div_rows.append(
            {
                "Group": group,
                "N_haplotype_copies": n,
                "N_haplotypes": int((sub["Count"] > 0).sum()),
                "Haplotype_diversity_Hd": hd,
            }
        )
    pd.DataFrame(div_rows).to_csv(div_path, sep="\t", index=False)

    presence = counts.pivot_table(
        index="Haplotype", columns="Group", values="Count", fill_value=0, aggfunc="sum"
    )
    private_rows = []
    # Private haplotypes are only meaningful when two or more biological groups are supplied.
    if len(presence.columns) >= 2:
        for hap, row in presence.iterrows():
            groups = [g for g, v in row.items() if v > 0]
            if len(groups) == 1:
                g = groups[0]
                sub = counts[
                    (counts["Group"] == g) & (counts["Haplotype"] == hap)
                ].iloc[0]
                private_rows.append(
                    {
                        "Group": g,
                        "Haplotype": hap,
                        "Count": int(sub["Count"]),
                        "Frequency": float(sub["Frequency"]),
                    }
                )
    pd.DataFrame(
        private_rows, columns=["Group", "Haplotype", "Count", "Frequency"]
    ).to_csv(private_path, sep="\t", index=False)
    return freq_path, div_path, private_path


def write_trait_association(
    prefix: str,
    hap_group_df: pd.DataFrame,
    trait_cols: Optional[Sequence[str]] = None,
    class_col: str = "Hap",
    min_class_n: int = 2,
) -> Tuple[Optional[str], Optional[str], Optional[str]]:
    """Report only classes supported by the plotting threshold for each trait."""
    if not isinstance(min_class_n, int) or isinstance(min_class_n, bool) or min_class_n < 1:
        raise ValueError("min_class_n must be a positive integer")
    # Singleton means are not reported as candidates or used in inference.
    min_class_n = max(2, min_class_n)
    if hap_group_df.empty or class_col not in hap_group_df.columns:
        return None, None, None
    available = [
        c
        for c in hap_group_df.columns
        if c not in {"Hap", "ClusterID", "Accession", "Type"}
    ]
    if trait_cols:
        wanted = set(trait_cols)
        available = [c for c in available if c in wanted]
    if not available:
        return None, None, None

    try:
        from scipy.stats import kruskal, mannwhitneyu
    except Exception as exc:
        raise RuntimeError("Trait association requires scipy>=1.10") from exc

    summary_rows: List[dict] = []
    test_rows: List[dict] = []
    superior_rows: List[dict] = []

    for trait in available:
        tmp = hap_group_df[[class_col, trait]].copy()
        tmp[trait] = pd.to_numeric(tmp[trait], errors="coerce")
        tmp = tmp[np.isfinite(tmp[trait])]
        counts = tmp[class_col].value_counts()
        eligible = counts[counts >= min_class_n].index
        tmp = tmp[tmp[class_col].isin(eligible)]
        if tmp.empty:
            continue
        overall_mean = float(tmp[trait].mean())
        class_stats = []
        for cls, sub in tmp.groupby(class_col):
            vals = sub[trait].astype(float)
            mean = float(vals.mean())
            sd = float(vals.std(ddof=1)) if len(vals) > 1 else float("nan")
            effect = mean - overall_mean
            rel = (
                (effect / abs(overall_mean) * 100.0)
                if overall_mean != 0
                else float("nan")
            )
            row = {
                "Trait": trait,
                "Class": cls,
                "N": int(len(vals)),
                "Mean": mean,
                "SD": sd,
                "Overall_mean": overall_mean,
                "Effect": effect,
                "Relative_effect_percent": rel,
            }
            summary_rows.append(row)
            class_stats.append((cls, vals, row))

        if len(class_stats) >= 2:
            try:
                p_overall = float(kruskal(*[x[1] for x in class_stats]).pvalue)
            except ValueError:
                p_overall = 1.0 if tmp[trait].nunique() == 1 else float("nan")
        else:
            p_overall = float("nan")
        if len(class_stats) >= 2:
            test_rows.append(
                {
                    "Trait": trait,
                    "Comparison": "overall",
                    "Class1": "ALL",
                    "Class2": "",
                    "Test": "Kruskal-Wallis",
                    "N1": len(tmp),
                    "N2": 0,
                    "N_classes": len(class_stats),
                    "pvalue": p_overall,
                }
            )

        pair_tmp = []
        for i in range(len(class_stats)):
            for j in range(i + 1, len(class_stats)):
                c1, x, _ = class_stats[i]
                c2, y, _ = class_stats[j]
                p = float(
                    mannwhitneyu(x, y, alternative="two-sided", method="auto").pvalue
                )
                pair_tmp.append(
                    {
                        "Trait": trait,
                        "Comparison": "pairwise",
                        "Class1": c1,
                        "Class2": c2,
                        "Test": "Mann-Whitney U",
                        "N1": len(x),
                        "N2": len(y),
                        "N_classes": 2,
                        "pvalue": p,
                    }
                )
        if pair_tmp:
            qs = bh_adjust([r["pvalue"] for r in pair_tmp])
            for r, q in zip(pair_tmp, qs):
                r["padj_BH"] = q
            test_rows.extend(pair_tmp)

        # Superior haplotypes are descriptive candidates, not causal claims.
        for cls, _vals, row in class_stats:
            superior_rows.append(
                {
                    **row,
                    "Candidate": True,
                    "Candidate_min_N": min_class_n,
                    "Direction": "higher"
                    if row["Effect"] > 0
                    else ("lower" if row["Effect"] < 0 else "neutral"),
                    "Overall_Kruskal_p": p_overall,
                }
            )

    summary_path = prefix + ".TraitHaplotypeSummary.tsv"
    tests_path = prefix + ".TraitAssociationTests.tsv"
    superior_path = prefix + ".SuperiorHaplotypeCandidates.tsv"
    summary_columns = ["Trait", "Class", "N", "Mean", "SD", "Overall_mean",
                       "Effect", "Relative_effect_percent"]
    pd.DataFrame(summary_rows, columns=summary_columns).to_csv(summary_path, sep="\t", index=False)
    tests_df = pd.DataFrame(test_rows, columns=["Trait", "Comparison", "Class1", "Class2",
                                               "Test", "N1", "N2", "N_classes", "pvalue", "padj_BH"])
    overall_idx = tests_df["Comparison"].eq("overall")
    tests_df.loc[overall_idx, "padj_BH"] = bh_adjust(
        tests_df.loc[overall_idx, "pvalue"].tolist()
    )
    tests_df.to_csv(tests_path, sep="\t", index=False)
    sup = pd.DataFrame(superior_rows, columns=summary_columns + ["Candidate", "Candidate_min_N",
                                                               "Direction", "Overall_Kruskal_p"])
    sup["Abs_effect"] = sup["Effect"].abs()
    sup.sort_values(["Trait", "Abs_effect", "N"], ascending=[True, False, False]).drop(
        columns=["Abs_effect"]
    ).to_csv(superior_path, sep="\t", index=False)
    return summary_path, tests_path, superior_path


def _alt_dosage(call, sample, allele=None):
    gt = call.genotypes.get(sample, ())
    if not gt or any(a is None or a < 0 for a in gt):
        return float("nan")
    return sum((a > 0 if allele is None else a == allele) for a in gt) / len(gt)


def calculate_ld_r2(calls, samples, absence_mask=None, allele_specific=False):
    items = [
        (i, c, a)
        for i, c in enumerate(calls)
        for a in (range(1, len(c.alts) + 1) if allele_specific else [None])
    ]
    labels, seen = [], {}
    for i, c, a in items:
        label = f"{c.chrom}:{c.pos}" + (f":ALT{a}" if a is not None else "")
        seen[label] = seen.get(label, 0) + 1
        labels.append(label if seen[label] == 1 else f"{label}#{seen[label]}")
    n = len(items)
    mat = np.full((n, n), np.nan)
    pairs = []
    dosage = np.array(
        [
            [
                np.nan if (absence_mask or {}).get((s, i)) else _alt_dosage(c, s, a)
                for s in samples
            ]
            for i, c, a in items
        ]
    )
    for i in range(n):
        x = dosage[i]
        valid = np.isfinite(x)
        if valid.sum() >= 3 and np.std(x[valid]) > 0:
            mat[i, i] = 1.0
        for j in range(i + 1, n):
            y = dosage[j]
            mask = np.isfinite(x) & np.isfinite(y)
            count = int(mask.sum())
            r2 = np.nan
            if count >= 3 and np.std(x[mask]) > 0 and np.std(y[mask]) > 0:
                r2 = float(np.clip(np.corrcoef(x[mask], y[mask])[0, 1] ** 2, 0, 1))
            mat[i, j] = mat[j, i] = r2
            c1, c2 = items[i][1], items[j][1]
            pairs.append(
                {
                    "Variant1": labels[i],
                    "Variant2": labels[j],
                    "Distance_bp": abs(c1.pos - c2.pos)
                    if c1.chrom == c2.chrom
                    else np.nan,
                    "N_called": count,
                    "r2": r2,
                    "Dosage_definition": "allele_specific"
                    if allele_specific
                    else "any_ALT_fraction",
                }
            )
    return pd.DataFrame(mat, index=labels, columns=labels), pd.DataFrame(pairs)


def write_ld_outputs(prefix, calls, samples, absence_mask=None):
    if len(calls) < 2:
        return None
    matrix, pairs = calculate_ld_r2(calls, samples, absence_mask)
    matrix_path = prefix + ".LD_r2_matrix.tsv"
    matrix.to_csv(matrix_path, sep="\t")
    pairs.to_csv(prefix + ".LD_pairs.tsv", sep="\t", index=False)
    matrix, pairs = calculate_ld_r2(calls, samples, absence_mask, allele_specific=True)
    matrix.to_csv(prefix + ".LD_allele_r2_matrix.tsv", sep="\t")
    pairs.to_csv(prefix + ".LD_allele_pairs.tsv", sep="\t", index=False)
    return matrix_path
