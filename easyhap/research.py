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


def _sample_hap_records(sample_haps: Dict[str, List[str]], group_map: Dict[str, str], mode: str) -> pd.DataFrame:
    rows = []
    for sample, haps in sample_haps.items():
        group = group_map.get(sample, "All")
        if mode == "inbred":
            hap_list = haps[:1]
        else:
            hap_list = list(haps)
        for copy_i, hap in enumerate(hap_list, 1):
            rows.append({"Accession": sample, "Group": group, "Copy": copy_i, "Haplotype": hap})
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
        pd.DataFrame(columns=["Group", "Haplotype", "Count", "Frequency"]).to_csv(freq_path, sep="\t", index=False)
        pd.DataFrame(columns=["Group", "N_haplotype_copies", "N_haplotypes", "Haplotype_diversity_Hd"]).to_csv(div_path, sep="\t", index=False)
        pd.DataFrame(columns=["Group", "Haplotype", "Count", "Frequency"]).to_csv(private_path, sep="\t", index=False)
        return freq_path, div_path, private_path

    counts = rec.groupby(["Group", "Haplotype"]).size().rename("Count").reset_index()
    totals = counts.groupby("Group")["Count"].transform("sum")
    counts["Frequency"] = counts["Count"] / totals
    counts.sort_values(["Group", "Count", "Haplotype"], ascending=[True, False, True]).to_csv(freq_path, sep="\t", index=False)

    div_rows = []
    for group, sub in counts.groupby("Group"):
        n = int(sub["Count"].sum())
        p = sub["Count"].astype(float) / max(n, 1)
        hd = float("nan") if n <= 1 else (n / (n - 1.0)) * (1.0 - float(np.sum(p * p)))
        div_rows.append({
            "Group": group,
            "N_haplotype_copies": n,
            "N_haplotypes": int((sub["Count"] > 0).sum()),
            "Haplotype_diversity_Hd": hd,
        })
    pd.DataFrame(div_rows).to_csv(div_path, sep="\t", index=False)

    presence = counts.pivot_table(index="Haplotype", columns="Group", values="Count", fill_value=0, aggfunc="sum")
    private_rows = []
    # Private haplotypes are only meaningful when two or more biological groups are supplied.
    if len(presence.columns) >= 2:
        for hap, row in presence.iterrows():
            groups = [g for g, v in row.items() if v > 0]
            if len(groups) == 1:
                g = groups[0]
                sub = counts[(counts["Group"] == g) & (counts["Haplotype"] == hap)].iloc[0]
                private_rows.append({"Group": g, "Haplotype": hap, "Count": int(sub["Count"]), "Frequency": float(sub["Frequency"])})
    pd.DataFrame(private_rows, columns=["Group", "Haplotype", "Count", "Frequency"]).to_csv(private_path, sep="\t", index=False)
    return freq_path, div_path, private_path


def write_trait_association(
    prefix: str,
    hap_group_df: pd.DataFrame,
    trait_cols: Optional[Sequence[str]] = None,
    class_col: str = "Hap",
    min_class_n: int = 2,
) -> Tuple[Optional[str], Optional[str], Optional[str]]:
    if hap_group_df.empty or class_col not in hap_group_df.columns:
        return None, None, None
    available = [c for c in hap_group_df.columns if c not in {"Hap", "ClusterID", "Accession", "Type"}]
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
        tmp = tmp.dropna(subset=[trait])
        if tmp.empty:
            continue
        overall_mean = float(tmp[trait].mean())
        class_stats = []
        all_stats = []
        for cls, sub in tmp.groupby(class_col):
            vals = sub[trait].astype(float)
            mean = float(vals.mean())
            sd = float(vals.std(ddof=1)) if len(vals) > 1 else 0.0
            effect = mean - overall_mean
            rel = (effect / abs(overall_mean) * 100.0) if overall_mean != 0 else float("nan")
            row = {"Trait": trait, "Class": cls, "N": int(len(vals)), "Mean": mean, "SD": sd,
                   "Overall_mean": overall_mean, "Effect": effect, "Relative_effect_percent": rel}
            summary_rows.append(row)
            all_stats.append((cls, vals, row))
            if len(vals) >= min_class_n:
                class_stats.append((cls, vals, row))

        if len(class_stats) >= 2:
            p_overall = float(kruskal(*[x[1] for x in class_stats]).pvalue)
        else:
            p_overall = float("nan")
        test_rows.append({"Trait": trait, "Comparison": "overall", "Class1": "ALL", "Class2": "", "Test": "Kruskal-Wallis", "pvalue": p_overall})

        pair_tmp = []
        for i in range(len(class_stats)):
            for j in range(i + 1, len(class_stats)):
                c1, x, _ = class_stats[i]
                c2, y, _ = class_stats[j]
                p = float(mannwhitneyu(x, y, alternative="two-sided", method="auto").pvalue)
                pair_tmp.append({"Trait": trait, "Comparison": "pairwise", "Class1": c1, "Class2": c2,
                                 "Test": "Mann-Whitney U", "pvalue": p})
        if pair_tmp:
            qs = bh_adjust([r["pvalue"] for r in pair_tmp])
            for r, q in zip(pair_tmp, qs):
                r["padj_BH"] = q
            test_rows.extend(pair_tmp)

        # Superior haplotypes are descriptive candidates, not causal claims.
        for cls, _vals, row in all_stats:
            superior_rows.append({**row, "Direction": "higher" if row["Effect"] > 0 else ("lower" if row["Effect"] < 0 else "neutral"),
                                    "Overall_Kruskal_p": p_overall})

    if not summary_rows:
        return None, None, None

    summary_path = prefix + ".TraitHaplotypeSummary.tsv"
    tests_path = prefix + ".TraitAssociationTests.tsv"
    superior_path = prefix + ".SuperiorHaplotypeCandidates.tsv"
    pd.DataFrame(summary_rows).to_csv(summary_path, sep="\t", index=False)
    tests_df = pd.DataFrame(test_rows)
    if "padj_BH" not in tests_df.columns:
        tests_df["padj_BH"] = np.nan
    tests_df.to_csv(tests_path, sep="\t", index=False)
    sup = pd.DataFrame(superior_rows)
    sup["Abs_effect"] = sup["Effect"].abs()
    sup.sort_values(["Trait", "Abs_effect", "N"], ascending=[True, False, False]).drop(columns=["Abs_effect"]).to_csv(superior_path, sep="\t", index=False)
    return summary_path, tests_path, superior_path


def _alt_dosage(call: VariantCall, sample: str) -> float:
    gt = call.genotypes.get(sample, tuple())
    called = [a for a in gt if a is not None and a >= 0]
    if not called:
        return float("nan")
    return float(sum(1 for a in called if a > 0)) / float(len(called))


def calculate_ld_r2(calls: Sequence[VariantCall], samples: Sequence[str]) -> Tuple[pd.DataFrame, pd.DataFrame]:
    labels = [f"{c.chrom}:{c.pos}" for c in calls]
    n = len(calls)
    mat = np.full((n, n), np.nan, dtype=float)
    pairs = []
    dosage = np.array([[_alt_dosage(c, s) for s in samples] for c in calls], dtype=float)
    for i in range(n):
        mat[i, i] = 1.0
        for j in range(i + 1, n):
            x, y = dosage[i], dosage[j]
            mask = np.isfinite(x) & np.isfinite(y)
            n_called = int(mask.sum())
            r2 = float("nan")
            if n_called >= 3 and np.nanstd(x[mask]) > 0 and np.nanstd(y[mask]) > 0:
                r = float(np.corrcoef(x[mask], y[mask])[0, 1])
                r2 = max(0.0, min(1.0, r * r)) if np.isfinite(r) else float("nan")
            mat[i, j] = mat[j, i] = r2
            pairs.append({"Variant1": labels[i], "Variant2": labels[j], "Distance_bp": abs(calls[j].pos - calls[i].pos), "N_called": n_called, "r2": r2})
    return pd.DataFrame(mat, index=labels, columns=labels), pd.DataFrame(pairs)


def write_ld_outputs(prefix: str, calls: Sequence[VariantCall], samples: Sequence[str]) -> Optional[str]:
    if len(calls) < 2:
        return None
    matrix, _ = calculate_ld_r2(calls, samples)
    matrix_path = prefix + ".LD_r2_matrix.tsv"
    matrix.to_csv(matrix_path, sep="\t")
    return matrix_path
