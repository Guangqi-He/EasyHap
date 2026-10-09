"""Optional, explicitly exploratory multilocus analyses with no covariates."""

from __future__ import annotations

import os
import numpy as np
import pandas as pd
from scipy import stats
from .stats import bh_adjust


def link_known_phase(merged, region_paths, outdir):
    """Connect existing copy order only when chromosome, input and PS agree.

    This never estimates phase from unphased genotypes. Homozygous regions do
    not constrain orientation. Explicit and implicit PS are not mixed.
    """
    tables = []
    for path in region_paths:
        sidecar = path[: -len(".HapGroup.tsv")] + ".HaplotypeCopies.tsv"
        if not os.path.isfile(sidecar):
            tables.append(None)
            continue
        df = pd.read_csv(sidecar, sep="\t", dtype=str, keep_default_na=False)
        required = {
            "Accession",
            "Copy",
            "Hap",
            "Chrom",
            "PhaseSet",
            "InputSHA256",
            "Mode",
        }
        if not required.issubset(df.columns):
            raise ValueError(f"Invalid copy metadata: {sidecar}")
        if df.duplicated(["Accession", "Copy"]).any():
            raise ValueError(f"Duplicate sample/copy in {sidecar}")
        tables.append(df)
    rows, copy_rows = [], []
    for sample in merged["Accession"]:
        sets = [
            None if df is None else df[df.Accession == sample].copy() for df in tables
        ]
        status, reason, linked = "UNRESOLVED", "", ""
        if any(df is None or df.empty for df in sets):
            reason = "missing_copy_metadata_or_genotype_mode"
        else:
            for df in sets:
                df["Copy"] = pd.to_numeric(df["Copy"], errors="raise").astype(int)
                df.sort_values("Copy", inplace=True)
            sizes = {len(df) for df in sets}
            hashes = {h for df in sets for h in df.InputSHA256}
            chroms = {c for df in sets for c in df.Chrom}
            anchors = {p for df in sets for p in df.PhaseSet if p != "UNINFORMATIVE"}
            if any(set(df.Mode) != {"copy"} for df in sets):
                reason = "not_copy_mode"
            elif len(sizes) != 1 or any(
                df.Copy.tolist() != list(range(1, len(df) + 1)) for df in sets
            ):
                reason = "inconsistent_copy_counts"
            elif len(chroms) != 1:
                reason = "different_chromosomes"
            elif len(hashes) != 1 or "" in hashes:
                reason = "different_or_unknown_inputs"
            elif len(anchors) > 1 or anchors & {"MULTIPLE", "", "UNRESOLVED"}:
                reason = "different_or_unknown_phase_sets"
            else:
                # Check that the sidecar corresponds to this HapGroup classification.
                sample_row = merged.loc[merged.Accession == sample].iloc[0]
                hap_cols = [
                    c for c in merged if c.startswith("Hap_") and c != "Hap_Combination"
                ]
                if any(
                    "_".join(sorted(df.Hap)) != sample_row[col]
                    for df, col in zip(sets, hap_cols)
                ):
                    reason = "copy_metadata_does_not_match_hapgroup"
                else:
                    copies = [
                        "+".join(df.iloc[i].Hap for df in sets)
                        for i in range(len(sets[0]))
                    ]
                    linked = "|".join(sorted(copies))
                    status, reason = (
                        "RESOLVED",
                        "all_regions_uninformative"
                        if not anchors
                        else "known_phase_block",
                    )
                    for i, combination in enumerate(copies, 1):
                        copy_rows.append(
                            {
                                "Accession": sample,
                                "Copy": i,
                                "Phased_Combination": combination,
                            }
                        )
        rows.append(
            {
                "Accession": sample,
                "Status": status,
                "Reason": reason,
                "Phased_Diplotype": linked,
            }
        )
    path = os.path.join(outdir, "Aggregate_KnownPhase.tsv")
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
    pd.DataFrame(copy_rows, columns=["Accession", "Copy", "Phased_Combination"]).to_csv(
        os.path.join(outdir, "Aggregate_PhasedCopies.tsv"), sep="\t", index=False
    )
    return path


def fit_categorical_interaction(df, factor1, factor2, traits, min_cell_n=2):
    """Nested OLS F test: y ~ A+B versus y ~ A+B+A:B.

    All observed factor levels must form a complete factorial table with
    replication; no partial interaction is presented as a full interaction.
    Contrasts use treatment coding and are differences of differences.
    """
    results, contrasts = [], []
    for trait in traits:
        tmp = df[[factor1, factor2, trait]].copy()
        tmp[trait] = pd.to_numeric(tmp[trait], errors="coerce")
        tmp = tmp[np.isfinite(tmp[trait])]
        a_levels, b_levels = (
            sorted(tmp[factor1].unique()),
            sorted(tmp[factor2].unique()),
        )
        counts = tmp.groupby([factor1, factor2]).size()
        row = {
            "Trait": trait,
            "Model": "categorical_OLS_A+B_vs_A*B",
            "N": len(tmp),
            "Factor1": factor1,
            "Factor2": factor2,
            "Levels1": len(a_levels),
            "Levels2": len(b_levels),
            "Cells_observed": len(counts),
            "Cells_required": len(a_levels) * len(b_levels),
            "Min_cell_n": int(counts.min()) if len(counts) else 0,
            "Status": "NOT_ESTIMABLE",
            "Reason": "",
            "F": np.nan,
            "pvalue": np.nan,
            "DF_interaction": np.nan,
            "DF_residual": np.nan,
            "R2_additive": np.nan,
            "R2_full": np.nan,
            "Delta_R2": np.nan,
            "Partial_eta2": np.nan,
        }
        if len(a_levels) < 2 or len(b_levels) < 2:
            row["Reason"] = "each_factor_requires_two_levels"
        elif len(counts) != row["Cells_required"]:
            row["Reason"] = "missing_factorial_cells"
        elif counts.min() < min_cell_n:
            row["Reason"] = "insufficient_cell_replication"
        else:
            a = np.column_stack(
                [
                    (tmp[factor1] == level).to_numpy(dtype=float)
                    for level in a_levels[1:]
                ]
            )
            b = np.column_stack(
                [
                    (tmp[factor2] == level).to_numpy(dtype=float)
                    for level in b_levels[1:]
                ]
            )
            terms = np.column_stack(
                [
                    a[:, i] * b[:, j]
                    for i in range(a.shape[1])
                    for j in range(b.shape[1])
                ]
            )
            reduced = np.column_stack([np.ones(len(tmp)), a, b])
            full = np.column_stack([reduced, terms])
            y = tmp[trait].to_numpy(dtype=float)
            beta0, _, rank0, _ = np.linalg.lstsq(reduced, y, rcond=None)
            beta, _, rank, _ = np.linalg.lstsq(full, y, rcond=None)
            df_num, df_res = int(rank - rank0), int(len(y) - rank)
            row["DF_interaction"], row["DF_residual"] = df_num, df_res
            rss0 = float(np.sum((y - reduced @ beta0) ** 2))
            rss = float(np.sum((y - full @ beta) ** 2))
            tss = float(np.sum((y - y.mean()) ** 2))
            if rank < full.shape[1] or df_num < 1 or df_res < 1:
                row["Reason"] = "rank_deficient_or_no_residual_degrees_of_freedom"
            elif tss <= 0 or rss <= np.finfo(float).eps * tss * 10:
                row["Reason"] = "zero_response_or_residual_variance"
            else:
                change = max(0.0, rss0 - rss)
                f = (change / df_num) / (rss / df_res)
                row.update(
                    Status="OK",
                    F=f,
                    pvalue=float(stats.f.sf(f, df_num, df_res)),
                    R2_additive=1 - rss0 / tss,
                    R2_full=1 - rss / tss,
                    Delta_R2=change / tss,
                    Partial_eta2=change / rss0 if rss0 > 0 else np.nan,
                )
                covariance = np.linalg.inv(full.T @ full) * (rss / df_res)
                for k, (i, j) in enumerate(
                    (i, j) for i in range(a.shape[1]) for j in range(b.shape[1])
                ):
                    index = reduced.shape[1] + k
                    effect, se = (
                        float(beta[index]),
                        float(np.sqrt(covariance[index, index])),
                    )
                    q = float(stats.t.ppf(0.975, df_res))
                    contrasts.append(
                        {
                            "Trait": trait,
                            "Factor1_level": a_levels[i + 1],
                            "Factor1_reference": a_levels[0],
                            "Factor2_level": b_levels[j + 1],
                            "Factor2_reference": b_levels[0],
                            "Interaction_effect": effect,
                            "SE": se,
                            "CI95_lower": effect - q * se,
                            "CI95_upper": effect + q * se,
                            "pvalue": float(2 * stats.t.sf(abs(effect / se), df_res))
                            if se > 0
                            else np.nan,
                        }
                    )
        results.append(row)
    res = pd.DataFrame(
        results,
        columns=list(row) if results else ["Trait", "Status", "Reason", "pvalue"],
    )
    con = pd.DataFrame(
        contrasts,
        columns=[
            "Trait",
            "Factor1_level",
            "Factor1_reference",
            "Factor2_level",
            "Factor2_reference",
            "Interaction_effect",
            "SE",
            "CI95_lower",
            "CI95_upper",
            "pvalue",
        ],
    )
    return res, con
