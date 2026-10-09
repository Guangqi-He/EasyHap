from __future__ import annotations

from math import comb
from typing import Iterable, List, Sequence, Tuple


def fisher_exact_2x2(a, b, c, d):
    """Stable two-sided exact test, including extremely small tail probabilities."""
    from scipy.stats import fisher_exact

    return float(fisher_exact([[a, b], [c, d]], alternative="two-sided").pvalue)


def bh_adjust(p_values: Sequence[float]) -> List[float]:
    """Benjamini-Hochberg adjusted p values, robust to NaN and extreme values."""
    import math

    n = len(p_values)
    adj = [float("nan")] * n
    valid_indices = [
        i for i, p in enumerate(p_values) if p is not None and not math.isnan(float(p))
    ]
    if not valid_indices:
        return adj

    m = len(valid_indices)
    sorted_indices = sorted(valid_indices, key=lambda i: float(p_values[i]))

    running = 1.0
    for rank, idx in enumerate(reversed(sorted_indices), 1):
        denom_rank = m - rank + 1
        val = float(p_values[idx]) * m / max(denom_rank, 1)
        running = min(running, val)
        adj[idx] = min(1.0, max(0.0, running))

    return adj


def comparable_differences(a, b, missing="N", absence="ABS"):
    if len(a) != len(b):
        raise ValueError("State sequences must have equal lengths")
    compared = diff = 0
    for x, y in zip(a, b):
        parts = set(str(x).replace("|", "/").split("/")) | set(
            str(y).replace("|", "/").split("/")
        )
        if missing in parts or absence in parts:
            continue
        compared += 1
        diff += x != y
    return int(diff), compared


def hamming_distance(a, b, missing="N", absence="ABS"):
    diff, compared = comparable_differences(a, b, missing, absence)
    return diff / compared if compared else float("nan")


def connected_component_clusters(
    seqs: Sequence[Sequence[str]], threshold: float = 0.15
) -> List[str]:
    """Cluster haplotypes by pairwise Hamming distance using connected components."""
    n = len(seqs)
    parent = list(range(n))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: int, b: int) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    for i in range(n):
        for j in range(i + 1, n):
            if hamming_distance(seqs[i], seqs[j]) <= threshold:
                union(i, j)

    roots = {}
    labels: List[str] = []
    for i in range(n):
        r = find(i)
        if r not in roots:
            roots[r] = f"HapC{len(roots) + 1}"
        labels.append(roots[r])
    return labels
