from __future__ import annotations

from math import comb
from typing import Iterable, List, Sequence, Tuple


def fisher_exact_2x2(a: int, b: int, c: int, d: int) -> float:
    """Two-sided Fisher exact test for [[a,b],[c,d]].

    Uses the conventional probability-summing definition. This avoids requiring
    SciPy for the filtering step.
    """
    n = a + b + c + d
    row1 = a + b
    row2 = c + d
    col1 = a + c

    def prob(x: int) -> float:
        return (comb(row1, x) * comb(row2, col1 - x)) / comb(n, col1)

    lo = max(0, col1 - row2)
    hi = min(row1, col1)
    p_obs = prob(a)
    p = 0.0
    eps = 1e-12
    for x in range(lo, hi + 1):
        px = prob(x)
        if px <= p_obs + eps:
            p += px
    return min(1.0, p)


def bh_adjust(p_values: Sequence[float]) -> List[float]:
    """Benjamini-Hochberg adjusted p values."""
    n = len(p_values)
    order = sorted(range(n), key=lambda i: p_values[i])
    adj = [1.0] * n
    running = 1.0
    for rank, idx in enumerate(reversed(order), 1):
        # rank from largest to smallest; BH denominator is n-rank+1
        denom_rank = n - rank + 1
        val = p_values[idx] * n / max(denom_rank, 1)
        running = min(running, val)
        adj[idx] = min(1.0, running)
    return adj


def hamming_distance(a: Sequence[str], b: Sequence[str], missing: str = "N") -> float:
    compared = 0
    diff = 0
    for x, y in zip(a, b):
        if x == missing or y == missing:
            continue
        compared += 1
        if x != y:
            diff += 1
    if compared == 0:
        return 1.0
    return diff / compared


def connected_component_clusters(seqs: Sequence[Sequence[str]], threshold: float = 0.15) -> List[str]:
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
