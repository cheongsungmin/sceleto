from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List, Literal, Optional, Tuple

import numpy as np
import pandas as pd

SigmaMethod = Literal["sd", "iqr", "mad"]


@dataclass(frozen=True)
class MarkerLabels:
    """Label matrices for genes across groups.

    Notes
    -----
    - Matrices are (n_groups, n_genes_used) aligned to ctx.groups and `genes`.
    - Values are -1/0/+1 for Low/Undetermined/High.
    """
    genes: List[str]
    level1: np.ndarray
    level2: np.ndarray
    level3: np.ndarray
    params: Dict[str, Any]


def robust_sigma(x: np.ndarray, method: SigmaMethod = "sd") -> float:
    """Robust scale estimator."""
    x = np.asarray(x, dtype=float)
    if x.size < 2:
        return np.nan

    if method == "sd":
        return float(x.std(ddof=1))
    if method == "iqr":
        q1, q3 = np.percentile(x, [25, 75])
        return float((q3 - q1) / 1.349)
    if method == "mad":
        med = np.median(x)
        mad = np.median(np.abs(x - med))
        return float(1.4826 * mad)

    raise ValueError("method must be 'sd', 'iqr', or 'mad'")


def _level1_seed_from_edges(
    edge_gene_df: pd.DataFrame,
    *,
    groups: List[str],
    group_to_idx: Dict[str, int],
    genes: Optional[List[str]] = None,
    fc_cutoff: float = 3.0,
) -> Tuple[np.ndarray, List[str]]:
    required = {"start", "end", "gene", "fc"}
    missing = required - set(edge_gene_df.columns)
    if missing:
        raise ValueError(f"edge_gene_df missing columns: {sorted(missing)}")

    df = edge_gene_df.loc[edge_gene_df["fc"] >= fc_cutoff, ["start", "end", "gene"]].copy()
    df["start"] = df["start"].astype(str)
    df["end"] = df["end"].astype(str)
    df["gene"] = df["gene"].astype(str)

    if genes is None:
        genes_used = sorted(df["gene"].unique().tolist())
    else:
        genes_used = list(map(str, genes))
        df = df[df["gene"].isin(set(genes_used))]

    C = len(groups)
    P = len(genes_used)
    labels1 = np.zeros((C, P), dtype=np.int8)

    for j, gene in enumerate(genes_used):
        sub = df[df["gene"] == gene]
        if sub.empty:
            continue

        low_set = set(sub["start"].unique().tolist())
        high_set = set(sub["end"].unique().tolist())
        inter = low_set & high_set

        for g in low_set:
            i = group_to_idx.get(g)
            if i is not None:
                labels1[i, j] = -1
        for g in high_set:
            i = group_to_idx.get(g)
            if i is not None:
                labels1[i, j] = 1
        for g in inter:
            i = group_to_idx.get(g)
            if i is not None:
                labels1[i, j] = 0

    return labels1, genes_used


def _level2_push_by_mean(mean: np.ndarray, labels1: np.ndarray) -> np.ndarray:
    if mean.shape != labels1.shape:
        raise ValueError("mean and labels1 must have the same shape.")

    _, P = mean.shape
    labels2 = np.zeros_like(labels1, dtype=np.int8)

    for j in range(P):
        low_mask = labels1[:, j] == -1
        high_mask = labels1[:, j] == 1

        if low_mask.any():
            max_low = float(np.max(mean[low_mask, j]))
            labels2[mean[:, j] <= max_low, j] = -1

        if high_mask.any():
            min_high = float(np.min(mean[high_mask, j]))
            labels2[mean[:, j] >= min_high, j] = 1

    return labels2


def _level3_expand_undetermined(
    mean: np.ndarray,
    labels2: np.ndarray,
    *,
    k: float = 2.0,
    sigma_method: SigmaMethod = "sd",
    min_gap: float = 0.2,      # legacy criterion (mu_high - mu_low)
    min_margin: float = 0.0,   # set 0.0 to match legacy behavior
) -> np.ndarray:
    labels3 = labels2.copy().astype(np.int8)
    _, P = mean.shape

    for j in range(P):
        low_vals = mean[labels2[:, j] == -1, j]
        high_vals = mean[labels2[:, j] == 1, j]
        if low_vals.size == 0 or high_vals.size == 0:
            continue

        mu_low = float(np.mean(low_vals))
        mu_high = float(np.mean(high_vals))
        gap = mu_high - mu_low

        sig_low = robust_sigma(low_vals, method=sigma_method)
        sig_high = robust_sigma(high_vals, method=sigma_method)
        if np.isnan(sig_low) or np.isnan(sig_high):
            continue

        low_upper = mu_low + k * sig_low
        high_lower = mu_high - k * sig_high
        margin = high_lower - low_upper

        # match old logic: avoid overlap OR too small gap
        if (low_upper >= high_lower) or (gap < min_gap):
            continue
        # extra safety if you still want margin constraint
        if margin < min_margin:
            continue

        und = labels3[:, j] == 0
        labels3[und & (mean[:, j] <= low_upper), j] = -1
        labels3[und & (mean[:, j] >= high_lower), j] = 1

    return labels3


def label_levels(
    ctx,
    edge_gene_df: pd.DataFrame,
    *,
    fc_cutoff: float = 3.0,
    genes: Optional[List[str]] = None,
    k: float = 2.0,
    sigma_method: SigmaMethod = "sd",
    min_margin: float = 0.2,
) -> MarkerLabels:
    """Run level1/2/3 labeling pipeline.

    Parameters
    ----------
    ctx
        MarkerContext from build_context().
    edge_gene_df
        DataFrame from compute_fc_delta(ctx, ...), must include start/end/gene/fc.
    genes
        If provided, restrict to these genes (still must exist in ctx.genes).
    """
    labels1, genes_used = _level1_seed_from_edges(
        edge_gene_df,
        groups=ctx.groups,
        group_to_idx=ctx.group_to_idx,
        genes=genes,
        fc_cutoff=fc_cutoff,
    )

    if len(genes_used) == 0:
        raise ValueError("No genes selected for labeling. Check fc_cutoff or input edge_gene_df.")

    # Subset ctx.mean_norm to genes_used (and align labels1 columns)
    gene_to_idx = {g: i for i, g in enumerate(ctx.genes.astype(str))}

    kept_cols = [c for c, g in enumerate(genes_used) if str(g) in gene_to_idx]
    if len(kept_cols) == 0:
        raise ValueError("None of the selected genes were found in ctx.genes.")

    labels1 = labels1[:, kept_cols]
    genes_used = [str(genes_used[c]) for c in kept_cols]
    idx = [gene_to_idx[g] for g in genes_used]

    mean_sub = ctx.mean_norm[:, idx]

    labels2 = _level2_push_by_mean(mean_sub, labels1)
    labels3 = _level3_expand_undetermined(
        mean_sub,
        labels2,
        k=k,
        sigma_method=sigma_method,
        min_margin=min_margin,
    )

    return MarkerLabels(
        genes=list(genes_used),
        level1=labels1,
        level2=labels2,
        level3=labels3,
        params={
            "fc_cutoff": fc_cutoff,
            "k": k,
            "sigma_method": sigma_method,
            "min_margin": min_margin,
        },
    )


def labels_to_note_df(
    ctx,
    labels: MarkerLabels,
    *,
    level: Literal[1, 2, 3] = 3,
    dtype=np.int8,
) -> pd.DataFrame:
    """Convert MarkerLabels to a note_df (groups×genes) DataFrame.

    Parameters
    ----------
    ctx
        MarkerContext (used for group order).
    labels
        Output of label_levels().
    level
        Which level to use: 1, 2, or 3.
    dtype
        Output dtype, default int8.

    Returns
    -------
    note_df
        DataFrame indexed by ctx.groups and columns labels.genes with values -1/0/+1.
    """
    if level == 1:
        mat = labels.level1
    elif level == 2:
        mat = labels.level2
    elif level == 3:
        mat = labels.level3
    else:
        raise ValueError("level must be one of {1,2,3}.")

    if mat.shape[0] != len(ctx.groups):
        raise ValueError(
            f"Group dimension mismatch: labels has {mat.shape[0]} rows, ctx.groups has {len(ctx.groups)}"
        )

    return pd.DataFrame(mat.astype(dtype, copy=False), index=ctx.groups, columns=labels.genes)
