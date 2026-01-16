from __future__ import annotations

from typing import Callable, Dict, Iterable, List, Optional, Literal

import numpy as np
import pandas as pd

from ._features import compute_gene_features

SpecificPreset = Literal[1, 2, 3, 4]
SpecificScoreFn = Callable[[pd.DataFrame], np.ndarray]


def compute_dst_gene_max_fc_delta(
    ctx,
    *,
    genes: Optional[Iterable[str]] = None,
    eps: float = 1e-3,
    drop_zeros: bool = False,
) -> pd.DataFrame:
    """Compute per-(dst,gene) maximum FC and delta over all PAGA edges (no filtering).

    Notes
    -----
    - Direction is defined per gene: low -> high using ctx.mean_norm.
    - Aggregation: for each dst and gene, take max over all incoming edges.
    - No biological thresholds/filters are applied (by design).
    """
    if ctx.undirected_edges is None:
        raise ValueError("ctx.undirected_edges is None. Build context with PAGA first.")

    all_genes = np.asarray(ctx.genes, dtype=str)
    if genes is None:
        gene_idx = np.arange(all_genes.size, dtype=int)
        use_genes = all_genes
    else:
        gene_set = set(map(str, genes))
        gene_idx = np.array([i for i, g in enumerate(all_genes) if g in gene_set], dtype=int)
        use_genes = all_genes[gene_idx]

    groups = list(map(str, ctx.groups))
    C = len(groups)
    P = len(gene_idx)

    max_fc = np.zeros((C, P), dtype=np.float32)
    max_delta = np.zeros((C, P), dtype=np.float32)

    for a, b in ctx.undirected_edges:
        ia = ctx.get_group_idx(a)
        ib = ctx.get_group_idx(b)

        mA = ctx.mean_norm[ia][gene_idx]
        mB = ctx.mean_norm[ib][gene_idx]

        high_is_A = mA >= mB
        m_high = np.where(high_is_A, mA, mB)
        m_low = np.where(high_is_A, mB, mA)

        fc = (m_high + eps) / (m_low + eps)
        delta = m_high - m_low

        idxA = high_is_A
        if np.any(idxA):
            max_fc[ia, idxA] = np.maximum(max_fc[ia, idxA], fc[idxA].astype(np.float32))
            max_delta[ia, idxA] = np.maximum(
                max_delta[ia, idxA], delta[idxA].astype(np.float32)
            )

        idxB = ~high_is_A
        if np.any(idxB):
            max_fc[ib, idxB] = np.maximum(max_fc[ib, idxB], fc[idxB].astype(np.float32))
            max_delta[ib, idxB] = np.maximum(
                max_delta[ib, idxB], delta[idxB].astype(np.float32)
            )

    wide_fc = pd.DataFrame(max_fc, index=groups, columns=use_genes)
    wide_delta = pd.DataFrame(max_delta, index=groups, columns=use_genes)

    out = (
        wide_fc.stack().rename("max_fc").to_frame()
        .join(wide_delta.stack().rename("max_delta"))
        .reset_index()
        .rename(columns={"level_0": "group", "level_1": "gene"})
    )

    if drop_zeros:
        out = out[(out["max_fc"] > 0) | (out["max_delta"] > 0)].reset_index(drop=True)

    return out


def compute_coverage_mats(
    ctx,
    *,
    genes: Optional[Iterable[str]] = None,
    use_weighted_mean: bool = True,
    fallback_to_rest_if_isolated: bool = True,
) -> Dict[str, pd.DataFrame]:
    """Compute coverage matrices for selected genes.

    Returns
    -------
    Dict[str, pd.DataFrame]
        DataFrames (clusters x genes):
        - "coverage": frac_expr within cluster
        - "coverage_rest": weighted mean across all other clusters
        - "coverage_neighbor": weighted mean across neighbor clusters (PAGA undirected)
    """
    groups = [str(g) for g in ctx.groups]
    group_to_idx = {g: i for i, g in enumerate(groups)}

    all_genes = [str(g) for g in list(ctx.genes)]
    gene_to_idx = {g: j for j, g in enumerate(all_genes)}

    if genes is None:
        genes_used = all_genes
        sel_idx = np.arange(len(all_genes), dtype=int)
    else:
        genes_list = [str(g) for g in genes]
        sel_idx = np.array([gene_to_idx[g] for g in genes_list if g in gene_to_idx], dtype=int)
        genes_used = [g for g in genes_list if g in gene_to_idx]

    frac = np.asarray(ctx.frac_expr)[:, sel_idx]  # (C,P)
    cov_df = pd.DataFrame(frac, index=groups, columns=genes_used)

    n_cells = np.asarray(ctx.n_cells, dtype=float).reshape(-1)  # (C,)

    # coverage_rest (vectorized)
    if use_weighted_mean:
        w = n_cells
        sum_w = float(w.sum())
        ws_all = (w[:, None] * frac).sum(axis=0)  # (P,)
        denom = (sum_w - w)[:, None]              # (C,1)
        numer = (ws_all[None, :] - (w[:, None] * frac))
        cov_rest = np.divide(numer, denom, out=np.zeros_like(numer), where=denom > 0)
    else:
        C = frac.shape[0]
        cov_rest = (frac.sum(axis=0, keepdims=True) - frac) / max(C - 1, 1)

    cov_rest_df = pd.DataFrame(cov_rest, index=groups, columns=genes_used)

    # coverage_neighbor (vectorized via adjacency)
    edges = getattr(ctx, "undirected_edges", None) or []
    C = len(groups)
    A = np.zeros((C, C), dtype=np.float32)
    for a, b in edges:
        ia = group_to_idx.get(str(a))
        ib = group_to_idx.get(str(b))
        if ia is None or ib is None:
            continue
        A[ia, ib] = 1.0
        A[ib, ia] = 1.0
    np.fill_diagonal(A, 0.0)

    if use_weighted_mean:
        Aw = A * n_cells[None, :]                 # (C,C)
        numer = Aw @ frac                         # (C,P)
        denom = Aw.sum(axis=1, keepdims=True)     # (C,1)
        cov_nbr = np.divide(numer, denom, out=np.zeros_like(numer), where=denom > 0)
    else:
        numer = A @ frac
        denom = A.sum(axis=1, keepdims=True)
        cov_nbr = np.divide(numer, denom, out=np.zeros_like(numer), where=denom > 0)

    if fallback_to_rest_if_isolated:
        iso = (A.sum(axis=1) == 0)
        if np.any(iso):
            cov_nbr[iso, :] = cov_rest[iso, :]

    cov_nbr_df = pd.DataFrame(cov_nbr, index=groups, columns=genes_used)

    return {
        "coverage": cov_df,
        "coverage_rest": cov_rest_df,
        "coverage_neighbor": cov_nbr_df,
    }


def local_score(df: pd.DataFrame, *, eps: float = 1e-9) -> np.ndarray:
    cov_in = df["coverage_one"].to_numpy(dtype=float)
    cov_out = (
        0.5 * df["coverage_rest"].to_numpy(dtype=float)
        + 0.5 * df["coverage_neighbor"].to_numpy(dtype=float)
    )

    # Specificity: in vs out (dominant term for local markers)
    spec = np.log((cov_in + eps) / (cov_out + eps))
    spec = np.maximum(0.0, spec)

    # Saturating reward for in-cluster expression fraction
    cov_term = np.sqrt(cov_in)

    # Strength (stabilize heavy tails)
    strength = np.log1p(df["max_fc"].to_numpy(dtype=float)) * np.log1p(
        df["max_delta"].to_numpy(dtype=float)
    )

    return strength * cov_term * spec


def global_score(df: pd.DataFrame, *, eps: float = 1e-9) -> np.ndarray:
    gap_term = 1.0 + df["gap"].to_numpy(dtype=float)
    low_term = np.log1p(df["n_low"].to_numpy(dtype=float))
    high_pen = 1.0 + df["n_high"].to_numpy(dtype=float)
    grey_pen = 1.0 + 0.5 * df["n_grey"].to_numpy(dtype=float)

    return (gap_term * (1.0 + low_term)) / (high_pen * grey_pen + eps)


def weight_local_prioritized(
    df: pd.DataFrame,
    *,
    A: float = 1.0,
    B: float = 0.5,
    eps: float = 1e-9,
) -> np.ndarray:
    L = local_score(df, eps=eps)
    G = global_score(df, eps=eps)
    return np.power(L, A) * np.power(G, B)

# ================================
# ========== Experiment ==========
# ================================
def local_score_soft(df: pd.DataFrame, *, eps: float = 1e-9, tau: float = 0.25) -> np.ndarray:
    # Soft specificity: softplus(log-ratio) for smooth behavior near zero
    cov_in = df["coverage_one"].to_numpy(dtype=float)
    cov_out = (
        0.5 * df["coverage_rest"].to_numpy(dtype=float)
        + 0.5 * df["coverage_neighbor"].to_numpy(dtype=float)
    )
    x = np.log((cov_in + eps) / (cov_out + eps))
    spec = np.log1p(np.exp(x / tau)) * tau

    cov_term = np.sqrt(cov_in)
    strength = np.log1p(df["max_fc"].to_numpy(dtype=float)) * np.log1p(df["max_delta"].to_numpy(dtype=float))
    return strength * cov_term * spec


def local_score_worst_out(df: pd.DataFrame, *, eps: float = 1e-9) -> np.ndarray:
    # Conservative specificity: compare to worst(out) = max(rest, neighbor)
    cov_in = df["coverage_one"].to_numpy(dtype=float)
    cov_out = np.maximum(
        df["coverage_rest"].to_numpy(dtype=float),
        df["coverage_neighbor"].to_numpy(dtype=float),
    )
    spec = np.log((cov_in + eps) / (cov_out + eps))
    spec = np.maximum(0.0, spec)

    cov_term = np.sqrt(cov_in)
    strength = np.log1p(df["max_fc"].to_numpy(dtype=float)) * np.log1p(df["max_delta"].to_numpy(dtype=float))
    return strength * cov_term * spec


def global_score_simple(df: pd.DataFrame, *, eps: float = 1e-9) -> np.ndarray:
    # Reward decisiveness (low grey), penalize ubiquity (high)
    decisiveness = 1.0 / (1.0 + df["n_grey"].to_numpy(dtype=float))
    ubiq_pen = 1.0 / (1.0 + df["n_high"].to_numpy(dtype=float))
    return (1.0 + df["gap"].to_numpy(dtype=float)) * decisiveness * ubiq_pen + eps


def make_specific_score_fn(
    preset: SpecificPreset,
    *,
    A: float = 1.0,
    B: float = 0.5,
    eps: float = 1e-9,
    tau: float = 0.25,
) -> SpecificScoreFn:
    """Return a score_fn(df)->np.ndarray for specific marker ranking."""

    if preset == 1:
        def _fn(df: pd.DataFrame) -> np.ndarray:
            L = local_score_soft(df, eps=eps, tau=tau)
            G = global_score(df, eps=eps)
            return np.power(L, A) * np.power(G, B)
        return _fn

    if preset == 2:
        def _fn(df: pd.DataFrame) -> np.ndarray:
            L = local_score_worst_out(df, eps=eps)
            G = global_score(df, eps=eps)
            return np.power(L, A) * np.power(G, B)
        return _fn

    if preset == 3:
        def _fn(df: pd.DataFrame) -> np.ndarray:
            L = local_score(df, eps=eps)
            G = global_score_simple(df, eps=eps)
            return np.power(L, A) * np.power(G, B)
        return _fn

    if preset == 4:
        def _fn(df: pd.DataFrame) -> np.ndarray:
            L = local_score_worst_out(df, eps=eps)
            G = global_score_simple(df, eps=eps)
            return np.power(L, A) * np.power(G, B)
        return _fn

    raise ValueError("preset must be one of {1,2,3,4}")
# ================================
# ========== Experiment ==========
# ================================


def build_local_marker_inputs(
    *,
    ctx,
    labels,
    note_df: pd.DataFrame,
    edge_fc: pd.DataFrame,
    edge_delta: pd.DataFrame,
    only_high_markers: bool = True,
) -> pd.DataFrame:
    """Build one merged long-form table for local marker ranking."""
    # Gene-level features (canonical source used in prioritization)
    gene_features = compute_gene_features(
        note_df=note_df,
        edge_fc=edge_fc,
        edge_delta=edge_delta,
        mean_norm=ctx.to_mean_norm_df(),
    )
    gene_features = gene_features[["n_low", "n_grey", "n_high", "n_edges", "gap"]]

    # Per-(group,gene) max fc/delta over incoming edges
    group_gene_max = compute_dst_gene_max_fc_delta(ctx, genes=labels.genes)

    # Coverage mats (vectorized)
    cov = compute_coverage_mats(ctx, genes=labels.genes, use_weighted_mean=True)
    coverage_df = cov["coverage"]
    coverage_rest_df = cov["coverage_rest"]
    coverage_neighbor_df = cov["coverage_neighbor"]

    # Select (group,gene) pairs to rank
    if only_high_markers:
        mask = note_df[labels.genes].eq(1)
        pairs = mask.stack()
        pairs = pairs[pairs].reset_index()
        pairs.columns = ["group", "gene", "_flag"]
        pairs = pairs.drop(columns=["_flag"])
    else:
        pairs = note_df[labels.genes].stack().reset_index()
        pairs.columns = ["group", "gene", "label"]

    # Coverage long
    cov_long = (
        coverage_df.stack().rename("coverage_one").to_frame()
        .join(coverage_rest_df.stack().rename("coverage_rest"))
        .join(coverage_neighbor_df.stack().rename("coverage_neighbor"))
        .reset_index()
        .rename(columns={"level_0": "group", "level_1": "gene"})
    )

    merged = (
        pairs
        .merge(group_gene_max, on=["group", "gene"], how="left")
        .merge(cov_long, on=["group", "gene"], how="left")
        .merge(gene_features.reset_index(), on="gene", how="left")
    )

    return merged


def rank_local_markers(
    merged_long: pd.DataFrame,
    *,
    A: float = 1.0,
    B: float = 0.5,
    score_fn: Optional[Callable[[pd.DataFrame], np.ndarray]] = None,
    score_col: str = "local_weight",
) -> Dict[str, List[str]]:
    """Rank genes within each group.

    - If score_fn is None: uses weight_local_prioritized(df, A=A, B=B).
    - Else: uses the provided score_fn(df).
    """
    df = merged_long.copy()

    if score_fn is None:
        df[score_col] = weight_local_prioritized(df, A=A, B=B)
    else:
        df[score_col] = score_fn(df)

    df = df.sort_values(["group", score_col], ascending=[True, False])

    out: Dict[str, List[str]] = {}
    for g, sdf in df.groupby("group", sort=False):
        out[str(g)] = sdf["gene"].astype(str).tolist()

    return out
