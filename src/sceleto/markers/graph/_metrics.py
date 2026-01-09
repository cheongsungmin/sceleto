from __future__ import annotations

from typing import Callable, Dict, Literal, Optional, Tuple, TYPE_CHECKING

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    import networkx as nx

Edge = Tuple[object, object]


def compute_fc_delta(
    ctx,
    *,
    thres_fc: float = 3.0,
    eps: float = 1e-3,
    min_mean_any: float = 0.05,
    min_mean_high: float = 0.5,
    min_frac_high: float = 0.2,
    max_mean_low: float = 0.2,
    min_nexpr_any: int = 0,
) -> pd.DataFrame:
    """Compute directional edge-gene FC/delta over PAGA undirected edges.

    Returns
    -------
    pd.DataFrame with columns:
      - start, end : directional (low -> high)
      - gene
      - fc
      - delta
    """
    if ctx.undirected_edges is None:
        raise ValueError("ctx.undirected_edges is None. Build context with PAGA first.")

    rows = []
    genes = ctx.genes

    for a, b in ctx.undirected_edges:
        ia = ctx.get_group_idx(a)
        ib = ctx.get_group_idx(b)

        mA = ctx.mean_norm[ia]
        mB = ctx.mean_norm[ib]

        fA = ctx.frac_expr[ia]
        fB = ctx.frac_expr[ib]

        nA = ctx.n_expr[ia]
        nB = ctx.n_expr[ib]

        # gene-wise decide high/low
        high_is_A = mA >= mB  # (G,) bool
        m_high = np.where(high_is_A, mA, mB)
        m_low = np.where(high_is_A, mB, mA)

        f_high = np.where(high_is_A, fA, fB)
        n_any = np.maximum(nA, nB)

        # filters
        keep = np.ones_like(m_high, dtype=bool)  # shape=(len(genes),)

        if min_nexpr_any > 0:
            keep &= (n_any >= min_nexpr_any)

        keep &= ~((mA < min_mean_any) & (mB < min_mean_any))
        keep &= (m_high >= min_mean_high)
        keep &= (f_high >= min_frac_high)
        keep &= (m_low <= max_mean_low)

        fc = (m_high + eps) / (m_low + eps)
        delta = m_high - m_low

        keep &= (fc >= thres_fc)

        idx = np.where(keep)[0]
        if idx.size == 0:
            continue

        # start/end = low->high (directional per gene)
        start = np.where(high_is_A, b, a)  # if A is high, low is B
        end = np.where(high_is_A, a, b)

        for j in idx:
            rows.append((str(start[j]), str(end[j]), str(genes[j]), float(fc[j]), float(delta[j])))

    return pd.DataFrame(rows, columns=["start", "end", "gene", "fc", "delta"])


def edge_gene_df_to_matrices(
    edge_gene_df: pd.DataFrame,
    *,
    agg: Literal["max", "mean", "first"] = "max",
    fill_value: float = 0.0,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Convert long-form edge-gene table to edge×gene matrices.

    Parameters
    ----------
    edge_gene_df
        Must contain columns: ["start","end","gene","fc","delta"].
        Each row corresponds to a directional edge (start->end) and a gene.
    agg
        How to aggregate if duplicate (edge, gene) entries exist.
        Usually "max" is safe.
    fill_value
        Fill value for missing (edge, gene) pairs.

    Returns
    -------
    edge_fc
        DataFrame indexed by edge string "start->end", columns are genes, values are fc.
    edge_delta
        Same shape as edge_fc, values are delta.
    """
    required = {"start", "end", "gene", "fc", "delta"}
    missing = required - set(edge_gene_df.columns)
    if missing:
        raise ValueError(f"edge_gene_df missing columns: {sorted(missing)}")

    df = edge_gene_df.copy()
    df["start"] = df["start"].astype(str)
    df["end"] = df["end"].astype(str)
    df["gene"] = df["gene"].astype(str)

    # canonical edge id used throughout the pipeline
    df["edge"] = df["start"] + "->" + df["end"]

    edge_fc = (
        df.pivot_table(index="edge", columns="gene", values="fc", aggfunc=agg)
        .fillna(fill_value)
        .astype(float)
    )
    edge_delta = (
        df.pivot_table(index="edge", columns="gene", values="delta", aggfunc=agg)
        .fillna(fill_value)
        .astype(float)
    )

    # ensure identical column order
    edge_delta = edge_delta.reindex(columns=edge_fc.columns, fill_value=fill_value)

    return edge_fc, edge_delta


def build_gene_edge_fc_from_edge_gene_df(
    edge_gene_fc_df: pd.DataFrame,
    *,
    G: Optional["nx.Graph"] = None,
    strip_prefix: bool = True,
    node_cast: Optional[Callable[[str], object]] = None,
) -> Dict[str, Dict[Edge, float]]:
    """Build {gene: {(u,v): fc}} from a DataFrame shaped (edge_str x gene).

    Parameters
    ----------
    edge_gene_fc_df
        DataFrame indexed by edge string (e.g., "start->end"), columns are genes, values are fc.
    G
        Optional networkx graph used to infer node casting.
    strip_prefix
        If True, helper will strip node prefixes like "cluster@33" -> "33" (depends on _as_edge_tuple).
    node_cast
        Optional callable to cast node strings into desired dtype (e.g., int).

    Returns
    -------
    Dict mapping gene -> dict of (u, v) edge tuple -> fc value
    """
    # Lazy import to avoid hard dependency / circular imports at import-time.
    from ._viz import _infer_node_cast, _as_edge_tuple  # type: ignore

    if node_cast is None and G is not None:
        node_cast = _infer_node_cast(G)

    out: Dict[str, Dict[Edge, float]] = {}
    for gene in edge_gene_fc_df.columns:
        m: Dict[Edge, float] = {}
        s = edge_gene_fc_df[gene]
        for edge_str, val in s.items():
            if val is None or float(val) == 0.0:
                continue
            u, v = _as_edge_tuple(edge_str, strip_prefix=strip_prefix, node_cast=node_cast)
            m[(u, v)] = float(val)
        out[str(gene)] = m

    return out
