from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from pandas.api.types import CategoricalDtype
from scipy import sparse


@dataclass(frozen=True)
class MarkerContext:
    """Cached matrices for marker pipeline.

    Notes
    -----
    - Matrices are (n_groups, n_genes) aligned to `groups` and `genes`.
    - `mean` is unnormalized group-wise mean expression.
    - `mean_norm` is per-gene normalized mean (0..1) by max across groups.
    - `n_expr` is number of expressing cells per gene per group (nnz).
    - `frac_expr` is fraction of expressing cells per gene per group (n_expr / n_cells_in_group).
    """
    groupby: str
    groups: List[str]          # length C
    genes: np.ndarray          # shape (G,)

    group_to_idx: Dict[str, int]

    mean: np.ndarray
    mean_norm: np.ndarray
    n_expr: np.ndarray
    frac_expr: np.ndarray

    # number of cells per group (aligned to `groups`)
    n_cells: np.ndarray        # shape (C,)

    # Graph-related caches (phase 1: only undirected PAGA edges)
    undirected_edges: Optional[List[Tuple[str, str]]] = None
    pos_df: Optional[pd.DataFrame] = None

    def get_group_idx(self, group: str) -> int:
        """Return row index for a group label."""
        return self.group_to_idx[str(group)]

    def to_mean_df(self) -> pd.DataFrame:
        """Convenience: convert `mean` matrix to DataFrame."""
        return pd.DataFrame(self.mean, index=self.groups, columns=self.genes)

    def to_mean_norm_df(self) -> pd.DataFrame:
        """Convenience: convert `mean_norm` matrix to DataFrame."""
        return pd.DataFrame(self.mean_norm, index=self.groups, columns=self.genes)

    def to_n_cells_series(self) -> pd.Series:
        """Return cell counts as a Series aligned to ctx.groups."""
        # Series index is cluster label, value is cell count
        return pd.Series(self.n_cells, index=self.groups, name="n_cells")


def get_groups(adata, groupby: str, exclude: Optional[List[str]] = None) -> List[str]:
    """Get group labels with a stable order."""
    if groupby not in adata.obs:
        raise KeyError(f"`groupby`='{groupby}' not found in adata.obs.")

    s = adata.obs[groupby]
    if isinstance(s.dtype, CategoricalDtype):
        groups = [str(x) for x in s.cat.categories]
    else:
        groups = [str(x) for x in pd.unique(s.astype(str))]

    if exclude:
        ex = set(map(str, exclude))
        groups = [g for g in groups if g not in ex]

    return groups


def top_k_edges_symmetric(con_df: pd.DataFrame, k: int) -> pd.DataFrame:
    """Keep top-k edges per node and symmetrize."""
    if con_df.shape[0] != con_df.shape[1]:
        raise ValueError("con_df must be square.")
    if not con_df.index.equals(con_df.columns):
        raise ValueError("con_df index/columns must match.")

    m = con_df.copy()
    np.fill_diagonal(m.values, 0.0)

    keep = np.zeros_like(m.values, dtype=bool)
    if k > 0:
        for i in range(m.shape[0]):
            row = m.values[i]
            idx = np.argsort(row)[::-1][:k]
            keep[i, idx] = row[idx] > 0

    keep = keep | keep.T
    return m.where(keep, other=0.0)


def build_context(
    adata,
    groupby: str,
    *,
    use_raw: bool = True,
    exclude: Optional[List[str]] = None,
    min_cells_per_group: int = 0,
    min_expr_cells_per_gene: int = 0,
    dtype=np.float32,
    k: int = 5,
) -> MarkerContext:
    """Build MarkerContext with expression stats + trimmed PAGA undirected edges."""
    # context with expression values
    ctx = _build_expression_context(
        adata,
        groupby,
        use_raw=use_raw,
        exclude=exclude,
        min_cells_per_group=min_cells_per_group,
        min_expr_cells_per_gene=min_expr_cells_per_gene,
        dtype=dtype,
    )

    # context with graph structures
    undirected_edges, pos_df = _extract_paga_undirected_edges(
        adata,
        groups=ctx.groups,  # enforce exact group order
        k=k,
    )

    return replace(ctx, undirected_edges=undirected_edges, pos_df=pos_df)


def _build_expression_context(
    adata,
    groupby: str,
    *,
    use_raw: bool,
    exclude: Optional[List[str]],
    min_cells_per_group: int,
    min_expr_cells_per_gene: int,
    dtype,
) -> MarkerContext:
    """Internal: compute expression summary matrices."""
    if use_raw and adata.raw is None:
        raise ValueError(
            "adata.raw is None but use_raw=True. Set use_raw=False or populate adata.raw."
        )

    groups = get_groups(adata, groupby, exclude=exclude)
    group_to_idx = {g: i for i, g in enumerate(groups)}

    if use_raw:
        X = adata.raw.X
        genes = adata.raw.var_names.to_numpy()
    else:
        X = adata.X
        genes = adata.var_names.to_numpy()

    if not sparse.issparse(X):
        X = sparse.csr_matrix(X)
    else:
        X = X.tocsr()

    C = len(groups)
    G = X.shape[1]
    if len(genes) != G:
        raise ValueError(f"Gene dimension mismatch: len(var_names)={len(genes)} vs X.shape[1]={G}")

    mean = np.zeros((C, G), dtype=dtype)        # mean expression
    n_expr = np.zeros((C, G), dtype=dtype)      # number of expressing cells
    frac_expr = np.zeros((C, G), dtype=dtype)   # fraction of expressing cells

    n_cells_arr = np.zeros(C, dtype=np.int32)   # cell count

    obs_groups = adata.obs[groupby].astype(str).to_numpy()

    for i, g in enumerate(groups):
        mask = (obs_groups == g)
        n_cells = int(mask.sum())
        n_cells_arr[i] = n_cells

        if n_cells < min_cells_per_group:
            continue

        Xg = X[mask]

        expr_cells = Xg.getnnz(axis=0).astype(dtype, copy=False)
        n_expr[i] = expr_cells
        frac_expr[i] = expr_cells / max(n_cells, 1)

        m = np.asarray(Xg.mean(axis=0)).ravel().astype(dtype, copy=False)

        if min_expr_cells_per_gene > 0:
            m = m.copy()
            m[expr_cells < min_expr_cells_per_gene] = 0

        mean[i] = m

    max_per_gene = mean.max(axis=0)
    mean_norm = mean.copy()
    nz = max_per_gene > 0
    mean_norm[:, nz] /= max_per_gene[nz]

    return MarkerContext(
        groupby=groupby,
        groups=groups,
        genes=genes,
        group_to_idx=group_to_idx,
        mean=mean,
        mean_norm=mean_norm,
        n_expr=n_expr,
        frac_expr=frac_expr,
        n_cells=n_cells_arr,
        undirected_edges=None,
        pos_df=None,
    )


def _extract_paga_undirected_edges(
    adata,
    *,
    groups: List[str],
    k: int,
) -> Tuple[List[Tuple[str, str]], pd.DataFrame]:
    """Extract trimmed PAGA undirected edges aligned to provided `groups` order.

    Requires:
    - adata.uns['paga']['connectivities']
    - adata.uns['paga']['pos']
    """
    paga = adata.uns.get("paga", None)
    if paga is None or "connectivities" not in paga:
        raise ValueError("PAGA connectivities not found. Run sc.tl.paga(adata, ...) first.")

    if "pos" not in paga:
        raise ValueError("PAGA pos not found. Run sc.pl.paga(adata, show=False) after sc.tl.paga(...) to populate paga['pos'].")

    con = paga["connectivities"]
    con = con.toarray() if sparse.issparse(con) else np.asarray(con)

    con_df = pd.DataFrame(con, index=groups, columns=groups)

    trimmed_con_df = top_k_edges_symmetric(con_df, k=k)
    mat = trimmed_con_df.to_numpy()

    pos_df = pd.DataFrame(paga["pos"], index=groups)

    rows, cols = np.where((mat > 0) & (~np.eye(mat.shape[0], dtype=bool)))

    undirected_edges = sorted({
        tuple(sorted((str(trimmed_con_df.index[r]), str(trimmed_con_df.columns[c]))))
        for r, c in zip(rows, cols)
    })

    return undirected_edges, pos_df
