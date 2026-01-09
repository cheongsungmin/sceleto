from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Dict, List, Optional, Tuple

import pandas as pd

from ._features import WeightFunc, add_default_weight, compute_gene_features

WeightFn = WeightFunc  # alias for readability


@dataclass
class PrioritizationState:
    """Mutable state for iterative prioritization."""
    edge_fc: pd.DataFrame            # edges x genes
    edge_delta: pd.DataFrame         # edges x genes
    note_df: pd.DataFrame            # clusters x genes (-1/0/+1), typically level3
    mean_norm: pd.DataFrame          # clusters x genes (0..1)

    corr_df: Optional[pd.DataFrame] = None
    marker_log: Optional[Dict[str, List[str]]] = None
    iter_count: int = 0


def _build_gene_corr(mean_norm: pd.DataFrame, genes: List[str]) -> pd.DataFrame:
    """Compute gene-gene correlation on mean_norm (cluster means)."""
    sub = mean_norm[genes]
    return sub.corr()


def _edges_to_keep_after_corr_filter(
    edge_fc: pd.DataFrame,
    corr_df: pd.DataFrame,
    seed_gene: str,
    *,
    corr_cutoff: float,
) -> List[str]:
    """Keep edges NOT incident to genes highly correlated with seed_gene."""
    if seed_gene not in corr_df.columns:
        raise KeyError(f"{seed_gene} not found in corr_df columns.")

    corr_genes = corr_df[seed_gene][corr_df[seed_gene] > corr_cutoff].index.astype(str).tolist()

    # exclude edges where any correlated gene is present (fc>0)
    exc_edges: List[str] = []
    for g in corr_genes:
        if g not in edge_fc.columns:
            continue
        exc_edges.extend(edge_fc.index[(edge_fc[g] > 0)].astype(str).tolist())

    exc_set = set(exc_edges)
    return [e for e in edge_fc.index.astype(str).tolist() if e not in exc_set]


def _append_gene_to_high_clusters(marker_log: Dict[str, List[str]], note_df: pd.DataFrame, gene: str) -> None:
    """Append gene to clusters where note_df[gene] == 1."""
    high_clusters = note_df.index[note_df[gene] == 1].astype(str).tolist()
    for cls in high_clusters:
        marker_log.setdefault(cls, []).append(gene)


def run_iterative_prioritization(
    state: PrioritizationState,
    *,
    weight_col: str = "weight",
    weight_fn: Optional[WeightFn] = None,
    corr_cutoff: float = 0.9,
    stop_if_unique: bool = True,
    max_iters: Optional[int] = None,
) -> Tuple[Dict[str, List[str]], List[pd.DataFrame]]:
    """Iterate: features -> compute weight -> pick top gene -> log -> prune edges."""
    if state.marker_log is None:
        state.marker_log = {str(c): [] for c in state.note_df.index.astype(str)}

    history: List[pd.DataFrame] = []

    while True:
        if max_iters is not None and state.iter_count >= max_iters:
            break
        if state.edge_fc.shape[0] == 0:
            break

        feats = compute_gene_features(state.note_df, state.edge_fc, state.edge_delta, state.mean_norm)
        feats = add_default_weight(feats, out_col="weight", inplace=False)

        if weight_fn is not None:
            w = weight_fn(feats)
            if not isinstance(w, pd.Series):
                w = pd.Series(w, index=feats.index)
            feats[weight_col] = w

        if weight_col not in feats.columns:
            raise KeyError(
                f"weight_col='{weight_col}' not found in features. "
                f"Available columns: {list(feats.columns)}"
            )

        feats = feats.sort_values(weight_col, ascending=False)
        history.append(feats)

        top_gene = str(feats.index[0])
        _append_gene_to_high_clusters(state.marker_log, state.note_df, top_gene)

        # Build correlation matrix once (based on mean_norm; stable across iterations)
        if state.corr_df is None:
            gene_list = feats.index.astype(str).tolist()
            state.corr_df = _build_gene_corr(state.mean_norm, gene_list)

        kept_edges = _edges_to_keep_after_corr_filter(
            state.edge_fc,
            state.corr_df,
            top_gene,
            corr_cutoff=corr_cutoff,
        )

        prev_n_edges = state.edge_fc.shape[0]
        state.edge_fc = state.edge_fc.loc[kept_edges]
        state.edge_delta = state.edge_delta.loc[kept_edges]
        state.iter_count += 1

        # stop if pruning makes no progress (avoid infinite loop)
        if state.edge_fc.shape[0] == prev_n_edges:
            break

        if stop_if_unique:
            sigs = [",".join(v) for _, v in sorted(state.marker_log.items(), key=lambda x: x[0])]
            if pd.Series(sigs).value_counts().max() == 1:
                break

    return state.marker_log, history
