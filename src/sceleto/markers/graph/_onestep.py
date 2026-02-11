from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable, Dict, List, Optional, Tuple, Union

import pandas as pd

from .._gene_filter import GeneFilter


@dataclass
class MarkerGraphRun:
    """Container for one-step marker-graph pipeline results.

    Notes
    -----
    - Keeps intermediate artifacts for debugging and inspection.
    """
    # Core artifacts
    ctx: Any
    edge_gene_df: pd.DataFrame
    edge_fc: pd.DataFrame
    edge_delta: pd.DataFrame
    labels: Any
    note_df: pd.DataFrame

    # Graph + viz
    G: Any
    pos: Any
    gene_edge_fc: Dict[str, Dict[Tuple[object, object], float]]
    gene_to_edges: Dict[str, List[str]]
    viz: Any

    # Specific marker outputs
    specific_ranking_df: pd.DataFrame
    specific_marker_log: Dict[str, List[str]]

    def plot_gene_edges_fc(self, gene: str, **kwargs):
        return self.viz.plot_gene_edges_fc(gene, **kwargs)

    def plot_gene_levels_with_edges(self, gene: str, **kwargs):
        return self.viz.plot_gene_levels_with_edges(gene, **kwargs)

    def plot_highlight_edges(self, edges, **kwargs):
        return self.viz.plot_highlight_edges(edges, **kwargs)

    def top_markers(
        self,
        group: str,
        n: int = 10,
        gene_filter: Optional[GeneFilter] = None,
    ) -> List[str]:
        """Return top-*n* markers for *group*, optionally filtered.

        Parameters
        ----------
        group
            Group key (must exist in ``specific_marker_log``).
        n
            Number of markers to return.
        gene_filter
            Optional :class:`GeneFilter`.  Excluded genes are skipped and
            the next-ranked gene fills the slot.
        """
        genes = self.specific_marker_log[group]
        if gene_filter is not None:
            genes = gene_filter.filter(genes)
        return genes[:n]

    def top_markers_dict(
        self,
        n: int = 10,
        gene_filter: Optional[GeneFilter] = None,
    ) -> Dict[str, List[str]]:
        """Return :meth:`top_markers` for every group as a dict."""
        return {
            g: self.top_markers(g, n=n, gene_filter=gene_filter)
            for g in self.specific_marker_log
        }


def run_marker_graph(
    adata: Any,
    *,
    groupby: str,
    thres_fc: float,
    # Specific ranking params
    specific_A: float = 1.0,
    specific_B: float = 0.5,
    specific_only_high_markers: bool = True,
    specific_score_col: str = "specific_weight",
    specific_score_preset: Optional[int] = None,
    specific_score_fn: Optional[Callable[[pd.DataFrame], object]] = None,
    # Context defaults
    use_raw: bool = True,
    k: int = 5,
    exclude: Optional[List[str]] = None,
    min_cells_per_group: int = 0,
    min_expr_cells_per_gene: int = 0,
    # FC/delta defaults
    eps: float = 1e-3,
    min_mean_any: float = 0.05,
    min_mean_high: float = 0.5,
    min_frac_high: float = 0.2,
    max_mean_low: float = 0.2,
    min_nexpr_any: int = 0,
    # Labeling defaults
    fc_cutoff: Optional[float] = None,
    label_k: float = 2.0,
    sigma_method: str = "sd",
    min_gap: float = 0.2,
    min_margin: float = 0.0,
    level: int = 3,
    # Graph/Viz defaults
    bidirectional: bool = True,
    node_size_scale: float = 10.0,
) -> MarkerGraphRun:
    """One-step wrapper: context -> edge metrics -> labels -> viz -> specific marker discovery"""
    from . import (
        build_context,
        compute_fc_delta,
        edge_gene_df_to_matrices,
        label_levels,
        labels_to_note_df,
        build_graph_and_pos_from_ctx,
        build_gene_edge_fc_from_edge_gene_df,
        GraphVizContext,
    )
    from ._local import (
        build_local_marker_inputs,
        weight_local_prioritized,
        make_specific_score_fn,
    )

    import scanpy as sc
    import matplotlib.pyplot as plt

    # --- Ensure PAGA exists and matches the current groupby ---
    paga = getattr(adata, "uns", {}).get("paga", None)
    need_paga = paga is None or "connectivities" not in paga

    if not need_paga:
        # Re-run if connectivities shape doesn't match the groupby categories
        n_groups = adata.obs[groupby].nunique()
        if paga["connectivities"].shape[0] != n_groups:
            need_paga = True

    if need_paga:
        # PAGA requires neighbors graph
        if "neighbors" not in getattr(adata, "uns", {}):
            sc.pp.neighbors(adata)
        sc.tl.paga(adata, groups=groupby)

    # --- Ensure PAGA positions exist; populate if missing or stale ---
    paga = adata.uns.get("paga", {})
    if "pos" not in paga or need_paga:
        # Populate paga['pos'] without showing any output
        # Prefer paga_compare; fallback to paga if something fails (e.g., missing embeddings)
        try:
            sc.pl.paga_compare(adata, show=False)
        except Exception:
            sc.pl.paga(adata, show=False)

        # Close any figures created internally
        plt.close("all")

    if fc_cutoff is None:
        fc_cutoff = thres_fc

    ctx = build_context(
        adata,
        groupby=groupby,
        use_raw=use_raw,
        exclude=exclude,
        min_cells_per_group=min_cells_per_group,
        min_expr_cells_per_gene=min_expr_cells_per_gene,
        k=k,
    )

    edge_gene_df = compute_fc_delta(
        ctx,
        thres_fc=thres_fc,
        eps=eps,
        min_mean_any=min_mean_any,
        min_mean_high=min_mean_high,
        min_frac_high=min_frac_high,
        max_mean_low=max_mean_low,
        min_nexpr_any=min_nexpr_any,
    )
    edge_fc, edge_delta = edge_gene_df_to_matrices(edge_gene_df)

    labels = label_levels(
        ctx,
        edge_gene_df,
        fc_cutoff=float(fc_cutoff),
        k=label_k,
        sigma_method=sigma_method,  # type: ignore[arg-type]
        min_gap=min_gap,
        min_margin=min_margin,
    )
    note_df = labels_to_note_df(ctx, labels, level=level)  # type: ignore[arg-type]

    G, pos = build_graph_and_pos_from_ctx(ctx, bidirectional=bidirectional)
    gene_edge_fc = build_gene_edge_fc_from_edge_gene_df(edge_fc, G=G)

    sub = edge_gene_df[edge_gene_df["fc"] >= float(thres_fc)]
    gene_to_edges: Dict[str, List[str]] = {}
    if len(sub) > 0:
        for g, sdf in sub.groupby("gene"):
            gene_to_edges[str(g)] = (
                sdf["start"].astype(str) + "->" + sdf["end"].astype(str)
            ).tolist()

    viz = GraphVizContext(
        G=G,
        ctx=ctx,
        note_df=note_df,
        gene_edge_fc=gene_edge_fc,
        gene_to_edges=gene_to_edges,
        node_size_scale=node_size_scale,
    )

    specific_ranking_df: pd.DataFrame
    specific_marker_log: Dict[str, List[str]]

    specific_inputs_df = build_local_marker_inputs(
        ctx=ctx,
        labels=labels,
        note_df=note_df,
        edge_fc=edge_fc,
        edge_delta=edge_delta,
        only_high_markers=specific_only_high_markers,
    )

    specific_ranking_df = specific_inputs_df.copy()

    if specific_score_fn is not None:
        # User-provided function wins
        specific_ranking_df[specific_score_col] = specific_score_fn(specific_ranking_df)
    elif specific_score_preset is not None:
        fn = make_specific_score_fn(
            int(specific_score_preset), A=specific_A, B=specific_B, eps=1e-9
        )
        specific_ranking_df[specific_score_col] = fn(specific_ranking_df)
    else:
        # Default behavior (current)
        specific_ranking_df[specific_score_col] = weight_local_prioritized(
            specific_ranking_df, A=specific_A, B=specific_B
        )

    specific_ranking_df = specific_ranking_df.sort_values(
        ["group", specific_score_col], ascending=[True, False]
    )

    specific_marker_log = {}
    for g, sdf in specific_ranking_df.groupby("group", sort=False):
        specific_marker_log[str(g)] = sdf["gene"].astype(str).tolist()

    return MarkerGraphRun(
        ctx=ctx,
        edge_gene_df=edge_gene_df,
        edge_fc=edge_fc,
        edge_delta=edge_delta,
        labels=labels,
        note_df=note_df,
        G=G,
        pos=pos,
        gene_edge_fc=gene_edge_fc,
        gene_to_edges=gene_to_edges,
        viz=viz,
        specific_ranking_df=specific_ranking_df,
        specific_marker_log=specific_marker_log,
    )