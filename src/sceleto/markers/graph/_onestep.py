from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable, Dict, List, Optional, Tuple

import pandas as pd


@dataclass
class MarkerGraphRun:
    """Container for one-step marker-graph pipeline results.

    Notes
    -----
    - Hierarchical/specific outputs can be toggled independently.
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

    # Hierarchical prioritization (optional)
    hierarchical_state: Optional[Any] = None
    hierarchical_marker_log: Optional[Dict[str, List[str]]] = None
    hierarchical_history: Optional[List[pd.DataFrame]] = None

    # Specific marker ranking (optional)
    specific_ranking_df: Optional[pd.DataFrame] = None
    specific_marker_log: Optional[Dict[str, List[str]]] = None

    def plot_gene_edges_fc(self, gene: str, **kwargs):
        return self.viz.plot_gene_edges_fc(gene, **kwargs)

    def plot_gene_levels_with_edges(self, gene: str, **kwargs):
        return self.viz.plot_gene_levels_with_edges(gene, **kwargs)

    def plot_highlight_edges(self, edges, **kwargs):
        return self.viz.plot_highlight_edges(edges, **kwargs)


def run_marker_graph(
    adata: Any,
    *,
    groupby: str,
    thres_fc: float,
    # Mode switches
    hierarchical_markers: bool = False,
    specific_markers: bool = True,
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
    # Hierarchical prioritization params
    weight_fn: Optional[Callable[[pd.DataFrame], Any]] = None,
    weight_col: str = "w1",
    stop_if_unique: bool = True,
    corr_cutoff: float = 0.9,
    max_iters: Optional[int] = None,
) -> MarkerGraphRun:
    """One-step wrapper: context -> edge metrics -> labels -> viz -> (optional) hierarchical/specific outputs."""
    from . import (
        build_context,
        compute_fc_delta,
        edge_gene_df_to_matrices,
        label_levels,
        labels_to_note_df,
        build_graph_and_pos_from_ctx,
        build_gene_edge_fc_from_edge_gene_df,
        GraphVizContext,
        PrioritizationState,
        run_iterative_prioritization,
    )
    from ._features import add_default_weight
    from ._local import build_local_marker_inputs, weight_local_prioritized, make_specific_score_fn

    if (not hierarchical_markers) and (not specific_markers):
        raise ValueError(
            "At least one of hierarchical_markers=True or specific_markers=True must be set."
        )

    paga = getattr(adata, "uns", {}).get("paga", None)
    if paga is None or "connectivities" not in paga:
        raise ValueError("PAGA not found. Run sc.tl.paga(adata, groups=groupby) first.")
    if "pos" not in paga:
        raise ValueError(
            "PAGA positions (adata.uns['paga']['pos']) not found. "
            "Run sc.pl.paga(adata, show=False) after sc.tl.paga(...) to populate paga['pos']."
        )

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

    specific_ranking_df: Optional[pd.DataFrame] = None
    specific_marker_log: Optional[Dict[str, List[str]]] = None

    if specific_markers:
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

    hierarchical_state: Optional[Any] = None
    hierarchical_marker_log: Optional[Dict[str, List[str]]] = None
    hierarchical_history: Optional[List[pd.DataFrame]] = None

    if hierarchical_markers:
        mean_norm = ctx.to_mean_norm_df()
        hierarchical_state = PrioritizationState(
            edge_fc=edge_fc,
            edge_delta=edge_delta,
            note_df=note_df,
            mean_norm=mean_norm,
        )

        if weight_fn is None and weight_col != "weight":

            def _default_weight_to_col(df: pd.DataFrame) -> pd.Series:
                tmp = add_default_weight(df, out_col=weight_col, inplace=False)
                return tmp[weight_col]

            weight_fn_used = _default_weight_to_col
        else:
            weight_fn_used = weight_fn

        hierarchical_marker_log, hierarchical_history = run_iterative_prioritization(
            hierarchical_state,
            weight_col=weight_col,
            weight_fn=weight_fn_used,
            corr_cutoff=corr_cutoff,
            stop_if_unique=stop_if_unique,
            max_iters=max_iters,
        )

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
        hierarchical_state=hierarchical_state,
        hierarchical_marker_log=hierarchical_marker_log,
        hierarchical_history=hierarchical_history,
        specific_ranking_df=specific_ranking_df,
        specific_marker_log=specific_marker_log,
    )