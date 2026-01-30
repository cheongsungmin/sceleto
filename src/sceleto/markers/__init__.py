from __future__ import annotations

from typing import Any, Optional

from ._classic import MarkersClassic


def classic(adata, groupby: str, **kwargs) -> MarkersClassic:
    """Factory for classic (cluster-level) marker workflow."""
    return MarkersClassic(adata, groupby, **kwargs)


def marker(
    adata,
    groupby: str,
    *,
    k: int = 5,
    thres_fc: float = 3.0,
    **kwargs: Any,
):
    """Graph-based marker workflow (one-word entry point).

    Parameters
    ----------
    adata
        AnnData object.
    groupby
        Key in `adata.obs` that contains cluster labels (e.g., "leiden").
    k
        Trim PAGA graph to top-k neighbors per node.
    thres_fc
        Fold-change threshold used to consider marker candidates.
    **kwargs
        Passed through to `sceleto.markers.graph.run_marker_graph`.
    """
    from .graph import run_marker_graph  # Lazy import to keep namespace clean
    return run_marker_graph(
        adata,
        groupby=groupby,
        k=k,
        thres_fc=thres_fc,
        **kwargs,
    )


def __dir__():
    return ["classic", "marker"]


__all__ = ["classic", "marker"]
