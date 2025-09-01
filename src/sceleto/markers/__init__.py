from ._classic import MarkersClassic
from ._graph import MarkersGraph

def classic(adata, groupby: str, **kwargs) -> MarkersClassic:
    """
    Factory for classic (cluster-level) marker workflow.
    Returns a MarkersClassic object (skeleton).
    """
    return MarkersClassic(adata, groupby, **kwargs)

def graph(adata, groupby: str, **kwargs) -> MarkersGraph:
    """
    Factory for graph (edge-level) marker workflow.
    Returns a MarkersGraph object (skeleton).
    """
    return MarkersGraph(adata, groupby, **kwargs)

__all__ = ["classic", "graph", "MarkersClassic", "MarkersGraph"]

