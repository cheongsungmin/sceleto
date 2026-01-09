from __future__ import annotations

from ._classic import MarkersClassic


def classic(adata, groupby: str, **kwargs) -> MarkersClassic:
    """Factory for classic (cluster-level) marker workflow."""
    return MarkersClassic(adata, groupby, **kwargs)


# Optional: make `from sceleto.markers import graph` work (module export)
from . import graph as graph  # noqa: F401  # Expose subpackage as attribute


__all__ = ["classic", "MarkersClassic", "graph"]