from __future__ import annotations
from ._base import MarkersBase, NotComputedError

from ._classic_v1compat import (
    find_markers as _find_markers,
    show_marker as _show_marker,
    show_marker_plot as _show_marker_plot,
)

class MarkersClassic(MarkersBase):
    """
    A thin wrapper that directly calls the original classic marker logic.
    - __init__: run the original find_markers() → store in self._mks
    - show_marker: wrap the original show_marker(...) as-is
    - plot_marker: wrap the original show_marker_plot(...) as-is
    """
    def __init__(self, adata, groupby: str, **kwargs):
        super().__init__(adata, groupby, **kwargs)
        self._mks = _find_markers(adata, groupby, **kwargs)
        self.mks = self._mks  # (compatibility) keep the same name self.mks as in the original class

    def show_marker(self, celltype=None, toshow: int = 40, result: bool = False, **kw):
        # Keep the original behavior
        return _show_marker(self._mks, celltype=celltype, toshow=toshow, result=result, **kw)

    def plot_marker(self, **kw):
        # Keep the original behavior (keyword-only args: toshow=, T=, etc.)
        return _show_marker_plot(self.adata, self.groupby, self._mks, **kw)

