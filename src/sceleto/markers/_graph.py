from ._base import MarkersBase, MissingPAGAError, NotComputedError

class MarkersGraph(MarkersBase):
    """
    Graph(edge)-level marker workflow (skeleton).
    """

    def build_graph(self, top_k: int = 5, **kwargs):
        """
        Placeholder for building/deriving the working graph.
        """
        raise NotImplementedError("MarkersGraph.build_graph is not implemented yet.")

    def compute_delta(self, **kwargs):
        """
        Placeholder for computing edge-wise delta (gene jumps).
        """
        raise NotImplementedError("MarkersGraph.compute_delta is not implemented yet.")

    def set_threshold(self, threshold: float):
        """
        Store threshold only (skeleton).
        """
        self._threshold = float(threshold)
        return self

    def to_edge_table(self, *args, **kwargs):
        """
        Placeholder for exporting edge×gene table.
        """
        raise NotImplementedError("MarkersGraph.to_edge_table is not implemented yet.")

    class _Plot:
        """
        Minimal plot namespace (skeleton).
        """
        def threshold_curve(self, *args, **kwargs):
            raise NotImplementedError("MarkersGraph.plot.threshold_curve is not implemented yet.")

        def gene(self, *args, **kwargs):
            raise NotImplementedError("MarkersGraph.plot.gene is not implemented yet.")

    @property
    def plot(self):
        return self._Plot()

