from ._base import MarkersBase, NotComputedError

class MarkersClassic(MarkersBase):
    """
    Classic cluster-level marker workflow (skeleton).
    """

    def top_genes(self, group: str, n: int = 20):
        """
        Placeholder. Implement later.
        """
        raise NotImplementedError("MarkersClassic.top_genes is not implemented yet.")

    class _Plot:
        """
        Minimal plot namespace (skeleton).
        """
        def dot(self, *args, **kwargs):
            raise NotImplementedError("MarkersClassic.plot.dot is not implemented yet.")

    @property
    def plot(self):
        return self._Plot()

