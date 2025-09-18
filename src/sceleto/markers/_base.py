from __future__ import annotations
from dataclasses import dataclass
from typing import Any

# ---- Minimal exceptions (names만 잡아둠) ----
class SceletoError(Exception):
    """Base error for sceleto."""
    pass

class MissingPAGAError(SceletoError):
    """Raised when PAGA info is required but missing."""
    pass

class GroupKeyError(SceletoError):
    """Raised when the given 'groupby' key is not found in adata.obs."""
    pass

class NotComputedError(SceletoError):
    """Raised when a result is requested before compute step."""
    pass

# ---- Minimal config/data holders (필요시 확장) ----
@dataclass
class MarkerConfig:
    groupby: str

# ---- Minimal base class ----
class MarkersBase:
    """
    Very small base for marker workflows. Only stores inputs.
    """
    def __init__(self, adata: Any, groupby: str, **kwargs) -> None:
        self.adata = adata
        self.groupby = groupby
        self.config = MarkerConfig(groupby=groupby)
        # NOTE: No validation here (skeleton). Add later.

    def summary(self) -> str:
        n_obs = getattr(self.adata, "n_obs", "?")
        n_vars = getattr(self.adata, "n_vars", "?")
        has_graph = hasattr(self, "_graph")
        return (
            f"{self.__class__.__name__}(groupby='{self.groupby}', "
            f"n_obs={n_obs}, n_vars={n_vars}, built_graph={has_graph})"
        )

    def __repr__(self) -> str:
        return self.summary()

