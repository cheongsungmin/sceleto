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
        """
        Return a minimal string summary (skeleton).
        """
        return f"{self.__class__.__name__}(groupby='{self.groupby}')"

