"""Composable gene-name filters for marker display.

Built-in exclude patterns cover common non-informative gene classes
(ncRNAs, ribosomal, mitochondrial).  The ``include`` parameter supports
future database-backed gene sets (e.g. transcription factors, surface
proteins).

Examples
--------
>>> gf = GeneFilter(exclude=["ncrna", "ribosomal", "mito"])
>>> gf("ACTB")
True
>>> gf("AC012345.1")
False
>>> gf("RPL11")
False

Combine with an include set (only keep genes in the set AND not excluded):

>>> tf_genes = {"PAX6", "SOX2", "FOXP2"}
>>> gf = GeneFilter(exclude=["ncrna"], include=tf_genes)
>>> gf("PAX6")
True
>>> gf("ACTB")       # not in include set
False
"""

from __future__ import annotations

import re
from typing import (
    Collection,
    Dict,
    List,
    Optional,
    Sequence,
    Union,
)

# ---------------------------------------------------------------------------
# Built-in exclude patterns (regex, applied with re.search)
# ---------------------------------------------------------------------------
EXCLUDE_PATTERNS: Dict[str, str] = {
    "ncrna": r"\.",              # AC012345.1, ENSG…, etc.
    "linc": r"^LINC",            # long intergenic non-coding RNA
    "ribosomal": r"^RP[SL]\d",   # ribosomal proteins (RPS6, RPL11, …)
    "mito": r"^MT-",             # mitochondrial genes
}


class GeneFilter:
    """Composable gene-name filter.

    Parameters
    ----------
    exclude
        Names of built-in exclude patterns (keys of ``EXCLUDE_PATTERNS``)
        and/or raw regex strings.  A gene is dropped if it matches **any**
        of the patterns.
    include
        Optional gene-name collection.  When provided, only genes present
        in this set are kept (after exclude filtering).
    """

    def __init__(
        self,
        exclude: Optional[Sequence[str]] = None,
        include: Optional[Collection[str]] = None,
    ) -> None:
        self._exclude_names: List[str] = []
        self._exclude_regexes: List[re.Pattern] = []

        for pat in (exclude or []):
            if pat in EXCLUDE_PATTERNS:
                self._exclude_names.append(pat)
                self._exclude_regexes.append(re.compile(EXCLUDE_PATTERNS[pat]))
            else:
                # Treat as raw regex
                self._exclude_names.append(pat)
                self._exclude_regexes.append(re.compile(pat))

        self._include: Optional[frozenset[str]] = (
            frozenset(include) if include is not None else None
        )

    # ----- core predicate -----

    def __call__(self, gene: str) -> bool:
        """Return ``True`` if *gene* should be **kept**."""
        for rx in self._exclude_regexes:
            if rx.search(gene):
                return False
        if self._include is not None and gene not in self._include:
            return False
        return True

    # ----- convenience: filter a list -----

    def filter(self, genes: Sequence[str]) -> List[str]:
        """Return the sublist of *genes* that pass the filter."""
        return [g for g in genes if self(g)]

    # ----- repr -----

    def __repr__(self) -> str:
        parts: List[str] = []
        if self._exclude_names:
            parts.append(f"exclude={self._exclude_names!r}")
        if self._include is not None:
            n = len(self._include)
            parts.append(f"include=<{n} genes>")
        return f"GeneFilter({', '.join(parts)})"
