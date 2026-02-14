"""Composable gene-name filters for marker display.

Uses curated gene lists from ``gene_categories.json`` to exclude or include
genes by biological category (e.g. mitochondrial, lncRNA, transcription
factors, surface proteins).

Examples
--------
>>> gf = GeneFilter(exclude=["Mito_RNA", "lncRNA"])
>>> gf("MT-CO1")
False
>>> gf("ACTB")
True

Keep only transcription factors:

>>> gf = GeneFilter(include=["Transcription Factor"])
>>> gf("PAX6")
True
>>> gf("ACTB")
False

Exclude by name pattern (prefix matching, "." means contains-dot):

>>> gf = GeneFilter(exclude=["lncRNA"], name_exclude=["LINC", "LOC", "."])
>>> gf("LINC00115")
False
>>> gf("AL035401.1")
False
"""

from __future__ import annotations

import json
from importlib import resources
from typing import (
    FrozenSet,
    List,
    Optional,
    Sequence,
)


def _load_categories() -> dict[str, list[str]]:
    """Load gene categories from the bundled JSON file."""
    ref = resources.files("sceleto.data").joinpath("gene_categories.json")
    with resources.as_file(ref) as path:
        with open(path) as f:
            return json.load(f)


def available_categories() -> list[str]:
    """Return the list of available gene category names."""
    return list(_load_categories().keys())


class GeneFilter:
    """Composable gene-name filter.

    Parameters
    ----------
    exclude
        Category names from ``gene_categories.json``.  A gene is dropped
        if it belongs to **any** of the listed categories.
    include
        Category names from ``gene_categories.json``.  When provided, only
        genes present in the union of these categories are kept (after
        exclude filtering).
    name_exclude
        List of name patterns to exclude.  Multi-character entries use
        **prefix** matching (``startswith``).  The special entry ``"."``
        excludes any gene whose name contains a dot (e.g. ``AL035401.1``).
    dot_filter
        Deprecated.  Use ``name_exclude=["."]`` instead.
    """

    def __init__(
        self,
        exclude: Optional[Sequence[str]] = None,
        include: Optional[Sequence[str]] = None,
        name_exclude: Optional[Sequence[str]] = None,
        dot_filter: bool = False,
    ) -> None:
        categories = _load_categories()

        self._exclude_names: List[str] = []
        self._exclude_genes: FrozenSet[str] = frozenset()
        if exclude:
            genes: set[str] = set()
            for cat in exclude:
                if cat not in categories:
                    raise ValueError(
                        f"Unknown category {cat!r}. "
                        f"Available: {list(categories.keys())}"
                    )
                self._exclude_names.append(cat)
                genes.update(categories[cat])
            self._exclude_genes = frozenset(genes)

        self._include_names: List[str] = []
        self._include_genes: Optional[FrozenSet[str]] = None
        if include:
            genes = set()
            for cat in include:
                if cat not in categories:
                    raise ValueError(
                        f"Unknown category {cat!r}. "
                        f"Available: {list(categories.keys())}"
                    )
                self._include_names.append(cat)
                genes.update(categories[cat])
            self._include_genes = frozenset(genes)

        # Build name_exclude lists
        raw = list(name_exclude) if name_exclude else []
        if dot_filter and "." not in raw:
            raw.append(".")
        self._name_prefixes: List[str] = [p for p in raw if p != "."]
        self._name_dot: bool = "." in raw
        self._name_exclude_raw: List[str] = raw

    # ----- core predicate -----

    def __call__(self, gene: str) -> bool:
        """Return ``True`` if *gene* should be **kept**."""
        if self._name_dot and "." in gene:
            return False
        if self._name_prefixes and any(gene.startswith(p) for p in self._name_prefixes):
            return False
        if gene in self._exclude_genes:
            return False
        if self._include_genes is not None and gene not in self._include_genes:
            return False
        return True

    # ----- convenience: filter a list -----

    def filter(self, genes: Sequence[str]) -> List[str]:
        """Return the sublist of *genes* that pass the filter."""
        return [g for g in genes if self(g)]

    # ----- repr -----

    def __repr__(self) -> str:
        parts: List[str] = []
        if self._name_exclude_raw:
            parts.append(f"name_exclude={self._name_exclude_raw!r}")
        if self._exclude_names:
            parts.append(f"exclude={self._exclude_names!r}")
        if self._include_names:
            parts.append(f"include={self._include_names!r}")
        return f"GeneFilter({', '.join(parts)})"
