from sceleto.markers._gene_filter import _load_categories


def available_gene_categories() -> list[str]:
    """Return the list of available gene category names."""
    return list(_load_categories().keys())


__all__ = ["available_gene_categories"]
