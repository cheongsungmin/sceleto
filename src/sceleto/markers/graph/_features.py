from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Sequence

import numpy as np
import pandas as pd


def compute_gene_features(
    note_df: pd.DataFrame,       # clusters x genes, -1/0/+1
    edge_fc: pd.DataFrame,       # edges x genes (0 if not present)
    edge_delta: pd.DataFrame,    # edges x genes (0 if not present)
    mean_norm: pd.DataFrame,     # clusters x genes (0..1)
) -> pd.DataFrame:
    """Compute base features per gene."""
    common_genes = [
        g for g in note_df.columns
        if (g in edge_fc.columns) and (g in edge_delta.columns) and (g in mean_norm.columns)
    ]
    if len(common_genes) == 0:
        raise ValueError("No common genes among note_df / edge_fc / edge_delta / mean_norm.")

    # Align cluster axis explicitly
    mean_aligned = mean_norm.reindex(note_df.index)

    note = note_df[common_genes].to_numpy(dtype=np.int8, copy=False)
    m = mean_aligned[common_genes].to_numpy(dtype=float, copy=False)

    fc = edge_fc[common_genes].to_numpy(dtype=float, copy=False)
    de = edge_delta[common_genes].to_numpy(dtype=float, copy=False)

    # counts from labels
    n_high = (note == 1).sum(axis=0)
    n_low = (note == -1).sum(axis=0)
    n_grey = (note == 0).sum(axis=0)

    # edge stats (positive only)
    fc_pos = np.where(fc > 0, fc, np.nan)
    de_pos = np.where(de > 0, de, np.nan)

    max_fc = np.nanmax(fc_pos, axis=0)
    mean_fc = np.nanmean(fc_pos, axis=0)
    max_de = np.nanmax(de_pos, axis=0)
    mean_de = np.nanmean(de_pos, axis=0)
    n_edges = (fc > 0).sum(axis=0)

    # gap = min(high) - max(low)
    high_mask = (note == 1)
    low_mask = (note == -1)

    any_high = high_mask.any(axis=0)
    any_low = low_mask.any(axis=0)

    high_min = np.min(np.where(high_mask, m, np.inf), axis=0)
    low_max = np.max(np.where(low_mask, m, -np.inf), axis=0)

    gap = high_min - low_max
    gap[~(any_high & any_low)] = np.nan

    out = pd.DataFrame(
        {
            "n_low": n_low,
            "n_grey": n_grey,
            "n_high": n_high,
            "max_fc": max_fc,
            "mean_fc": mean_fc,
            "max_delta": max_de,
            "mean_delta": mean_de,
            "n_edges": n_edges,
            "gap": gap,
        },
        index=pd.Index(common_genes, name="gene"),
    )
    return out
