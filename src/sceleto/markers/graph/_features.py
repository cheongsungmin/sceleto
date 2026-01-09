from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Sequence

import numpy as np
import pandas as pd

WeightFunc = Callable[[pd.DataFrame], pd.Series]


def add_default_weight(
    features: pd.DataFrame,
    *,
    out_col: str = "weight",
    eps: float = 1e-9,
    inplace: bool = False,
) -> pd.DataFrame:
    """Add default weight column.

    Default
    -------
    weight = n_high * mean_delta * n_edges

    Notes
    -----
    - `eps` is kept for API stability; currently not used by default formula.
    """
    df = features if inplace else features.copy()

    df[out_col] = (
        df["n_high"]
        * df["mean_delta"].to_numpy(dtype=float)
        * df["n_edges"].to_numpy(dtype=float)
    )
    return df


@dataclass(frozen=True)
class WeightRule:
    """A derived score column computed from a features DataFrame."""
    name: str
    func: WeightFunc
    description: str = ""


def expr_rule(name: str, expr: str, *, description: str = "") -> WeightRule:
    """Create a WeightRule from a pandas eval expression.

    Notes
    -----
    - Expression can use feature columns and a small safe namespace (np).
    - Uses python engine to avoid numexpr surprises.
    """
    def _fn(df: pd.DataFrame) -> pd.Series:
        local_dict = {c: df[c] for c in df.columns}
        local_dict["np"] = np
        out = pd.eval(expr, local_dict=local_dict, engine="python")
        return out if isinstance(out, pd.Series) else pd.Series(out, index=df.index)

    return WeightRule(name=name, func=_fn, description=description)


def add_score_columns(
    features: pd.DataFrame,
    rules: Sequence[WeightRule],
    *,
    inplace: bool = False,
) -> pd.DataFrame:
    """Add derived score columns to a features DataFrame."""
    df = features if inplace else features.copy()
    for rule in rules:
        s = rule.func(df)
        if not isinstance(s, pd.Series):
            s = pd.Series(s, index=df.index)
        df[rule.name] = s
    return df


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
