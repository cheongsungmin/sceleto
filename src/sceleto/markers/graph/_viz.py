from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Dict, Iterable, List, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
import matplotlib.patches as mpatches

import networkx as nx

Edge = Tuple[object, object]
EdgeLike = Union[str, Edge]


# ----------------------------
# Internal helpers (private)
# ----------------------------
def _infer_node_cast(G: nx.Graph) -> Callable[[str], object]:
    """Infer node casting rule from graph node types.

    If all nodes are int -> cast tokens to int; otherwise keep str.
    """
    node_types = {type(n) for n in G.nodes()}
    if len(node_types) == 1 and int in node_types:
        return lambda x: int(x)
    return lambda x: x


def _as_edge_tuple(
    e: EdgeLike,
    *,
    strip_prefix: bool = True,
    node_cast: Optional[Callable[[str], object]] = None,
) -> Edge:
    """Convert edge representations into (src, dst)."""
    if isinstance(e, tuple) and len(e) == 2:
        return (e[0], e[1])

    if not isinstance(e, str) or "->" not in e:
        raise ValueError(f"Unsupported edge format: {e}")

    a, b = e.split("->", 1)
    a = a.split("@")[-1] if strip_prefix else a
    b = b.split("@")[-1] if strip_prefix else b

    if node_cast is None:
        return (a, b)

    return (node_cast(a), node_cast(b))


def _normalize_edges(
    edges: Iterable[EdgeLike],
    *,
    strip_prefix: bool = True,
    node_cast: Optional[Callable[[str], object]] = None,
) -> List[Edge]:
    """Normalize edges into list of (u, v)."""
    return [_as_edge_tuple(e, strip_prefix=strip_prefix, node_cast=node_cast) for e in edges]


def _resolve_node_sizes(
    G: nx.Graph,
    node_sizes: Optional[Union[int, float, Dict[object, float], Sequence[float]]],
    default: float = 300.0,
) -> Union[float, List[float]]:
    """Resolve node_sizes into networkx-compatible format."""
    if node_sizes is None:
        return float(default)
    if isinstance(node_sizes, (int, float)):
        return float(node_sizes)
    if isinstance(node_sizes, dict):
        # keep node keys as-is
        return [float(node_sizes.get(n, default)) for n in G.nodes()]
    return list(node_sizes)


def _draw_node_labels(ax, pos_dict: Dict[object, Tuple[float, float]], fontsize: int = 8) -> None:
    """Draw node labels with white outline for readability."""
    for node, (x, y) in pos_dict.items():
        txt = ax.text(x, y, str(node), fontsize=fontsize, ha="center", va="center", color="black")
        txt.set_path_effects(
            [
                path_effects.Stroke(linewidth=2, foreground="white"),
                path_effects.Normal(),
            ]
        )


def _safe_vmin_vmax(values: List[float]) -> Tuple[float, float]:
    """Avoid vmin==vmax causing colorbar issues."""
    vmin = float(np.min(values)) if len(values) else 0.0
    vmax = float(np.max(values)) if len(values) else 1.0
    if vmax <= vmin:
        vmax = vmin + 1e-9
    return vmin, vmax


def _pos_df_to_pos_dict(pos_df: pd.DataFrame, groups: Sequence[str]) -> Dict[str, Tuple[float, float]]:
    """Convert pos_df (index=groups) to pos_dict."""
    if pos_df is None:
        raise ValueError("pos_df is None")

    candidates = [
        ("x", "y"),
        ("X", "Y"),
        ("pos_x", "pos_y"),
        ("umap1", "umap2"),
        ("PAGA1", "PAGA2"),
    ]
    cols = None
    for c1, c2 in candidates:
        if c1 in pos_df.columns and c2 in pos_df.columns:
            cols = (c1, c2)
            break

    if cols is None:
        # fallback to the first two numeric columns
        num_cols = [c for c in pos_df.columns if np.issubdtype(pos_df[c].dtype, np.number)]
        if len(num_cols) < 2:
            raise ValueError(f"pos_df must have >=2 numeric columns. Got: {list(pos_df.columns)}")
        cols = (num_cols[0], num_cols[1])

    c1, c2 = cols

    missing = [g for g in groups if g not in pos_df.index]
    if missing:
        raise ValueError(f"pos_df is missing groups: {missing[:5]} ... (n={len(missing)})")

    return {g: (float(pos_df.loc[g, c1]), float(pos_df.loc[g, c2])) for g in groups}


# ----------------------------
# Public API
# ----------------------------
@dataclass
class GraphVizContext:
    G: nx.Graph
    ctx: Optional["MarkerContext"] = None

    # Allow None -> auto-filled from ctx.pos_df
    pos_dict: Optional[Dict[object, Tuple[float, float]]] = None
    # Allow None -> auto-filled from ctx.n_cells with sqrt scaling
    node_sizes: Optional[Union[int, float, Dict[object, float], Sequence[float]]] = None
    node_size_scale: float = 10.0

    # Optional data slots (set once, reuse many times)
    mean_mat: Optional[pd.DataFrame] = None
    note_df: Optional[pd.DataFrame] = None
    gene_to_edges: Optional[Dict[str, List[str]]] = None
    gene_edge_fc: Optional[Dict[str, Dict[Edge, float]]] = None

    # Parsing
    strip_prefix: bool = True
    node_cast: Optional[Callable[[str], object]] = None

    # Style defaults
    connectionstyle: str = "arc3,rad=0.05"
    arrowsize: int = 12
    min_target_margin: int = 10
    label_fontsize: int = 8

    def __post_init__(self) -> None:
        if self.node_cast is None:
            self.node_cast = _infer_node_cast(self.G)

        # Auto-fill mean_mat from ctx (if not provided)
        if self.mean_mat is None and self.ctx is not None:
            # default to mean expression for node coloring
            self.mean_mat = self.ctx.to_mean_df()

        # Auto-fill pos_dict from ctx.pos_df (if not provided)
        if self.pos_dict is None:
            if self.ctx is None or self.ctx.pos_df is None:
                raise ValueError("pos_dict is None and ctx.pos_df is None. Provide pos_dict or ensure ctx has pos_df.")
            # build positions mapping matching node types via node_cast
            self.pos_dict = {
                self.node_cast(str(idx)): (float(row.iloc[0]), float(row.iloc[1]))
                for idx, row in self.ctx.pos_df.iterrows()
            }

        # Auto-fill node sizes from ctx.n_cells (if not provided)
        if self.node_sizes is None:
            if self.ctx is None or getattr(self.ctx, "n_cells", None) is None:
                raise ValueError("node_sizes is None and ctx.n_cells is missing. Provide node_sizes or ensure ctx has n_cells.")
            # sqrt(cell_count) scaling with user-provided scale factor
            n_cells_s = pd.Series(self.ctx.n_cells, index=self.ctx.groups)
            self.node_sizes = [
                float(np.sqrt(n_cells_s.get(str(node), 1))) * float(self.node_size_scale)
                for node in self.G.nodes()
            ]

    def plot_gene_edges_fc(
        self,
        gene: str,
        fc_bins: Tuple[float, ...] = (0, 3.0, 4.0, 5.0, 6.0, np.inf),
        fc_colors: Tuple[str, ...] = ("lightgrey", "#9fdab8", "#57b8d0", "#1d7eb7", "#084081"),
        background_edge_color: str = "lightgrey",
        background_edge_alpha: float = 0.25,
        edge_width: float = 1.8,
        cmap: str = "Greys",
        figsize: Tuple[float, float] = (6, 6),
        ax=None,
        title: Optional[str] = None,
    ):
        """Nodes: expression; Edges: FC bins."""
        if self.mean_mat is None:
            raise ValueError("mean_mat is required for expression-colored nodes.")
        if self.gene_edge_fc is None:
            raise ValueError("gene_edge_fc is required for FC-colored edges.")
        if len(fc_bins) - 1 != len(fc_colors):
            raise ValueError("len(fc_colors) must be len(fc_bins)-1")

        own_fig = False
        if ax is None:
            own_fig = True
            _, ax = plt.subplots(figsize=figsize)

        # Node colors (expression)
        node_values: List[float] = []
        for n in self.G.nodes():
            val = 0.0
            if (n in self.mean_mat.index) and (gene in self.mean_mat.columns):
                val = float(self.mean_mat.loc[n, gene])
            node_values.append(val)
        vmin, vmax = _safe_vmin_vmax(node_values)

        ns = _resolve_node_sizes(self.G, self.node_sizes, default=300.0)
        nodes = nx.draw_networkx_nodes(
            self.G,
            pos=self.pos_dict,
            node_color=node_values,
            node_size=ns,
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            linewidths=0.0,
            ax=ax,
        )
        _draw_node_labels(ax, self.pos_dict, fontsize=self.label_fontsize)

        def _bin_index(fc: float) -> int:
            # return bin index 0..n_bins-1
            idx = int(np.digitize([fc], fc_bins)[0] - 1)
            return max(0, min(idx, len(fc_bins) - 2))

        edge_fc_map = self.gene_edge_fc.get(gene, {})

        # Background edges
        bg_edges = [(u, v) for (u, v) in self.G.edges() if (u, v) not in edge_fc_map]
        if bg_edges:
            nx.draw_networkx_edges(
                self.G,
                pos=self.pos_dict,
                edgelist=bg_edges,
                edge_color=background_edge_color,
                alpha=background_edge_alpha,
                width=edge_width,
                arrows=True,
                arrowsize=self.arrowsize,
                connectionstyle=self.connectionstyle,
                min_target_margin=self.min_target_margin,
                ax=ax,
            )

        # FC edges grouped by bins
        bin_to_edges: Dict[int, List[Edge]] = {i: [] for i in range(len(fc_bins) - 1)}
        g_edges = set(self.G.edges())
        for (u, v), fc in edge_fc_map.items():
            if (u, v) in g_edges:
                bin_to_edges[_bin_index(float(fc))].append((u, v))

        for i, edges_i in bin_to_edges.items():
            if not edges_i:
                continue
            nx.draw_networkx_edges(
                self.G,
                pos=self.pos_dict,
                edgelist=edges_i,
                edge_color=fc_colors[i],
                alpha=1.0,
                width=edge_width,
                arrows=True,
                arrowsize=self.arrowsize,
                connectionstyle=self.connectionstyle,
                min_target_margin=self.min_target_margin,
                ax=ax,
            )

        # FC legend
        patches = []
        for i in range(len(fc_bins) - 1):
            left, right = fc_bins[i], fc_bins[i + 1]
            label = f"≥ {left:.0f}" if np.isinf(right) else f"{left:.0f}–{right:.0f}"
            patches.append(mpatches.Patch(color=fc_colors[i], label=label))
        ax.legend(
            handles=patches,
            title="FC bin",
            loc="upper left",
            bbox_to_anchor=(1.30, 0.60),
            frameon=False,
            fontsize=10,
            title_fontsize=10,
        )

        cbar = ax.figure.colorbar(nodes, ax=ax, shrink=0.7, pad=0.02)
        cbar.set_label(f"{gene} mean expression", fontsize=10)

        ax.set_axis_off()
        ax.set_title(title or f"FC on edges + expression nodes: {gene}", fontsize=11)

        if own_fig:
            plt.tight_layout()
            plt.show()
        return ax

    def plot_gene_levels_with_edges(
        self,
        gene: str,
        category_colors: Optional[Dict[int, str]] = None,
        highlight_edge_color: str = "crimson",
        background_edge_color: str = "lightgrey",
        background_edge_alpha: float = 0.25,
        edge_width: float = 1.4,
        figsize: Tuple[float, float] = (6, 6),
        ax=None,
        title: Optional[str] = None,
    ):
        """Nodes: categorical level (-1/0/1); Edges: highlight gene edges."""
        if self.note_df is None:
            raise ValueError("note_df is required for level-colored nodes.")
        if self.gene_to_edges is None:
            raise ValueError("gene_to_edges is required to highlight gene edges.")

        own_fig = False
        if ax is None:
            own_fig = True
            _, ax = plt.subplots(figsize=figsize)

        if category_colors is None:
            category_colors = {-1: "#3b82f6", 0: "#d1d5db", 1: "#ef4444"}  # blue/gray/red

        node_colors: List[str] = []
        for n in self.G.nodes():
            val = 0
            if (gene in self.note_df.columns) and (n in self.note_df.index):
                try:
                    val = int(self.note_df.loc[n, gene])
                except Exception:
                    val = 0
            node_colors.append(category_colors.get(val, category_colors[0]))

        ns = _resolve_node_sizes(self.G, self.node_sizes, default=300.0)
        nx.draw_networkx_nodes(
            self.G,
            pos=self.pos_dict,
            node_color=node_colors,
            node_size=ns,
            linewidths=0.0,
            ax=ax,
        )
        _draw_node_labels(ax, self.pos_dict, fontsize=self.label_fontsize)

        gene_edges_raw = self.gene_to_edges.get(gene, [])
        hl_edges = set(
            _normalize_edges(
                gene_edges_raw,
                strip_prefix=self.strip_prefix,
                node_cast=self.node_cast,
            )
        )

        bg_edges = [(u, v) for (u, v) in self.G.edges() if (u, v) not in hl_edges]
        if bg_edges:
            nx.draw_networkx_edges(
                self.G,
                pos=self.pos_dict,
                edgelist=bg_edges,
                edge_color=background_edge_color,
                alpha=background_edge_alpha,
                width=edge_width,
                arrows=True,
                arrowsize=self.arrowsize,
                connectionstyle=self.connectionstyle,
                min_target_margin=self.min_target_margin,
                ax=ax,
            )

        hl_edges_inG = [(u, v) for (u, v) in self.G.edges() if (u, v) in hl_edges]
        if hl_edges_inG:
            nx.draw_networkx_edges(
                self.G,
                pos=self.pos_dict,
                edgelist=hl_edges_inG,
                edge_color=highlight_edge_color,
                alpha=1.0,
                width=edge_width * 1.8,
                arrows=True,
                arrowsize=self.arrowsize,
                connectionstyle=self.connectionstyle,
                min_target_margin=self.min_target_margin,
                ax=ax,
            )

        handles = [
            mpatches.Patch(color=category_colors[-1], label="Low (-1)"),
            mpatches.Patch(color=category_colors[0], label="Undetermined (0)"),
            mpatches.Patch(color=category_colors[1], label="High (1)"),
        ]
        ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.02, 1.0), frameon=False)

        ax.set_axis_off()
        ax.set_title(title or f"Levels + gene edges: {gene}", fontsize=11)

        if own_fig:
            plt.tight_layout()
            plt.show()
        return ax

    def plot_highlight_edges(
        self,
        highlight_edges: Iterable[EdgeLike],
        node_color: str = "white",
        node_edgecolor: str = "black",
        background_edge_color: str = "lightgrey",
        background_edge_alpha: float = 0.5,
        highlight_edge_color: str = "crimson",
        highlight_edge_width: float = 2.8,
        background_edge_width: float = 1.2,
        figsize: Tuple[float, float] = (6, 6),
        ax=None,
        title: Optional[str] = None,
    ):
        """Highlight specified edges only (no gene dependency)."""
        own_fig = False
        if ax is None:
            own_fig = True
            _, ax = plt.subplots(figsize=figsize)

        hl_set = set(
            _normalize_edges(
                highlight_edges,
                strip_prefix=self.strip_prefix,
                node_cast=self.node_cast,
            )
        )

        ns = _resolve_node_sizes(self.G, self.node_sizes, default=300.0)
        nx.draw_networkx_nodes(
            self.G,
            pos=self.pos_dict,
            node_size=ns,
            node_color=node_color,
            edgecolors=node_edgecolor,
            linewidths=1.0,
            ax=ax,
        )
        _draw_node_labels(ax, self.pos_dict, fontsize=self.label_fontsize)

        bg_edges = [(u, v) for (u, v) in self.G.edges() if (u, v) not in hl_set]
        if bg_edges:
            nx.draw_networkx_edges(
                self.G,
                pos=self.pos_dict,
                edgelist=bg_edges,
                edge_color=background_edge_color,
                alpha=background_edge_alpha,
                width=background_edge_width,
                arrows=True,
                arrowsize=self.arrowsize,
                connectionstyle=self.connectionstyle,
                min_target_margin=self.min_target_margin,
                ax=ax,
            )

        hl_edges_inG = [(u, v) for (u, v) in self.G.edges() if (u, v) in hl_set]
        if hl_edges_inG:
            nx.draw_networkx_edges(
                self.G,
                pos=self.pos_dict,
                edgelist=hl_edges_inG,
                edge_color=highlight_edge_color,
                alpha=1.0,
                width=highlight_edge_width,
                arrows=True,
                arrowsize=self.arrowsize,
                connectionstyle=self.connectionstyle,
                min_target_margin=self.min_target_margin,
                ax=ax,
            )

        ax.set_axis_off()
        ax.set_title(title or "Highlighted edges", fontsize=11)

        if own_fig:
            plt.tight_layout()
            plt.show()
        return ax


def build_graph_and_pos_from_ctx(
    ctx,
    *,
    bidirectional: bool = True,
) -> Tuple[nx.DiGraph, Dict[str, Tuple[float, float]]]:
    """Build directed graph + pos_dict from MarkerContext.

    Notes
    -----
    - ctx.undirected_edges are undirected by nature.
    - If bidirectional=True, add both u->v and v->u so arrows are always visible.
    """
    if ctx.undirected_edges is None:
        raise ValueError("ctx.undirected_edges is None")
    if ctx.pos_df is None:
        raise ValueError("ctx.pos_df is None")

    G = nx.DiGraph()
    G.add_nodes_from(ctx.groups)

    if bidirectional:
        edges = []
        for u, v in ctx.undirected_edges:
            edges.append((u, v))
            edges.append((v, u))
        G.add_edges_from(edges)
    else:
        G.add_edges_from(ctx.undirected_edges)

    pos_dict = _pos_df_to_pos_dict(ctx.pos_df, ctx.groups)
    return G, pos_dict


__all__ = [
    "GraphVizContext",
    "build_graph_and_pos_from_ctx",
    "_infer_node_cast",
    "_as_edge_tuple",
]
