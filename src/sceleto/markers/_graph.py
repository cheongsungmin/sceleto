from __future__ import annotations
from typing import Dict, List, Tuple, Any
import numpy as np
import pandas as pd
from collections import Counter, defaultdict
from ._base import MarkersBase, MissingPAGAError
import scanpy as sc
import time

class MarkersGraph(MarkersBase):
    """
    Edge-level marker workflow (skeleton).

    Usage (minimal):
        MG = sceleto.markers.graph(adata, 'leiden')
        MG.build_graph(k=5, thres=0.0, min_mean_cut=0.2, min_drop_cut=0.2, min_cnt_cut=0)

    After build_graph(), the following are prepared in self._graph:
        - 'con'                : full PAGA connectivity (pd.DataFrame)
        - 'con_trim'           : k-trimmed symmetric connectivity (pd.DataFrame)
        - 'pos'                : PAGA 2D positions (pd.DataFrame)
        - 'unique_pairs'       : undirected connected pairs (List[Tuple[str, str]])
        - 'edges_df'           : directed edge list DataFrame with cols ['from','to']
        - 'cluster_stats'      : {'mean':dict, 'cnt':dict, 'drop':dict} (per cluster -> np.array)
                                 (mean values are max-normalized per gene across clusters)
        - 'jumping_genes'      : dict(edge_name->List[(gene, delta)])
                                 (edge_name format: '{groupby}@A->{groupby}@B')
        - 'edges_matrix'       : edge x gene matrix (pd.DataFrame, values = Δ≥0 else 0)
        - 'mean_mat_df'        : cluster x gene mean (non-log) DataFrame (raw.to_adata().to_df() then groupby)
        - 'node_sizes'         : dict(cluster -> count)

    The stored pieces are ready for later:
        - plot_threshold() : scan thresholds over 'edges_matrix' and count survivors
        - set_threshold()  : store chosen threshold in self._graph['threshold']
        - get_markers()    : make ef (edge→markers list) & aggregate/rank

    Notes:
        * We reuse raw-layer values if present (and recommend setting adata.raw beforehand).
        * We do a simple per-gene max normalization across clusters (for Δ calculation stability).
        * No figures are drawn here; this function only prepares data structures.
    """

    # --- public API ---
    def build_graph(
        self,
        k: int = 5,
#        thres: float = 0.0,
        min_mean_cut: float = 0.2,
        min_drop_cut: float = 0.2,
        min_cnt_cut: int = 0,
        use_existing_paga: bool = True,
		verbose: bool = True
    ) -> "MarkersGraph":
        adata = self.adata
        groupby = self.groupby

        t0 = time.perf_counter()
        if verbose: print(f"[sceleto] build_graph: start (groupby='{groupby}')")

        # 0) Make sure `groupby` is categorical and PAGA is available (or compute)
        if groupby not in adata.obs:
            raise KeyError(f"'{groupby}' not found in adata.obs")
        if not pd.api.types.is_categorical_dtype(adata.obs[groupby]):
            adata.obs[groupby] = adata.obs[groupby].astype("category")

        if (not use_existing_paga) or ("paga" not in adata.uns):
            if verbose: print("[sceleto]   computing PAGA...")
            if "neighbors" not in adata.uns:
                if verbose: print("[sceleto]   neighbors missing → sc.pp.neighbors(...)")
                sc.pp.neighbors(adata)
            sc.tl.paga(adata, groups=groupby)
        else:
            if verbose: print("[sceleto]   using existing PAGA in adata.uns['paga']")

        # 1) PAGA connectivity / position → DataFrames
        con = pd.DataFrame(
            adata.uns["paga"]["connectivities"].A,
            index=adata.obs[groupby].cat.categories,
            columns=adata.obs[groupby].cat.categories,
        )
        pos = pd.DataFrame(
            adata.uns["paga"]["pos"],
            index=adata.obs[groupby].cat.categories,
            columns=["x", "y"],
        )

        con_trim = self._top_k_edges_symmetric(con, k=k)
        if verbose:
            n_nodes = con.shape[0]
            n_undir_edges = int(np.count_nonzero(np.triu(con_trim.values, 1)))
            print(f"[sceleto]   PAGA nodes={n_nodes}, trimmed undirected edges={n_undir_edges} (k={k})")

        # 2) Directed edge list (from k-trimmed connectivity)
        rows, cols = np.where((con_trim.values > 0) & (np.eye(con_trim.shape[0]) == 0))
        edges_df = pd.DataFrame({"from": con_trim.index[rows], "to": con_trim.columns[cols]})
        edges_df["from"] = [f"{groupby}@{x}" for x in edges_df["from"]]
        edges_df["to"] = [f"{groupby}@{x}" for x in edges_df["to"]]

        # undirected unique pairs for looping once per pair
        nz_pairs = np.transpose(np.nonzero(con_trim.values))
        unique_pairs = sorted({
            tuple(sorted((con_trim.index[i], con_trim.columns[j])))
            for i, j in nz_pairs
            if i != j
        })

        # 3) Per-cluster stats (cnt/drop/mean) on raw; then per-gene max-normalize the means
        if verbose: print("[sceleto]   computing per-cluster stats (cnt/drop/mean)...")
        if adata.raw is None:
            # minimal guard: mirror current to raw so downstream works
            adata.raw = adata

        clusters = list(adata.obs[groupby].cat.categories)
        mtx = adata.raw
        n_genes = mtx.shape[1]

        cnt_dt: Dict[str, np.ndarray] = {}
        drop_dt: Dict[str, np.ndarray] = {}
        mean_dt: Dict[str, np.ndarray] = {}

        for ct in clusters:
            mask = np.array(adata.obs[groupby] == ct)
            X_ct = mtx.X[mask]
            n_ct = int(mask.sum())

            # counts per gene (non-zero cells)
            pos_count = Counter(X_ct.nonzero()[1])
            cnt_arr = np.zeros(n_genes, dtype=float)
            for gidx, c in pos_count.items():
                cnt_arr[gidx] = float(c)

            drop_dt[ct] = cnt_arr / max(n_ct, 1)
            cnt_dt[ct] = cnt_arr

            # mean per gene (dense 1D)
            mean_arr = np.asarray(X_ct.mean(axis=0)).ravel()
            # (optional) mask very-low support genes
            if min_cnt_cut > 0:
                mean_arr[cnt_arr < min_cnt_cut] = 0.0
            mean_dt[ct] = mean_arr

        # per-gene max-normalization across clusters
        mean_mat = np.stack([mean_dt[ct] for ct in clusters], axis=0)  # (C, G)
        max_per_gene = mean_mat.max(axis=0)
        nz = max_per_gene > 0
        mean_mat[:, nz] = mean_mat[:, nz] / max_per_gene[nz]
        for i, ct in enumerate(clusters):
            mean_dt[ct] = mean_mat[i]
        var_names = list(adata.raw.var_names)
        normalized_mean_mat_df = pd.DataFrame(mean_mat, index=clusters, columns=var_names)

        # 4) delta per edge and gene ("jumping genes")
        if verbose: print("[sceleto]   scanning delta per edge/gene...")
        edge_keys = [f"{src}->{dst}" for src, dst in zip(edges_df["from"], edges_df["to"])]
        jumping_genes: Dict[str, List[Tuple[str, float]]] = {k: [] for k in edge_keys}

        var_names = list(adata.raw.var_names)
        for a, b in unique_pairs:
            ct_A, ct_B = str(a), str(b)
            # iterate genes once per pair; push into the appropriate directed edge
            for gi, gname in enumerate(var_names):
                # quick skip by min_cnt_cut
                if (cnt_dt[ct_A][gi] < min_cnt_cut) and (cnt_dt[ct_B][gi] < min_cnt_cut):
                    continue

                mA, mB = mean_dt[ct_A][gi], mean_dt[ct_B][gi]
                dA, dB = drop_dt[ct_A][gi], drop_dt[ct_B][gi]
                delta = float(mA - mB)

                if (delta > 0) and (mA > min_mean_cut) and (dA > min_drop_cut):
                    jumping_genes[f"{groupby}@{ct_B}->{groupby}@{ct_A}"].append((gname, delta))
                elif (-delta > 0) and (mB > min_mean_cut) and (dB > min_drop_cut):
                    jumping_genes[f"{groupby}@{ct_A}->{groupby}@{ct_B}"].append((gname, -delta))

        # 5) Edge x gene matrix (pivot)
        flat_rows: List[Tuple[str, str, float]] = []
        for edge, lst in jumping_genes.items():
            for gene, score in lst:
                flat_rows.append((edge, gene, float(score)))

        if len(flat_rows):
            df = pd.DataFrame(flat_rows, columns=["edge", "gene", "score"])
            edges_matrix = df.pivot(index="edge", columns="gene", values="score").fillna(0.0)
        else:
            edges_matrix = pd.DataFrame(index=edge_keys, columns=[], dtype=float)

        # 6) extra handy bits kept for later plotters/aggregators
        # cluster-level mean (not normalized) for coloring nodes on graphs later
        mean_mat_df = adata.raw.to_adata().to_df()
        mean_mat_df[groupby] = adata.obs[groupby].values
        mean_mat_df = mean_mat_df.groupby(groupby).mean()
        node_sizes = adata.obs[groupby].value_counts().to_dict()

        # 7) stash everything
        self._graph: Dict[str, Any] = {
            "params": {
                "k": k,
                # "thres": thres,
                "min_mean_cut": min_mean_cut,
                "min_drop_cut": min_drop_cut,
                "min_cnt_cut": min_cnt_cut,
            },
            "con": con,
            "con_trim": con_trim,
            "pos": pos,
            "edges_df": edges_df,
            "unique_pairs": unique_pairs,
            "cluster_stats": {"mean": mean_dt, "cnt": cnt_dt, "drop": drop_dt},
            "jumping_genes": jumping_genes,
            "edges_matrix": edges_matrix,
            "mean_mat_df": mean_mat_df,
			"normalized_mean_mat": normalized_mean_mat_df,
			"clusters": clusters,
			"genes": var_names,
            "node_sizes": node_sizes,
        }
        t1 = time.perf_counter()
        if verbose: print(f"[sceleto] build_graph: done in {t1 - t0:.2f}s")

        return self

    # --- helpers ---
    @staticmethod
    def _top_k_edges_symmetric(connectivity_matrix: pd.DataFrame, k: int = 5) -> pd.DataFrame:
        """Keep top-k per row/col, then symmetrize by max."""
        A = np.asarray(connectivity_matrix, dtype=float).copy()
        new_A = np.zeros_like(A)
        # top-k by row
        for i in range(A.shape[0]):
            topk_i = A[i].argsort()[-k:]
            new_A[i, topk_i] = A[i, topk_i]
        # top-k by col
        for j in range(A.shape[1]):
            topk_j = A[:, j].argsort()[-k:]
            new_A[topk_j, j] = A[topk_j, j]
        # symmetrize
        new_A = np.maximum(new_A, new_A.T)
        return pd.DataFrame(new_A, index=connectivity_matrix.index, columns=connectivity_matrix.columns)

    # --- Threshold scanning / set ---
    def plot_threshold(self, thresholds=None, ax=None, target_lines=(1500, 1000)):
        """
        Δ(threshold)별 '살아남는 고유 유전자 수'를 시각화.
        - thresholds: iterable of floats in [0,1]. None이면 np.arange(0,1.0,0.1)
        - target_lines: 수평선 표시값 (e.g., (1500, 1000))
        """
        import numpy as np
        import matplotlib.pyplot as plt
        G = self._graph
        ef = G["edges_matrix"]
        if ef is None or ef.shape[1] == 0:
            raise RuntimeError("edges_matrix is empty. Run build_graph() first.")

        if thresholds is None:
            thresholds = np.round(np.arange(0.0, 1.0, 0.1), 2)

        thr_gene_dict = {}
        for thr in thresholds:
            tmp = ef.copy()
            tmp[tmp <= thr] = 0.0
            survived = (tmp.max(axis=0) > 0).sum()
            thr_gene_dict[float(thr)] = int(survived)

        s = pd.Series(thr_gene_dict).sort_index()

        own = False
        if ax is None:
            own = True
            fig, ax = plt.subplots(figsize=(5,5))
        s.plot(kind="bar", width=0.9, ax=ax)
        for y in (target_lines or []):
            ax.axhline(y=y, linestyle="--", linewidth=1.0, color="red")
        ax.set_xlabel("Δ threshold")
        ax.set_ylabel("# of unique genes")
        ax.set_title("Survived genes per Δ threshold")
        ax.grid(True, axis="y")
        if own:
            plt.tight_layout()
            plt.show()
        return s  # (threshold -> survived_genes) 시리즈 반환

    def set_threshold(self, delta_threshold: float):
        """
        엣지 마커 선별용 Δ 임계값만 저장.
        유전자별 low/high 경계는 get_markers()에서 edge 기반으로 동적으로 계산.
        """
        if not hasattr(self, "_graph"):
            raise RuntimeError("Run build_graph() first.")
        self._graph["threshold"] = float(delta_threshold)
        return self

    # --- Marker aggregation / stats ---
    def get_markers(self, verbose: bool = True):
        import numpy as np
        import pandas as pd
        from tqdm import tqdm
        import time

        t0 = time.perf_counter()
        G = self._graph
        groupby = self.groupby
        thr = G.get("threshold", None)
        if thr is None:
            raise RuntimeError("No Δ threshold set. Call set_threshold(delta_threshold) first.")

        edges_matrix = G["edges_matrix"]       # edge x gene (Δ>0만 담겨 있음; build에서 0 기준)
        con_trim = G["con_trim"]
        pos_df = G["pos"]
        mean_mat_df = G["mean_mat_df"]         # 비정규화 평균
        norm_mean = G.get("normalized_mean_mat", None)

        if verbose:
            n_edges, n_genes = edges_matrix.shape[0], edges_matrix.shape[1]
            n_cells_pos = int((edges_matrix.values > 0).sum())
            print(f"[sceleto] get_markers: start (Δ-thr={thr})")
            print(f"[sceleto]   edges_matrix: edges={n_edges}, genes={n_genes}, positive cells(Δ>0)={n_cells_pos}")

        # Δ>threshold로 melt 필터링
        eg_melt = edges_matrix.reset_index().melt(id_vars="edge", var_name="gene", value_name="value")
        eg_melt = eg_melt[eg_melt["value"] > thr]   # ★ 여기서 확실히 '현재 threshold' 적용
        if verbose:
            print(f"[sceleto]   after Δ>{thr}: pairs={len(eg_melt)} (edge,gene) kept")

        # src/dst 파싱
        eg_melt["src"] = eg_melt["edge"].str.split("->").str[0].str.split("@").str[1]
        eg_melt["dst"] = eg_melt["edge"].str.split("->").str[1].str.split("@").str[1]

        # edge→markers 테이블
        rows = []
        for edge, sub in eg_melt.groupby("edge"):
            kept = sub.sort_values("value", ascending=False)["gene"].tolist()
            src, dst = edge.split("->")
            rows.append((edge, src, dst, ",".join(kept)))
        ef = pd.DataFrame(rows, columns=["edge", "from", "to", "markers"]).set_index("edge")
        ef["n_markers"] = np.where(ef["markers"] == "", 0, ef["markers"].str.split(",").apply(len))
        ef = ef[ef["n_markers"] > 0].copy()
        if verbose:
            print(f"[sceleto]   ef: edges with ≥1 marker = {ef.shape[0]}")

        # (gene,dst)별 max Δ
        if len(eg_melt):
            agg_df = eg_melt.groupby(["gene", "dst"])["value"].max().reset_index()
            mdf = agg_df.pivot(index="dst", columns="gene", values="value").fillna(0.0)
        else:
            mdf = pd.DataFrame(index=[], columns=[], dtype=float)

        # coverage 계산 (단순 스켈레톤)
        mat = self.adata.to_df()
        mat[groupby] = self.adata.obs[groupby].values
        cov_df = (mat.drop(columns=[groupby]) > 0).groupby(self.adata.obs[groupby]).mean()

        # note_df(-1/0/1) + gene 요약
        search_genes = sorted(eg_melt["gene"].unique().tolist()) if len(eg_melt) else []
        note_df = pd.DataFrame(0, index=mean_mat_df.index.astype(str), columns=search_genes, dtype=np.int8)
        recs = {}

        def _levels_by_edges_for_gene(mean_mat_df, eg_sub, gene):
            lows  = sorted(set(eg_sub["src"].tolist()))
            highs = sorted(set(eg_sub["dst"].tolist()))
            if len(highs) == 0:
                return [], [], 0.0, 0.0
            max_for_low  = mean_mat_df.loc[lows,  gene].max() if len(lows)  else mean_mat_df[gene].min()
            min_for_high = mean_mat_df.loc[highs, gene].min()
            low_clusters  = mean_mat_df.index[(mean_mat_df[gene] <= max_for_low)].astype(str).tolist()
            high_clusters = mean_mat_df.index[(mean_mat_df[gene] >= min_for_high)].astype(str).tolist()
            return low_clusters, high_clusters, float(max_for_low), float(min_for_high)

        for g in tqdm(search_genes, disable=not verbose):
            sub = eg_melt[eg_melt["gene"] == g]
            if len(sub) == 0:
                recs[g] = {"n_high": 0, "n_low": 0, "gap": 0.0}
                continue
            low_clusters, high_clusters, max_for_low, min_for_high = _levels_by_edges_for_gene(mean_mat_df, sub, g)
            note_df.loc[note_df.index.isin(low_clusters),  g] = -1
            note_df.loc[note_df.index.isin(high_clusters), g] =  1
            recs[g] = {
                "n_high": len(high_clusters),
                "n_low":  len(low_clusters),
                "gap":    float(min_for_high - max_for_low),
            }

        rf = pd.DataFrame(recs).T

        # celltype별 간단 랭킹
        mdict = {}
        sub_cov = cov_df.reindex(columns=search_genes).copy() if len(search_genes) else pd.DataFrame(index=cov_df.index)
        sub_mdf = mdf.reindex(columns=search_genes).copy() if len(search_genes) else pd.DataFrame(index=mdf.index)
        for celltype in mean_mat_df.index.astype(str):
            marker_genes = note_df.columns[note_df.loc[celltype] > 0].tolist() if len(note_df.columns) else []
            if len(marker_genes) == 0:
                mdict[celltype] = []
                continue
            resdf = pd.concat(
                [rf.loc[marker_genes], sub_cov.loc[celltype, marker_genes], sub_mdf.loc[celltype, marker_genes]],
                axis=1
            )
            resdf.columns = ["n_high", "n_low", "gap", "coverage", "max_delta"]
            resdf['n_high'] = -resdf['n_high']
            resdf["score_mean"] = resdf[["n_high", "n_low", "gap", "coverage", "max_delta"]].mean(axis=1)
            mdict[celltype] = list(resdf.sort_values("score_mean", ascending=False).index)

        # 저장
        G["ef"] = ef
        G["rf"] = rf
        G["note_df"] = note_df
        G["coverage_df"] = cov_df
        G["maxdelta_df"] = mdf
        G["markers_by_celltype"] = mdict

        t1 = time.perf_counter()
        if verbose:
            n_hi = int((note_df.values == 1).sum()) if len(note_df.columns) else 0
            n_lo = int((note_df.values == -1).sum()) if len(note_df.columns) else 0
            print(f"[sceleto]   note_df: genes={len(search_genes)}, highs={n_hi}, lows={n_lo}")
            print(f"[sceleto] get_markers: done in {t1 - t0:.2f}s")

        return {"ef": ef, "rf": rf, "note_df": note_df, "coverage": cov_df, "max_delta": mdf, "markers": mdict}

    # --- Plots: edge / level ---
    def _nx_graph_and_pos(self):
        """내부용: trimmed connectivity로 DiGraph와 좌표 딕셔너리 구성."""
        import networkx as nx
        con_trim = self._graph["con_trim"]
        pos = self._graph["pos"]
        Gnx = nx.DiGraph()
        for i in con_trim.index:
            for j in con_trim.columns:
                if i != j and con_trim.loc[i, j] > 0:
                    Gnx.add_edge(str(i), str(j))
        pos_dict = {str(ix): (float(x), float(y)) for ix, (x, y) in pos[["x", "y"]].iterrows()}
        return Gnx, pos_dict

    def plot_edge(self, gene_name: str, cmap="Reds", highlight_color="red", background_color="lightgrey",
                  figsize=(6,6), title_prefix="Edges with", ax=None):
        """
        특정 gene이 marker로 잡힌 edge를 강조하고, 노드 컬러는 (비정규화) 평균 발현으로 칠함.
        (참조한 plot_gene_marker_edges_with_expression의 클래스 버전)
        """
        import numpy as np
        import matplotlib.pyplot as plt
        import matplotlib.patheffects as path_effects
        import networkx as nx

        if "ef" not in self._graph:
            raise RuntimeError("Run get_markers() first to build edge→markers table (ef).")

        # gene → edges 매핑
        ef = self._graph["ef"]
        gene_to_edges = {}
        for edge, row in ef.iterrows():
            genes = row["markers"].split(",") if row["markers"] else []
            for g in genes:
                if g:
                    gene_to_edges.setdefault(g, []).append(edge)

        Gnx, pos_dict = self._nx_graph_and_pos()
        mean_mat = self._graph["mean_mat_df"]
        node_sizes_dict = self._graph.get("node_sizes", {})
        node_sizes = [np.sqrt(node_sizes_dict.get(n, 1)) * 10 for n in Gnx.nodes()]

        own = False
        if ax is None:
            own = True
            fig, ax = plt.subplots(figsize=figsize)

        # 강조 엣지
        hl_edges = [(a.split("@")[-1], b.split("@")[-1]) for a, b in (e.split("->") for e in gene_to_edges.get(gene_name, []))]

        # edge 스타일
        for (u, v) in Gnx.edges():
            is_hi = (u, v) in hl_edges
            nx.draw_networkx_edges(
                Gnx,
                pos=pos_dict,
                edgelist=[(u, v)],
                edge_color=highlight_color if is_hi else background_color,
                alpha=1.0 if is_hi else 0.5,
                arrows=True,
                arrowsize=12,
                width=1.2,
                connectionstyle="arc3,rad=0.05",
                min_target_margin=10,
                ax=ax,
            )

        # 노드 컬러는 평균 발현
        node_vals = []
        for n in Gnx.nodes():
            try:
                node_vals.append(float(mean_mat.loc[n, gene_name]))
            except Exception:
                node_vals.append(0.0)

        nodes = nx.draw_networkx_nodes(
            Gnx, pos=pos_dict, node_color=node_vals, node_size=node_sizes, cmap=cmap,
            vmin=min(node_vals), vmax=max(node_vals) if max(node_vals)>0 else 1.0, ax=ax
        )

        for node, (x, y) in pos_dict.items():
            txt = ax.text(x, y, node, fontsize=8, ha="center", va="center", color="black")
            txt.set_path_effects([path_effects.Stroke(linewidth=2, foreground="white"), path_effects.Normal()])

        ax.set_title(f"{title_prefix} '{gene_name}'")
        ax.axis("off")
        cbar = ax.figure.colorbar(nodes, ax=ax, shrink=0.7, pad=0.02)
        cbar.set_label(f"{gene_name} mean expression", fontsize=9)
        if own:
            plt.tight_layout()
            plt.show()
        return ax

    def plot_level(self, gene_name: str, figsize=(6,6), category_colors=None,
                   highlight_color="red", background_edge_color="lightgrey", edge_width=1.4,
                   show_labels=True, title_prefix="Levels for", ax=None):
        """
        gene의 low/intermediate/high(-1/0/1)를 노드 색으로, 해당 gene이 잡힌 edge를 강조해서 표시.
        (참조한 plot_gene_levels_with_edges의 클래스 버전)
        """
        import matplotlib.pyplot as plt
        import matplotlib.patheffects as path_effects
        import networkx as nx
        from matplotlib.patches import Patch

        if "note_df" not in self._graph:
            raise RuntimeError("Run get_markers() first (builds note_df).")
        note_df = self._graph["note_df"]

        # gene → edges
        ef = self._graph["ef"]
        gene_to_edges = {}
        for edge, row in ef.iterrows():
            genes = row["markers"].split(",") if row["markers"] else []
            for g in genes:
                if g:
                    gene_to_edges.setdefault(g, []).append(edge)

        Gnx, pos_dict = self._nx_graph_and_pos()
        node_sizes_dict = self._graph.get("node_sizes", {})
        if category_colors is None:
            category_colors = {-1: "#3b82f6", 0: "#d1d5db", 1: "#ef4444"}

        own = False
        if ax is None:
            own = True
            fig, ax = plt.subplots(figsize=figsize)

        # 강조 엣지
        hl_edges = {(a.split("@")[-1], b.split("@")[-1]) for a, b in (e.split("->") for e in gene_to_edges.get(gene_name, []))}
        # 노드 색
        levels = {n: int(note_df.loc[n, gene_name]) if (n in note_df.index and gene_name in note_df.columns) else 0 for n in Gnx.nodes()}
        node_colors = [category_colors.get(levels[n], category_colors[0]) for n in Gnx.nodes()]
        node_sizes = [np.sqrt(node_sizes_dict.get(n, 1)) * 10 for n in Gnx.nodes()]

        nx.draw_networkx_nodes(Gnx, pos=pos_dict, node_color=node_colors, node_size=node_sizes, linewidths=0, edgecolors="white", ax=ax)

        if show_labels:
            for node, (x, y) in pos_dict.items():
                txt = ax.text(x, y, node, fontsize=8, ha="center", va="center", color="black")
                txt.set_path_effects([path_effects.Stroke(linewidth=2, foreground="white"), path_effects.Normal()])

        for (u, v) in Gnx.edges():
            is_hi = (u, v) in hl_edges
            nx.draw_networkx_edges(
                Gnx, pos=pos_dict, edgelist=[(u, v)],
                edge_color=highlight_color if is_hi else background_edge_color,
                alpha=1.0 if is_hi else 0.5,
                arrows=True, arrowsize=12, width=edge_width,
                connectionstyle="arc3,rad=0.05", min_target_margin=10, ax=ax
            )

        ax.set_title(f"{title_prefix} '{gene_name}'")
        ax.axis("off")
        legend_handles = [Patch(color=category_colors[-1], label="Low (-1)"),
                          Patch(color=category_colors[0],  label="Intermediate (0)"),
                          Patch(color=category_colors[1],  label="High (1)")]
        ax.legend(handles=legend_handles, loc="upper left", bbox_to_anchor=(1.02, 1))
        if own:
            plt.tight_layout()
            plt.show()
        return ax

