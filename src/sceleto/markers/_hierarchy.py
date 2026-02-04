from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List, Sequence, Tuple, Optional

import numpy as np
import pandas as pd


@dataclass
class HierarchyRun:
    # Inputs / meta
    levels: List[str]
    params: Dict[str, Any]

    # Key artifacts you already build
    icls: pd.Series
    path: pd.Series
    icls_full_dict: Dict[str, str]
    icls_path_df: pd.DataFrame
    marker_rank_df: pd.DataFrame
    icls_gene_presence_df: pd.DataFrame
    gene_freq_df: pd.DataFrame
    score_df: pd.DataFrame

    # Tree artifacts
    tree_root: Dict[str, Any]
    icls_to_path: Dict[int, List[str]]

    # Output from your printing traversal (markers per branching)
    markers: Any


def hierarchy(
    adata: Any,
    markers_list: Sequence[Any],
    *,
    min_cells_for_path: int = 500,
    n_top_markers: int = 10,
) -> HierarchyRun:
    """
    Run the user's hierarchy pipeline as-is and return all intermediate tables.

    Notes
    -----
    - This function intentionally preserves the user's original logic.
    - It writes `adata.obs["path"]` and `adata.obs["icls"]` as side effects (same as the script).
    """
    # ----------------------------
    # (User code) 그대로 시작
    # ----------------------------
    markers_list = list(markers_list)

    g0, g1, g2 = [marker_output.ctx.groupby for marker_output in markers_list]

    # 1) Ensure categorical dtype and preserve desired category order
    #    (If already categorical with correct order, this is harmless.)
    for g in (g0, g1, g2):
        if not pd.api.types.is_categorical_dtype(adata.obs[g]):
            adata.obs[g] = adata.obs[g].astype("category")
        # Keep current categories order explicitly (optional but makes intent clear)
        adata.obs[g] = adata.obs[g].cat.set_categories(adata.obs[g].cat.categories, ordered=True)

    # 2) Build per-cell path strings (cell order is preserved automatically)
    adata.obs["path"] = (
        f"{g0}@" + adata.obs[g0].astype(str)
        + "|" + f"{g1}@" + adata.obs[g1].astype(str)
        + "|" + f"{g2}@" + adata.obs[g2].astype(str)
    )

    # 3) (Optional, but usually what people mean by "keep the order" after concatenation)
    #    Make `path` categorical with categories ordered by the cartesian product of cat.categories
    path_categories = [
        f"{g0}@{x}|{g1}@{y}|{g2}@{z}"
        for x in adata.obs[g0].cat.categories
        for y in adata.obs[g1].cat.categories
        for z in adata.obs[g2].cat.categories
    ]

    adata.obs["path"] = pd.Categorical(adata.obs["path"], categories=path_categories, ordered=True)

    # excluding small
    path_idx = adata.obs["path"].value_counts() < min_cells_for_path

    # Compute small paths (keep your logic)
    small_path_idx = adata.obs["path"].value_counts().loc[lambda s: s < min_cells_for_path].index

    # Mark small as missing (do NOT recast)
    mask_small = adata.obs["path"].isin(small_path_idx)
    adata.obs.loc[mask_small, "path"] = pd.NA

    # Rebuild categories in the original intended order
    # - Keep only categories that still appear (non-missing)
    present = set(adata.obs["path"].dropna().unique())

    # Keep order defined by the original cartesian-product list
    new_categories = [c for c in path_categories if c in present]

    # Apply categories explicitly (this preserves your intended order)
    adata.obs["path"] = adata.obs["path"].cat.set_categories(new_categories, ordered=True)

    icls_full_dict: Dict[str, str] = {}
    for i, path in enumerate(adata.obs["path"].cat.categories):
        icls_full_dict[str(i)] = path

    path_to_key = {v: k for k, v in icls_full_dict.items()}
    adata.obs["icls"] = adata.obs["path"].map(path_to_key).astype("string")

    df_icls_path = pd.DataFrame(pd.Series(icls_full_dict), columns=["icls_full"])
    df_icls_path[g0] = [x.split("|")[0] for x in df_icls_path["icls_full"]]
    df_icls_path[g1] = [x.split("|")[1] for x in df_icls_path["icls_full"]]
    df_icls_path[g2] = [x.split("|")[2] for x in df_icls_path["icls_full"]]
    df_icls_path["root"] = df_icls_path[g0]
    df_icls_path = df_icls_path.reset_index(names="icls")

    rows: List[List[Any]] = []
    for k, v in markers_list[0].specific_marker_log.items():
        for i, gene in enumerate(v[:n_top_markers]):
            rows.append([f"{g0}", f"{g0}@{k}", gene, i + 1])

    for k, v in markers_list[1].specific_marker_log.items():
        for i, gene in enumerate(v[:n_top_markers]):
            rows.append([f"{g1}", f"{g1}@{k}", gene, i + 1])

    for k, v in markers_list[2].specific_marker_log.items():
        for i, gene in enumerate(v[:n_top_markers]):
            rows.append([f"{g2}", f"{g2}@{k}", gene, i + 1])

    df_marker_rank = pd.DataFrame(rows, columns=["resolution", "leiden", "gene", "rank"])

    temp: List[pd.DataFrame] = []

    for k, v in icls_full_dict.items():
        l0 = v.split("|")[0]
        l1 = v.split("|")[1]
        l2 = v.split("|")[2]

        # (A, rank) Pivot to gene x resolution table
        piv = df_marker_rank[df_marker_rank["leiden"].isin([l0, l1, l2])].pivot(
            index="gene", columns="resolution", values="rank"
        )
        # (B, binary) Convert presence/absence to 1/0
        df_binary = piv.notna().astype("int8")  # or "int"

        # merge (A) and (B)
        df = pd.merge(piv, df_binary, left_index=True, right_index=True, how="left")
        df.columns = ["rank_0", "rank_1", "rank_2", "present_0", "present_1", "present_2"]
        df["n_levels"] = df["present_0"] + df["present_1"] + df["present_2"]
        df = pd.merge(
            pd.Series([k for _ in range(df.shape[0])], name="icls"),
            df.reset_index(),
            how="left",
            left_index=True,
            right_index=True,
        )
        df = df.sort_values("n_levels", ascending=False)

        temp.append(df)

    icls_gene_presence = pd.concat(temp, axis=0).reset_index(drop=True)

    # --- Robustly define "gene is present in this icls" ---
    # If present_* columns exist, use them; otherwise assume each (icls, gene) row is already a union member.
    present_cols = [c for c in icls_gene_presence.columns if c.startswith("present_")]
    if present_cols:
        present_any = icls_gene_presence[present_cols].fillna(False).astype(bool).any(axis=1)
        df_use = icls_gene_presence.loc[present_any, ["icls", "gene"]].copy()
    else:
        df_use = icls_gene_presence.loc[:, ["icls", "gene"]].copy()

    # --- Global DF/IDF across icls (icls-level document frequency) ---
    N_icls = df_use["icls"].nunique()

    gene_freq = (
        df_use.groupby("gene")["icls"]
        .nunique()
        .rename("df_global_icls")
        .reset_index()
    )

    gene_freq["frac_icls"] = gene_freq["df_global_icls"] / N_icls

    # Smooth IDF (always >= 1)
    gene_freq["idf_global_icls"] = np.log((N_icls + 1) / (gene_freq["df_global_icls"] + 1)) + 1.0

    # Optional: index by gene and sort
    gene_freq = (
        gene_freq.set_index("gene")
        .sort_values(["df_global_icls", "idf_global_icls"], ascending=[False, True])
    )

    gene_freq.columns = ["n_icls", "frac_icls", "idf_icls"]

    score_df = pd.merge(icls_gene_presence, gene_freq.reset_index(), how="left", on="gene")
    score_df["icls"] = score_df["icls"].astype("int")

    # Example input: list of "path" strings
    paths_str = icls_full_dict.values()

    # Parse paths
    paths = [s.split("|") for s in paths_str]

    icls_tree: Dict[str, Any] = {}
    for p in paths:
        cur = icls_tree
        for node in p:
            cur = cur.setdefault(node, {})

    # Hierarchy Dictionary
    hierarchy_map = icls_full_dict

    # 2. Build the Tree Structure
    # Tree node structure: { 'name': 'leiden_X.0@Y', 'children': {}, 'icls_indices': [] }
    tree_root: Dict[str, Any] = {}

    # icls to leaf node map (나중에 데이터 조회용)
    icls_to_path: Dict[int, List[str]] = {}

    for icls_idx, path_str in hierarchy_map.items():
        parts = path_str.split("|")
        current_level = tree_root

        path_list: List[str] = []
        for part in parts:
            path_list.append(part)
            if part not in current_level:
                current_level[part] = {"children": {}, "icls_indices": [], "level_name": part}

            # 해당 노드 하위에 속하는 모든 icls index 저장
            current_level[part]["icls_indices"].append(int(icls_idx))
            current_level = current_level[part]["children"]

        icls_to_path[int(icls_idx)] = path_list

    # 4. Scoring Logic (The Core)
    def find_branching_markers(parent_node_name, children_dict, score_df, target_level_suffix):
        """
        Sibling Competition을 수행하여 각 자식 노드를 대표하는 마커 선정
        target_level_suffix: '2.0' or '4.0' (비교에 사용할 데이터 컬럼 접미사)
        """
        if not children_dict:
            return {}

        rank_col = f"rank_{target_level_suffix}"
        present_col = f"present_{target_level_suffix}"

        # 1. 모든 자식 노드의 데이터 준비
        children_stats = {}
        all_genes = set()

        for child_name, child_node in children_dict.items():
            stats = get_node_stats(child_node["icls_indices"], score_df, rank_col, present_col)
            children_stats[child_name] = stats
            all_genes.update(stats.index.tolist())

        results = {}

        # 2. 각 자식 노드별로 스코어 계산
        for target_child, target_stats in children_stats.items():
            scores = []

            # 형제 노드들의 리스트
            siblings = [name for name in children_stats if name != target_child]

            for gene in target_stats.index:
                # A. Basic Info
                my_present = target_stats.loc[gene, present_col]
                if my_present == 0:
                    continue  # 내가 안 가지고 있으면 마커 후보 탈락 (Positive Marker Only)

                my_rank = target_stats.loc[gene, rank_col]
                idf = target_stats.loc[gene, "idf_icls"]

                # B. Exclusivity Calculation
                sibling_present_vals = []
                for sib in siblings:
                    if gene in children_stats[sib].index:
                        sibling_present_vals.append(children_stats[sib].loc[gene, present_col])
                    else:
                        sibling_present_vals.append(0)

                # 형제들의 평균 발현 여부 (0~1)
                avg_sibling_present = np.mean(sibling_present_vals) if siblings else 0

                # Exclusivity: 나한테는 있는데(1), 남들은 없을수록(0) -> 1.0에 가까워짐
                exclusivity = my_present - avg_sibling_present

                if exclusivity <= 0.1:
                    continue  # 변별력이 너무 낮으면 스킵

                # C. Intensity (Rank Score)
                # Rank 1 -> 1.0, Rank 10 -> 0.1, Rank 100 -> 0.01
                rank_score = 1.0 / (my_rank + 1.0)

                # D. Priority Boost
                # removed

                # E. Final Score
                # idf는 보통 log scale이므로 그대로 곱하거나, 너무 흔한걸 죽이기 위해 사용
                final_score = exclusivity * (1 + rank_score) * np.log1p(idf)

                scores.append((gene, final_score, my_rank, exclusivity))

            # Sort by score desc
            scores.sort(key=lambda x: x[1], reverse=True)
            results[target_child] = scores[:5]  # Top 5 markers

        return results

    # 3. Helper Function: Get aggregated stats for a node
    def get_node_stats(icls_indices, score_df, target_rank_col, target_present_col):
        """
        해당 노드(여러 icls의 집합)에 대한 유전자 통계를 구함.
        여기서는 해당 노드에 속한 icls 중 하나라도 present면 1, rank는 min(best) rank를 사용.
        """
        # Filter DF for relevant icls
        subset = score_df[score_df["icls"].isin(icls_indices)].copy()

        # rank 컬럼이 NaN이면 매우 큰 값으로 대체
        subset[target_rank_col] = subset[target_rank_col].fillna(100)
        subset[target_present_col] = subset[target_present_col].fillna(0)

        stats = subset.groupby("gene").agg(
            {
                target_rank_col: "min",
                target_present_col: "max",  # max: any member has it
                "idf_icls": "first",
            }
        )
        return stats

    # 5. Recursive Execution & Printing
    def print_tree(node, level_prefix=g0, depth=0):
        indent = "    " * depth

        if depth == 0:  # Root (Virtual or 1.0 container)
            print("Hierarchical Marker Tree")
            print("=" * 30)

        children = node
        if not children:
            return

        # Determine stats level based on child names
        first_child = next(iter(children))
        if g1 in first_child:
            target_suffix = "1"
        elif g2 in first_child:
            target_suffix = "2"
        else:
            target_suffix = "0"  # Fallback

        # Calculate Markers for this branching
        markers = find_branching_markers("parent", children, score_df, target_suffix)

        for child_name, child_node in children.items():
            # Format Markers
            marker_str = ""
            if child_name in markers:
                top_genes = [f"{m[0]}" for m in markers[child_name]]  # (Gene, Score, Rank, Excl)
                marker_str = f" :: Markers: {', '.join(top_genes)}"

            # ICLS info (if leaf)
            icls_info = ""
            if not child_node["children"]:
                icls_list = child_node["icls_indices"]
                icls_info = f" (icls {icls_list})"

            print(f"{indent}├── {child_name}{icls_info}{marker_str}")

            # Recursive Call
            if child_node["children"]:
                print_tree(child_node["children"], level_prefix=child_name, depth=depth + 1)

        return markers

    # Execute
    markers = print_tree(tree_root)

    # ----------------------------
    # (User code) 그대로 끝
    # ----------------------------

    return HierarchyRun(
        levels=[str(g0), str(g1), str(g2)],
        params={"min_cells_for_path": int(min_cells_for_path), "n_top_markers": int(n_top_markers)},
        icls=adata.obs["icls"],
        path=adata.obs["path"],
        icls_full_dict=icls_full_dict,
        icls_path_df=df_icls_path,
        marker_rank_df=df_marker_rank,
        icls_gene_presence_df=icls_gene_presence,
        gene_freq_df=gene_freq.reset_index(),
        score_df=score_df,
        tree_root=tree_root,
        icls_to_path=icls_to_path,
        markers=markers,
    )
