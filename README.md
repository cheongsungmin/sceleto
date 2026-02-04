# 🏗️ sceleto

This repository is currently focused on updating the marker functionality in the sceleto package (https://github.com/scmgl-kaist/sceleto). Only the marker-related components are implemented here for now; additional modules will be integrated later.

## Installation

- Since this is a development version, creating a fresh environment is recommended.

```bash
pip install --user git+https://github.com/cheongsungmin/sceleto.git
```

## Graph-based marker API

### Prerequisites

- Your AnnData object should have log1p-transformed expression values stored in adata.raw.
- To use the graph-based marker workflow, you must provide:

1. A PAGA graph stored in `adata.uns["paga"]`. (We recommend running PAGA separately and saving it into `adata`.)
2. Node position coordinates (`pos`) for plotting. Running `sc.pl.paga_compare(adata)` will store cluster node positions in `adata.uns["paga"]["pos"]`.

### Usage

```python
import sceleto as scl
```

Choose the options depending on your goal. The defaults are `hierarchical_markers=False` and `specific_markers=True`
```python
MG = scl.markers.marker(
    adata,
    groupby="leiden",            # cluster labels stored in adata.obs
    k=5,                         # keep top-k neighbors per node in the PAGA graph
    thres_fc=3.0,                # genes with FC >= 3.0 are treated as marker candidates
)
```

```python
# Marker genes are stored in a dictionary format
MG.specific_marker_log

show_genes = []
for g in MG.ctx.groups:
    show_genes += MG.specific_marker_log[str(g)][:5]
    
sc.pl.dotplot(
    adata,
    show_genes,
    groupby=leiden,
    standard_scale='var',
    dot_max=0.5
)
```

**Recovering hierarchical marker genes**

You can trace marker genes across multiple clustering resolutions to build a hierarchical lineage.

> Note: Currently, this feature supports exactly three levels of grouping. You must compute the markers for each resolution (e.g., MG0, MG1, MG2) before building the hierarchy.

```python
# 1. Compute markers for different resolutions
MG0 = scl.markers.marker(adata, groupby="leiden_1.0")
MG1 = scl.markers.marker(adata, groupby="leiden_2.0")
MG2 = scl.markers.marker(adata, groupby="leiden_4.0")

# 2. Build the hierarchy
res = scl.markers.hierarchy(
    adata, 
    [MG0, MG1, MG2], 
    min_cells_for_path=500, 
    n_top_markers=50
)
```

The hierarchical structure is stored in `res.icls` and `res.path`. When printed, it displays the marker tree as follows:
```
Hierarchical Marker Tree
==============================
├── leiden_0.5@0 :: Markers: ACTN1, FHIT, AIF1, BEX3, ARMH1
    ├── leiden_1.0@0 :: Markers: TSHZ2, FHIT, SYPL1, SOCS3, TRAC
        ├── leiden_1.5@0 (icls [0]) :: Markers: FHIT, ACTN1, TSHZ2, AIF1, OXNAD1
        ├── leiden_1.5@3 (icls [1]) :: Markers: NSG1, CTSH, CISH, TIMP1, PTGER2
    ...
```

### Graph visualization

```python
MG.plot_gene_edges_fc("CD3D", figsize=(7, 5))
MG.plot_gene_levels_with_edges("CD3D", figsize=(7, 5))
```

## Note

- The legacy marker workflow (`sceleto.markers.marker`) is still available via `scl.markers.classic()`:
```python
import sceleto as scl
MC = scl.markers.classic(adata, 'leiden')
MC.plot_marker()  # Plotting dot plot
MC.mks            # Dictionary of marker genes for each cluster
```
- The legacy `sceleto.us` function is also available now.

- The function runs PAGA when no PAGA information is available. Otherwise, it uses the precomputed PAGA results. Ensure that PAGA information for the group of interest is present before computing markers.

## Dependencies

```
| Package    | Version |
| ---------- | ------- |
| sceleto    | 0.0.1   |
| scanpy     | 1.11.5  |
| matplotlib | 3.10.7  |
| numpy      | 2.3.4   |
| pandas     | 2.3.3   |
| seaborn    | 0.13.2  |

| Dependency         | Version     |
| ------------------ | ----------- |
| six                | 1.17.0      |
| psutil             | 7.2.1       |
| PyYAML             | 6.0.3       |
| scikit-learn       | 1.7.2       |
| parso              | 0.8.5       |
| kiwisolver         | 1.4.9       |
| natsort            | 8.4.0       |
| h5py               | 3.15.1      |
| decorator          | 5.2.1       |
| cycler             | 0.12.1      |
| joblib             | 1.5.2       |
| pure_eval          | 0.2.3       |
| igraph             | 1.0.0       |
| pillow             | 12.0.0      |
| matplotlib-inline  | 0.2.1       |
| anndata            | 0.12.6      |
| statsmodels        | 0.14.5      |
| donfig             | 0.8.1.post1 |
| texttable          | 1.7.0       |
| wcwidth            | 0.2.14      |
| prompt_toolkit     | 3.0.52      |
| crc32c             | 2.8         |
| typing_extensions  | 4.15.0      |
| numcodecs          | 0.16.3      |
| platformdirs       | 4.5.1       |
| comm               | 0.2.3       |
| patsy              | 1.0.2       |
| jupyter_client     | 8.8.0       |
| asttokens          | 3.0.1       |
| tornado            | 6.5.4       |
| executing          | 2.2.1       |
| traitlets          | 5.14.3      |
| numba              | 0.62.1      |
| packaging          | 25.0        |
| scipy              | 1.16.3      |
| jupyter_core       | 5.9.1       |
| zarr               | 3.1.3       |
| charset-normalizer | 3.4.4       |
| python-dateutil    | 2.9.0.post0 |
| threadpoolctl      | 3.6.0       |
| ipykernel          | 7.1.0       |
| jedi               | 0.19.2      |
| pyparsing          | 3.2.5       |
| ipython            | 9.9.0       |
| setuptools         | 80.9.0      |
| pyzmq              | 27.1.0      |
| networkx           | 3.5         |
| Pygments           | 2.19.2      |
| legacy-api-wrap    | 1.5         |
| llvmlite           | 0.45.1      |
| pytz               | 2025.2      |
| leidenalg          | 0.11.0      |
| stack-data         | 0.6.3       |
| session-info2      | 0.2.3       |
| debugpy            | 1.8.19      |

| Component | Info                                                                 |
| --------- | -------------------------------------------------------------------- |
| Python    | 3.11.13 (main, Jun  5 2025, 13:12:00) [GCC 11.2.0]                   |
| OS        | Linux-5.4.0-150-generic-x86_64-with-glibc2.27                        |
| CPU       | 96 logical CPU cores, x86_64                                         |
| GPU       | ID: 0, NVIDIA GeForce RTX 3090, Driver: 530.41.03, Memory: 24576 MiB |
| Updated   | 2026-01-12 05:58                                                     |
```

