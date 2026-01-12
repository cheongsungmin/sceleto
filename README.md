# sceleto

This repository is currently focused on updating the marker functionality in the sceleto package. Only the marker-related components are implemented here for now; additional modules will be integrated later.

## Installation

- Since this is a development version, creating a fresh environment is recommended.
- If you already have a local folder named sceleto, cloning will overwrite it — please be careful.
- `git clone https://github.com/cheongsungmin/sceleto.git`
- `cd sceleto`
- `pip install .`


## Graph-based marker API (`sceleto.markers.graph`)

### Prerequisites

- Your AnnData object should have log1p-transformed expression values stored in adata.raw.

- To use the graph-based marker workflow, you must provide:

1. A **PAGA graph** stored in `adata.uns["paga"]`
   (We recommend running PAGA separately and saving it into `adata`.)
2. Node position coordinates (`pos`) for plotting
   Running `sc.pl.paga_compare(adata)` will store cluster node positions in `adata.uns["paga"]["pos"]`.

### Usage

```python
MG = mg.run_marker_graph(
    adata,
    groupby=leiden,
    thres_fc=3.0,
)

MG.marker_log  # Dictionary of marker genes for each cluster
```

### Visualization

```python
MG.plot_gene_edges_fc("CD3D", figsize=(7, 5))
MG.plot_gene_levels_with_edges("CD3D", figsize=(7, 5))
```

## Note

- The legacy marker workflow (sceleto.markers.marker) is still available via scl.markers.classic():
```python
import sceleto as scl
M = scl.markers.classic(adata, 'Level2')
M.plot_marker()
```

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

