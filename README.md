# sceleto (version 3)

- 아래 패키지 버전에서 개발 중입니다.
```python
scanpy==1.9.6 anndata==0.10.3 umap==0.5.5 numpy==1.26.4 scipy==1.11.4 pandas==2.1.3 scikit-learn==1.3.2 statsmodels==0.14.0 pynndescent==0.5.11
```

- 설치 방법
- 개발 버전이니까 새로운 환경을 만드시는 것을 추천드립니다.
- sceleto라는 폴더가 이미 있다면 덮어써집니다.
- `git clone https://github.com/cheongsungmin/sceleto.git`
- `cd sceleto`
- `pip install .`

- adata.raw에는 log1p-transformed value가 저장되어있어야 합니다.
- 기존 버전의 sceleto.markers.marker는 아래와 같이 `scl.markers.classic()` 함수로 사용 가능합니다.
```python
import sceleto as scl
M = scl.markers.classic(adata, 'Level2')
M.plot_marker()
```

- graph 버전을 사용하시려면 PAGA graph가 있어야 합니다. 현재는 별도로 수행해서 넣는 것을 추천드립니다. `adata.uns['paga']` 가 존재해야 합니다.
```python
MG = scl.markers.graph(adata, 'Level2')
MG.build_graph(k=5, use_existing_paga=True, verbose=True)
MG.plot_threshold()
MG.set_threshold(0.6)
out = MG.get_markers()
```

- visualization
```python
MG.plot_edge("CD3D", figsize=(6,5))
MG.plot_level("CD8A", figsize=(7,5))
```

- dot plot
```python
show_genes = []
for v in out['markers'].values():
    show_genes += v[:5][::-1]

sc.pl.dotplot(adata, show_genes, groupby='Level2', standard_scale='var')
```
