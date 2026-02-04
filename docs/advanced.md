### Advanced: Specific marker scoring (presets & custom)

By default, `scl.markers.marker` ranks specific markers using an internal scoring function that combines local evidence (in-vs-out specificity, effect size, within-cluster coverage) and global evidence (how decisively the gene separates clusters across the graph).
Advanced users can either (1) switch among built-in presets, or (2) provide a fully custom scoring function.

**Option A) Built-in score presets (1–4)**

You can switch the specific-marker scoring behavior via `specific_score_preset`.

* **Preset 1**: Smooth / forgiving around borderline genes (less brittle thresholds).
* **Preset 2**: More conservative specificity (compares against the “worst” out-group; stricter).
* **Preset 3**: Simpler global term (more interpretable, lighter regularization).
* **Preset 4**: Most conservative + simple global term (strict and stable).

```python
import sceleto as scl

MG = scl.markers.marker(
    adata,
    groupby="leiden",
    specific_score_preset=2,   # 1, 2, 3, or 4
)
```

**Option B) Custom scoring function (advanced users)**

You can also provide your own scoring function via `specific_score_fn`.
If you’re unsure what columns are available, run the default pipeline once and inspect
`MG.specific_ranking_df` (it contains the computed per-gene features used for ranking).

```python
import sceleto as scl

def my_score(df):
    # Example: reward strong edge separation and in-cluster coverage,
    # while penalizing genes that are high in many clusters.
    return (df["max_delta"] * df["coverage_one"]) / (1.0 + df["n_high"])

MG = scl.markers.marker(
    adata,
    groupby="leiden",
    specific_score_fn=my_score,
)

# Explore available features/columns:
MG.specific_ranking_df.head()
```