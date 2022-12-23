# doppelganger

An R package to analyze correlations and select a subset of non-redundant variables.

For details, see the post on Towards Data Science:

[Fighting doppelgängers: How to rid data of evil twins reducing the feature space](https://medium.com/towards-data-science/fighting-doppelg%C3%A4ngers-2fc28762e169)

### Quick usage

```
dg <- doppelganger(data, priority={...}, threshold={...})
```

Argument | Description
:--- | :---
data | Data frame containing numerical variables.
variables | Columns of `data` to consider in the analysis (default: all).
priority | Ranking method to prioritize variables (`"centrality"`, `"peripherality"`, `"raw_order"`).
threshold | Correlation cut-off (absolute value).

### Output

Variabiles to keep / drop:

```
dg$keep
dg$drop
```