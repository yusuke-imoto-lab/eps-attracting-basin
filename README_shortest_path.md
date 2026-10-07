# Shortest-path revision for `epsbasin`

This revision **does not replace** the existing `epsbasin.py`.  It adds a new
module, `epsbasin/shortest_path.py`, and exports its functions from
`epsbasin/__init__.py`.

## Assumption

For every sequence (`adata.obs['seq_id']`), the terminal observation is the last
observation in that sequence.  Every terminal observation must be labelled
`good` or `bad` in `adata.obs['cluster']` (the labels can be changed through
function arguments).  Intermediate observations may keep labels such as
`other`.

The directed graph has one vertex per observation:

- same sequence, forward direction: weight `0`;
- same sequence, backward direction: weight `inf`;
- different sequences: pairwise cost (`adata.uns['cost_matrix']`, or Euclidean
  distance computed from `adata.X`).

By default all forward pairs in a sequence have zero weight.  Set
`sequence_edge_mode='adjacent'` to keep only observed one-step zero-cost edges.

## Two path costs

### Ordinary epsilon

`eps_attracting_basin_shortest_path` uses the minimax (bottleneck) path value

`min_path max(edge_cost)`.

This corresponds to the smallest uniform per-step control bound.

### epsilon_Sigma

`eps_sum_attracting_basin_shortest_path` uses the ordinary additive path value

`min_path sum(edge_cost)`.

This corresponds to the smallest total control cost.

## Signed debut functions from the Good/Bad terminal partition

Let `d_G(y)` and `d_B(y)` be the graph distances from observation `y` to the
Good and Bad terminal sets, using either the bottleneck or additive path cost.
The uncontrolled terminal label of an observation is the label of the terminal
state in its own sequence.

If the natural terminal is Good:

- debut to Good = `-d_B(y)`;
- debut to Bad  = `+d_B(y)`.

If the natural terminal is Bad:

- debut to Good = `+d_G(y)`;
- debut to Bad  = `-d_G(y)`.

The landscape column is the pointwise minimum of the two debut functions,
matching the current `plot_landscape` implementation.

## Usage

```python
import epsbasin

# epsilon (maximum edge / bottleneck)
epsbasin.eps_attracting_basin_shortest_path(
    adata,
    output_key="eps_attracting_basin_sp",
)

# epsilon_Sigma (sum of edge costs)
epsbasin.eps_sum_attracting_basin_shortest_path(
    adata,
    output_key="eps_sum_attracting_basin_sp",
)

# Existing plotting code can be reused.
epsbasin.plot_debut(
    adata,
    eps_key="eps_attracting_basin_sp",
    target_cluster_key="good",
)

epsbasin.plot_landscape(
    adata,
    eps_key="eps_attracting_basin_sp",
)
```

If the terminal labels are `Good` and `Bad` rather than lowercase:

```python
epsbasin.eps_attracting_basin_shortest_path(
    adata,
    good_cluster_key="Good",
    bad_cluster_key="Bad",
)
```

If within-sequence time order is not the current row order, provide the column
used for ordering:

```python
epsbasin.eps_attracting_basin_shortest_path(
    adata,
    time_key="time",
)
```

## Applying to the current repository

Copy `epsbasin/shortest_path.py` into the repository and add

```python
from .shortest_path import *
```

to `epsbasin/__init__.py`, or apply the included patch from the repository root:

```bash
git apply shortest_path_revision.patch
```
