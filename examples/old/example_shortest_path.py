"""Minimal usage example for the terminal-state shortest-path implementation."""

import anndata
import numpy as np
import pandas as pd

import epsbasin


adata = anndata.AnnData(
    X=np.array(
        [
            [1.0, 6.0],
            [5.0, 6.0],
            [9.0, 9.0],
            [1.0, 5.0],
            [5.0, 5.0],
            [9.0, 2.0],
        ]
    ),
    obs=pd.DataFrame(
        {
            "seq_id": [0, 0, 0, 1, 1, 1],
            "cluster": ["other", "other", "good", "other", "other", "bad"],
        }
    ),
)

# Ordinary epsilon: minimax / bottleneck shortest path.
epsbasin.eps_attracting_basin_shortest_path(
    adata,
    output_key="eps_attracting_basin_sp",
    distance_key="bottleneck_distance",
    terminal_class_key="terminal_class",
)

# epsilon_Sigma: additive shortest path.
epsbasin.eps_sum_attracting_basin_shortest_path(
    adata,
    output_key="eps_sum_attracting_basin_sp",
    distance_key="sum_distance",
)

print(
    adata.obs[
        [
            "seq_id",
            "cluster",
            "terminal_class",
            "eps_attracting_basin_sp_good",
            "eps_attracting_basin_sp_bad",
            "eps_attracting_basin_sp_landscape",
            "eps_sum_attracting_basin_sp_good",
            "eps_sum_attracting_basin_sp_bad",
            "eps_sum_attracting_basin_sp_landscape",
        ]
    ]
)

# Existing plotting functions can be reused by setting eps_key to the new prefix:
# epsbasin.plot_debut(adata, eps_key="eps_attracting_basin_sp", target_cluster_key="good")
# epsbasin.plot_landscape(adata, eps_key="eps_attracting_basin_sp")
