import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.tri import Triangulation
from sklearn.metrics import pairwise_distances
import matplotlib.patches as mpatches
from scipy.spatial import Delaunay
from collections import OrderedDict
from matplotlib.lines import Line2D
import os

def eps_attracting_basin(
    adata, 
    cluster_key="cluster", 
    target_cluster_key="target_cluster", 
    cost_matrix_key="cost_matrix", 
    output_key="eps_attracting_basin",
    epsilon_delta=0.01
):
    """
    Calculates the ε-attracting basin for a specified cluster.

    Parameters
    ----------
    adata : anndata.AnnData
        AnnData object whose obs attribute contains information such as 'lon', 'lat', 'cluster', 'seq_id', etc.
    cluster_key : str, default "cluster"
        Column name in obs attribute where cluster information is stored.
    target_cluster_key : str
        The cluster name or label for which to calculate the ε-attracting basin.
    cost_matrix_key : str, default "cost_matrix"
        Key name in uns attribute where the cost matrix is stored. If not present, it will be computed automatically.
    epsilon_delta : float, default 0.01
        Increment of ε (value to increase per iteration).

    Returns
    -------
    None
        The result is stored in adata.obs[f'eps_attracting_basin_{target_cluster_key}'].
    """

    # Check if cluster_key exists in adata.obs
    if cluster_key not in adata.obs.columns:
        raise KeyError(f"'{cluster_key}' not found in adata.obs columns")
    
    # Check if cost matrix exists in adata.uns, otherwise compute from adata.X
    if cost_matrix_key in adata.uns:
        cost_matrix = adata.uns[cost_matrix_key]
    else:
        cost_matrix = pairwise_distances(adata.X)
        adata.uns[cost_matrix_key] = cost_matrix

    # Initialize epsilons with infinity
    epsilons = np.full(adata.shape[0], np.inf, dtype=float)
    matching_cluster_indices = np.arange(adata.shape[0])[adata.obs[cluster_key] == target_cluster_key]
    if len(matching_cluster_indices) == 0:
        raise ValueError("No points found for the specified target cluster.")
    matching_rows = []
    
    for idx in matching_cluster_indices:
        group_id = adata.obs['seq_id'].values[idx]
        group_indices = np.where(adata.obs['seq_id'].values == group_id)[0]
        valid_indices = group_indices[group_indices <= idx]
        matching_rows.extend(adata.obs.index[valid_indices])
    
    idx_0 = np.array(np.unique(matching_rows), dtype=int)
    epsilons[idx_0] = 0
    cost_mat = cost_matrix.copy()
    cost_mat[:, matching_cluster_indices] = np.inf  

    # Phase 1: Increase epsilon for points with cost below epsilon
    epsilon = 0
    dist_min = np.full(adata.shape[0], np.inf, dtype=float)
    while all(epsilons < epsilon) == False:
        remain_indices = np.where(epsilons > epsilon)[0]
        asign_indices = np.where(epsilons <= epsilon)[0]
        cost_submatrix = cost_mat[np.ix_(remain_indices, asign_indices)]
        dist_min[remain_indices] = np.min(cost_submatrix, axis=1)
        for gid in np.unique(adata.obs['seq_id'].values[remain_indices]):
            group_mask = np.where(adata.obs['seq_id'] == gid)[0]
            group_epsilon = np.array([np.min(dist_min[group_mask][i:]) for i in range(len(group_mask))])
            epsilons[group_mask[(group_epsilon < epsilon) & (epsilons[group_mask] > epsilon)]] = epsilon

        if len(remain_indices) == len(np.where(epsilons > epsilon)[0]):
            epsilon = epsilon + epsilon_delta
    epsilons[idx_0] = 0
    
    # Phase 2: Decrease epsilon for points still at 0 based on their minimum costs
    epsilon = 0
    dist_min = np.min(cost_matrix, axis=1)
    while any(epsilons == 0):
        remain_indices = np.where(epsilons == 0)[0]
        asign_indices = np.where(epsilons != 0)[0]
        cost_submatrix = cost_matrix[np.ix_(remain_indices, asign_indices)]
        # cost_submatrix = cost_mat[np.ix_(remain_indices, asign_indices)]
        dist_min[remain_indices] = np.min(cost_submatrix, axis=1)
        for gid in np.unique(adata.obs['seq_id'].values[remain_indices]):
            group_mask = np.where(adata.obs['seq_id'] == gid)[0]
            group_epsilon = np.array([np.min(dist_min[group_mask][i:]) for i in range(len(group_mask))])
            epsilons[group_mask[(group_epsilon < epsilon) & (epsilons[group_mask] == 0)]] = -epsilon

        if len(remain_indices) == len(np.where(epsilons == 0)[0]):
            epsilon = epsilon + epsilon_delta
    
    adata.obs[f'{output_key}_{target_cluster_key}'] = epsilons


def eps_sum_attracting_basin(
    adata, 
    cluster_key="cluster", 
    target_cluster_key="target_cluster", 
    cost_matrix_key="cost_matrix", 
    output_key="eps_sum_attracting_basin",
    epsilon_delta=0.01
):
    """
    Calculates the $\varepsilon_\Sigma$ attracting basin for a specified cluster.

    Parameters
    ----------
    adata : anndata.AnnData
        AnnData object whose obs attribute contains information such as 'lon', 'lat', 'cluster', 'seq_id', etc.
    cluster_key : str, default "cluster"
        Column name in obs attribute where cluster information is stored.
    target_cluster_key : str
        The cluster name or label for which to calculate the ε-attracting basin.
    cost_matrix_key : str, default "cost_matrix"
        Key name in uns attribute where the cost matrix is stored. If not present, it will be computed automatically.
    output_key : str, default "eps_sum_attracting_basin"
        Prefix for the output column name in adata.obs.
    epsilon_delta : float, default 0.01
        Increment of ε (value to increase per iteration).

    Returns
    -------
    None
        The result is stored in adata.obs[f'{output_key}_{target_cluster_key}'].
    """

    # Check if cluster_key exists in adata.obs
    if cluster_key not in adata.obs.columns:
        raise KeyError(f"'{cluster_key}' not found in adata.obs columns")
    
    # Check if cost matrix exists in adata.uns, otherwise compute from adata.X
    if cost_matrix_key in adata.uns:
        cost_matrix = adata.uns[cost_matrix_key]
    else:
        cost_matrix = pairwise_distances(adata.X)
        adata.uns[cost_matrix_key] = cost_matrix

    # Initialize epsilons with infinity
    sig_epsilons = np.full(adata.shape[0], np.inf, dtype=float)
    matching_cluster_indices = np.arange(adata.shape[0])[adata.obs[cluster_key] == target_cluster_key]
    if len(matching_cluster_indices) == 0:
        raise ValueError("No points found for the specified target cluster.")
    matching_rows = []
    
    for idx in matching_cluster_indices:
        group_id = adata.obs['seq_id'].values[idx]
        group_indices = np.where(adata.obs['seq_id'].values == group_id)[0]
        valid_indices = group_indices[group_indices <= idx]
        matching_rows.extend(adata.obs.index[valid_indices])
    idx_0 = np.array(np.unique(matching_rows), dtype=int)
    sig_epsilons[idx_0] = 0
    cost_mat = cost_matrix.copy()
    cost_mat[:, matching_cluster_indices] = np.inf  

    # Phase 1: Calculate ε > 0 by incrementing ε based on the minimum cost for each point
    epsilon = 0
    dist_min = np.full(adata.shape[0], np.inf, dtype=float)
    while not np.all(sig_epsilons < epsilon):
        remain_indices = np.where(sig_epsilons > epsilon)[0]
        asign_indices = np.where(sig_epsilons <= epsilon)[0]
        if len(remain_indices) == 0 or len(asign_indices) == 0:
            break
        cost_submatrix = cost_mat[np.ix_(remain_indices, asign_indices)]
        dist_min[remain_indices] = np.min(cost_submatrix, axis=1)
        for gid in np.unique(adata.obs['seq_id'].values[remain_indices]):
            group_mask = np.where(adata.obs['seq_id'] == gid)[0]
            group_epsilon = np.array([np.min(dist_min[group_mask][i:]) for i in range(len(group_mask))])
            idx = sig_epsilons[group_mask] > epsilon
            # Update sig_epsilons with the minimum of current and (group_epsilon + epsilon)
            sig_epsilons[group_mask[idx]] = np.minimum(sig_epsilons[group_mask[idx]], group_epsilon[idx] + epsilon)
        epsilon += epsilon_delta

    # Phase 2: Readjust points with ε value of 0 based on cost information
    epsilon = 0
    dist_min = np.full(adata.shape[0], np.inf, dtype=float)
    sig_epsilons_tmp = sig_epsilons.copy()
    sig_epsilons_tmp[idx_0] = -np.inf
    while np.any(sig_epsilons_tmp < -epsilon):
        remain_indices = np.where(sig_epsilons_tmp < -epsilon)[0]
        asign_indices = np.where(sig_epsilons_tmp >= -epsilon)[0]
        if len(remain_indices) == 0 or len(asign_indices) == 0:
            break
        cost_submatrix = cost_matrix[np.ix_(remain_indices, asign_indices)]
        # cost_submatrix = cost_mat[np.ix_(remain_indices, asign_indices)]
        dist_min[remain_indices] = np.min(cost_submatrix, axis=1)
        for gid in np.unique(adata.obs['seq_id'].values[remain_indices]):
            group_mask = np.where(adata.obs['seq_id'] == gid)[0]
            group_epsilon = np.array([np.min(dist_min[group_mask][i:]) for i in range(len(group_mask))])
            idx = sig_epsilons_tmp[group_mask] < -epsilon
            sig_epsilons_tmp[group_mask[idx]] = np.maximum(sig_epsilons_tmp[group_mask[idx]], -group_epsilon[idx] - epsilon)
        epsilon += epsilon_delta
    sig_epsilons[idx_0] = sig_epsilons_tmp[idx_0]

    adata.obs[f'{output_key}_{target_cluster_key}'] = sig_epsilons

def sublevel_set_visualization(
    adata,
    plot_key=None,
    eps_key="eps_attracting_basin",
    target_cluster_key="target_cluster",
    linewidth=0.5,
    color_name="$\\varepsilon$",
    vmin=None,
    vmax=None,
    xlim=None,
    ylim=None,
    show_xticks=True,
    show_yticks=True,
    show_xlabel=True,
    show_ylabel=True,
    show_colorbar=True,
    show_colorbar_label=True,
    show_xticklabels=True,
    show_yticklabels=True,
    figsize=(10, 8),
    area_percentile=99.5,
    edge_percentile=99.5,
    save_fig=False,
    save_fig_dir=".",
    save_fig_name="sublevel_set_visualization"
):
    """
    Sublevel sets visualization of $\varepsilon$- or $\varepsilon_\Sigma$-attracting basin.

    Parameters:
        ...
        area_percentile (float): Percentile threshold for triangle area outlier removal.
        edge_percentile (float): Percentile threshold for edge length outlier removal.
    """
    # Get coordinates
    if plot_key is None:
        points = adata.X[:, :2]
    else:
        points = adata.obsm[plot_key][:, :2]
    # Get values for coloring
    colors = adata.obs[f"{eps_key}_{target_cluster_key}"].values

    # Delaunay triangulation
    tri = Delaunay(points)
    simplices = tri.simplices

    # Remove outliers by triangle area
    tri_pts = points[simplices]
    vec1 = np.hstack([tri_pts[:, 1] - tri_pts[:, 0], np.zeros((tri_pts.shape[0], 1))])
    vec2 = np.hstack([tri_pts[:, 2] - tri_pts[:, 0], np.zeros((tri_pts.shape[0], 1))])
    cross = np.cross(vec1, vec2)
    areas = 0.5 * np.abs(cross[:, 2])
    area_threshold = np.percentile(areas, area_percentile)
    mask_area = areas < area_threshold

    # Remove outliers by maximum edge length
    d01 = np.linalg.norm(tri_pts[:, 0] - tri_pts[:, 1], axis=1)
    d12 = np.linalg.norm(tri_pts[:, 1] - tri_pts[:, 2], axis=1)
    d20 = np.linalg.norm(tri_pts[:, 2] - tri_pts[:, 0], axis=1)
    max_edge = np.maximum.reduce([d01, d12, d20])
    edge_threshold = np.percentile(max_edge, edge_percentile)
    mask_edge = max_edge < edge_threshold

    # Combine masks
    mask = mask_area & mask_edge
    filtered_simplices = simplices[mask]

    triang = Triangulation(points[:, 0], points[:, 1], triangles=filtered_simplices)

    # Set color normalization
    if vmin is None:
        vmin = np.percentile(colors, 5)
    if vmax is None:
        vmax = np.percentile(colors, 95)
    
    plt.figure(figsize=figsize)
    norm = TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)
    sf = plt.tripcolor(triang, colors, cmap='bwr', norm=norm, alpha=1, lw=0, zorder=10)
    plt.tricontour(points[:, 0], points[:, 1], filtered_simplices, colors,
                   levels=10, colors='k', linewidths=linewidth, zorder=20)
    
    cbar = None
    if show_colorbar:
        cbar = plt.colorbar(sf)
        cbar.ax.tick_params(labelsize=20)
        if show_colorbar_label:
            cbar.set_label(color_name, fontsize=22)
        else:
            cbar.set_label("")
    xticks = np.arange(np.floor(points[:, 0].min()/5) * 5, np.ceil(points[:, 0].max()/5) * 5 + 1, 5).astype(int)
    yticks = np.arange(np.floor(points[:, 1].min()/5) * 5, np.ceil(points[:, 1].max()/5) * 5 + 1, 5).astype(int)
    
    if show_xticks:
        plt.xticks(xticks)
        if not show_xticklabels:
            plt.gca().set_xticklabels([])
    else:
        plt.xticks([])
    if show_yticks:
        plt.yticks(yticks)
        if not show_yticklabels:
            plt.gca().set_yticklabels([])
    else:
        plt.yticks([])
    if xlim is not None:
        plt.xlim(xlim)
    if ylim is not None:
        plt.ylim(ylim)
    plt.grid(ls="--", zorder=0)
    plt.tick_params(axis='both', which='major', labelsize=20)
    if show_xlabel:
        plt.xlabel("Longitude", fontsize=20)
    if show_ylabel:
        plt.ylabel("Latitude", fontsize=20)
    if save_fig:
        os.makedirs(save_fig_dir, exist_ok=True)
        plt.savefig(f"{save_fig_dir}/{save_fig_name}.png", bbox_inches="tight")
    plt.show()


def plot_attracting_basin(
    adata, 
    plot_key=None,
    eps_key="eps_attracting_basin", 
    cluster_key="cluster", 
    good_cluster_key="good", 
    bad_cluster_key="bad", 
    eps_threshold=0, 
    background=None, 
    fontsize = 14,
    pointsize = 20,
    lon_lim=None, 
    lat_lim=None,
    show_legend=True,
    show_label=True,
    show_ticks=True,
    show_title=True,
    title=None,
    label_clusters = ["$G$", "$B$"],
    label_basins = ["$G_{F, \\varepsilon}$", "$B_{F, \\varepsilon}$"],
    save_fig=False,
    save_fig_dir=".",
    save_fig_name="attracting_basin"
):
    if plot_key is None:
        points = adata.X[:, :2]
    else:
        points = adata.obsm[plot_key][:, :2]

    fig, ax = plt.subplots(figsize=(10, 10))
    # Plot background map using the global 'world' GeoDataFrame
    if background is not None:
        background.plot(ax=ax, color='lightgray', zorder=0)
    ax.set_facecolor('white')

    # Use a normalization for the scatter plot color scale
    other = ax.scatter(points[:,0], points[:,1], color="lightgray", s=10,
               alpha=0.8, zorder=10)

    # Highlight 'good' cluster points where epsilon is below the threshold
    good_mask = adata.obs[f'{eps_key}_{good_cluster_key}'] <= eps_threshold
    good_basin = ax.scatter(points[good_mask][:,0], points[good_mask][:,1], 
               color="red", s=pointsize, linewidths=1.5, zorder=20)

    # Highlight 'bad' cluster points where epsilon is below the negative threshold
    bad_mask = adata.obs[f'{eps_key}_{bad_cluster_key}'] <= -eps_threshold
    bad_basin = ax.scatter(points[bad_mask][:,0], points[bad_mask][:,1], 
               color="blue", s=pointsize, linewidths=1.5, zorder=20)

    # Highlight 'good' cluster points where epsilon is below the threshold
    good_mask2 = adata.obs[cluster_key] == good_cluster_key
    good_cluster = ax.scatter(points[good_mask2][:,0], points[good_mask2][:,1], 
               color="red", s=pointsize, linewidths=1.5, zorder=20, label=rf"{label_basins[0]}")

    # Highlight 'bad' cluster points where epsilon is below the negative threshold
    bad_mask2 = adata.obs[cluster_key] == bad_cluster_key
    bad_cluster = ax.scatter(points[bad_mask2][:,0], points[bad_mask2][:,1], 
               color="blue", s=pointsize, linewidths=1.5, zorder=20, label=rf"{label_basins[1]}")
    # Draw dashed lines between forecast points for each ensemble member
    members = sorted(adata.obs['Member'].unique())
    for member in members:
        member_group = adata.obs[adata.obs['Member'] == member].sort_values('ForecastTime')
        indices = member_group.index.values
        ax.plot(points[indices, 0], points[indices, 1], marker='o', markersize=5,
                color="lightgray", ls="--", lw=1, zorder=15)

    # Add circles to mark the last forecast point for each cluster based on cluster_key values
    good_last = adata.obs[adata.obs[cluster_key] == good_cluster_key].groupby('Member').tail(1)
    bad_last = adata.obs[adata.obs[cluster_key] == bad_cluster_key].groupby('Member').tail(1)
    circ_handles = []
    circ_handles.append(Line2D([0], [0], marker='o', color='red', markerfacecolor='none', 
                               markeredgewidth=2.5, markersize=10, linestyle='None', label=rf"{label_clusters[0]}"))
    circ_handles.append(Line2D([0], [0], marker='o', color='blue', markerfacecolor='none', 
                               markeredgewidth=2.5, markersize=10, linestyle='None', label=rf"{label_clusters[1]}"))
    for idx in good_last.index:
        circ = mpatches.Circle((points[idx, 0], points[idx, 1]), 0.2, edgecolor='red', facecolor='none', 
                          linewidth=2.5, zorder=20)
        ax.add_patch(circ)
    for idx in bad_last.index:
        circ = mpatches.Circle((points[idx, 0], points[idx, 1]), 0.2, edgecolor='blue', facecolor='none', 
                          linewidth=2.5, zorder=20)
        ax.add_patch(circ)

    # Set axis limits and ticks if not provided
    if lon_lim is None:
        lon_lim = ax.get_xlim()
    if lat_lim is None:
        lat_lim = ax.get_ylim()
    xticks = np.arange(np.floor(lon_lim[0] / 2) * 2, np.ceil(lon_lim[1] / 2) * 2 + 1, 2).astype(int)
    yticks = np.arange(np.floor(lat_lim[0] / 2) * 2, np.ceil(lat_lim[1] / 2) * 2 + 1, 2).astype(int)
    ax.set_xlim(lon_lim)
    ax.set_ylim(lat_lim)

    if show_ticks:
        ax.set_xticks(xticks)
        ax.set_yticks(yticks)
        ax.tick_params(axis='both', labelsize=fontsize)
    else:
        ax.set_xticks([])
        ax.set_yticks([])

    ax.grid(ls="--", zorder=0)

    if show_label:
        ax.set_xlabel("Longitude", fontsize=fontsize)
        ax.set_ylabel("Latitude", fontsize=fontsize)
    else:
        ax.set_xticklabels([])
        ax.set_yticklabels([])

    if show_legend:
        # Remove duplicate labels
        handles, labels = ax.get_legend_handles_labels()
        # Remove duplicate scatter labels as well
        label_dict = OrderedDict()
        for h, l in zip(handles, labels):
            if l not in label_dict:
                label_dict[l] = h
        # Merge
        for h in circ_handles:
            label_dict[h.get_label()] = h
        ax.legend(list(label_dict.values()), list(label_dict.keys()), loc="best", fontsize=fontsize)

    if show_title:
        if title is not None:
            ax.set_title(title, fontsize=fontsize)
        else:
            ax.set_title(rf"$\varepsilon$-attracting basin ($\varepsilon$={eps_threshold})", fontsize=fontsize)

    if save_fig:
        os.makedirs(save_fig_dir, exist_ok=True)
        plt.savefig(f"{save_fig_dir}/{save_fig_name}_{eps_threshold:.2f}.png", bbox_inches="tight")

    plt.show()