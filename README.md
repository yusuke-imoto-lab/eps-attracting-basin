# $\varepsilon$-Attracting Basin

> Tools for computing and visualizing **$\varepsilon$-attracting basin** as introduced in 
> *Filtrations Indexed by Attracting Levels and their Applications*. Y. Imoto, T. Mitsui, K. Tokuda, and T. Yokoyama. *in preparation*. 

---

## Table of Contents

* [Overview](#overview)
* [Installation](#installation)
* [Usage](#usage)
* [API Reference](#api-reference)
* [Requirements](#requirements)
* [Exmaples](#exmaples)
* [License](#license)
* [Contact](#contact)

---

## Overview

- The $\varepsilon$-attracting basin $A_{F, \varepsilon}$ denotes the set of which the system underlined by the dynamical system $F$ can be driven into the cluster $A$ if we continue to control it with energy $\varepsilon>0$ *in every time step*.
- The $-\varepsilon$-attracting basin $A_{F, -\varepsilon}$ denotes the set of which the system underlined by the dynamical system $F$ never escape from the cluster $A$ even if we continue to control it with energy $\varepsilon>0$ *in every time step*.
- The $\varepsilon_\Sigma$-attracting basins $A_{F, \varepsilon_\Sigma}$ and $A_{F, -\varepsilon_\Sigma}$ denote the similar set if we control the system with energy $\varepsilon>0$ *in total*. 

<div style="text-align:left"><img style="width:100%; height: auto" src="https://github.com/yusuke-imoto-lab/eps-attracting-basin/blob/main/images/image_eps_attracting_basin.jpg"/></div>

For more theoretical background, see
*Filtrations Indexed by Attracting Levels and their Applications*. Yusuke Imoto, Takahito Mitsui, Keita Tokuda, and Tomoo Yokoyama. arXiv**** (2025)

---

## Installation

```bash
git clone https://github.com/yusuke-imoto-lab/eps-attracting-basin.git
cd eps-attracting-basin
pip install -r requirements.txt
```

---

## Usage

Put epsbasin.py in the same directory as the executable file and import it:

```python
import epsbasin
```

### Input Data Structure

The functions expect an [AnnData](https://anndata.readthedocs.io/en/stable/) object with the following structure:

- **`.X`**  
  A numeric data matrix (observations × features) used for computing pairwise costs if no precomputed cost matrix exists.

- **`.obs`**  
  A pandas DataFrame containing at least:
  - **`'cluster'`**: Categorical labels indicating cluster membership for each observation.  
  - **`'seq_id'`**: Identifiers grouping observations into sequences or trajectories.  

- **`.uns`** (optional)
  - **`'cost_matrix'`**: Optional precomputed pairwise cost matrix (array of shape n_obs × n_obs). If absent, functions will compute it from `.X` and store it here.

- **`.obsm`** (optional)
  - **`'plot_data'`**: 2-dimensional data for 2D visualiztion plotting instead of `.X`.


---

## API Reference
### $\varepsilon$-attracting basin
- **Function:** `eps_attracting_basin(adata, cluster_key="cluster", target_cluster_key="target_cluster", cost_matrix_key="cost_matrix", output_key="eps_attracting_basin", epsilon_delta=0.01)`
- **Description:** Compute the debut function $\underline{\varepsilon}(\ast;A) = \inf\{\varepsilon\in\mathbb{R} \mid\,\ast\in A_{F, \varepsilon} \}$ for the target cluster $A$ (indexed by `target_cluster_key`) for input data `adata.X`. If the cost matrix has not been precomputed (no `adata.uns[cost_matrix_key]`), the Euclidean distance matrix is ​​calculated as the cost matrix.

### $\varepsilon_\Sigma$-attracting basin
- **Function:** `eps_sum_attracting_basin(adata, cluster_key="cluster", target_cluster_key="target_cluster", cost_matrix_key="cost_matrix", output_key="eps_sum_attracting_basin", epsilon_delta=0.01)`
- **Description:** Compute the debut function $\underline{\varepsilon}_{\Sigma}(\ast;A) = \inf\{\varepsilon\in\mathbb{R} \mid\,\ast\in A_{F, \varepsilon_{\Sigma}} \}$ for the target cluster $A$ (indexed by `target_cluster_key`) for input data `adata.X`. If the cost matrix has not been precomputed (no `adata.uns[cost_matrix_key]`), the Euclidean distance matrix is ​​calculated as the cost matrix.

###  Sublevel set visualization
- **Function:** `sublevel_set_visualization(adata, plot_key=None, eps_key="eps_attracting_basin", target_cluster_key="target_cluster", linewidth=0.5, color_name="$\\varepsilon$", vmin=None, vmax=None, xlim=None, ylim=None, show_xticks=True, show_yticks=True, show_xlabel=True, show_ylabel=True, show_colorbar=True, show_colorbar_label=True, show_xticklabels=True, show_yticklabels=True, figsize=(10,8), area_percentile=99.5, edge_percentile=99.5, save_fig=False, save_fig_dir=".", save_fig_name="sublevel_set_visualization")`
- **Description:** Visualize the $\varepsilon-$ or $\varepsilon_\Sigma$-attracting basin (tha values on points are the debut function $\underline{\varepsilon}(\ast;A)$ or $\underline{\varepsilon}_{\Sigma}(\ast;A)$) for the target cluster $A$ (indexed by `target_cluster_key`) on a 2D Delaunay triangulation generated from `adata.X` (or `adata.obsm[plot_key]`). Triangles with too large sides or areas can be removed by adjusting the parameters `area_percentile` and `edge_percentile`.

### Plot attracting basin
- **Function:**: `plot_attracting_basin(adata, plot_key=None, eps_key="eps_attracting_basin", cluster_key="cluster", good_cluster_key="good", bad_cluster_key="bad", eps_threshold=0, background=None, fontsize=14, pointsize=20, lon_lim=None, lat_lim=None, show_legend=True, show_label=True, show_ticks=True, show_title=True, title=None, label_clusters=["$G$", "$B$"], label_basins=["$G_{F, \\varepsilon}$", "$B_{F, \\varepsilon}$"], save_fig=False, save_fig_dir=".", save_fig_name="attracting_basin")`
- **Description:**: Plot the $\varepsilon-$ or $\varepsilon_\Sigma$-attracting basin of good and bad clusters (indexed by `good_cluster_key` anc `bad_cluster_key`) by mapping 2D coordinates from `adata.X` or `adata.obsm[plot_key]`.  

*For detailed usage, see* `espbasin.py`.

---

## Examples

* [**example_Dorphin.ipynb**]("https://github.com/yusuke-imoto-lab/eps-attracting-basin/blob/main/example_Dorphin.ipynb"): An application exmaple pf ensemble weather forecast data for tropical cyclone (typhoon) ``Drophin'', generated by the Meso-scale Ensemble Prediction System (MEPS), provided by the Japan Meteorological Agency (JSA).

<div style="text-align:left"><img style="width:100%; height: auto" src="https://github.com/yusuke-imoto-lab/eps-attracting-basin/blob/main/images/example_Dorphin.jpg"/></div>

**a**, Best‐track trajectory (real orbit, blue line) overlaid with all forecast center positions (gray dots) of Drophin. **b**, Selected ensemble forecast tracks at representative initial times. **c**, Clustering and assigning good and bad clusters. **d**, Distribution of debut functions $\underline{\varepsilon}$ (top) and $\underline{\varepsilon}_{\Sigma}$ (bottom) for clusters $G$ (reft) and $B$ (right). 


---

## Requirements

* Python ≥ 3.8
* numpy
* pandas
* scipy
* scikit-learn
* matplotlib

Install all dependencies via:

```bash
pip install -r requirements.txt
```

---

## License

MIT © 2025 Yusuke Imoto

---

## Contact

* **Yusuke Imoto**
* Email: [imoto.yusuke.4e@kyoto-u.ac.jp](mailto:imoto.yusuke.4e@kyoto-u.ac.jp)
* GitHub: [yusuke-imoto-lab/eps-attracting-basin](https://github.com/yusuke-imoto-lab/eps-attracting-basin)

