# ε-Attracting Basin
<!--
> Tools for computing and visualizing **ε-attracting basin** as introduced in 
> *Filtrations Indexed by Attracting Levels and their Applications*.

---

## Table of Contents

* [Overview](#overview)
* [Installation](#installation)
* [Dataset Structure](#dataset-structure)
* [Usage](#usage)

  * [Compute ε and εΣ debut series](#compute-ε-and-εΣ-debut-series)
  * [Plot time-series of debut functions](#plot-time-series-of-debut-functions)
* [API Reference](#api-reference)
* [Requirements](#requirements)
* [Contributing](#contributing)
* [License](#license)
* [Contact](#contact)

---

## Overview

This repository provides:

1. **Computation** of per-step (`ε`) and total‑energy (`εΣ`) debut functions for an attractor cluster in ensemble forecast data.
2. **Visualization** routines to plot how these debut values evolve over forecast time.

Under the hood, you build a distance matrix among ensemble tracks, identify “good” vs “bad” clusters, then for each time step compute the smallest perturbation needed to steer the system into (or keep it in) the target cluster.

For more theoretical background, see
*Filtrations Indexed by Attracting Levels and their Applications*
Yusuke Imoto, Takahito Mitsui, Keita Tokuda, and Tomoo Yokoyama (2025)

---

## Installation

```bash
git clone https://github.com/yusuke-imoto-lab/eps-attracting-basin.git
cd eps-attracting-basin
pip install -r requirements.txt
```

---

## Dataset Structure

Place your typhoon track CSVs under:

```
data/
└── Typhoon/
    └── Tracks.2020-09-Dolphin.csv
```

Each CSV should contain the following columns:

* `ForecastTime` (e.g. `"2020-09-01T00:00:00"`)
* `Latitude`
* `Longitude`
* `MemberID`

---

## Usage

Open and run the Jupyter notebook to step through the full workflow:

```bash
jupyter lab esp_attracting_basin.ipynb
```

### Compute ε and εΣ debut series

```python
from esp_attracting_basin import (
    load_typhoon_csv,
    identify_clusters,
    compute_epsilon_debut_series,
    compute_epsilon_sum_debut_series,
)

# 1) Load the CSV
df = load_typhoon_csv("data/Typhoon/Tracks.2020-09-Dolphin.csv")

# 2) Identify “good” vs “bad” clusters at final time
good_cluster, bad_cluster = identify_clusters(df, time_column="ForecastTime")

# 3) Compute per‑step debut (ε) series
eps_series = compute_epsilon_debut_series(df, cluster_id=good_cluster)

# 4) Compute total‑energy debut (εΣ) series
eps_sum_series = compute_epsilon_sum_debut_series(df, cluster_id=good_cluster)
```

### Plot time‑series of debut functions

```python
from esp_attracting_basin import plot_eps_time_series

# Plot ε‑debut over forecast lead time
plot_eps_time_series(df, eps_series, cluster_id=good_cluster, title="ε‑Debut Series")

# Plot εΣ‑debut (total-energy) over forecast lead time
plot_eps_time_series(df, eps_sum_series, cluster_id=good_cluster, title="εΣ‑Debut Series")
```

---

## API Reference

| Function                                              | Description                                                     |
| ----------------------------------------------------- | --------------------------------------------------------------- |
| `load_typhoon_csv(path: str) → DataFrame`             | Load a forecast CSV into a pandas DataFrame                     |
| `identify_clusters(df, time_column: str)`             | Return `(good_cluster, bad_cluster)` IDs at the final timestamp |
| `compute_epsilon_debut_series(df, cluster_id)`        | Compute per‑step ε‑debut for each forecast lead time            |
| `compute_epsilon_sum_debut_series(df, cluster_id)`    | Compute total‑energy εΣ‑debut over the entire time series       |
| `plot_eps_time_series(df, series, cluster_id, title)` | Plot a Matplotlib time‑series of ε or εΣ debut values           |

*For detailed usage, see* `esp_attracting_basin.ipynb`.

---

## Requirements

* Python ≥ 3.8
* numpy
* pandas
* scipy
* scikit-learn
* matplotlib
* geopandas (optional, for geographic plotting)

Install all dependencies via:

```bash
pip install -r requirements.txt
```

---

## Contributing

1. Fork the repository
2. Create a feature branch:

   ```bash
   git checkout -b feature/YourFeature
   ```
3. Commit your changes:

   ```bash
   git commit -m "Add new feature"
   ```
4. Push and open a pull request:

   ```bash
   git push origin feature/YourFeature
   ```

We welcome bug reports and feature requests via GitHub Issues.

---

## License

MIT © 2025 Yusuke Imoto

---

## Contact

* **Yusuke Imoto**
* Email: [imoto.yusuke.4e@kyoto-u.ac.jp](mailto:imoto.yusuke.4e@kyoto-u.ac.jp)
* GitHub: [yusuke-imoto-lab/eps-attracting-basin](https://github.com/yusuke-imoto-lab/eps-attracting-basin)
-->
