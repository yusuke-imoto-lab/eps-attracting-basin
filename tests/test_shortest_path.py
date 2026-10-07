import numpy as np
import pandas as pd

from epsbasin.shortest_path import (
    build_terminal_state_graph,
    eps_attracting_basin_shortest_path,
    eps_sum_attracting_basin_shortest_path,
)


class MiniAnnData:
    """Small AnnData-like object sufficient for algorithmic unit tests."""

    def __init__(self, X, obs):
        self.X = np.asarray(X, dtype=float)
        self.obs = obs.copy()
        self.uns = {}
        self.obsm = {}

    @property
    def shape(self):
        return self.X.shape


def make_two_track_adata():
    # Two parallel trajectories.  The upper one terminates in Bad and the lower
    # one in Good.  Cross-track Euclidean distance is exactly 2 everywhere.
    x = np.array(
        [
            [0.0, 0.0],
            [1.0, 0.0],
            [2.0, 0.0],
            [0.0, 2.0],
            [1.0, 2.0],
            [2.0, 2.0],
        ]
    )
    obs = pd.DataFrame(
        {
            "seq_id": [0, 0, 0, 1, 1, 1],
            "cluster": ["other", "other", "good", "other", "other", "bad"],
            "time": [0, 1, 2, 0, 1, 2],
        }
    )
    return MiniAnnData(X=x, obs=obs)


def test_graph_orientation_and_weights():
    adata = make_two_track_adata()
    graph = build_terminal_state_graph(adata, time_key="time")
    w = graph.weights

    assert w[0, 1] == 0.0
    assert w[0, 2] == 0.0
    assert np.isinf(w[1, 0])
    assert np.isinf(w[2, 0])

    # Different sequences retain the metric cost in both directions.
    assert np.isclose(w[0, 3], 2.0)
    assert np.isclose(w[3, 0], 2.0)


def test_signed_debut_sum_parallel_tracks():
    adata = make_two_track_adata()
    eps_sum_attracting_basin_shortest_path(
        adata,
        time_key="time",
        distance_key="d_sum",
        terminal_class_key="terminal_class",
    )

    good = adata.obs["eps_sum_attracting_basin_sp_good"].to_numpy()
    bad = adata.obs["eps_sum_attracting_basin_sp_bad"].to_numpy()
    landscape = adata.obs["eps_sum_attracting_basin_sp_landscape"].to_numpy()

    assert np.allclose(good[:3], -2.0)
    assert np.allclose(bad[:3], 2.0)
    assert np.allclose(good[3:], 2.0)
    assert np.allclose(bad[3:], -2.0)
    assert np.allclose(landscape, -2.0)


def test_max_distance_never_exceeds_sum_distance():
    adata = make_two_track_adata()
    eps_attracting_basin_shortest_path(
        adata,
        time_key="time",
        distance_key="d_max",
    )
    eps_sum_attracting_basin_shortest_path(
        adata,
        time_key="time",
        distance_key="d_sum",
    )

    for label in ["good", "bad"]:
        d_max = adata.obs[f"d_max_{label}"].to_numpy()
        d_sum = adata.obs[f"d_sum_{label}"].to_numpy()
        assert np.all(d_max <= d_sum + 1e-12)


def test_terminal_labels_are_required():
    adata = make_two_track_adata()
    adata.obs.loc[adata.obs.index[-1], "cluster"] = "other"

    try:
        build_terminal_state_graph(adata, time_key="time")
    except ValueError as exc:
        assert "Every terminal state" in str(exc)
    else:
        raise AssertionError("Expected ValueError for an unlabelled terminal state")


def test_bottleneck_and_sum_paths_are_distinct_when_expected():
    # Three two-point sequences with a custom cross-sequence cost matrix.
    # Bad sequence 0 can reach Good sequence 2 either directly at cost 3,
    # or through sequence 1 using two jumps of cost 1.  Therefore:
    # bottleneck distance = 1, additive distance = 2.
    X = np.zeros((6, 1))
    obs = pd.DataFrame(
        {
            "seq_id": [0, 0, 1, 1, 2, 2],
            "cluster": ["other", "bad", "other", "bad", "other", "good"],
        }
    )
    adata = MiniAnnData(X, obs)

    cost = np.full((6, 6), 10.0)
    np.fill_diagonal(cost, 0.0)
    # Cross-sequence links used by the cheaper two-jump path.
    cost[1, 2] = 1.0
    cost[3, 4] = 1.0
    # Direct Bad -> Good jump.
    cost[1, 4] = 3.0
    # Keep reverse cross links finite as well.
    cost[2, 1] = 1.0
    cost[4, 3] = 1.0
    cost[4, 1] = 3.0
    adata.uns["cost_matrix"] = cost

    eps_attracting_basin_shortest_path(adata, distance_key="d_max")
    eps_sum_attracting_basin_shortest_path(adata, distance_key="d_sum")

    assert np.isclose(adata.obs["d_max_good"].iloc[0], 1.0)
    assert np.isclose(adata.obs["d_sum_good"].iloc[0], 2.0)
