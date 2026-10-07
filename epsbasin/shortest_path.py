"""Shortest-path and time-layered implementations for epsilon-attracting basins.

Two transition models are available.

``transition_mode='unrestricted'`` preserves the original graph construction:
forward motion along an observed sequence is free, backward motion is forbidden,
and cross-sequence motion uses the supplied pairwise cost irrespective of time.

``transition_mode='time_layered'`` is intended for non-autonomous systems such as
ensemble weather forecasts.  For a non-final time ``t_k``, a control may switch
from the current state to another ensemble member only at the same time ``t_k``;
the chosen member then advances naturally to ``t_{k+1}``.  Thus at most one
control/switch is allowed per forecast step and transitions cannot jump to a
different time layer.  No control or transition is allowed from the final time
layer.  The resulting costs are computed by backward dynamic programming.
"""

from __future__ import annotations

from dataclasses import dataclass
import heapq
from typing import Literal, Sequence

import numpy as np
from sklearn.metrics import pairwise_distances


PathCost = Literal["max", "sum"]
SequenceEdgeMode = Literal["all_forward", "adjacent"]
TransitionMode = Literal["unrestricted", "time_layered"]

__all__ = [
    "TerminalStateGraph",
    "ExternalStateEpsilon",
    "PathStep",
    "OptimalPath",
    "ExternalStatePathResult",
    "build_terminal_state_graph",
    "terminal_state_shortest_path_debut",
    "eps_attracting_basin_shortest_path",
    "eps_sum_attracting_basin_shortest_path",
    "minimum_epsilon_to_good_bad",
    "minimum_epsilon_and_paths_to_good_bad",
    "reconstruct_terminal_path",
]


@dataclass(frozen=True)
class TerminalStateGraph:
    """Container for a directed graph and terminal-state metadata."""

    weights: np.ndarray
    terminal_indices: np.ndarray
    terminal_labels: np.ndarray
    terminal_label_by_observation: np.ndarray
    sequence_ids: np.ndarray


@dataclass(frozen=True)
class ExternalStateEpsilon:
    """Good/Bad reach costs for an external query state."""

    epsilon_good: float
    epsilon_bad: float
    entry_index_good: int
    entry_index_bad: int
    path_cost: str
    good_cluster_key: str
    bad_cluster_key: str

    @property
    def score_good(self) -> float:
        """Signed score C_B - C_G; positive values favor Good."""
        return float(self.epsilon_bad - self.epsilon_good)


@dataclass(frozen=True)
class PathStep:
    """One transition in a path realizing an external-state reach cost.

    ``source_index=None`` denotes the external query state ``x``.
    ``kind`` is one of ``entry``, ``switch``, ``advance``, or ``graph``.
    ``advance`` steps have zero control cost in the time-layered model.
    """

    kind: str
    source_index: int | None
    destination_index: int
    cost: float


@dataclass(frozen=True)
class OptimalPath:
    """One optimal path from an external query to a target outcome."""

    target_label: str
    path_cost: str
    value: float
    steps: tuple[PathStep, ...]
    terminal_index: int

    @property
    def control_costs(self) -> tuple[float, ...]:
        return tuple(float(step.cost) for step in self.steps)

    @property
    def realized_value(self) -> float:
        costs = self.control_costs
        if not costs:
            return np.inf if np.isinf(self.value) else 0.0
        if self.path_cost == "sum":
            return float(np.sum(costs))
        if self.path_cost == "max":
            return float(np.max(costs))
        raise ValueError(f"Unknown path_cost={self.path_cost!r}.")

    @property
    def observation_indices(self) -> tuple[int, ...]:
        """Destination observations visited after the external query state."""
        return tuple(int(step.destination_index) for step in self.steps)


@dataclass(frozen=True)
class ExternalStatePathResult:
    """Good/Bad reach costs and one realizing path for each target class."""

    epsilon_good: float
    epsilon_bad: float
    entry_index_good: int
    entry_index_bad: int
    path_cost: str
    good_cluster_key: str
    bad_cluster_key: str
    path_good: OptimalPath
    path_bad: OptimalPath

    @property
    def score_good(self) -> float:
        """Signed score C_B - C_G; positive values favor Good."""
        return float(self.epsilon_bad - self.epsilon_good)


def _validate_obs_columns(adata, columns: Sequence[str]) -> None:
    missing = [column for column in columns if column not in adata.obs.columns]
    if missing:
        missing_text = ", ".join(repr(column) for column in missing)
        raise KeyError(f"Missing required adata.obs column(s): {missing_text}")


def _sequence_order(adata, seq_key: str, time_key: str | None) -> list[np.ndarray]:
    seq_values = np.asarray(adata.obs[seq_key])
    result: list[np.ndarray] = []

    for seq_id in adata.obs[seq_key].unique():
        indices = np.flatnonzero(seq_values == seq_id)
        if time_key is not None:
            times = np.asarray(adata.obs[time_key])[indices]
            order = np.argsort(times, kind="stable")
            indices = indices[order]
        result.append(indices.astype(int, copy=False))
    return result


def _get_cost_matrix(adata, cost_matrix_key: str) -> np.ndarray:
    if cost_matrix_key in adata.uns:
        raw_cost_matrix = adata.uns[cost_matrix_key]
        if hasattr(raw_cost_matrix, "toarray"):
            raw_cost_matrix = raw_cost_matrix.toarray()
        cost_matrix = np.asarray(raw_cost_matrix, dtype=float)
    else:
        cost_matrix = np.asarray(pairwise_distances(adata.X), dtype=float)
        adata.uns[cost_matrix_key] = cost_matrix.copy()

    n_obs = int(adata.shape[0])
    if cost_matrix.shape != (n_obs, n_obs):
        raise ValueError(
            f"adata.uns[{cost_matrix_key!r}] must have shape {(n_obs, n_obs)}, "
            f"got {cost_matrix.shape}."
        )
    if np.isnan(cost_matrix).any():
        raise ValueError("The cost matrix contains NaN values.")
    if np.isneginf(cost_matrix).any():
        raise ValueError("The cost matrix may contain +inf, but not -inf.")
    finite = np.isfinite(cost_matrix)
    if np.any(cost_matrix[finite] < 0):
        raise ValueError("The cost matrix must be non-negative (or +inf).")
    return cost_matrix


def _terminal_metadata(
    adata,
    *,
    seq_key: str,
    cluster_key: str,
    good_cluster_key: str,
    bad_cluster_key: str,
    time_key: str | None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[np.ndarray]]:
    ordered_sequences = _sequence_order(adata, seq_key=seq_key, time_key=time_key)
    seq_values = np.asarray(adata.obs[seq_key])
    cluster_values = np.asarray(adata.obs[cluster_key], dtype=object)

    terminal_indices: list[int] = []
    terminal_labels: list[object] = []
    terminal_label_by_observation = np.empty(adata.shape[0], dtype=object)

    for indices in ordered_sequences:
        if len(indices) == 0:
            continue
        terminal_idx = int(indices[-1])
        terminal_label = cluster_values[terminal_idx]
        if terminal_label not in {good_cluster_key, bad_cluster_key}:
            seq_id = seq_values[terminal_idx]
            raise ValueError(
                "Every terminal state must be labelled good or bad. "
                f"Sequence {seq_id!r} ends at observation {terminal_idx}, whose "
                f"{cluster_key!r} value is {terminal_label!r}."
            )
        terminal_indices.append(terminal_idx)
        terminal_labels.append(terminal_label)
        terminal_label_by_observation[indices] = terminal_label

    terminal_indices_array = np.asarray(terminal_indices, dtype=int)
    terminal_labels_array = np.asarray(terminal_labels, dtype=object)
    if terminal_indices_array.size == 0:
        raise ValueError("No sequences / terminal states were found.")
    if not np.any(terminal_labels_array == good_cluster_key):
        raise ValueError(f"No terminal state is labelled {good_cluster_key!r}.")
    if not np.any(terminal_labels_array == bad_cluster_key):
        raise ValueError(f"No terminal state is labelled {bad_cluster_key!r}.")

    return (
        terminal_indices_array,
        terminal_labels_array,
        terminal_label_by_observation,
        ordered_sequences,
    )


def _time_layer_layout(
    adata,
    *,
    seq_key: str,
    time_key: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (sequence_ids, time_values, observation_index_matrix).

    Time-layered dynamics requires a rectangular ensemble: each sequence must
    contain exactly one state for every time value.
    """
    _validate_obs_columns(adata, [seq_key, time_key])
    seq_values = np.asarray(adata.obs[seq_key])
    time_values_obs = np.asarray(adata.obs[time_key])

    sequence_ids = np.asarray(list(adata.obs[seq_key].unique()), dtype=object)
    time_values = np.unique(time_values_obs)
    try:
        time_values = np.sort(time_values)
    except TypeError as exc:
        raise ValueError("time_key values must have a common sortable type.") from exc

    index_matrix = np.full(
        (sequence_ids.size, time_values.size),
        -1,
        dtype=int,
    )

    for s_pos, seq_id in enumerate(sequence_ids):
        seq_idx = np.flatnonzero(seq_values == seq_id)
        seq_times = time_values_obs[seq_idx]
        if seq_idx.size != time_values.size:
            raise ValueError(
                "time_layered mode requires every sequence to have exactly one "
                f"observation at each time. Sequence {seq_id!r} has {seq_idx.size} "
                f"observations but the common time grid has {time_values.size}."
            )
        for t_pos, time_value in enumerate(time_values):
            matches = seq_idx[seq_times == time_value]
            if matches.size != 1:
                raise ValueError(
                    "time_layered mode requires exactly one observation per "
                    f"(sequence, time). Sequence {seq_id!r}, time {time_value!r} "
                    f"has {matches.size} observations."
                )
            index_matrix[s_pos, t_pos] = int(matches[0])

    return sequence_ids, time_values, index_matrix


def build_terminal_state_graph(
    adata,
    *,
    seq_key: str = "seq_id",
    cluster_key: str = "cluster",
    good_cluster_key: str = "good",
    bad_cluster_key: str = "bad",
    cost_matrix_key: str = "cost_matrix",
    time_key: str | None = None,
    sequence_edge_mode: SequenceEdgeMode = "all_forward",
    transition_mode: TransitionMode = "unrestricted",
    graph_key: str | None = None,
) -> TerminalStateGraph:
    """Build a directed graph for terminal-state analysis.

    In ``unrestricted`` mode this preserves the original graph: same-sequence
    forward motion is free, backward motion is forbidden, and cross-sequence
    edges use the pairwise cost irrespective of time.

    In ``time_layered`` mode, the inspection graph follows the literal
    non-autonomous rule: cross-sequence edges are finite only when the two
    observations have the same ``time_key`` value and are not in the final
    time layer; same-sequence adjacent forward motion has weight zero; all other
    inter-time jumps are infinite. No control or transition is allowed from the
    final time layer. The actual time-layered reach costs are *not* obtained by
    Dijkstra on this inspection graph, because that would allow repeated
    same-time switches. They are computed by backward dynamic programming in
    :func:`terminal_state_shortest_path_debut`, which permits at most one
    same-time switch per non-final forecast step.
    """
    required = [seq_key, cluster_key]
    if time_key is not None:
        required.append(time_key)
    _validate_obs_columns(adata, required)

    if good_cluster_key == bad_cluster_key:
        raise ValueError("good_cluster_key and bad_cluster_key must be different.")
    if sequence_edge_mode not in {"all_forward", "adjacent"}:
        raise ValueError(
            "sequence_edge_mode must be either 'all_forward' or 'adjacent'."
        )
    if transition_mode not in {"unrestricted", "time_layered"}:
        raise ValueError(
            "transition_mode must be either 'unrestricted' or 'time_layered'."
        )
    if transition_mode == "time_layered" and time_key is None:
        raise ValueError("time_key is required when transition_mode='time_layered'.")

    cost_matrix = _get_cost_matrix(adata, cost_matrix_key)
    (
        terminal_indices,
        terminal_labels,
        terminal_label_by_observation,
        ordered_sequences,
    ) = _terminal_metadata(
        adata,
        seq_key=seq_key,
        cluster_key=cluster_key,
        good_cluster_key=good_cluster_key,
        bad_cluster_key=bad_cluster_key,
        time_key=time_key,
    )

    if transition_mode == "unrestricted":
        weights = cost_matrix.copy()
        np.fill_diagonal(weights, 0.0)
        for indices in ordered_sequences:
            weights[np.ix_(indices, indices)] = np.inf
            weights[indices, indices] = 0.0
            if sequence_edge_mode == "all_forward":
                for position, source in enumerate(indices[:-1]):
                    destinations = indices[position + 1 :]
                    weights[source, destinations] = 0.0
            elif len(indices) > 1:
                weights[indices[:-1], indices[1:]] = 0.0
    else:
        _, _, index_matrix = _time_layer_layout(
            adata,
            seq_key=seq_key,
            time_key=time_key,
        )
        weights = np.full_like(cost_matrix, np.inf, dtype=float)
        np.fill_diagonal(weights, 0.0)

        # Cross-sequence switching is allowed only within the same non-final
        # time layer. No switch or other outgoing transition is allowed from the
        # final time layer.
        for t_pos in range(index_matrix.shape[1] - 1):
            current = index_matrix[:, t_pos]
            weights[np.ix_(current, current)] = cost_matrix[
                np.ix_(current, current)
            ]

        # Natural forward dynamics of each member is free and moves exactly one
        # time layer forward. Repeated same-time switches are possible in this
        # inspection matrix, which is why the actual reach costs use the DP
        # recursion below rather than Dijkstra on this matrix.
        for indices in ordered_sequences:
            if len(indices) > 1:
                weights[indices[:-1], indices[1:]] = 0.0

    if graph_key is not None:
        adata.uns[graph_key] = weights.copy()

    return TerminalStateGraph(
        weights=weights,
        terminal_indices=terminal_indices,
        terminal_labels=terminal_labels,
        terminal_label_by_observation=terminal_label_by_observation,
        sequence_ids=np.asarray(adata.obs[seq_key]).copy(),
    )


def _multi_source_distance_to_targets(
    weights: np.ndarray,
    target_indices: np.ndarray,
    *,
    path_cost: PathCost,
) -> tuple[np.ndarray, np.ndarray]:
    """Distances from every vertex to the closest target in a directed graph."""
    if path_cost not in {"sum", "max"}:
        raise ValueError("path_cost must be either 'sum' or 'max'.")

    n = weights.shape[0]
    distance = np.full(n, np.inf, dtype=float)
    next_hop = np.full(n, -1, dtype=int)
    heap: list[tuple[float, int]] = []

    for target in np.unique(np.asarray(target_indices, dtype=int)):
        distance[target] = 0.0
        heapq.heappush(heap, (0.0, int(target)))

    while heap:
        current_distance, vertex = heapq.heappop(heap)
        if current_distance != distance[vertex]:
            continue
        incoming = weights[:, vertex]
        predecessors = np.flatnonzero(np.isfinite(incoming))
        for predecessor in predecessors:
            edge_weight = float(incoming[predecessor])
            if path_cost == "sum":
                candidate = current_distance + edge_weight
            else:
                candidate = max(current_distance, edge_weight)
            if candidate < distance[predecessor]:
                distance[predecessor] = candidate
                next_hop[predecessor] = vertex
                heapq.heappush(heap, (candidate, int(predecessor)))

    return distance, next_hop


def _time_layered_distance_to_target(
    adata,
    *,
    target_label: str,
    path_cost: PathCost,
    seq_key: str,
    time_key: str,
    cluster_key: str,
    cost_matrix_key: str,
) -> tuple[np.ndarray, np.ndarray]:
    """Backward dynamic programming for one target terminal class."""
    if path_cost not in {"sum", "max"}:
        raise ValueError("path_cost must be either 'sum' or 'max'.")

    _, _, index_matrix = _time_layer_layout(
        adata,
        seq_key=seq_key,
        time_key=time_key,
    )
    cost_matrix = _get_cost_matrix(adata, cost_matrix_key)
    cluster_values = np.asarray(adata.obs[cluster_key], dtype=object)

    n_seq, n_times = index_matrix.shape
    distance = np.full(adata.shape[0], np.inf, dtype=float)
    next_hop = np.full(adata.shape[0], -1, dtype=int)

    final_indices = index_matrix[:, -1]
    final_labels = cluster_values[final_indices]
    if not np.any(final_labels == target_label):
        raise ValueError(f"No terminal state is labelled {target_label!r}.")

    # No perturbation or time evolution is allowed from the final time layer.
    # Therefore a terminal state has zero cost to its own class and infinite
    # cost to the other class. next_hop remains -1 at the final time.
    distance[final_indices] = np.where(
        final_labels == target_label,
        0.0,
        np.inf,
    )

    # Backward recursion. At t_k choose one sequence n at the same time, pay
    # d(x_m^k, x_n^k), then follow n to its next-time state x_n^{k+1}.
    for t_pos in range(n_times - 2, -1, -1):
        current = index_matrix[:, t_pos]
        following = index_matrix[:, t_pos + 1]
        switching_cost = cost_matrix[np.ix_(current, current)]
        future = distance[following]

        if path_cost == "max":
            candidates = np.maximum(switching_cost, future[None, :])
        else:
            candidates = switching_cost + future[None, :]

        best_seq = np.argmin(candidates, axis=1)
        distance[current] = candidates[np.arange(n_seq), best_seq]
        next_hop[current] = following[best_seq]

    return distance, next_hop


def _signed_debut(
    *,
    distance_to_good: np.ndarray,
    distance_to_bad: np.ndarray,
    terminal_label_by_observation: np.ndarray,
    good_cluster_key: str,
    bad_cluster_key: str,
) -> tuple[np.ndarray, np.ndarray]:
    """Convert target distances into complementary signed debut functions."""
    natural_good = terminal_label_by_observation == good_cluster_key
    natural_bad = terminal_label_by_observation == bad_cluster_key

    good_debut = np.empty_like(distance_to_good)
    bad_debut = np.empty_like(distance_to_bad)

    good_debut[natural_good] = -distance_to_bad[natural_good]
    bad_debut[natural_good] = distance_to_bad[natural_good]
    good_debut[natural_bad] = distance_to_good[natural_bad]
    bad_debut[natural_bad] = -distance_to_good[natural_bad]
    return good_debut, bad_debut


def terminal_state_shortest_path_debut(
    adata,
    *,
    path_cost: PathCost = "sum",
    seq_key: str = "seq_id",
    cluster_key: str = "cluster",
    good_cluster_key: str = "good",
    bad_cluster_key: str = "bad",
    cost_matrix_key: str = "cost_matrix",
    time_key: str | None = None,
    sequence_edge_mode: SequenceEdgeMode = "all_forward",
    transition_mode: TransitionMode = "unrestricted",
    output_key: str = "terminal_shortest_path",
    landscape_key: str | None = None,
    distance_key: str | None = None,
    terminal_class_key: str | None = None,
    graph_key: str | None = None,
    store_next_hop: bool = False,
) -> None:
    """Compute Good/Bad debut functions from terminal-state reach costs.

    ``transition_mode='time_layered'`` requires ``time_key`` and implements a
    non-autonomous ensemble model: at each non-final forecast time, switching is
    allowed only between states at the same time, at most once per forecast
    step, followed by natural advance to the next time layer. No perturbation or
    transition is allowed from the final time layer.
    """
    required = [seq_key, cluster_key]
    if time_key is not None:
        required.append(time_key)
    _validate_obs_columns(adata, required)
    if transition_mode not in {"unrestricted", "time_layered"}:
        raise ValueError(
            "transition_mode must be either 'unrestricted' or 'time_layered'."
        )
    if transition_mode == "time_layered" and time_key is None:
        raise ValueError("time_key is required when transition_mode='time_layered'.")

    (
        terminal_indices,
        terminal_labels,
        terminal_label_by_observation,
        _,
    ) = _terminal_metadata(
        adata,
        seq_key=seq_key,
        cluster_key=cluster_key,
        good_cluster_key=good_cluster_key,
        bad_cluster_key=bad_cluster_key,
        time_key=time_key,
    )

    if transition_mode == "unrestricted":
        graph = build_terminal_state_graph(
            adata,
            seq_key=seq_key,
            cluster_key=cluster_key,
            good_cluster_key=good_cluster_key,
            bad_cluster_key=bad_cluster_key,
            cost_matrix_key=cost_matrix_key,
            time_key=time_key,
            sequence_edge_mode=sequence_edge_mode,
            transition_mode="unrestricted",
            graph_key=graph_key,
        )
        good_terminals = graph.terminal_indices[
            graph.terminal_labels == good_cluster_key
        ]
        bad_terminals = graph.terminal_indices[
            graph.terminal_labels == bad_cluster_key
        ]
        distance_to_good, next_hop_good = _multi_source_distance_to_targets(
            graph.weights, good_terminals, path_cost=path_cost
        )
        distance_to_bad, next_hop_bad = _multi_source_distance_to_targets(
            graph.weights, bad_terminals, path_cost=path_cost
        )
    else:
        # Validate/store the inspection graph only if requested.
        if graph_key is not None:
            build_terminal_state_graph(
                adata,
                seq_key=seq_key,
                cluster_key=cluster_key,
                good_cluster_key=good_cluster_key,
                bad_cluster_key=bad_cluster_key,
                cost_matrix_key=cost_matrix_key,
                time_key=time_key,
                sequence_edge_mode=sequence_edge_mode,
                transition_mode="time_layered",
                graph_key=graph_key,
            )
        distance_to_good, next_hop_good = _time_layered_distance_to_target(
            adata,
            target_label=good_cluster_key,
            path_cost=path_cost,
            seq_key=seq_key,
            time_key=time_key,
            cluster_key=cluster_key,
            cost_matrix_key=cost_matrix_key,
        )
        distance_to_bad, next_hop_bad = _time_layered_distance_to_target(
            adata,
            target_label=bad_cluster_key,
            path_cost=path_cost,
            seq_key=seq_key,
            time_key=time_key,
            cluster_key=cluster_key,
            cost_matrix_key=cost_matrix_key,
        )

    good_debut, bad_debut = _signed_debut(
        distance_to_good=distance_to_good,
        distance_to_bad=distance_to_bad,
        terminal_label_by_observation=terminal_label_by_observation,
        good_cluster_key=good_cluster_key,
        bad_cluster_key=bad_cluster_key,
    )

    adata.obs[f"{output_key}_{good_cluster_key}"] = good_debut
    adata.obs[f"{output_key}_{bad_cluster_key}"] = bad_debut

    if landscape_key is None:
        landscape_key = f"{output_key}_landscape"
    adata.obs[landscape_key] = np.minimum(good_debut, bad_debut)

    if distance_key is not None:
        adata.obs[f"{distance_key}_{good_cluster_key}"] = distance_to_good
        adata.obs[f"{distance_key}_{bad_cluster_key}"] = distance_to_bad

    if terminal_class_key is not None:
        adata.obs[terminal_class_key] = terminal_label_by_observation

    params = {
        "method": (
            "terminal_state_time_layered_dp"
            if transition_mode == "time_layered"
            else "terminal_state_shortest_path"
        ),
        "transition_mode": transition_mode,
        "path_cost": path_cost,
        "seq_key": seq_key,
        "cluster_key": cluster_key,
        "good_cluster_key": good_cluster_key,
        "bad_cluster_key": bad_cluster_key,
        "cost_matrix_key": cost_matrix_key,
        "sequence_edge_mode": sequence_edge_mode,
        "terminal_indices": terminal_indices.copy(),
        "terminal_labels": terminal_labels.copy(),
    }
    if time_key is not None:
        params["time_key"] = time_key
    if distance_key is not None:
        params["distance_key"] = distance_key
    adata.uns[f"{output_key}_params"] = params

    if store_next_hop:
        adata.uns[f"{output_key}_next_hop_{good_cluster_key}"] = next_hop_good
        adata.uns[f"{output_key}_next_hop_{bad_cluster_key}"] = next_hop_bad


def eps_attracting_basin_shortest_path(
    adata,
    *,
    cluster_key: str = "cluster",
    good_cluster_key: str = "good",
    bad_cluster_key: str = "bad",
    cost_matrix_key: str = "cost_matrix",
    seq_key: str = "seq_id",
    time_key: str | None = None,
    sequence_edge_mode: SequenceEdgeMode = "all_forward",
    transition_mode: TransitionMode = "unrestricted",
    output_key: str = "eps_attracting_basin_sp",
    landscape_key: str | None = None,
    distance_key: str | None = None,
    terminal_class_key: str | None = None,
    graph_key: str | None = None,
    store_next_hop: bool = False,
) -> None:
    """Compute the minimax/bottleneck epsilon debut function."""
    terminal_state_shortest_path_debut(
        adata,
        path_cost="max",
        seq_key=seq_key,
        cluster_key=cluster_key,
        good_cluster_key=good_cluster_key,
        bad_cluster_key=bad_cluster_key,
        cost_matrix_key=cost_matrix_key,
        time_key=time_key,
        sequence_edge_mode=sequence_edge_mode,
        transition_mode=transition_mode,
        output_key=output_key,
        landscape_key=landscape_key,
        distance_key=distance_key,
        terminal_class_key=terminal_class_key,
        graph_key=graph_key,
        store_next_hop=store_next_hop,
    )


def eps_sum_attracting_basin_shortest_path(
    adata,
    *,
    cluster_key: str = "cluster",
    good_cluster_key: str = "good",
    bad_cluster_key: str = "bad",
    cost_matrix_key: str = "cost_matrix",
    seq_key: str = "seq_id",
    time_key: str | None = None,
    sequence_edge_mode: SequenceEdgeMode = "all_forward",
    transition_mode: TransitionMode = "unrestricted",
    output_key: str = "eps_sum_attracting_basin_sp",
    landscape_key: str | None = None,
    distance_key: str | None = None,
    terminal_class_key: str | None = None,
    graph_key: str | None = None,
    store_next_hop: bool = False,
) -> None:
    """Compute the additive epsilon_Sigma debut function."""
    terminal_state_shortest_path_debut(
        adata,
        path_cost="sum",
        seq_key=seq_key,
        cluster_key=cluster_key,
        good_cluster_key=good_cluster_key,
        bad_cluster_key=bad_cluster_key,
        cost_matrix_key=cost_matrix_key,
        time_key=time_key,
        sequence_edge_mode=sequence_edge_mode,
        transition_mode=transition_mode,
        output_key=output_key,
        landscape_key=landscape_key,
        distance_key=distance_key,
        terminal_class_key=terminal_class_key,
        graph_key=graph_key,
        store_next_hop=store_next_hop,
    )


def _external_time_layer_indices(
    adata,
    *,
    seq_key: str,
    time_key: str,
    time_value,
) -> tuple[np.ndarray, np.ndarray | None, bool]:
    _, times, index_matrix = _time_layer_layout(
        adata,
        seq_key=seq_key,
        time_key=time_key,
    )
    matches = np.flatnonzero(times == time_value)
    if matches.size != 1:
        raise ValueError(
            f"time_value={time_value!r} is not a unique value in adata.obs[{time_key!r}]."
        )
    t_pos = int(matches[0])
    current = index_matrix[:, t_pos]
    if t_pos == index_matrix.shape[1] - 1:
        return current, None, True
    return current, index_matrix[:, t_pos + 1], False


def minimum_epsilon_to_good_bad(
    adata,
    x,
    costs,
    *,
    output_key: str = "eps_attracting_basin_sp",
    time_value=None,
) -> ExternalStateEpsilon:
    """Evaluate Good/Bad reach costs for an external state.

    For ``transition_mode='unrestricted'`` this uses the usual one-entry
    extension to the precomputed graph distances.

    For ``transition_mode='time_layered'``, ``time_value`` is required. Only
    training states at that same time are eligible entry states. After the
    switch, the chosen member advances to the next time layer, so the future
    cost is read from that member's next-time state. External queries at the
    final time layer are not defined because no perturbation or further
    evolution is allowed there.
    """
    params_key = f"{output_key}_params"
    if params_key not in adata.uns:
        raise KeyError(
            f"{params_key!r} is missing. Compute the shortest-path debut "
            "functions first."
        )
    params = adata.uns[params_key]
    path_cost = params["path_cost"]
    transition_mode = params.get("transition_mode", "unrestricted")
    good_key = params.get("good_cluster_key", "good")
    bad_key = params.get("bad_cluster_key", "bad")

    x_array = np.asarray(x, dtype=float).reshape(-1)
    cost_array = np.asarray(costs, dtype=float).reshape(-1)
    if x_array.shape[0] != adata.X.shape[1]:
        raise ValueError("x has a different dimension from adata.X.")
    if cost_array.shape[0] != adata.n_obs:
        raise ValueError(f"costs must have length {adata.n_obs}.")
    if np.isnan(cost_array).any():
        raise ValueError("costs may not contain NaN.")
    if np.isneginf(cost_array).any():
        raise ValueError("costs may contain +inf, but not -inf.")
    finite = np.isfinite(cost_array)
    if np.any(cost_array[finite] < 0):
        raise ValueError("costs must be non-negative (or +inf).")

    d_good = np.maximum(
        np.asarray(adata.obs[f"{output_key}_{good_key}"], dtype=float),
        0.0,
    )
    d_bad = np.maximum(
        np.asarray(adata.obs[f"{output_key}_{bad_key}"], dtype=float),
        0.0,
    )

    if transition_mode == "unrestricted":
        if path_cost == "max":
            candidate_good = np.maximum(cost_array, d_good)
            candidate_bad = np.maximum(cost_array, d_bad)
        else:
            candidate_good = cost_array + d_good
            candidate_bad = cost_array + d_bad
        eps_good = float(np.min(candidate_good))
        eps_bad = float(np.min(candidate_bad))
        entry_good = -1 if np.isinf(eps_good) else int(np.argmin(candidate_good))
        entry_bad = -1 if np.isinf(eps_bad) else int(np.argmin(candidate_bad))
    else:
        time_key = params.get("time_key")
        if time_key is None:
            raise KeyError("time_layered output is missing time_key metadata.")
        if time_value is None:
            raise ValueError(
                "time_value is required for external queries in time_layered mode."
            )
        seq_key = params.get("seq_key", "seq_id")
        current, following, is_terminal = _external_time_layer_indices(
            adata,
            seq_key=seq_key,
            time_key=time_key,
            time_value=time_value,
        )

        if is_terminal:
            raise ValueError(
                "External queries are not defined at the final time layer in "
                "time_layered mode because no perturbation or further evolution "
                "is allowed."
            )
        else:
            entry_cost = cost_array[current]
            future_good = d_good[following]
            future_bad = d_bad[following]
            if path_cost == "max":
                candidate_good = np.maximum(entry_cost, future_good)
                candidate_bad = np.maximum(entry_cost, future_bad)
            else:
                candidate_good = entry_cost + future_good
                candidate_bad = entry_cost + future_bad
            eps_good = float(np.min(candidate_good))
            eps_bad = float(np.min(candidate_bad))
            entry_good = (
                -1
                if np.isinf(eps_good)
                else int(current[int(np.argmin(candidate_good))])
            )
            entry_bad = (
                -1
                if np.isinf(eps_bad)
                else int(current[int(np.argmin(candidate_bad))])
            )

    return ExternalStateEpsilon(
        epsilon_good=eps_good,
        epsilon_bad=eps_bad,
        entry_index_good=entry_good,
        entry_index_bad=entry_bad,
        path_cost=path_cost,
        good_cluster_key=good_key,
        bad_cluster_key=bad_key,
    )



def _target_terminal_indices_from_params(
    params,
    target_label: str,
) -> np.ndarray:
    terminal_indices = np.asarray(params["terminal_indices"], dtype=int)
    terminal_labels = np.asarray(params["terminal_labels"], dtype=object)
    result = terminal_indices[terminal_labels == target_label]
    if result.size == 0:
        raise ValueError(f"No terminal state is labelled {target_label!r}.")
    return result


def _ensure_next_hop_to_target(
    adata,
    *,
    output_key: str,
    target_label: str,
) -> np.ndarray:
    """Return/store next-hop data for one learned target cost."""
    params_key = f"{output_key}_params"
    if params_key not in adata.uns:
        raise KeyError(
            f"{params_key!r} is missing. Compute the shortest-path debut "
            "functions first."
        )
    params = adata.uns[params_key]
    next_key = f"{output_key}_next_hop_{target_label}"
    if next_key in adata.uns:
        next_hop = np.asarray(adata.uns[next_key], dtype=int)
        if next_hop.shape != (adata.n_obs,):
            raise ValueError(
                f"adata.uns[{next_key!r}] has shape {next_hop.shape}; "
                f"expected {(adata.n_obs,)}."
            )
        return next_hop

    path_cost = params["path_cost"]
    transition_mode = params.get("transition_mode", "unrestricted")
    seq_key = params.get("seq_key", "seq_id")
    cluster_key = params.get("cluster_key", "cluster")
    cost_matrix_key = params.get("cost_matrix_key", "cost_matrix")

    if transition_mode == "time_layered":
        time_key = params.get("time_key")
        if time_key is None:
            raise KeyError("time_layered output is missing time_key metadata.")
        _, next_hop = _time_layered_distance_to_target(
            adata,
            target_label=target_label,
            path_cost=path_cost,
            seq_key=seq_key,
            time_key=time_key,
            cluster_key=cluster_key,
            cost_matrix_key=cost_matrix_key,
        )
    elif transition_mode == "unrestricted":
        graph = build_terminal_state_graph(
            adata,
            seq_key=seq_key,
            cluster_key=cluster_key,
            good_cluster_key=params.get("good_cluster_key", "good"),
            bad_cluster_key=params.get("bad_cluster_key", "bad"),
            cost_matrix_key=cost_matrix_key,
            time_key=params.get("time_key"),
            sequence_edge_mode=params.get("sequence_edge_mode", "all_forward"),
            transition_mode="unrestricted",
        )
        targets = _target_terminal_indices_from_params(params, target_label)
        _, next_hop = _multi_source_distance_to_targets(
            graph.weights,
            targets,
            path_cost=path_cost,
        )
    else:
        raise ValueError(
            f"Unsupported transition_mode={transition_mode!r}."
        )

    next_hop = np.asarray(next_hop, dtype=int)
    adata.uns[next_key] = next_hop.copy()
    return next_hop


def _time_layer_positions(
    adata,
    *,
    seq_key: str,
    time_key: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return layout plus observation -> sequence/time-position maps."""
    sequence_ids, times, index_matrix = _time_layer_layout(
        adata,
        seq_key=seq_key,
        time_key=time_key,
    )
    obs_to_seq = np.full(adata.n_obs, -1, dtype=int)
    obs_to_time = np.full(adata.n_obs, -1, dtype=int)
    for seq_pos in range(index_matrix.shape[0]):
        for time_pos in range(index_matrix.shape[1]):
            obs_index = int(index_matrix[seq_pos, time_pos])
            obs_to_seq[obs_index] = seq_pos
            obs_to_time[obs_index] = time_pos
    if np.any(obs_to_seq < 0) or np.any(obs_to_time < 0):
        raise RuntimeError("Could not map all observations into time layers.")
    return sequence_ids, times, index_matrix, obs_to_seq, obs_to_time


def _time_layered_suffix_steps(
    adata,
    *,
    start_index: int,
    target_label: str,
    output_key: str,
    next_hop: np.ndarray,
) -> tuple[list[PathStep], int]:
    """Expand stored DP next hops into explicit switch/advance steps."""
    params = adata.uns[f"{output_key}_params"]
    seq_key = params.get("seq_key", "seq_id")
    time_key = params.get("time_key")
    cluster_key = params.get("cluster_key", "cluster")
    cost_matrix_key = params.get("cost_matrix_key", "cost_matrix")
    if time_key is None:
        raise KeyError("time_layered output is missing time_key metadata.")

    _, _, index_matrix, obs_to_seq, obs_to_time = _time_layer_positions(
        adata,
        seq_key=seq_key,
        time_key=time_key,
    )
    cost_matrix = _get_cost_matrix(adata, cost_matrix_key)
    labels = np.asarray(adata.obs[cluster_key], dtype=object)

    current = int(start_index)
    steps: list[PathStep] = []
    seen: set[int] = set()

    while True:
        if current in seen:
            raise RuntimeError("Cycle detected while reconstructing time-layered path.")
        seen.add(current)

        t_pos = int(obs_to_time[current])
        if t_pos == index_matrix.shape[1] - 1:
            if labels[current] == target_label:
                return steps, current
            # A finite optimal path cannot change class after reaching the final
            # time layer, because no perturbation or transition is allowed there.
            return [], -1

        future = int(next_hop[current])
        if future == -1:
            return [], -1
        if int(obs_to_time[future]) != t_pos + 1:
            raise RuntimeError(
                "A time-layered next_hop did not advance exactly one time layer."
            )

        selected_seq_pos = int(obs_to_seq[future])
        switched_current = int(index_matrix[selected_seq_pos, t_pos])
        if switched_current != current:
            steps.append(
                PathStep(
                    kind="switch",
                    source_index=current,
                    destination_index=switched_current,
                    cost=float(cost_matrix[current, switched_current]),
                )
            )

        steps.append(
            PathStep(
                kind="advance",
                source_index=switched_current,
                destination_index=future,
                cost=0.0,
            )
        )
        current = future


def _unrestricted_suffix_steps(
    adata,
    *,
    start_index: int,
    target_label: str,
    output_key: str,
    next_hop: np.ndarray,
) -> tuple[list[PathStep], int]:
    params = adata.uns[f"{output_key}_params"]
    targets = _target_terminal_indices_from_params(params, target_label)
    vertex_path = reconstruct_terminal_path(start_index, next_hop, targets)
    if not vertex_path:
        return [], -1

    graph = build_terminal_state_graph(
        adata,
        seq_key=params.get("seq_key", "seq_id"),
        cluster_key=params.get("cluster_key", "cluster"),
        good_cluster_key=params.get("good_cluster_key", "good"),
        bad_cluster_key=params.get("bad_cluster_key", "bad"),
        cost_matrix_key=params.get("cost_matrix_key", "cost_matrix"),
        time_key=params.get("time_key"),
        sequence_edge_mode=params.get("sequence_edge_mode", "all_forward"),
        transition_mode="unrestricted",
    )
    seq_key = params.get("seq_key", "seq_id")
    seq_values = np.asarray(adata.obs[seq_key])

    steps: list[PathStep] = []
    for source, destination in zip(vertex_path[:-1], vertex_path[1:]):
        weight = float(graph.weights[source, destination])
        kind = (
            "advance"
            if seq_values[source] == seq_values[destination] and weight == 0.0
            else "graph"
        )
        steps.append(
            PathStep(
                kind=kind,
                source_index=int(source),
                destination_index=int(destination),
                cost=weight,
            )
        )
    return steps, int(vertex_path[-1])


def _aggregate_path_cost(steps: Sequence[PathStep], path_cost: PathCost) -> float:
    costs = [float(step.cost) for step in steps]
    if not costs:
        return 0.0
    if path_cost == "sum":
        return float(np.sum(costs))
    if path_cost == "max":
        return float(np.max(costs))
    raise ValueError("path_cost must be either 'sum' or 'max'.")


def _build_external_optimal_path(
    adata,
    *,
    costs: np.ndarray,
    output_key: str,
    time_value,
    target_label: str,
    epsilon: float,
    entry_index: int,
    atol: float,
) -> OptimalPath:
    params = adata.uns[f"{output_key}_params"]
    path_cost: PathCost = params["path_cost"]
    transition_mode = params.get("transition_mode", "unrestricted")

    if np.isinf(epsilon) or entry_index < 0:
        return OptimalPath(
            target_label=target_label,
            path_cost=path_cost,
            value=float(epsilon),
            steps=tuple(),
            terminal_index=-1,
        )

    next_hop = _ensure_next_hop_to_target(
        adata,
        output_key=output_key,
        target_label=target_label,
    )
    steps: list[PathStep] = [
        PathStep(
            kind="entry",
            source_index=None,
            destination_index=int(entry_index),
            cost=float(costs[int(entry_index)]),
        )
    ]

    if transition_mode == "time_layered":
        time_key = params.get("time_key")
        seq_key = params.get("seq_key", "seq_id")
        if time_key is None:
            raise KeyError("time_layered output is missing time_key metadata.")
        if time_value is None:
            raise ValueError(
                "time_value is required for external queries in time_layered mode."
            )
        current, following, is_terminal = _external_time_layer_indices(
            adata,
            seq_key=seq_key,
            time_key=time_key,
            time_value=time_value,
        )
        if int(entry_index) not in set(int(v) for v in current):
            raise RuntimeError(
                "The selected external entry is not in the requested time layer."
            )

        if is_terminal:
            raise ValueError(
                "External path reconstruction is not defined at the final time "
                "layer in time_layered mode because no perturbation or further "
                "evolution is allowed."
            )
        else:
            entry_pos = int(np.flatnonzero(current == int(entry_index))[0])
            first_future = int(following[entry_pos])
            steps.append(
                PathStep(
                    kind="advance",
                    source_index=int(entry_index),
                    destination_index=first_future,
                    cost=0.0,
                )
            )
            suffix, terminal_index = _time_layered_suffix_steps(
                adata,
                start_index=first_future,
                target_label=target_label,
                output_key=output_key,
                next_hop=next_hop,
            )
            if terminal_index < 0:
                raise RuntimeError(
                    f"Could not reconstruct a finite path to {target_label!r}."
                )
            steps.extend(suffix)
    elif transition_mode == "unrestricted":
        suffix, terminal_index = _unrestricted_suffix_steps(
            adata,
            start_index=int(entry_index),
            target_label=target_label,
            output_key=output_key,
            next_hop=next_hop,
        )
        if terminal_index < 0:
            raise RuntimeError(
                f"Could not reconstruct a finite path to {target_label!r}."
            )
        steps.extend(suffix)
    else:
        raise ValueError(
            f"Unsupported transition_mode={transition_mode!r}."
        )

    realized = _aggregate_path_cost(steps, path_cost)
    if not np.isclose(realized, epsilon, atol=atol, rtol=1e-9):
        raise RuntimeError(
            f"Reconstructed path to {target_label!r} has value {realized}, "
            f"but the external reach cost is {epsilon}."
        )

    return OptimalPath(
        target_label=target_label,
        path_cost=path_cost,
        value=float(epsilon),
        steps=tuple(steps),
        terminal_index=int(terminal_index),
    )


def minimum_epsilon_and_paths_to_good_bad(
    adata,
    x,
    costs,
    *,
    output_key: str = "eps_attracting_basin_sp",
    time_value=None,
    atol: float = 1e-8,
) -> ExternalStatePathResult:
    """Return C_G,p(x), C_B,p(x), and one optimal path realizing each cost.

    The debut/reach functions must already have been learned with
    :func:`eps_attracting_basin_shortest_path` (p=infinity) or
    :func:`eps_sum_attracting_basin_shortest_path` (p=1).

    Parameters
    ----------
    adata
        Ensemble trajectories used to learn the debut functions.
    x
        External state. It is used for dimension validation; the supplied
        ``costs`` define the actual entry costs to ensemble observations.
    costs
        Length-``adata.n_obs`` array of costs from ``x`` to each ensemble
        observation. In ``time_layered`` mode only the values at ``time_value``
        are used; other entries may be ``+inf``.
    output_key
        Output key used when the debut function was learned.
    time_value
        Required for ``transition_mode='time_layered'``.
    atol
        Absolute tolerance for checking that each reconstructed path realizes
        the reported optimum.

    Notes
    -----
    In time-layered mode an optimal path is expanded explicitly into
    ``entry`` / optional same-time ``switch`` / zero-cost ``advance`` steps.
    Thus the returned path directly represents the Bellman recursion used to
    compute C_G,p and C_B,p.
    """
    base = minimum_epsilon_to_good_bad(
        adata,
        x,
        costs,
        output_key=output_key,
        time_value=time_value,
    )
    cost_array = np.asarray(costs, dtype=float).reshape(-1)

    path_good = _build_external_optimal_path(
        adata,
        costs=cost_array,
        output_key=output_key,
        time_value=time_value,
        target_label=base.good_cluster_key,
        epsilon=base.epsilon_good,
        entry_index=base.entry_index_good,
        atol=atol,
    )
    path_bad = _build_external_optimal_path(
        adata,
        costs=cost_array,
        output_key=output_key,
        time_value=time_value,
        target_label=base.bad_cluster_key,
        epsilon=base.epsilon_bad,
        entry_index=base.entry_index_bad,
        atol=atol,
    )

    return ExternalStatePathResult(
        epsilon_good=base.epsilon_good,
        epsilon_bad=base.epsilon_bad,
        entry_index_good=base.entry_index_good,
        entry_index_bad=base.entry_index_bad,
        path_cost=base.path_cost,
        good_cluster_key=base.good_cluster_key,
        bad_cluster_key=base.bad_cluster_key,
        path_good=path_good,
        path_bad=path_bad,
    )


def reconstruct_terminal_path(
    start_index: int,
    next_hop: np.ndarray,
    target_indices: Sequence[int],
) -> list[int]:
    """Reconstruct one stored optimal path from a vertex to a target terminal."""
    start_index = int(start_index)
    if start_index < 0 or start_index >= len(next_hop):
        raise IndexError("start_index is outside the graph.")

    targets = {int(index) for index in target_indices}
    path = [start_index]
    current = start_index
    seen = {current}

    while current not in targets:
        following = int(next_hop[current])
        if following == -1:
            return []
        current = following
        if current in seen:
            raise RuntimeError("Cycle detected in next_hop data.")
        path.append(current)
        seen.add(current)
    return path
