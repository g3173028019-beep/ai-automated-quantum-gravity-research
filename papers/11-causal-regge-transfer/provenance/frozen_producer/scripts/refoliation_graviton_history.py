"""Exact refoliation certificate and Euclidean boundary response for PCRG.

The same four-simplex complex can be read in two ways:

* as one 4-dimensional tent move with a single internal edge; or
* as a 3-dimensional Pachner history in which a 2->3 move creates that edge,
  a 1->4 and a distinct 2->3 change its neighbourhood, and a later 3->2
  removes it.

The combinatorial part below is exact.  The response calculation is a finite-
difference calibration of the Euclidean Regge Hamilton principal function.  It
is not a Lorentzian retarded propagator.
"""

from collections import Counter
from itertools import combinations

import numpy as np
from scipy.optimize import root_scalar

try:
    from scripts.first_order_angle_measure import (
        initial_first_order_action_at_geometric_saddle,
    )
    from scripts.regge_four_two_action import (
        INITIAL_SIMPLICES,
        canonical_radon_configuration,
        regge_action_with_boundary,
        squared_lengths,
    )
    from scripts.sourced_regge_dynamics import vacuum_regge_equation
except ModuleNotFoundError:  # Allow ``python scripts/refoliation_graviton_history.py``.
    from first_order_angle_measure import initial_first_order_action_at_geometric_saddle
    from regge_four_two_action import (
        INITIAL_SIMPLICES,
        canonical_radon_configuration,
        regge_action_with_boundary,
        squared_lengths,
    )
    from sourced_regge_dynamics import vacuum_regge_equation


Tetrahedron = tuple[int, int, int, int]
Complex3 = set[Tetrahedron]

INITIAL_SLICE: Complex3 = {(0, 1, 2, 3), (0, 1, 2, 4)}
GRAVITON_EDGE = (3, 4)
LATER_EDGE = (2, 5)

# Map the old A7 labels to the labels of the Pachner history.  In particular,
# the old internal/tent edge (01) becomes the history edge (34).
A7_TO_HISTORY = {0: 3, 1: 4, 2: 0, 3: 1, 4: 2, 5: 5}


def normalize(simplex) -> tuple[int, ...]:
    return tuple(sorted(simplex))


def faces(simplex: tuple[int, ...], size: int) -> set[tuple[int, ...]]:
    return {normalize(face) for face in combinations(simplex, size)}


def complex_faces(complex_: Complex3, size: int) -> set[tuple[int, ...]]:
    return set().union(*(faces(tetrahedron, size) for tetrahedron in complex_))


def boundary_triangles(complex_: Complex3) -> set[tuple[int, int, int]]:
    counts = Counter(
        triangle for tetrahedron in complex_ for triangle in faces(tetrahedron, 3)
    )
    if any(count > 2 for count in counts.values()):
        raise ValueError("non-manifold triangle incidence")
    return {triangle for triangle, count in counts.items() if count == 1}


def euler_characteristic(complex_: Complex3) -> int:
    return sum(
        (-1) ** (size - 1) * len(complex_faces(complex_, size))
        for size in range(1, 5)
    )


def replace_patch(
    complex_: Complex3, old_patch: Complex3, new_patch: Complex3
) -> Complex3:
    if not old_patch <= complex_:
        raise ValueError("old Pachner patch is not present")
    if boundary_triangles(old_patch) != boundary_triangles(new_patch):
        raise ValueError("replacement does not preserve the patch boundary")
    unchanged = complex_ - old_patch
    if unchanged & new_patch:
        raise ValueError("replacement would duplicate a tetrahedron")
    result = unchanged | new_patch
    boundary_triangles(result)  # Also checks triangle manifold incidence.
    if euler_characteristic(result) != 1:
        raise ValueError("the local hypersurface ceased to be a three-ball")
    return result


def refoliated_history_certificate() -> dict:
    """Return all exact hypersurfaces in the non-inverse graviton history."""
    states = [set(INITIAL_SLICE)]

    # 2->3 on (0123),(0124): create edge (34).
    states.append(
        replace_patch(
            states[-1],
            {(0, 1, 2, 3), (0, 1, 2, 4)},
            {(0, 1, 3, 4), (0, 2, 3, 4), (1, 2, 3, 4)},
        )
    )

    # 1->4 on (0134): introduce vertex 5 (lapse/shift sector).
    states.append(
        replace_patch(
            states[-1],
            {(0, 1, 3, 4)},
            {
                (0, 1, 3, 5),
                (0, 1, 4, 5),
                (0, 3, 4, 5),
                (1, 3, 4, 5),
            },
        )
    )

    # A different 2->3 move creates (25) and changes the star of (34).
    states.append(
        replace_patch(
            states[-1],
            {(1, 2, 3, 4), (1, 3, 4, 5)},
            {(1, 2, 3, 5), (1, 2, 4, 5), (2, 3, 4, 5)},
        )
    )

    # 3->2 on the now-valence-three edge (34), not the inverse of the first move.
    states.append(
        replace_patch(
            states[-1],
            {(0, 2, 3, 4), (0, 3, 4, 5), (2, 3, 4, 5)},
            {(0, 2, 3, 5), (0, 2, 4, 5)},
        )
    )

    bulk_simplices = {
        (0, 1, 2, 3, 4),
        (0, 1, 3, 4, 5),
        (1, 2, 3, 4, 5),
        (0, 2, 3, 4, 5),
    }
    tetrahedron_counts = Counter(
        tetrahedron
        for simplex in bulk_simplices
        for tetrahedron in faces(simplex, 4)
    )
    bulk_boundary = {
        tetrahedron for tetrahedron, count in tetrahedron_counts.items() if count == 1
    }
    boundary_edges = complex_faces(bulk_boundary, 2)
    bulk_edges = complex_faces({normalize(s) for s in bulk_simplices}, 2)
    internal_edges = bulk_edges - boundary_edges
    internal_triangles = {
        triangle
        for triangle, count in Counter(
            triangle
            for simplex in bulk_simplices
            for triangle in faces(simplex, 3)
        ).items()
        if count == 3
    }

    relabeled_a7 = {
        normalize(A7_TO_HISTORY[vertex] for vertex in simplex)
        for simplex in INITIAL_SIMPLICES
    }
    return {
        "states": states,
        "state_edges": [complex_faces(state, 2) for state in states],
        "state_euler_characteristics": [euler_characteristic(state) for state in states],
        "state_boundary_triangles": [boundary_triangles(state) for state in states],
        "bulk_simplices": bulk_simplices,
        "bulk_boundary_tetrahedra": bulk_boundary,
        "bulk_internal_edges": internal_edges,
        "bulk_internal_triangles": internal_triangles,
        "relabeled_a7_simplices": relabeled_a7,
        "isomorphic_to_a7": relabeled_a7 == bulk_simplices,
    }


def _varied_lengths(
    background: np.ndarray, edges: tuple[tuple[int, int], ...], values: np.ndarray
) -> np.ndarray:
    varied = np.array(background, copy=True)
    for edge, value in zip(edges, values):
        varied[edge] = varied[edge[::-1]] = value
    return varied


def numerical_hessian(function, point: np.ndarray, relative_step: float) -> np.ndarray:
    """Symmetric second-order central Hessian with independent coordinate steps."""
    point = np.asarray(point, dtype=float)
    steps = relative_step * np.maximum(1.0, np.abs(point))
    hessian = np.zeros((len(point), len(point)))
    center = function(point)
    for first in range(len(point)):
        first_step = np.zeros_like(point)
        first_step[first] = steps[first]
        hessian[first, first] = (
            function(point + first_step)
            - 2.0 * center
            + function(point - first_step)
        ) / steps[first] ** 2
        for second in range(first):
            second_step = np.zeros_like(point)
            second_step[second] = steps[second]
            hessian[first, second] = hessian[second, first] = (
                function(point + first_step + second_step)
                - function(point + first_step - second_step)
                - function(point - first_step + second_step)
                + function(point - first_step - second_step)
            ) / (4.0 * steps[first] * steps[second])
    return hessian


def boundary_transfer_certificate(relative_step: float = 2e-4) -> dict:
    """Euclidean initial-to-later mixed response after eliminating the internal edge.

    History edges ``(34),(01),(25)`` map respectively to old A7 edges
    ``(01),(23),(45)``.  The first is internal, ``(01)`` is present on the
    initial slice, and ``(25)`` is created only by the second 2->3 move.  No
    four-simplex contains both boundary edges, so their fixed-internal mixed
    Hessian vanishes; the effective response is mediated by the internal edge.
    """
    background = squared_lengths(canonical_radon_configuration())
    old_edges = ((0, 1), (2, 3), (4, 5))
    point = np.asarray([background[edge] for edge in old_edges])

    def regge(values: np.ndarray) -> float:
        lengths = _varied_lengths(background, old_edges, values)
        return regge_action_with_boundary(lengths, values[0])

    def first_order(values: np.ndarray) -> float:
        lengths = _varied_lengths(background, old_edges, values)
        return initial_first_order_action_at_geometric_saddle(lengths)

    regge_hessian = numerical_hessian(regge, point, relative_step)
    first_order_hessian = numerical_hessian(first_order, point, relative_step)

    def eliminate_internal(hessian: np.ndarray) -> np.ndarray:
        return hessian[1:, 1:] - np.outer(hessian[1:, 0], hessian[0, 1:]) / hessian[0, 0]

    regge_effective = eliminate_internal(regge_hessian)
    first_order_effective = eliminate_internal(first_order_hessian)

    boundary_point = point[1:]
    boundary_steps = relative_step * np.maximum(1.0, np.abs(boundary_point))

    def on_shell_action(boundary_values: np.ndarray) -> float:
        lengths = _varied_lengths(background, old_edges[1:], boundary_values)
        solution = root_scalar(
            lambda internal: vacuum_regge_equation(lengths, internal),
            x0=point[0],
            x1=point[0] * (1.0 + 1e-5),
            method="secant",
            xtol=1e-13,
            maxiter=30,
        )
        if not solution.converged:
            raise RuntimeError("near-flat internal Regge equation did not converge")
        return regge_action_with_boundary(lengths, solution.root)

    first_step = np.array((boundary_steps[0], 0.0))
    second_step = np.array((0.0, boundary_steps[1]))
    on_shell_mixed_response = (
        on_shell_action(boundary_point + first_step + second_step)
        - on_shell_action(boundary_point + first_step - second_step)
        - on_shell_action(boundary_point - first_step + second_step)
        + on_shell_action(boundary_point - first_step - second_step)
    ) / (4.0 * boundary_steps[0] * boundary_steps[1])
    return {
        "history_edges": (GRAVITON_EDGE, (0, 1), LATER_EDGE),
        "a7_edges": old_edges,
        "point": point,
        "regge_hessian": regge_hessian,
        "first_order_hessian": first_order_hessian,
        "regge_effective_boundary_hessian": regge_effective,
        "first_order_effective_boundary_hessian": first_order_effective,
        "initial_to_later_mixed_response": float(regge_effective[0, 1]),
        "fixed_internal_direct_response": float(regge_hessian[1, 2]),
        "exact_separated_transfer": 3.0 / 3584.0,
        "on_shell_mixed_response": float(on_shell_mixed_response),
        "action_value_residual": abs(regge(point) - first_order(point)),
    }


if __name__ == "__main__":
    history = refoliated_history_certificate()
    response = boundary_transfer_certificate()
    print(f"state_tetrahedron_counts={[len(state) for state in history['states']]}")
    print(f"state_euler_characteristics={history['state_euler_characteristics']}")
    print(
        "common_patch_boundary="
        f"{all(boundary == history['state_boundary_triangles'][0] for boundary in history['state_boundary_triangles'])}"
    )
    print(f"edge_34_presence={[GRAVITON_EDGE in edges for edges in history['state_edges']]}")
    print(f"edge_25_presence={[LATER_EDGE in edges for edges in history['state_edges']]}")
    print(f"bulk_internal_edges={sorted(history['bulk_internal_edges'])}")
    print(f"bulk_internal_triangles={sorted(history['bulk_internal_triangles'])}")
    print(f"isomorphic_to_a7={history['isomorphic_to_a7']}")
    print(
        "initial_to_later_mixed_response="
        f"{response['initial_to_later_mixed_response']:.16g}"
    )
    print(
        "first_order_response="
        f"{response['first_order_effective_boundary_hessian'][0, 1]:.16g}"
    )
    print(f"on_shell_response={response['on_shell_mixed_response']:.16g}")
    print(f"action_value_residual={response['action_value_residual']:.3e}")
