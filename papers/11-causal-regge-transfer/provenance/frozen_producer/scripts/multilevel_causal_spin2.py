"""Nested longitudinal refinement test for the causal C6 spin-two carrier.

The coarse C5*C6 star has one moving q vertex on each spacelike cap.  A level
``m`` uniformly subdivides the time-reflected cap-edge pair into ``m`` pieces.
Thus m=2 reproduces the C7*C6 witness and m=4 gives its genuine nested C11*C6
refinement.  Exact D6 symmetry lets us differentiate only the physical radial
bulk directions and the spin-two boundary directions.  This is a symmetry
reduction of the Hessian, not a fitted truncation.
"""

from itertools import product

import numpy as np

try:
    from scripts.causal_join_star_transfer import time_reflection_symmetric_points
    from scripts.complex_lorentzian_regge import (
        MINKOWSKI_METRIC,
        richardson_complex_hessian,
        squared_intervals,
        varied_lengths,
    )
    from scripts.lorentzian_weyl_star import (
        _general_outward_normal,
        _topology_for,
        general_complex_action,
    )
except ModuleNotFoundError:
    from causal_join_star_transfer import time_reflection_symmetric_points
    from complex_lorentzian_regge import (
        MINKOWSKI_METRIC,
        richardson_complex_hessian,
        squared_intervals,
        varied_lengths,
    )
    from lorentzian_weyl_star import (
        _general_outward_normal,
        _topology_for,
        general_complex_action,
    )


def _path_edges(path: tuple[int, ...]) -> tuple[tuple[int, int], ...]:
    return tuple(
        tuple(sorted((path[index], path[index + 1])))
        for index in range(len(path) - 1)
    )


def _relative_edges(
    path: tuple[int, ...],
    interior_q: tuple[int, ...],
    r_vertices: tuple[int, ...],
) -> tuple[tuple[int, int], ...]:
    return tuple(
        sorted(
            set(
                _path_edges(path)
                + tuple(
                    tuple(sorted((q_vertex, r_vertex)))
                    for q_vertex in interior_q
                    for r_vertex in r_vertices
                )
            )
        )
    )


def refinement_geometry(segment_count: int) -> dict:
    """Build a nested cap-edge subdivision of the fixed coarse witness."""
    if segment_count < 1:
        raise ValueError("segment_count must be positive")
    coarse = time_reflection_symmetric_points()
    q0, q1, q2, q3, q4 = coarse[1:6]
    initial_insertions = [
        (1.0 - j / segment_count) * q0 + (j / segment_count) * q1
        for j in range(1, segment_count)
    ]
    final_insertions = [
        (1.0 - j / segment_count) * q4 + (j / segment_count) * q0
        for j in range(1, segment_count)
    ]
    q_points = np.asarray(
        [q0, *initial_insertions, q1, q2, q3, q4, *final_insertions]
    )
    points = np.vstack((np.zeros((1, 4)), q_points, coarse[6:]))
    q_vertices = tuple(range(1, len(q_points) + 1))
    r_vertices = tuple(range(len(q_points) + 1, len(q_points) + 7))
    q_cycle_edges = tuple(
        tuple(sorted((q_vertices[index], q_vertices[(index + 1) % len(q_vertices)])))
        for index in range(len(q_vertices))
    )
    r_cycle_edges = tuple(
        tuple(sorted((r_vertices[index], r_vertices[(index + 1) % 6])))
        for index in range(6)
    )
    boundary_edges = tuple(
        sorted(
            set(
                q_cycle_edges
                + r_cycle_edges
                + tuple(
                    tuple(sorted(edge)) for edge in product(q_vertices, r_vertices)
                )
            )
        )
    )
    internal_edges = tuple((0, vertex) for vertex in (*q_vertices, *r_vertices))
    simplices = tuple(
        tuple(
            sorted(
                (
                    0,
                    q_vertices[q_index],
                    q_vertices[(q_index + 1) % len(q_vertices)],
                    r_vertices[r_index],
                    r_vertices[(r_index + 1) % 6],
                )
            )
        )
        for q_index in range(len(q_vertices))
        for r_index in range(6)
    )
    q0_index = 1
    q1_index = segment_count + 1
    q2_index = segment_count + 2
    q3_index = segment_count + 3
    q4_index = segment_count + 4
    final_insertions_indices = tuple(range(segment_count + 5, 2 * segment_count + 4))
    initial_path = tuple(range(q0_index, q2_index + 1))
    final_path = (q3_index, q4_index, *final_insertions_indices, q0_index)
    initial_interior_q = tuple(range(2, q1_index + 1))
    final_interior_q = (*reversed(final_insertions_indices), q4_index)
    initial_relative_edges = _relative_edges(
        initial_path, initial_interior_q, r_vertices
    )
    final_relative_edges = _relative_edges(final_path, final_interior_q, r_vertices)
    return {
        "segment_count": segment_count,
        "points": points,
        "q_vertices": q_vertices,
        "r_vertices": r_vertices,
        "q_cycle_edges": q_cycle_edges,
        "r_cycle_edges": r_cycle_edges,
        "boundary_edges": boundary_edges,
        "internal_edges": internal_edges,
        "simplices": simplices,
        "initial_path": initial_path,
        "final_path": final_path,
        "initial_interior_q": initial_interior_q,
        "final_interior_q": final_interior_q,
        "initial_relative_edges": initial_relative_edges,
        "final_relative_edges": final_relative_edges,
    }


def _multi_vertex_normal_basis(
    points: np.ndarray,
    relative_edges: tuple[tuple[int, int], ...],
    moving_vertices: tuple[int, ...],
) -> tuple[np.ndarray, np.ndarray, int]:
    tangent = np.zeros((len(relative_edges), 4 * len(moving_vertices)))
    column_vertex = {vertex: index for index, vertex in enumerate(moving_vertices)}
    for row, edge in enumerate(relative_edges):
        difference = points[edge[0]] - points[edge[1]]
        covector = 2.0 * difference @ MINKOWSKI_METRIC
        for endpoint, sign in ((edge[0], 1.0), (edge[1], -1.0)):
            if endpoint in column_vertex:
                start = 4 * column_vertex[endpoint]
                tangent[row, start : start + 4] = sign * covector
    left, singular_values, _ = np.linalg.svd(tangent, full_matrices=True)
    rank = int(np.count_nonzero(singular_values > singular_values[0] * 1e-10))
    return left[:, rank:], singular_values, rank


def _spin_two_edge_basis(
    points: np.ndarray,
    relative_edges: tuple[tuple[int, int], ...],
    interior_q: tuple[int, ...],
    r_vertices: tuple[int, ...],
    normal_basis: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, float]:
    """Return orthonormal edge representatives and the raw-to-basis transform."""
    raw = np.zeros((len(relative_edges), 2 * len(interior_q)))
    q_position = {vertex: index for index, vertex in enumerate(interior_q)}
    for row, edge in enumerate(relative_edges):
        q_candidates = [vertex for vertex in edge if vertex in q_position]
        r_candidates = [vertex for vertex in edge if vertex in r_vertices]
        if len(q_candidates) != 1 or len(r_candidates) != 1:
            continue
        site = q_position[q_candidates[0]]
        r_vertex = r_candidates[0]
        angle = np.arctan2(points[r_vertex, 2], points[r_vertex, 1])
        raw[row, 2 * site : 2 * site + 2] = (
            np.cos(2.0 * angle),
            np.sin(2.0 * angle),
        )
    coordinates = normal_basis.T @ raw
    orthonormal, upper = np.linalg.qr(coordinates)
    width = raw.shape[1]
    if np.min(np.abs(np.diag(upper[:width, :width]))) < 1e-10:
        raise ArithmeticError("spin-two samples collapsed in the cap quotient")
    projected = normal_basis @ normal_basis.T @ raw
    residual = float(np.linalg.norm(raw - projected) / np.linalg.norm(raw))
    return normal_basis @ orthonormal[:, :width], upper[:width, :], residual


def _boundary_classification(geometry: dict) -> dict:
    points = geometry["points"]
    coordinates = {index: point for index, point in enumerate(points)}
    counts = {"past_spacelike": 0, "future_spacelike": 0, "timelike_side": 0}
    margins = []
    for q_edge in geometry["q_cycle_edges"]:
        for r_edge in geometry["r_cycle_edges"]:
            normal = _general_outward_normal(
                tuple(sorted(q_edge + r_edge)), 0, coordinates
            )
            norm_squared = float(normal @ MINKOWSKI_METRIC @ normal)
            margins.append(abs(norm_squared))
            if norm_squared < 0.0:
                kind = "future_spacelike" if normal[0] > 0.0 else "past_spacelike"
            else:
                kind = "timelike_side"
            counts[kind] += 1
    return {"counts": counts, "minimum_normal_margin": min(margins)}


def _coarse_injection(segment_count: int) -> np.ndarray:
    injection = np.zeros((2 * segment_count, 2))
    for site in range(segment_count):
        injection[2 * site : 2 * site + 2] = (
            (site + 1.0) / segment_count
        ) * np.eye(2)
    return injection


def multilevel_spin_two_certificate(
    segment_count: int = 4, hessian_step: float = 1.8e-4
) -> dict:
    geometry = refinement_geometry(segment_count)
    points = geometry["points"]
    simplices = geometry["simplices"]
    internal_edges = geometry["internal_edges"]
    initial_edges = geometry["initial_relative_edges"]
    final_edges = geometry["final_relative_edges"]
    topology = _topology_for(simplices)
    background = squared_intervals(points)

    initial_normal, initial_tangent_singular_values, initial_tangent_rank = (
        _multi_vertex_normal_basis(
            points, initial_edges, geometry["initial_interior_q"]
        )
    )
    final_normal, final_tangent_singular_values, final_tangent_rank = (
        _multi_vertex_normal_basis(points, final_edges, geometry["final_interior_q"])
    )
    initial_spin_edges, initial_raw_transform, initial_spin_residual = (
        _spin_two_edge_basis(
            points,
            initial_edges,
            geometry["initial_interior_q"],
            geometry["r_vertices"],
            initial_normal,
        )
    )
    final_spin_edges, final_raw_transform, final_spin_residual = (
        _spin_two_edge_basis(
            points,
            final_edges,
            geometry["final_interior_q"],
            geometry["r_vertices"],
            final_normal,
        )
    )

    central_tangent = np.asarray(
        [
            2.0 * (points[first] - points[second]) @ MINKOWSKI_METRIC
            for first, second in internal_edges
        ]
    )
    gauge_left, gauge_singular_values, _ = np.linalg.svd(
        central_tangent, full_matrices=True
    )
    physical_internal = gauge_left[:, 4:]
    physical_internal_count = physical_internal.shape[1]
    spin_count = 2 * segment_count
    origin = np.zeros(physical_internal_count + 2 * spin_count)
    internal_background = np.asarray([background[edge] for edge in internal_edges])
    initial_background = np.asarray([background[edge] for edge in initial_edges])
    final_background = np.asarray([background[edge] for edge in final_edges])

    def lengths_from_coordinates(coordinates: np.ndarray) -> np.ndarray:
        lengths = varied_lengths(
            background,
            internal_edges,
            internal_background
            + physical_internal @ coordinates[:physical_internal_count],
        )
        initial_start = physical_internal_count
        final_start = initial_start + spin_count
        lengths = varied_lengths(
            lengths,
            initial_edges,
            initial_background
            + initial_spin_edges @ coordinates[initial_start:final_start],
        )
        return varied_lengths(
            lengths,
            final_edges,
            final_background + final_spin_edges @ coordinates[final_start:],
        )

    def action(coordinates: np.ndarray) -> complex:
        return general_complex_action(
            lengths_from_coordinates(coordinates),
            simplices,
            topology["bulk_internal_triangles"],
        )

    hessian = richardson_complex_hessian(action, origin, hessian_step)
    internal_hessian = hessian[:physical_internal_count, :physical_internal_count]
    internal_boundary = hessian[:physical_internal_count, physical_internal_count:]
    effective = (
        hessian[physical_internal_count:, physical_internal_count:]
        - internal_boundary.T
        @ np.linalg.solve(internal_hessian, internal_boundary)
    )
    spin_initial = effective[:spin_count, :spin_count]
    spin_transfer = effective[:spin_count, spin_count:]
    spin_final = effective[spin_count:, spin_count:]
    transfer_left, transfer_singular_values, transfer_right_h = np.linalg.svd(
        spin_transfer, full_matrices=True
    )
    transfer_rank = int(
        np.count_nonzero(transfer_singular_values > transfer_singular_values[0] * 1e-5)
    )
    surviving_initial = transfer_left[:, :transfer_rank]
    surviving_final = transfer_right_h.T[:, :transfer_rank]
    constrained_initial = transfer_left[:, transfer_rank:]
    constrained_final = transfer_right_h.T[:, transfer_rank:]
    injection = _coarse_injection(segment_count)
    injected_initial = np.linalg.qr(initial_raw_transform @ injection)[0][:, :2]
    injected_final = np.linalg.qr(final_raw_transform @ injection)[0][:, :2]
    initial_overlap = np.linalg.svd(
        surviving_initial.T @ injected_initial, compute_uv=False
    )
    final_overlap = np.linalg.svd(
        surviving_final.T @ injected_final, compute_uv=False
    )

    surviving_initial_self = surviving_initial.T @ spin_initial @ surviving_initial
    surviving_transfer = surviving_initial.T @ spin_transfer @ surviving_final
    surviving_final_self = surviving_final.T @ spin_final @ surviving_final
    inverse_transfer = np.linalg.inv(surviving_transfer)
    canonical_map = np.block(
        [
            [-inverse_transfer @ surviving_initial_self, -inverse_transfer],
            [
                surviving_transfer.T
                - surviving_final_self
                @ inverse_transfer
                @ surviving_initial_self,
                -surviving_final_self @ inverse_transfer,
            ],
        ]
    )
    identity = np.eye(transfer_rank)
    omega = np.block(
        [[np.zeros_like(identity), identity], [-identity, np.zeros_like(identity)]]
    )
    eigenvalues = np.linalg.eigvals(canonical_map)
    positive_phases = np.sort(np.angle(eigenvalues[np.imag(eigenvalues) > 0.0]))
    simplex_eigenvalues = []
    for simplex in simplices:
        tangents = points[list(simplex[1:])] - points[simplex[0]]
        simplex_eigenvalues.append(
            np.linalg.eigvalsh(tangents @ MINKOWSKI_METRIC @ tangents.T)
        )
    return {
        "segment_count": segment_count,
        "q_cycle_vertex_count": len(geometry["q_vertices"]),
        "simplex_count": len(simplices),
        "boundary_classification": _boundary_classification(geometry),
        "cap_relative_edge_count": len(initial_edges),
        "cap_tangent_rank": initial_tangent_rank,
        "cap_physical_mode_count": initial_normal.shape[1],
        "spin_sample_count": spin_count,
        "initial_spin_projection_residual": initial_spin_residual,
        "final_spin_projection_residual": final_spin_residual,
        "initial_tangent_singular_values": initial_tangent_singular_values,
        "final_tangent_singular_values": final_tangent_singular_values,
        "central_gauge_singular_values": gauge_singular_values,
        "physical_internal_hessian_singular_values": np.linalg.svd(
            internal_hessian, compute_uv=False
        ),
        "spin_transfer_singular_values": transfer_singular_values,
        "spin_transfer_rank": transfer_rank,
        "coarse_injection_initial_overlap_singular_values": initial_overlap,
        "coarse_injection_final_overlap_singular_values": final_overlap,
        "initial_survivor_constraint_leakage": float(
            np.linalg.norm(surviving_initial.T @ spin_initial @ constrained_initial)
        ),
        "final_survivor_constraint_leakage": float(
            np.linalg.norm(surviving_final.T @ spin_final @ constrained_final)
        ),
        "surviving_initial_self": surviving_initial_self,
        "surviving_transfer": surviving_transfer,
        "surviving_final_self": surviving_final_self,
        "surviving_canonical_eigenvalues": eigenvalues,
        "surviving_positive_oscillation_phases": positive_phases,
        "surviving_unit_circle_residual": float(
            np.max(np.abs(np.abs(eigenvalues) - 1.0))
        ),
        "surviving_symplectic_residual": float(
            np.max(np.abs(canonical_map.T @ omega @ canonical_map - omega))
        ),
        "minimum_negative_simplex_eigenvalue_magnitude": min(
            abs(values[0]) for values in simplex_eigenvalues
        ),
        "minimum_positive_simplex_eigenvalue": min(
            values[1] for values in simplex_eigenvalues
        ),
    }


if __name__ == "__main__":
    for level in (2, 4):
        result = multilevel_spin_two_certificate(level)
        for key in (
            "segment_count",
            "q_cycle_vertex_count",
            "simplex_count",
            "boundary_classification",
            "cap_relative_edge_count",
            "cap_tangent_rank",
            "cap_physical_mode_count",
            "spin_sample_count",
            "initial_spin_projection_residual",
            "physical_internal_hessian_singular_values",
            "spin_transfer_singular_values",
            "spin_transfer_rank",
            "coarse_injection_initial_overlap_singular_values",
            "coarse_injection_final_overlap_singular_values",
            "initial_survivor_constraint_leakage",
            "final_survivor_constraint_leakage",
            "surviving_canonical_eigenvalues",
            "surviving_positive_oscillation_phases",
            "surviving_unit_circle_residual",
            "surviving_symplectic_residual",
            "minimum_negative_simplex_eigenvalue_magnitude",
            "minimum_positive_simplex_eigenvalue",
        ):
            print(f"m={level} {key}={result[key]}")
