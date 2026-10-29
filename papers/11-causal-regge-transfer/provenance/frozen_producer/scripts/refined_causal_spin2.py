"""First longitudinal refinement of the causal C6 spin-two carrier.

The C5*C6 cone has one interior q vertex on each spacelike cap.  This module
subdivides a time-reflected pair of cap edges, producing C7*C6 with two interior
q vertices per cap.  Exact transverse C6 symmetry isolates four configuration
directions: plus/cross samples at two longitudinal sites.  The resulting
canonical map tests whether the coarse stable oscillator persists as two
spin-two lattice bands rather than a one-cell accident.
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


def refined_points() -> np.ndarray:
    coarse = time_reflection_symmetric_points()
    coarse_q = coarse[1:6]
    q_points = np.asarray(
        (
            coarse_q[0],
            0.5 * (coarse_q[0] + coarse_q[1]),
            coarse_q[1],
            coarse_q[2],
            coarse_q[3],
            coarse_q[4],
            0.5 * (coarse_q[4] + coarse_q[0]),
        )
    )
    return np.vstack((np.zeros((1, 4)), q_points, coarse[6:]))


POINTS = refined_points()
Q_VERTICES = tuple(range(1, 8))
R_VERTICES = tuple(range(8, 14))
Q_CYCLE_EDGES = tuple(
    tuple(sorted((Q_VERTICES[index], Q_VERTICES[(index + 1) % 7])))
    for index in range(7)
)
R_CYCLE_EDGES = tuple(
    tuple(sorted((R_VERTICES[index], R_VERTICES[(index + 1) % 6])))
    for index in range(6)
)
BOUNDARY_EDGES = tuple(
    sorted(
        set(
            Q_CYCLE_EDGES
            + R_CYCLE_EDGES
            + tuple(
                tuple(sorted(edge)) for edge in product(Q_VERTICES, R_VERTICES)
            )
        )
    )
)
INTERNAL_EDGES = tuple((0, vertex) for vertex in range(1, 14))
SIMPLICES = tuple(
    tuple(
        sorted(
            (
                0,
                Q_VERTICES[q_index],
                Q_VERTICES[(q_index + 1) % 7],
                R_VERTICES[r_index],
                R_VERTICES[(r_index + 1) % 6],
            )
        )
    )
    for q_index in range(7)
    for r_index in range(6)
)
INITIAL_Q_PATH = (1, 2, 3, 4)
FINAL_Q_PATH = (5, 6, 7, 1)
INITIAL_INTERIOR_Q = (2, 3)
FINAL_INTERIOR_Q = (7, 6)  # time-reflection order matching (2,3)


def _path_edges(path: tuple[int, ...]) -> tuple[tuple[int, int], ...]:
    return tuple(
        tuple(sorted((path[index], path[index + 1])))
        for index in range(len(path) - 1)
    )


def _relative_edges(
    path: tuple[int, ...], interior_q: tuple[int, ...]
) -> tuple[tuple[int, int], ...]:
    return tuple(
        sorted(
            set(
                _path_edges(path)
                + tuple(
                    tuple(sorted((q_vertex, r_vertex)))
                    for q_vertex in interior_q
                    for r_vertex in R_VERTICES
                )
            )
        )
    )


INITIAL_RELATIVE_EDGES = _relative_edges(INITIAL_Q_PATH, INITIAL_INTERIOR_Q)
FINAL_RELATIVE_EDGES = _relative_edges(FINAL_Q_PATH, FINAL_INTERIOR_Q)


def _multi_vertex_normal_basis(
    relative_edges: tuple[tuple[int, int], ...], moving_vertices: tuple[int, ...]
) -> tuple[np.ndarray, np.ndarray, int]:
    tangent = np.zeros((len(relative_edges), 4 * len(moving_vertices)))
    column_vertex = {vertex: index for index, vertex in enumerate(moving_vertices)}
    for row, edge in enumerate(relative_edges):
        difference = POINTS[edge[0]] - POINTS[edge[1]]
        covector = 2.0 * difference @ MINKOWSKI_METRIC
        for endpoint, sign in ((edge[0], 1.0), (edge[1], -1.0)):
            if endpoint in column_vertex:
                start = 4 * column_vertex[endpoint]
                tangent[row, start : start + 4] = sign * covector
    left, singular_values, _ = np.linalg.svd(tangent, full_matrices=True)
    rank = int(np.count_nonzero(singular_values > singular_values[0] * 1e-10))
    return left[:, rank:], singular_values, rank


def _spin_two_coordinates(
    relative_edges: tuple[tuple[int, int], ...],
    interior_q: tuple[int, ...],
    normal_basis: np.ndarray,
) -> tuple[np.ndarray, float, np.ndarray]:
    raw = np.zeros((len(relative_edges), 2 * len(interior_q)))
    q_position = {vertex: index for index, vertex in enumerate(interior_q)}
    for row, edge in enumerate(relative_edges):
        q_candidates = [vertex for vertex in edge if vertex in q_position]
        r_candidates = [vertex for vertex in edge if vertex in R_VERTICES]
        if len(q_candidates) != 1 or len(r_candidates) != 1:
            continue
        q_index = q_position[q_candidates[0]]
        r_vertex = r_candidates[0]
        angle = np.arctan2(POINTS[r_vertex, 2], POINTS[r_vertex, 1])
        raw[row, 2 * q_index : 2 * q_index + 2] = (
            np.cos(2.0 * angle),
            np.sin(2.0 * angle),
        )
    coordinates = normal_basis.T @ raw
    orthonormal, upper = np.linalg.qr(coordinates)
    if np.min(np.abs(np.diag(upper))) < 1e-10:
        raise ArithmeticError("refined spin-two samples collapsed")
    projected = normal_basis @ normal_basis.T @ raw
    residual = float(np.linalg.norm(raw - projected) / np.linalg.norm(raw))
    return orthonormal[:, : raw.shape[1]], residual, upper[: raw.shape[1], :]


def _boundary_classification() -> dict:
    coordinates = {index: point for index, point in enumerate(POINTS)}
    counts = {"past_spacelike": 0, "future_spacelike": 0, "timelike_side": 0}
    margins = []
    for q_edge in Q_CYCLE_EDGES:
        for r_edge in R_CYCLE_EDGES:
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


def refined_spin_two_certificate(hessian_step: float = 1.5e-4) -> dict:
    topology = _topology_for(SIMPLICES)
    background = squared_intervals(POINTS)
    initial_basis, initial_tangent_singular_values, initial_tangent_rank = (
        _multi_vertex_normal_basis(INITIAL_RELATIVE_EDGES, INITIAL_INTERIOR_Q)
    )
    final_basis, final_tangent_singular_values, final_tangent_rank = (
        _multi_vertex_normal_basis(FINAL_RELATIVE_EDGES, FINAL_INTERIOR_Q)
    )
    initial_spin, initial_spin_projection_residual, initial_spin_raw_transform = _spin_two_coordinates(
        INITIAL_RELATIVE_EDGES, INITIAL_INTERIOR_Q, initial_basis
    )
    final_spin, final_spin_projection_residual, final_spin_raw_transform = _spin_two_coordinates(
        FINAL_RELATIVE_EDGES, FINAL_INTERIOR_Q, final_basis
    )
    internal_count = len(INTERNAL_EDGES)
    initial_count = initial_basis.shape[1]
    final_count = final_basis.shape[1]
    internal_background = np.asarray([background[edge] for edge in INTERNAL_EDGES])
    initial_background = np.asarray(
        [background[edge] for edge in INITIAL_RELATIVE_EDGES]
    )
    final_background = np.asarray([background[edge] for edge in FINAL_RELATIVE_EDGES])
    coordinate_origin = np.concatenate(
        (internal_background, np.zeros(initial_count + final_count))
    )

    def lengths_from_coordinates(coordinates: np.ndarray) -> np.ndarray:
        lengths = varied_lengths(
            background, INTERNAL_EDGES, coordinates[:internal_count]
        )
        initial_values = initial_background + initial_basis @ coordinates[
            internal_count : internal_count + initial_count
        ]
        lengths = varied_lengths(lengths, INITIAL_RELATIVE_EDGES, initial_values)
        final_values = final_background + final_basis @ coordinates[
            internal_count + initial_count :
        ]
        return varied_lengths(lengths, FINAL_RELATIVE_EDGES, final_values)

    def action(coordinates: np.ndarray) -> complex:
        return general_complex_action(
            lengths_from_coordinates(coordinates),
            SIMPLICES,
            topology["bulk_internal_triangles"],
        )

    hessian = richardson_complex_hessian(action, coordinate_origin, hessian_step)
    central_tangent = np.asarray(
        [
            2.0 * (POINTS[first] - POINTS[second]) @ MINKOWSKI_METRIC
            for first, second in INTERNAL_EDGES
        ]
    )
    gauge_left, gauge_singular_values, _ = np.linalg.svd(
        central_tangent, full_matrices=True
    )
    physical_internal_basis = gauge_left[:, 4:]
    internal_hessian = (
        physical_internal_basis.T
        @ hessian[:internal_count, :internal_count]
        @ physical_internal_basis
    )
    internal_boundary = (
        physical_internal_basis.T @ hessian[:internal_count, internal_count:]
    )
    effective = (
        hessian[internal_count:, internal_count:]
        - internal_boundary.T
        @ np.linalg.solve(internal_hessian, internal_boundary)
    )
    initial_self = effective[:initial_count, :initial_count]
    transfer = effective[:initial_count, initial_count:]
    final_self = effective[initial_count:, initial_count:]
    spin_initial = initial_spin.T @ initial_self @ initial_spin
    spin_transfer = initial_spin.T @ transfer @ final_spin
    spin_final = final_spin.T @ final_self @ final_spin
    transfer_left, transfer_singular_values, transfer_right_h = np.linalg.svd(
        spin_transfer, full_matrices=True
    )
    spin_transfer_rank = int(
        np.count_nonzero(transfer_singular_values > transfer_singular_values[0] * 1e-5)
    )
    surviving_initial = transfer_left[:, :spin_transfer_rank]
    surviving_final = transfer_right_h.T[:, :spin_transfer_rank]
    constrained_initial = transfer_left[:, spin_transfer_rank:]
    constrained_final = transfer_right_h.T[:, spin_transfer_rank:]
    coarse_site_injection = np.asarray(
        ((0.5, 0.0), (0.0, 0.5), (1.0, 0.0), (0.0, 1.0))
    )
    injected_initial = np.linalg.qr(
        initial_spin_raw_transform @ coarse_site_injection
    )[0][:, :2]
    injected_final = np.linalg.qr(
        final_spin_raw_transform @ coarse_site_injection
    )[0][:, :2]
    initial_injection_overlap = np.linalg.svd(
        surviving_initial.T @ injected_initial, compute_uv=False
    )
    final_injection_overlap = np.linalg.svd(
        surviving_final.T @ injected_final, compute_uv=False
    )
    surviving_initial_self = (
        surviving_initial.T @ spin_initial @ surviving_initial
    )
    surviving_transfer = (
        surviving_initial.T @ spin_transfer @ surviving_final
    )
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
    identity = np.eye(spin_transfer_rank)
    symplectic_form = np.block(
        [[np.zeros_like(identity), identity], [-identity, np.zeros_like(identity)]]
    )
    eigenvalues = np.linalg.eigvals(canonical_map)
    positive_phases = np.sort(np.angle(eigenvalues[np.imag(eigenvalues) > 0.0]))
    simplex_eigenvalues = []
    for simplex in SIMPLICES:
        tangents = POINTS[list(simplex[1:])] - POINTS[simplex[0]]
        simplex_eigenvalues.append(
            np.linalg.eigvalsh(tangents @ MINKOWSKI_METRIC @ tangents.T)
        )
    return {
        "topology": topology,
        "boundary_classification": _boundary_classification(),
        "initial_relative_edge_count": len(INITIAL_RELATIVE_EDGES),
        "final_relative_edge_count": len(FINAL_RELATIVE_EDGES),
        "initial_tangent_rank": initial_tangent_rank,
        "final_tangent_rank": final_tangent_rank,
        "initial_mode_count": initial_count,
        "final_mode_count": final_count,
        "initial_spin_projection_residual": initial_spin_projection_residual,
        "final_spin_projection_residual": final_spin_projection_residual,
        "initial_tangent_singular_values": initial_tangent_singular_values,
        "final_tangent_singular_values": final_tangent_singular_values,
        "central_gauge_singular_values": gauge_singular_values,
        "physical_internal_hessian_singular_values": np.linalg.svd(
            internal_hessian, compute_uv=False
        ),
        "spin_initial": spin_initial,
        "spin_transfer": spin_transfer,
        "spin_final": spin_final,
        "spin_transfer_singular_values": np.linalg.svd(
            spin_transfer, compute_uv=False
        ),
        "spin_transfer_rank": spin_transfer_rank,
        "surviving_initial_axes": surviving_initial,
        "surviving_final_axes": surviving_final,
        "coarse_injection_initial_overlap_singular_values": initial_injection_overlap,
        "coarse_injection_final_overlap_singular_values": final_injection_overlap,
        "surviving_initial_self": surviving_initial_self,
        "surviving_transfer": surviving_transfer,
        "surviving_final_self": surviving_final_self,
        "initial_survivor_constraint_leakage": float(
            np.linalg.norm(
                surviving_initial.T @ spin_initial @ constrained_initial
            )
        ),
        "final_survivor_constraint_leakage": float(
            np.linalg.norm(surviving_final.T @ spin_final @ constrained_final)
        ),
        "surviving_canonical_map": canonical_map,
        "surviving_canonical_eigenvalues": eigenvalues,
        "surviving_positive_oscillation_phases": positive_phases,
        "surviving_unit_circle_residual": float(
            np.max(np.abs(np.abs(eigenvalues) - 1.0))
        ),
        "surviving_symplectic_residual": float(
            np.max(
                np.abs(
                    canonical_map.T @ symplectic_form @ canonical_map
                    - symplectic_form
                )
            )
        ),
        "minimum_negative_simplex_eigenvalue_magnitude": min(
            abs(values[0]) for values in simplex_eigenvalues
        ),
        "minimum_positive_simplex_eigenvalue": min(
            values[1] for values in simplex_eigenvalues
        ),
    }


if __name__ == "__main__":
    result = refined_spin_two_certificate()
    for key in (
        "boundary_classification",
        "initial_relative_edge_count",
        "initial_tangent_rank",
        "initial_mode_count",
        "initial_spin_projection_residual",
        "physical_internal_hessian_singular_values",
        "spin_transfer",
        "spin_transfer_singular_values",
        "spin_transfer_rank",
        "surviving_transfer",
        "coarse_injection_initial_overlap_singular_values",
        "coarse_injection_final_overlap_singular_values",
        "initial_survivor_constraint_leakage",
        "final_survivor_constraint_leakage",
        "surviving_canonical_eigenvalues",
        "surviving_positive_oscillation_phases",
        "surviving_unit_circle_residual",
        "surviving_symplectic_residual",
    ):
        print(f"{key}={result[key]}")
