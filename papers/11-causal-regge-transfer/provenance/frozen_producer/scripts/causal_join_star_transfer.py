"""Causal-boundary and cap-to-cap transfer audit for the C5*C6 join-star.

The coordinate embedding used in A27 was previously treated only as a closed
boundary.  Its boundary tetrahedra actually split by outward-normal signature
into two past-spacelike q-strips, two future-spacelike q-strips, and one
timelike q-strip.  Fixing the endpoint cones and the timelike strip leaves
eight relative edge lengths on each spacelike cap.  This module removes their
embedded one-vertex tangent directions, eliminates the seven physical central
radial modes, and measures the cap-to-cap Hamilton transfer rank.

Normal classification is geometric.  The choice of which shared endpoint cone
is fixed is part of the declared local-cobordism polarization; no Weyl or wave
condition is used to choose it.
"""

import numpy as np

try:
    from scripts.complex_lorentzian_regge import (
        MINKOWSKI_METRIC,
        richardson_complex_hessian,
        squared_intervals,
        varied_lengths,
    )
    from scripts.lorentzian_weighted_weyl_star import (
        INTERNAL_EDGES,
        POINTS,
        Q_CYCLE_EDGES,
        R_CYCLE_EDGES,
        SIMPLICES,
    )
    from scripts.lorentzian_vertex_hybrid_audit import (
        hinge_vertex_area_share,
        signature_corrected_hinge_riemann,
    )
    from scripts.lorentzian_weyl_star import (
        _general_outward_normal,
        _symmetric_components,
        _topology_for,
        branch_curvature,
        general_complex_action,
    )
    from scripts.algebraic_regge_wave_target import electric_magnetic_weyl
except ModuleNotFoundError:
    from complex_lorentzian_regge import (
        MINKOWSKI_METRIC,
        richardson_complex_hessian,
        squared_intervals,
        varied_lengths,
    )
    from lorentzian_weighted_weyl_star import (
        INTERNAL_EDGES,
        POINTS,
        Q_CYCLE_EDGES,
        R_CYCLE_EDGES,
        SIMPLICES,
    )
    from lorentzian_vertex_hybrid_audit import (
        hinge_vertex_area_share,
        signature_corrected_hinge_riemann,
    )
    from lorentzian_weyl_star import (
        _general_outward_normal,
        _symmetric_components,
        _topology_for,
        branch_curvature,
        general_complex_action,
    )
    from algebraic_regge_wave_target import electric_magnetic_weyl


PAST_Q_EDGES = ((1, 2), (2, 3))
FUTURE_Q_EDGES = ((4, 5), (1, 5))
SIDE_Q_EDGE = (3, 4)
INITIAL_INTERIOR_Q_VERTEX = 2
FINAL_INTERIOR_Q_VERTEX = 5
INITIAL_RELATIVE_EDGES = tuple(
    sorted(
        PAST_Q_EDGES
        + tuple((INITIAL_INTERIOR_Q_VERTEX, r) for r in range(6, 12))
    )
)
FINAL_RELATIVE_EDGES = tuple(
    sorted(
        FUTURE_Q_EDGES
        + tuple(
            tuple(sorted((FINAL_INTERIOR_Q_VERTEX, r)))
            for r in range(6, 12)
        )
    )
)


def regular_transverse_points(phase_degrees: float = 0.0) -> np.ndarray:
    """Restore an exact C6 action in the transverse x-y plane."""
    points = np.array(POINTS, copy=True)
    angles = np.deg2rad(phase_degrees + 60.0 * np.arange(6))
    points[6:, 0] = 0.0
    points[6:, 1] = np.cos(angles)
    points[6:, 2] = np.sin(angles)
    points[6:, 3] = 0.0
    return points


def time_reflection_symmetric_points() -> np.ndarray:
    """Well-conditioned C6- and time-reflection-symmetric causal witness.

    Time reflection fixes q1 and exchanges q2<->q5 and q3<->q4.  The values
    were selected without reference to the action or transfer spectrum: they
    maximize elementary causal-normal and q-plane nondegeneracy margins within
    a bounded random geometry scan.
    """
    points = regular_transverse_points()
    points[1:6] = np.asarray(
        (
            (0.0, 0.0, 0.0, -1.5709397991101974),
            (0.3762757793700152, 0.0, 0.0, -0.42646145218382514),
            (0.7702937141192545, 0.0, 0.0, 1.1723475896580613),
            (-0.7702937141192545, 0.0, 0.0, 1.1723475896580613),
            (-0.3762757793700152, 0.0, 0.0, -0.42646145218382514),
        )
    )
    return points


def _relative_normal_basis(
    points: np.ndarray,
    relative_edges: tuple[tuple[int, int], ...],
    moving_vertex: int,
) -> tuple[np.ndarray, np.ndarray]:
    tangent = np.zeros((len(relative_edges), 4), dtype=float)
    for row, edge in enumerate(relative_edges):
        if moving_vertex not in edge:
            continue
        other = edge[0] if edge[1] == moving_vertex else edge[1]
        tangent[row] = 2.0 * (
            points[moving_vertex] - points[other]
        ) @ MINKOWSKI_METRIC
    left, singular_values, _ = np.linalg.svd(tangent, full_matrices=True)
    rank = int(np.count_nonzero(singular_values > singular_values[0] * 1e-10))
    return left[:, rank:], singular_values


def boundary_normal_classification(points: np.ndarray = POINTS) -> dict:
    coordinates = {index: point for index, point in enumerate(points)}
    records = []
    for q_edge in Q_CYCLE_EDGES:
        for r_edge in R_CYCLE_EDGES:
            face = tuple(sorted(q_edge + r_edge))
            normal = _general_outward_normal(face, 0, coordinates)
            norm_squared = float(normal @ MINKOWSKI_METRIC @ normal)
            if norm_squared < 0.0:
                kind = "future_spacelike" if normal[0] > 0.0 else "past_spacelike"
            else:
                kind = "timelike_side"
            records.append(
                {
                    "face": face,
                    "q_edge": q_edge,
                    "r_edge": r_edge,
                    "normal": normal,
                    "normal_norm_squared": norm_squared,
                    "kind": kind,
                }
            )
    counts = {
        kind: sum(record["kind"] == kind for record in records)
        for kind in ("past_spacelike", "future_spacelike", "timelike_side")
    }
    return {"records": records, "counts": counts}


def _spin_two_coordinates(
    points: np.ndarray,
    relative_edges: tuple[tuple[int, int], ...],
    normal_basis: np.ndarray,
    q_vertex: int,
) -> tuple[np.ndarray, float]:
    """Orthonormal normal-coordinate representatives of cos(2phi), sin(2phi)."""
    raw = np.zeros((len(relative_edges), 2), dtype=float)
    for row, edge in enumerate(relative_edges):
        if q_vertex not in edge:
            continue
        other = edge[0] if edge[1] == q_vertex else edge[1]
        if other < 6:
            continue
        angle = np.arctan2(points[other, 2], points[other, 1])
        raw[row] = (np.cos(2.0 * angle), np.sin(2.0 * angle))
    coordinates = normal_basis.T @ raw
    orthonormal, upper = np.linalg.qr(coordinates)
    if np.min(np.abs(np.diag(upper))) < 1e-10:
        raise ArithmeticError("spin-two transverse samples collapsed in the normal quotient")
    projected = normal_basis @ normal_basis.T @ raw
    projection_residual = float(
        np.linalg.norm(raw - projected) / np.linalg.norm(raw)
    )
    return orthonormal[:, :2], projection_residual


def causal_join_star_transfer_certificate(
    points: np.ndarray | None = None,
    hessian_step: float = 1.2e-4,
    curvature_step: float = 1e-5,
) -> dict:
    if points is None:
        points = POINTS
    points = np.asarray(points, dtype=float)
    topology = _topology_for(SIMPLICES)
    background = squared_intervals(points)
    initial_basis, initial_tangent_singular_values = _relative_normal_basis(
        points, INITIAL_RELATIVE_EDGES, INITIAL_INTERIOR_Q_VERTEX
    )
    final_basis, final_tangent_singular_values = _relative_normal_basis(
        points, FINAL_RELATIVE_EDGES, FINAL_INTERIOR_Q_VERTEX
    )
    internal_background = np.asarray([background[edge] for edge in INTERNAL_EDGES])
    initial_background = np.asarray(
        [background[edge] for edge in INITIAL_RELATIVE_EDGES]
    )
    final_background = np.asarray(
        [background[edge] for edge in FINAL_RELATIVE_EDGES]
    )
    internal_count = len(INTERNAL_EDGES)
    initial_count = initial_basis.shape[1]
    final_count = final_basis.shape[1]
    coordinate_origin = np.concatenate(
        (internal_background, np.zeros(initial_count + final_count))
    )

    def lengths_from_coordinates(coordinates: np.ndarray) -> np.ndarray:
        lengths = varied_lengths(
            background, INTERNAL_EDGES, coordinates[:internal_count]
        )
        initial_values = (
            initial_background
            + initial_basis
            @ coordinates[internal_count : internal_count + initial_count]
        )
        lengths = varied_lengths(lengths, INITIAL_RELATIVE_EDGES, initial_values)
        final_values = (
            final_background
            + final_basis
            @ coordinates[internal_count + initial_count :]
        )
        return varied_lengths(lengths, FINAL_RELATIVE_EDGES, final_values)

    def action_from_coordinates(coordinates: np.ndarray) -> complex:
        return general_complex_action(
            lengths_from_coordinates(coordinates),
            SIMPLICES,
            topology["bulk_internal_triangles"],
        )

    hessian = richardson_complex_hessian(
        action_from_coordinates, coordinate_origin, hessian_step
    )
    central_displacement_tangent = np.asarray(
        [
            2.0 * (points[first] - points[second]) @ MINKOWSKI_METRIC
            for first, second in INTERNAL_EDGES
        ]
    )
    gauge_left, gauge_singular_values, _ = np.linalg.svd(
        central_displacement_tangent, full_matrices=True
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
    boundary_hessian = hessian[internal_count:, internal_count:]
    internal_response = -np.linalg.solve(internal_hessian, internal_boundary)
    effective_boundary_hessian = (
        boundary_hessian + internal_boundary.T @ internal_response
    )
    transfer = effective_boundary_hessian[:initial_count, initial_count:]
    transfer_singular_values = np.linalg.svd(transfer, compute_uv=False)
    transfer_rank = int(
        np.count_nonzero(
            transfer_singular_values > transfer_singular_values[0] * 1e-7
        )
    )
    initial_spin_two, initial_spin_projection_residual = _spin_two_coordinates(
        points,
        INITIAL_RELATIVE_EDGES,
        initial_basis,
        INITIAL_INTERIOR_Q_VERTEX,
    )
    final_spin_two, final_spin_projection_residual = _spin_two_coordinates(
        points,
        FINAL_RELATIVE_EDGES,
        final_basis,
        FINAL_INTERIOR_Q_VERTEX,
    )
    initial_complement = np.linalg.svd(
        initial_spin_two.T, full_matrices=True
    )[2].T[:, 2:]
    final_complement = np.linalg.svd(
        final_spin_two.T, full_matrices=True
    )[2].T[:, 2:]
    spin_two_transfer = initial_spin_two.T @ transfer @ final_spin_two
    spin_to_complement = initial_spin_two.T @ transfer @ final_complement
    complement_to_spin = initial_complement.T @ transfer @ final_spin_two
    complement_transfer = initial_complement.T @ transfer @ final_complement
    initial_self = effective_boundary_hessian[:initial_count, :initial_count]
    final_self = effective_boundary_hessian[initial_count:, initial_count:]
    spin_two_initial_self = initial_spin_two.T @ initial_self @ initial_spin_two
    spin_two_final_self = final_spin_two.T @ final_self @ final_spin_two
    initial_spin_self_leakage = initial_spin_two.T @ initial_self @ initial_complement
    final_spin_self_leakage = final_spin_two.T @ final_self @ final_complement
    inverse_transfer = np.linalg.inv(spin_two_transfer)
    canonical_map = np.block(
        [
            [
                -inverse_transfer @ spin_two_initial_self,
                -inverse_transfer,
            ],
            [
                spin_two_transfer.T
                - spin_two_final_self
                @ inverse_transfer
                @ spin_two_initial_self,
                -spin_two_final_self @ inverse_transfer,
            ],
        ]
    )
    identity = np.eye(initial_spin_two.shape[1])
    symplectic_form = np.block(
        [[np.zeros_like(identity), identity], [-identity, np.zeros_like(identity)]]
    )
    canonical_eigenvalues = np.linalg.eigvals(canonical_map)
    c6_spin_rotation = np.asarray(
        (
            (np.cos(2.0 * np.pi / 3.0), -np.sin(2.0 * np.pi / 3.0)),
            (np.sin(2.0 * np.pi / 3.0), np.cos(2.0 * np.pi / 3.0)),
        )
    )
    phase_space_c6_rotation = np.block(
        [
            [c6_spin_rotation, np.zeros((2, 2))],
            [np.zeros((2, 2)), c6_spin_rotation],
        ]
    )
    positive_phase_eigenvalues = canonical_eigenvalues[
        np.imag(canonical_eigenvalues) > 0.0
    ]
    oscillation_phase = (
        float(np.mean(np.angle(positive_phase_eigenvalues)))
        if positive_phase_eigenvalues.size
        else float("nan")
    )

    internal_triangles = tuple(sorted(topology["bulk_internal_triangles"]))

    def curvatures(coordinates: np.ndarray) -> np.ndarray:
        lengths = lengths_from_coordinates(coordinates)
        return np.asarray(
            [
                branch_curvature(
                    triangle,
                    lengths,
                    topology["bulk_internal_triangles"],
                    SIMPLICES,
                )
                for triangle in internal_triangles
            ]
        )

    curvature_jacobian = np.zeros(
        (len(internal_triangles), len(coordinate_origin)), dtype=complex
    )
    for column in range(len(coordinate_origin)):
        step = curvature_step * max(1.0, abs(coordinate_origin[column]))
        displacement = np.zeros_like(coordinate_origin)
        displacement[column] = step
        curvature_jacobian[:, column] = (
            curvatures(coordinate_origin + displacement)
            - curvatures(coordinate_origin - displacement)
        ) / (2.0 * step)
    on_shell_curvature = (
        curvature_jacobian[:, internal_count:]
        + curvature_jacobian[:, :internal_count]
        @ physical_internal_basis
        @ internal_response
    )
    area_shares = np.asarray(
        [
            hinge_vertex_area_share(triangle, points)
            for triangle in internal_triangles
        ]
    )
    riemann_tensors = [
        signature_corrected_hinge_riemann(triangle, points)
        for triangle in internal_triangles
    ]
    riemann_map = np.column_stack(
        [riemann.reshape(-1) for riemann in riemann_tensors]
    )
    einstein_tensors = []
    for riemann in riemann_tensors:
        ricci = np.einsum("bd,abcd->ac", MINKOWSKI_METRIC, riemann)
        scalar = np.einsum("ac,ac", MINKOWSKI_METRIC, ricci)
        einstein_tensors.append(ricci - 0.5 * MINKOWSKI_METRIC * scalar)
    einstein_map = np.column_stack(
        [_symmetric_components(einstein) for einstein in einstein_tensors]
    )
    positive_frequency_multiplier = (
        complex(np.mean(positive_phase_eigenvalues))
        if positive_phase_eigenvalues.size
        else complex("nan")
    )
    wave_boundary_axes = np.vstack(
        (
            initial_spin_two,
            positive_frequency_multiplier * final_spin_two,
        )
    )
    wave_curvatures = on_shell_curvature @ wave_boundary_axes
    wave_riemann = [
        (
            riemann_map @ np.diag(area_shares) @ wave_curvatures[:, column]
        ).reshape(4, 4, 4, 4)
        for column in range(2)
    ]
    wave_einstein = (
        einstein_map @ np.diag(area_shares) @ wave_curvatures
    )
    wave_electric_magnetic = [
        electric_magnetic_weyl(riemann) for riemann in wave_riemann
    ]
    plus_electric, plus_magnetic = wave_electric_magnetic[0]
    cross_electric, cross_magnetic = wave_electric_magnetic[1]
    wave_scale = max(
        np.linalg.norm(plus_electric),
        np.linalg.norm(plus_magnetic),
        np.linalg.norm(cross_electric),
        np.linalg.norm(cross_magnetic),
    )
    right_moving_residual = float(
        max(
            np.linalg.norm(plus_magnetic - cross_electric),
            np.linalg.norm(cross_magnetic + plus_electric),
        )
        / wave_scale
    )
    left_moving_residual = float(
        max(
            np.linalg.norm(plus_magnetic + cross_electric),
            np.linalg.norm(cross_magnetic - plus_electric),
        )
        / wave_scale
    )
    impedance_target = np.concatenate(
        (cross_electric.reshape(-1), -plus_electric.reshape(-1))
    )
    impedance_magnetic = np.concatenate(
        (plus_magnetic.reshape(-1), cross_magnetic.reshape(-1))
    )
    fitted_impedance = complex(
        np.vdot(impedance_target, impedance_magnetic)
        / np.vdot(impedance_target, impedance_target)
    )
    fitted_impedance_residual = float(
        np.linalg.norm(
            impedance_magnetic - fitted_impedance * impedance_target
        )
        / np.linalg.norm(impedance_magnetic)
    )

    def scalar_matrix_residual(matrix: np.ndarray) -> float:
        scalar_part = np.trace(matrix) * np.eye(2) / 2.0
        return float(np.linalg.norm(matrix - scalar_part) / np.linalg.norm(matrix))
    classification = boundary_normal_classification(points)
    return {
        "classification": classification,
        "initial_relative_edges": INITIAL_RELATIVE_EDGES,
        "final_relative_edges": FINAL_RELATIVE_EDGES,
        "initial_tangent_singular_values": initial_tangent_singular_values,
        "final_tangent_singular_values": final_tangent_singular_values,
        "initial_mode_count": initial_count,
        "final_mode_count": final_count,
        "central_gauge_singular_values": gauge_singular_values,
        "physical_internal_hessian_singular_values": np.linalg.svd(
            internal_hessian, compute_uv=False
        ),
        "effective_boundary_hessian": effective_boundary_hessian,
        "transfer": transfer,
        "transfer_singular_values": transfer_singular_values,
        "transfer_rank": transfer_rank,
        "initial_spin_two_coordinates": initial_spin_two,
        "final_spin_two_coordinates": final_spin_two,
        "initial_spin_projection_residual": initial_spin_projection_residual,
        "final_spin_projection_residual": final_spin_projection_residual,
        "spin_two_transfer": spin_two_transfer,
        "spin_two_transfer_singular_values": np.linalg.svd(
            spin_two_transfer, compute_uv=False
        ),
        "spin_to_complement_norm": float(np.linalg.norm(spin_to_complement)),
        "complement_to_spin_norm": float(np.linalg.norm(complement_to_spin)),
        "complement_transfer_singular_values": np.linalg.svd(
            complement_transfer, compute_uv=False
        ),
        "spin_two_initial_self": spin_two_initial_self,
        "spin_two_final_self": spin_two_final_self,
        "initial_spin_self_leakage_norm": float(
            np.linalg.norm(initial_spin_self_leakage)
        ),
        "final_spin_self_leakage_norm": float(
            np.linalg.norm(final_spin_self_leakage)
        ),
        "spin_two_canonical_map": canonical_map,
        "spin_two_canonical_determinant": complex(np.linalg.det(canonical_map)),
        "spin_two_symplectic_residual": float(
            np.max(
                np.abs(
                    canonical_map.T @ symplectic_form @ canonical_map
                    - symplectic_form
                )
            )
        ),
        "spin_two_canonical_eigenvalues": canonical_eigenvalues,
        "spin_two_eigenvalue_unit_circle_residual": float(
            np.max(np.abs(np.abs(canonical_eigenvalues) - 1.0))
        ),
        "spin_two_oscillation_phase": oscillation_phase,
        "spin_two_transfer_isotropy_residual": scalar_matrix_residual(
            spin_two_transfer
        ),
        "spin_two_initial_isotropy_residual": scalar_matrix_residual(
            spin_two_initial_self
        ),
        "spin_two_final_isotropy_residual": scalar_matrix_residual(
            spin_two_final_self
        ),
        "c6_canonical_commutator_residual": float(
            np.max(
                np.abs(
                    canonical_map @ phase_space_c6_rotation
                    - phase_space_c6_rotation @ canonical_map
                )
            )
        ),
        "on_shell_curvature_singular_values": np.linalg.svd(
            on_shell_curvature, compute_uv=False
        ),
        "positive_frequency_multiplier": positive_frequency_multiplier,
        "wave_curvatures": wave_curvatures,
        "wave_riemann": wave_riemann,
        "wave_einstein_residual": float(np.max(np.abs(wave_einstein))),
        "wave_einstein_relative_residual": float(
            np.linalg.norm(wave_einstein)
            / np.linalg.norm(
                riemann_map @ np.diag(area_shares) @ wave_curvatures
            )
        ),
        "wave_electric_magnetic": wave_electric_magnetic,
        "right_moving_wave_residual": right_moving_residual,
        "left_moving_wave_residual": left_moving_residual,
        "fitted_wave_impedance": fitted_impedance,
        "fitted_wave_impedance_modulus": float(abs(fitted_impedance)),
        "fitted_wave_impedance_phase": float(np.angle(fitted_impedance)),
        "fitted_wave_impedance_residual": fitted_impedance_residual,
        "minimum_spacelike_normal_magnitude": min(
            abs(record["normal_norm_squared"])
            for record in classification["records"]
            if "spacelike" in record["kind"]
        ),
        "minimum_side_normal_magnitude": min(
            record["normal_norm_squared"]
            for record in classification["records"]
            if record["kind"] == "timelike_side"
        ),
    }


if __name__ == "__main__":
    result = causal_join_star_transfer_certificate()
    print(f"boundary_counts={result['classification']['counts']}")
    print(f"initial_relative_edges={result['initial_relative_edges']}")
    print(f"final_relative_edges={result['final_relative_edges']}")
    print(f"initial_mode_count={result['initial_mode_count']}")
    print(f"final_mode_count={result['final_mode_count']}")
    print(
        "physical_internal_hessian_singular_values="
        f"{result['physical_internal_hessian_singular_values']}"
    )
    print(f"transfer=\n{result['transfer']}")
    print(f"transfer_singular_values={result['transfer_singular_values']}")
    print(f"transfer_rank={result['transfer_rank']}")
    print(f"spin_two_transfer=\n{result['spin_two_transfer']}")
    print(
        "spin_two_transfer_singular_values="
        f"{result['spin_two_transfer_singular_values']}"
    )
    print(f"spin_to_complement_norm={result['spin_to_complement_norm']}")
    print(f"complement_to_spin_norm={result['complement_to_spin_norm']}")
    print(f"spin_two_initial_self=\n{result['spin_two_initial_self']}")
    print(f"spin_two_final_self=\n{result['spin_two_final_self']}")
    print(f"spin_two_canonical_map=\n{result['spin_two_canonical_map']}")
    print(
        "spin_two_canonical_eigenvalues="
        f"{result['spin_two_canonical_eigenvalues']}"
    )
    print(
        "spin_two_symplectic_residual="
        f"{result['spin_two_symplectic_residual']}"
    )
    print(f"wave_einstein_residual={result['wave_einstein_residual']}")
    print(
        "wave_einstein_relative_residual="
        f"{result['wave_einstein_relative_residual']}"
    )
    print(f"right_moving_wave_residual={result['right_moving_wave_residual']}")
    print(f"left_moving_wave_residual={result['left_moving_wave_residual']}")
    print(f"fitted_wave_impedance={result['fitted_wave_impedance']}")
    print(
        "fitted_wave_impedance_residual="
        f"{result['fitted_wave_impedance_residual']}"
    )
