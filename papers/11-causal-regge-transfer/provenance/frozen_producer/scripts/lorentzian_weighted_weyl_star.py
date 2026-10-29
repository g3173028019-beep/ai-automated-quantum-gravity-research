"""Dual-weight-robust C5*C6 Lorentzian Weyl-star certificate.

The C5*C5 star of A26 has five boundary non-embeddability modes; plausible
dual-area weights leave only one Einstein-null direction.  Adding one vertex
to the transverse cycle gives seven boundary modes.  After seven physical
radial equations, both the unsigned and incidence-oriented circumcentric
screens retain exactly two null directions.

The boundary remains a closed S^3 and the final signed Lorentzian hybrid-volume
convention remains open.  This is not a causal graviton propagator.
"""

from itertools import product

import numpy as np

try:
    from scripts.complex_lorentzian_regge import (
        MINKOWSKI_METRIC,
        richardson_complex_hessian,
        squared_intervals,
        varied_lengths,
    )
    from scripts.algebraic_regge_wave_target import electric_magnetic_weyl
    from scripts.lorentzian_weyl_star import (
        EX,
        EY,
        FIVE_Q_POINTS,
        _hinge_normal_riemann,
        _oriented_flag_dual_area_magnitude,
        _rank,
        _symmetric_components,
        _topology_for,
        _unsigned_circumcentric_dual_area,
        branch_curvature,
        general_complex_action,
        length_tangent,
    )
except ModuleNotFoundError:
    from complex_lorentzian_regge import (
        MINKOWSKI_METRIC,
        richardson_complex_hessian,
        squared_intervals,
        varied_lengths,
    )
    from algebraic_regge_wave_target import electric_magnetic_weyl
    from lorentzian_weyl_star import (
        EX,
        EY,
        FIVE_Q_POINTS,
        _hinge_normal_riemann,
        _oriented_flag_dual_area_magnitude,
        _rank,
        _symmetric_components,
        _topology_for,
        _unsigned_circumcentric_dual_area,
        branch_curvature,
        general_complex_action,
        length_tangent,
    )


Edge = tuple[int, int]
Triangle = tuple[int, int, int]
Simplex4 = tuple[int, int, int, int, int]

Q_VERTICES = (1, 2, 3, 4, 5)
R_VERTICES = (6, 7, 8, 9, 10, 11)
R_ANGLES = np.deg2rad((-135.0, -45.0, 0.0, 45.0, 90.0, 135.0))
R_POINTS = tuple(
    np.cos(angle) * EX + np.sin(angle) * EY for angle in R_ANGLES
)
POINTS = np.asarray((np.zeros(4),) + FIVE_Q_POINTS + R_POINTS)
Q_CYCLE_EDGES: tuple[Edge, ...] = tuple(
    tuple(sorted((Q_VERTICES[index], Q_VERTICES[(index + 1) % 5])))
    for index in range(5)
)
R_CYCLE_EDGES: tuple[Edge, ...] = tuple(
    tuple(sorted((R_VERTICES[index], R_VERTICES[(index + 1) % 6])))
    for index in range(6)
)
CROSS_EDGES: tuple[Edge, ...] = tuple(
    tuple(sorted(edge)) for edge in product(Q_VERTICES, R_VERTICES)
)
BOUNDARY_EDGES: tuple[Edge, ...] = tuple(
    sorted(set(Q_CYCLE_EDGES + R_CYCLE_EDGES + CROSS_EDGES))
)
INTERNAL_EDGES: tuple[Edge, ...] = tuple((0, vertex) for vertex in range(1, 12))
ALL_EDGES: tuple[Edge, ...] = BOUNDARY_EDGES + INTERNAL_EDGES
SIMPLICES: tuple[Simplex4, ...] = tuple(
    tuple(
        sorted(
            (
                0,
                Q_VERTICES[q_index],
                Q_VERTICES[(q_index + 1) % 5],
                R_VERTICES[r_index],
                R_VERTICES[(r_index + 1) % 6],
            )
        )
    )
    for q_index in range(5)
    for r_index in range(6)
)


def transverse_angle_deformed_points(angle_degrees: float = 50.0) -> np.ndarray:
    """Move only the 45-degree transverse vertex while preserving combinatorics."""
    points = np.array(POINTS, copy=True)
    angle = np.deg2rad(angle_degrees)
    points[9] = np.cos(angle) * EX + np.sin(angle) * EY
    return points


def _radiation_screen(riemann_basis: list[np.ndarray]) -> dict:
    electric_magnetic = [electric_magnetic_weyl(riemann) for riemann in riemann_basis]
    dimension = len(riemann_basis)
    invariant_one = np.zeros((dimension, dimension), dtype=complex)
    invariant_two = np.zeros_like(invariant_one)
    for first, (first_electric, first_magnetic) in enumerate(electric_magnetic):
        for second, (second_electric, second_magnetic) in enumerate(electric_magnetic):
            invariant_one[first, second] = np.sum(
                first_electric * second_electric
                - first_magnetic * second_magnetic
            )
            invariant_two[first, second] = 0.5 * np.sum(
                first_electric * second_magnetic
                + second_electric * first_magnetic
            )
    roots = np.asarray((), dtype=complex)
    null_invariants = []
    right_residuals = []
    left_residuals = []
    if dimension == 2:
        roots = np.roots(
            (
                invariant_one[1, 1],
                2.0 * invariant_one[0, 1],
                invariant_one[0, 0],
            )
        )
        for root in roots:
            riemann = riemann_basis[0] + root * riemann_basis[1]
            electric, magnetic = electric_magnetic_weyl(riemann)
            rotated_electric = np.asarray(
                (
                    (-electric[0, 1], electric[0, 0], 0.0),
                    (electric[0, 0], electric[0, 1], 0.0),
                    (0.0, 0.0, 0.0),
                ),
                dtype=complex,
            )
            scale = max(np.linalg.norm(electric), np.linalg.norm(magnetic))
            right_residuals.append(
                float(np.linalg.norm(magnetic - rotated_electric) / scale)
            )
            left_residuals.append(
                float(np.linalg.norm(magnetic + rotated_electric) / scale)
            )
            null_invariants.append(
                np.asarray(
                    (
                        np.sum(electric * electric - magnetic * magnetic),
                        np.sum(electric * magnetic),
                    )
                )
            )
    return {
        "electric_magnetic_basis": electric_magnetic,
        "invariant_one": invariant_one,
        "invariant_two": invariant_two,
        "null_roots": roots,
        "null_invariants": np.asarray(null_invariants),
        "right_residuals": np.asarray(right_residuals),
        "left_residuals": np.asarray(left_residuals),
    }


def weighted_weyl_star_certificate(
    hessian_step: float = 1.2e-4,
    curvature_step: float = 1e-5,
    points: np.ndarray | None = None,
) -> dict:
    if points is None:
        points = POINTS
    points = np.asarray(points, dtype=float)
    if points.shape != POINTS.shape:
        raise ValueError("deformed witness must retain the C5*C6 vertex set")
    topology = _topology_for(SIMPLICES)
    internal_triangles = tuple(sorted(topology["bulk_internal_triangles"]))
    background = squared_intervals(points)
    relative_boundary_edges = tuple(
        (first - 1, second - 1) for first, second in BOUNDARY_EDGES
    )
    boundary_tangent = length_tangent(points[1:], relative_boundary_edges)
    boundary_rank = _rank(boundary_tangent)
    boundary_u, _, _ = np.linalg.svd(boundary_tangent, full_matrices=True)
    boundary_normal_basis = boundary_u[:, boundary_rank:]
    background_boundary = np.asarray([background[edge] for edge in BOUNDARY_EDGES])
    point = np.concatenate(
        (
            np.asarray([background[edge] for edge in INTERNAL_EDGES]),
            np.zeros(boundary_normal_basis.shape[1]),
        )
    )

    def lengths_from_coordinates(coordinates: np.ndarray) -> np.ndarray:
        lengths = varied_lengths(
            background, INTERNAL_EDGES, coordinates[: len(INTERNAL_EDGES)]
        )
        boundary_values = (
            background_boundary
            + boundary_normal_basis @ coordinates[len(INTERNAL_EDGES) :]
        )
        return varied_lengths(lengths, BOUNDARY_EDGES, boundary_values)

    def action_from_coordinates(coordinates: np.ndarray) -> complex:
        return general_complex_action(
            lengths_from_coordinates(coordinates),
            SIMPLICES,
            topology["bulk_internal_triangles"],
        )

    hessian = richardson_complex_hessian(
        action_from_coordinates, point, hessian_step
    )
    internal_count = len(INTERNAL_EDGES)
    central_displacement_tangent = np.asarray(
        [
            2.0 * (points[first] - points[second]) @ MINKOWSKI_METRIC
            for first, second in INTERNAL_EDGES
        ]
    )
    gauge_u, _, _ = np.linalg.svd(central_displacement_tangent, full_matrices=True)
    physical_internal_basis = gauge_u[:, 4:]
    physical_internal_hessian = (
        physical_internal_basis.T
        @ hessian[:internal_count, :internal_count]
        @ physical_internal_basis
    )
    physical_boundary_cross = (
        physical_internal_basis.T @ hessian[:internal_count, internal_count:]
    )
    internal_response = -np.linalg.solve(
        physical_internal_hessian, physical_boundary_cross
    )

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

    def curvature_jacobian(step_scale: float) -> np.ndarray:
        result = np.zeros((len(internal_triangles), len(point)), dtype=complex)
        for column in range(len(point)):
            step = step_scale * max(1.0, abs(point[column]))
            displacement = np.zeros_like(point)
            displacement[column] = step
            result[:, column] = (
                curvatures(point + displacement)
                - curvatures(point - displacement)
            ) / (2.0 * step)
        return result

    curvature_jacobian_ = (
        4.0 * curvature_jacobian(curvature_step)
        - curvature_jacobian(2.0 * curvature_step)
    ) / 3.0
    on_shell_curvature = (
        curvature_jacobian_[:, internal_count:]
        + curvature_jacobian_[:, :internal_count]
        @ physical_internal_basis
        @ internal_response
    )

    hinge_riemann = [
        _hinge_normal_riemann(triangle, points)
        for triangle in internal_triangles
    ]
    hinge_einstein = []
    for riemann in hinge_riemann:
        ricci = np.einsum("bd,abcd->ac", MINKOWSKI_METRIC, riemann)
        scalar = np.einsum("ac,ac", MINKOWSKI_METRIC, ricci)
        hinge_einstein.append(ricci - 0.5 * MINKOWSKI_METRIC * scalar)
    einstein_map = np.column_stack(
        [_symmetric_components(tensor) for tensor in hinge_einstein]
    )
    riemann_map = np.column_stack(
        [riemann.reshape(-1) for riemann in hinge_riemann]
    )
    unsigned_areas = np.asarray(
        [
            _unsigned_circumcentric_dual_area(triangle, points, SIMPLICES)
            for triangle in internal_triangles
        ]
    )
    oriented_area_magnitudes = np.asarray(
        [
            _oriented_flag_dual_area_magnitude(triangle, points, SIMPLICES)
            for triangle in internal_triangles
        ]
    )

    screens = {}
    for name, weights in (
        ("unit", np.ones(len(internal_triangles))),
        ("unsigned", 1.0 / unsigned_areas),
        ("oriented_magnitude", 1.0 / oriented_area_magnitudes),
    ):
        weighted_einstein = einstein_map @ np.diag(weights) @ on_shell_curvature
        _, singular_values, right = np.linalg.svd(
            weighted_einstein, full_matrices=True
        )
        rank = int(
            np.count_nonzero(singular_values > singular_values[0] * 1e-5)
        )
        vacuum_axes = right.T[:, rank:]
        vacuum_riemann = [
            (
                riemann_map
                @ np.diag(weights)
                @ on_shell_curvature
                @ vacuum_axes[:, column]
            ).reshape(4, 4, 4, 4)
            for column in range(vacuum_axes.shape[1])
        ]
        screens[name] = {
            "weights": weights,
            "einstein_singular_values": singular_values,
            "einstein_rank": rank,
            "vacuum_dimension": vacuum_axes.shape[1],
            "vacuum_boundary_axes": vacuum_axes,
            "vacuum_riemann": vacuum_riemann,
            "radiation": _radiation_screen(vacuum_riemann),
        }

    simplex_eigenvalues = []
    for simplex in SIMPLICES:
        base = simplex[0]
        tangents = points[list(simplex[1:])] - points[base]
        simplex_eigenvalues.append(
            np.linalg.eigvalsh(tangents @ MINKOWSKI_METRIC @ tangents.T)
        )
    return {
        "points": points,
        "internal_triangles": internal_triangles,
        "boundary_normal_basis": boundary_normal_basis,
        "on_shell_curvature": on_shell_curvature,
        "hinge_riemann": hinge_riemann,
        "hinge_einstein_map": einstein_map,
        "hinge_riemann_map": riemann_map,
        "topology": topology,
        "bulk_f_vector": topology["bulk_f_vector"],
        "boundary_f_vector": topology["boundary_f_vector"],
        "boundary_length_tangent_rank": boundary_rank,
        "boundary_nonembeddability_dimension": boundary_normal_basis.shape[1],
        "full_length_tangent_rank": _rank(length_tangent(points, ALL_EDGES)),
        "full_nonembeddability_dimension": len(ALL_EDGES)
        - _rank(length_tangent(points, ALL_EDGES)),
        "maximum_flat_curvature_residual": float(
            max(
                abs(
                    branch_curvature(
                        triangle,
                        background,
                        topology["bulk_internal_triangles"],
                        SIMPLICES,
                    )
                )
                for triangle in internal_triangles
            )
        ),
        "minimum_negative_eigenvalue_magnitude": min(
            abs(eigenvalues[0]) for eigenvalues in simplex_eigenvalues
        ),
        "minimum_positive_eigenvalue": min(
            eigenvalues[1] for eigenvalues in simplex_eigenvalues
        ),
        "physical_internal_hessian_singular_values": np.linalg.svd(
            physical_internal_hessian, compute_uv=False
        ),
        "physical_boundary_cross_singular_values": np.linalg.svd(
            physical_boundary_cross, compute_uv=False
        ),
        "on_shell_curvature_singular_values": np.linalg.svd(
            on_shell_curvature, compute_uv=False
        ),
        "unsigned_dual_areas": unsigned_areas,
        "oriented_flag_dual_area_magnitudes": oriented_area_magnitudes,
        "screens": screens,
        "inverse_dual_screens_are_diagnostic_only": True,
        "vertex_hybrid_screen_implemented_elsewhere": True,
        "causal_boundary_decomposition_available": True,
        "closed_boundary_not_causal_cobordism": False,
        "deformed_witness": not np.array_equal(points, POINTS),
    }


if __name__ == "__main__":
    certificate = weighted_weyl_star_certificate()
    print(f"bulk_f_vector={certificate['bulk_f_vector']}")
    print(f"boundary_f_vector={certificate['boundary_f_vector']}")
    print(
        "boundary_nonembeddability_dimension="
        f"{certificate['boundary_nonembeddability_dimension']}"
    )
    print(
        "on_shell_curvature_singular_values="
        f"{certificate['on_shell_curvature_singular_values']}"
    )
    for name, screen in certificate["screens"].items():
        print(f"{name}_einstein_singular_values={screen['einstein_singular_values']}")
        print(f"{name}_vacuum_dimension={screen['vacuum_dimension']}")
        print(f"{name}_right_residuals={screen['radiation']['right_residuals']}")
