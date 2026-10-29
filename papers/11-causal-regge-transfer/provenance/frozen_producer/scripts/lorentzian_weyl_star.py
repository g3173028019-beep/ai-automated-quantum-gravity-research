"""A connected Lorentzian Regge carrier for the A25 hinge-wave target.

The boundary is the join C4*C4, a triangulated three-sphere with eight
vertices.  Coning it from one interior vertex gives a sixteen-simplex
four-ball.  A factorized flat Minkowski embedding realizes the nine A25
spacelike hinge planes among its twenty-four internal triangles and leaves
exactly two non-embeddable boundary-length directions.

The boundary is closed rather than a causal initial/final cobordism.  This
module is therefore a geometric and dynamical intermediate target, not a
retarded propagation certificate.
"""

from collections import Counter
from itertools import combinations, product

import numpy as np

try:
    from scripts.complex_lorentzian_regge import (
        MINKOWSKI_METRIC,
        complex_log_minus,
        complex_sqrt_plus,
        simplex_coordinates_from_lengths,
        simplex_length_gram,
        squared_intervals,
        triangle_area,
        triangles_of,
        varied_lengths,
        richardson_complex_hessian,
    )
    from scripts.refoliation_graviton_history import complex_faces
    from scripts.algebraic_regge_wave_target import (
        RAPIDITIES as TARGET_RAPIDITIES,
        TRANSVERSE_DIRECTIONS as TARGET_TRANSVERSE_DIRECTIONS,
        hinge_tensors as target_hinge_tensors,
        electric_magnetic_weyl,
    )
except ModuleNotFoundError:
    from complex_lorentzian_regge import (
        MINKOWSKI_METRIC,
        complex_log_minus,
        complex_sqrt_plus,
        simplex_coordinates_from_lengths,
        simplex_length_gram,
        squared_intervals,
        triangle_area,
        triangles_of,
        varied_lengths,
        richardson_complex_hessian,
    )
    from refoliation_graviton_history import complex_faces
    from algebraic_regge_wave_target import (
        RAPIDITIES as TARGET_RAPIDITIES,
        TRANSVERSE_DIRECTIONS as TARGET_TRANSVERSE_DIRECTIONS,
        hinge_tensors as target_hinge_tensors,
        electric_magnetic_weyl,
    )


Edge = tuple[int, int]
Triangle = tuple[int, int, int]
Simplex4 = tuple[int, int, int, int, int]

E0, EX, EY, EZ = np.eye(4)
EP = (EX + EY) / np.sqrt(2.0)
RP = (EX - EY) / np.sqrt(2.0)


def q_direction(rapidity: float) -> np.ndarray:
    """Unit spacelike vector in the t-z plane, orthogonal to t(rapidity)."""
    return np.sinh(rapidity) * E0 + np.cosh(rapidity) * EZ


# Convex cyclic orders in two complementary two-planes.  Each polygon contains
# the origin strictly in its interior.  Signs do not change a hinge plane.
Q_VERTICES = (1, 2, 3, 4)
R_VERTICES = (5, 6, 7, 8)
POINTS = np.asarray(
    (
        np.zeros(4),
        q_direction(0.0),
        -q_direction(1.0),
        -q_direction(-0.5),
        q_direction(0.5),
        EY,
        -EP,
        RP,
        EX,
    )
)

Q_CYCLE_EDGES: tuple[Edge, ...] = tuple(
    tuple(sorted((Q_VERTICES[index], Q_VERTICES[(index + 1) % 4])))
    for index in range(4)
)
R_CYCLE_EDGES: tuple[Edge, ...] = tuple(
    tuple(sorted((R_VERTICES[index], R_VERTICES[(index + 1) % 4])))
    for index in range(4)
)
CROSS_EDGES: tuple[Edge, ...] = tuple(
    tuple(sorted(edge)) for edge in product(Q_VERTICES, R_VERTICES)
)
BOUNDARY_EDGES: tuple[Edge, ...] = tuple(
    sorted(set(Q_CYCLE_EDGES + R_CYCLE_EDGES + CROSS_EDGES))
)
INTERNAL_EDGES: tuple[Edge, ...] = tuple((0, vertex) for vertex in range(1, 9))
ALL_EDGES: tuple[Edge, ...] = BOUNDARY_EDGES + INTERNAL_EDGES

SIMPLICES: tuple[Simplex4, ...] = tuple(
    tuple(
        sorted(
            (
                0,
                Q_VERTICES[q_index],
                Q_VERTICES[(q_index + 1) % 4],
                R_VERTICES[r_index],
                R_VERTICES[(r_index + 1) % 4],
            )
        )
    )
    for q_index in range(4)
    for r_index in range(4)
)

# Nine cross hinges reproduce the A25 rapidities (0,1,1/2) and transverse
# normals (x,y,p).  q vertex 3 is the additional -1/2 direction and r vertex
# 6 supplies the additional transverse direction required to close both cycles.
# This order exactly matches A25: rapidities (1/2,1,0), each followed by
# transverse normals (x,y,p).  Tangent r vertices (5,8,7) are respectively
# orthogonal to those transverse directions.
TARGET_Q_VERTICES = (4, 2, 1)
TARGET_R_VERTICES = (5, 8, 7)
TARGET_HINGES: tuple[Triangle, ...] = tuple(
    tuple(sorted((0, q_vertex, r_vertex)))
    for q_vertex in TARGET_Q_VERTICES
    for r_vertex in TARGET_R_VERTICES
)

# A larger closed star with five vertices in each cycle.  The non-unit radial
# scales make all five t-z rays extreme points of a convex polygon while
# retaining the A25 rapidity lines.  The cyclic orders below contain the origin
# strictly in their interiors.
FIVE_Q_VERTICES = (1, 2, 3, 4, 5)
FIVE_R_VERTICES = (6, 7, 8, 9, 10)
FIVE_Q_POINTS = (
    -0.77197331 * q_direction(0.5),
    -0.31528597 * q_direction(0.0),
    0.83915891 * q_direction(1.0),
    0.24140403 * q_direction(-0.5),
    -0.17554963 * q_direction(1.5),
)
FIVE_R_ANGLES = np.deg2rad((-135.0, -45.0, 0.0, 90.0, 135.0))
FIVE_R_POINTS = tuple(
    np.cos(angle) * EX + np.sin(angle) * EY for angle in FIVE_R_ANGLES
)
FIVE_POINTS = np.asarray((np.zeros(4),) + FIVE_Q_POINTS + FIVE_R_POINTS)
FIVE_Q_CYCLE_EDGES: tuple[Edge, ...] = tuple(
    tuple(
        sorted(
            (
                FIVE_Q_VERTICES[index],
                FIVE_Q_VERTICES[(index + 1) % 5],
            )
        )
    )
    for index in range(5)
)
FIVE_R_CYCLE_EDGES: tuple[Edge, ...] = tuple(
    tuple(
        sorted(
            (
                FIVE_R_VERTICES[index],
                FIVE_R_VERTICES[(index + 1) % 5],
            )
        )
    )
    for index in range(5)
)
FIVE_CROSS_EDGES: tuple[Edge, ...] = tuple(
    tuple(sorted(edge)) for edge in product(FIVE_Q_VERTICES, FIVE_R_VERTICES)
)
FIVE_BOUNDARY_EDGES: tuple[Edge, ...] = tuple(
    sorted(
        set(FIVE_Q_CYCLE_EDGES + FIVE_R_CYCLE_EDGES + FIVE_CROSS_EDGES)
    )
)
FIVE_INTERNAL_EDGES: tuple[Edge, ...] = tuple(
    (0, vertex) for vertex in range(1, 11)
)
FIVE_ALL_EDGES: tuple[Edge, ...] = FIVE_BOUNDARY_EDGES + FIVE_INTERNAL_EDGES
FIVE_SIMPLICES: tuple[Simplex4, ...] = tuple(
    tuple(
        sorted(
            (
                0,
                FIVE_Q_VERTICES[q_index],
                FIVE_Q_VERTICES[(q_index + 1) % 5],
                FIVE_R_VERTICES[r_index],
                FIVE_R_VERTICES[(r_index + 1) % 5],
            )
        )
    )
    for q_index in range(5)
    for r_index in range(5)
)
FIVE_TARGET_Q_VERTICES = (1, 3, 2)  # rapidities (1/2,1,0)
FIVE_TARGET_R_VERTICES = (9, 8, 7)  # transverse normals (x,y,p)
FIVE_TARGET_HINGES: tuple[Triangle, ...] = tuple(
    tuple(sorted((0, q_vertex, r_vertex)))
    for q_vertex in FIVE_TARGET_Q_VERTICES
    for r_vertex in FIVE_TARGET_R_VERTICES
)


def _topology_for(simplices: tuple[Simplex4, ...]) -> dict:
    tetrahedron_counts = Counter(
        tuple(sorted(face))
        for simplex in simplices
        for face in combinations(simplex, 4)
    )
    boundary_tetrahedra = {
        face for face, count in tetrahedron_counts.items() if count == 1
    }
    bulk_edges = complex_faces(set(simplices), 2)
    boundary_edges = complex_faces(boundary_tetrahedra, 2)
    bulk_triangles = complex_faces(set(simplices), 3)
    boundary_triangles = complex_faces(boundary_tetrahedra, 3)
    bulk_tetrahedra = complex_faces(set(simplices), 4)
    bulk_vertices = complex_faces(set(simplices), 1)
    return {
        "boundary_tetrahedra": boundary_tetrahedra,
        "bulk_internal_edges": bulk_edges - boundary_edges,
        "bulk_internal_triangles": bulk_triangles - boundary_triangles,
        "bulk_f_vector": (
            len(bulk_vertices),
            len(bulk_edges),
            len(bulk_triangles),
            len(bulk_tetrahedra),
            len(simplices),
        ),
        "boundary_f_vector": (
            len(complex_faces(boundary_tetrahedra, 1)),
            len(boundary_edges),
            len(boundary_triangles),
            len(boundary_tetrahedra),
        ),
        "bulk_euler_characteristic": (
            len(bulk_vertices)
            - len(bulk_edges)
            + len(bulk_triangles)
            - len(bulk_tetrahedra)
            + len(simplices)
        ),
        "boundary_euler_characteristic": (
            len(complex_faces(boundary_tetrahedra, 1))
            - len(boundary_edges)
            + len(boundary_triangles)
            - len(boundary_tetrahedra)
        ),
        "maximum_tetrahedron_incidence": max(tetrahedron_counts.values()),
    }


def _topology() -> dict:
    return _topology_for(SIMPLICES)


def length_tangent(points: np.ndarray, edges: tuple[Edge, ...]) -> np.ndarray:
    tangent = np.zeros((len(edges), 4 * len(points)), dtype=float)
    for row, (first, second) in enumerate(edges):
        difference = points[first] - points[second]
        covector = 2.0 * difference @ MINKOWSKI_METRIC
        tangent[row, 4 * first : 4 * first + 4] = covector
        tangent[row, 4 * second : 4 * second + 4] = -covector
    return tangent


def _rank(matrix: np.ndarray, tolerance: float = 1e-10) -> int:
    singular_values = np.linalg.svd(matrix, compute_uv=False)
    return int(np.count_nonzero(singular_values > singular_values[0] * tolerance))


def _hinge_normal_riemann(
    triangle: Triangle, points: np.ndarray = POINTS
) -> np.ndarray:
    origin = points[triangle[0]]
    tangents = np.asarray(
        [points[vertex] - origin for vertex in triangle[1:]], dtype=float
    )
    normal_plane = np.linalg.svd(tangents @ MINKOWSKI_METRIC)[2][-2:].T
    normal_gram = normal_plane.T @ MINKOWSKI_METRIC @ normal_plane
    eigenvalues, eigenvectors = np.linalg.eigh(normal_gram)
    if np.min(np.abs(eigenvalues)) < 1e-12:
        raise ValueError("hinge normal plane is null or degenerate")
    orthonormal = normal_plane @ eigenvectors @ np.diag(
        1.0 / np.sqrt(np.abs(eigenvalues))
    )
    bivector = np.outer(orthonormal[:, 0], orthonormal[:, 1]) - np.outer(
        orthonormal[:, 1], orthonormal[:, 0]
    )
    covariant = MINKOWSKI_METRIC @ bivector @ MINKOWSKI_METRIC
    return np.einsum("ab,cd->abcd", covariant, covariant)


def _minkowski_circumcenter(
    face: tuple[int, ...], points: np.ndarray
) -> np.ndarray:
    vertices = points[list(face)]
    tangents = vertices[1:] - vertices[0]
    gram = tangents @ MINKOWSKI_METRIC @ tangents.T
    coefficients = np.linalg.solve(gram, 0.5 * np.diag(gram))
    return vertices[0] + coefficients @ tangents


def _unsigned_circumcentric_dual_area(
    triangle: Triangle,
    points: np.ndarray,
    simplices: tuple[Simplex4, ...],
) -> float:
    """Sum magnitudes of flag triangles h<circ(tau)<circ(sigma).

    This is a diagnostic only.  The physical Lorentzian Ricci weight requires
    oriented signed flag areas and hybrid-volume conventions; absolute values
    are returned here so that failure cannot be hidden by an arbitrary sign.
    """
    triangle_center = _minkowski_circumcenter(triangle, points)
    area = 0.0
    for simplex in simplices:
        if not set(triangle) <= set(simplex):
            continue
        simplex_center = _minkowski_circumcenter(simplex, points)
        for extra in set(simplex) - set(triangle):
            tetrahedron = tuple(sorted((*triangle, extra)))
            tetrahedron_center = _minkowski_circumcenter(tetrahedron, points)
            first = tetrahedron_center - triangle_center
            second = simplex_center - triangle_center
            gram_determinant = (
                (first @ MINKOWSKI_METRIC @ first)
                * (second @ MINKOWSKI_METRIC @ second)
                - (first @ MINKOWSKI_METRIC @ second) ** 2
            )
            area += 0.5 * np.sqrt(abs(gram_determinant))
    return float(area)


def _oriented_flag_dual_area_magnitude(
    triangle: Triangle,
    points: np.ndarray,
    simplices: tuple[Simplex4, ...],
) -> float:
    """Norm of the incidence-oriented circumcentric flag bivector sum.

    This advances beyond the unsigned screen but still does not fix the final
    Lorentzian scalar sign/hybrid-volume convention.  It is kept as a second
    failure-oriented diagnostic, not advertised as the completed Ricci weight.
    """
    triangle_center = _minkowski_circumcenter(triangle, points)
    dual_bivector = np.zeros((4, 4), dtype=float)
    for simplex in simplices:
        if not set(triangle) <= set(simplex):
            continue
        determinant = np.linalg.det(
            (points[list(simplex[1:])] - points[simplex[0]]).T
        )
        simplex_orientation = 1 if determinant > 0.0 else -1
        simplex_center = _minkowski_circumcenter(simplex, points)
        for extra in set(simplex) - set(triangle):
            tetrahedron = tuple(sorted((*triangle, extra)))
            omitted = next(vertex for vertex in simplex if vertex not in tetrahedron)
            simplex_face_position = list(simplex).index(omitted)
            tetrahedron_face_position = list(tetrahedron).index(extra)
            incidence_sign = simplex_orientation * (-1) ** (
                simplex_face_position + tetrahedron_face_position
            )
            tetrahedron_center = _minkowski_circumcenter(
                tetrahedron, points
            )
            first = tetrahedron_center - triangle_center
            second = simplex_center - triangle_center
            dual_bivector += 0.5 * incidence_sign * (
                np.outer(first, second) - np.outer(second, first)
            )
    covariant = MINKOWSKI_METRIC @ dual_bivector @ MINKOWSKI_METRIC
    norm_squared = 0.5 * np.einsum("ab,ab", dual_bivector, covariant)
    return float(np.sqrt(abs(norm_squared)))


def _general_outward_normal(
    face: tuple[int, int, int, int],
    opposite: int,
    coordinates: dict[int, np.ndarray],
) -> np.ndarray:
    origin = coordinates[face[0]]
    tangents = np.asarray(
        [coordinates[vertex] - origin for vertex in face[1:]], dtype=float
    )
    normal = np.linalg.svd(tangents @ MINKOWSKI_METRIC)[2][-1]
    if normal @ MINKOWSKI_METRIC @ (coordinates[opposite] - origin) > 0.0:
        normal = -normal
    norm_squared = float(normal @ MINKOWSKI_METRIC @ normal)
    if abs(norm_squared) < 1e-12:
        raise ValueError("null tetrahedral facet is outside the current branch")
    return normal


def general_complex_dihedral_angle(
    simplex: Simplex4, hinge: Triangle, lengths: np.ndarray
) -> complex:
    """Theta-plus angle allowing either spacelike or timelike tetrahedral facets."""
    coordinates = simplex_coordinates_from_lengths(simplex, lengths)
    opposite = [vertex for vertex in simplex if vertex not in hinge]
    first_face = tuple(vertex for vertex in simplex if vertex != opposite[0])
    second_face = tuple(vertex for vertex in simplex if vertex != opposite[1])
    first = _general_outward_normal(first_face, opposite[0], coordinates)
    second = -_general_outward_normal(second_face, opposite[1], coordinates)
    first_norm = float(first @ MINKOWSKI_METRIC @ first)
    second_norm = float(second @ MINKOWSKI_METRIC @ second)
    inner = float(first @ MINKOWSKI_METRIC @ second)
    discriminant = first_norm * second_norm - inner**2
    # On this mixed-facet star the Schlaefli sheet uses opposite Hodge-dual
    # roots for timelike and spacelike hinges.  The choice is not guessed from
    # flat deficit cancellation: differentiating all ten angles of every
    # simplex gives sum_h A_h d theta_h=0 below 9e-10.  It is kept local to this
    # mixed causal sector; the all-spacelike-facet evaluator in A21 retains its
    # independently opposite-edge-certified branch.
    area = triangle_area(hinge, lengths)
    root_sign = 1.0 if abs(area.imag) > 1e-10 else -1.0
    argument = (
        inner + root_sign * 1j * complex_sqrt_plus(discriminant)
    ) / (complex_sqrt_plus(first_norm) * complex_sqrt_plus(second_norm))
    return complex(-1j * complex_log_minus(argument))


def angle_sum(
    triangle: Triangle,
    lengths: np.ndarray,
    simplices: tuple[Simplex4, ...] = SIMPLICES,
) -> complex:
    return sum(
        (
            general_complex_dihedral_angle(simplex, triangle, lengths)
            for simplex in simplices
            if set(triangle) <= set(simplex)
        ),
        0j,
    )


def branch_curvature(
    triangle: Triangle,
    lengths: np.ndarray,
    internal_triangles: set[Triangle],
    simplices: tuple[Simplex4, ...] = SIMPLICES,
) -> complex:
    """Curvature with the locally constant flat-sheet integer selected by incidence."""
    total = angle_sum(triangle, lengths, simplices)
    period = 2.0 * np.pi if triangle in internal_triangles else np.pi
    branch_integer = int(np.rint(-total.real / period))
    return branch_integer * period + total


def general_complex_action(
    lengths: np.ndarray,
    simplices: tuple[Simplex4, ...] = SIMPLICES,
    internal_triangles: set[Triangle] | None = None,
) -> complex:
    if internal_triangles is None:
        internal_triangles = _topology()["bulk_internal_triangles"]
    all_triangles = set().union(*(triangles_of(simplex) for simplex in simplices))
    i_action = sum(
        triangle_area(triangle, lengths)
        * branch_curvature(triangle, lengths, internal_triangles, simplices)
        for triangle in all_triangles
    )
    return complex(-1j * i_action)


def weyl_star_geometry_certificate() -> dict:
    topology = _topology()
    background = squared_intervals(POINTS)
    boundary_tangent = length_tangent(POINTS[1:], tuple(
        (first - 1, second - 1) for first, second in BOUNDARY_EDGES
    ))
    full_tangent = length_tangent(POINTS, ALL_EDGES)
    simplex_eigenvalues = [
        np.linalg.eigvalsh(simplex_length_gram(simplex, background))
        for simplex in SIMPLICES
    ]
    internal_curvatures = np.asarray(
        [
            branch_curvature(
                triangle, background, topology["bulk_internal_triangles"]
            )
            for triangle in sorted(topology["bulk_internal_triangles"])
        ]
    )
    realized_target_riemann = [
        _hinge_normal_riemann(triangle) for triangle in TARGET_HINGES
    ]
    expected_target_riemann = [
        target_hinge_tensors(rapidity, transverse)["riemann"]
        for rapidity in TARGET_RAPIDITIES
        for transverse in TARGET_TRANSVERSE_DIRECTIONS
    ]
    return {
        **topology,
        "background_lengths": background,
        "background_action": general_complex_action(background),
        "boundary_length_tangent_rank": _rank(boundary_tangent),
        "boundary_nonembeddability_dimension": len(BOUNDARY_EDGES)
        - _rank(boundary_tangent),
        "full_length_tangent_rank": _rank(full_tangent),
        "full_nonembeddability_dimension": len(ALL_EDGES) - _rank(full_tangent),
        "simplex_gram_eigenvalues": simplex_eigenvalues,
        "minimum_negative_eigenvalue_magnitude": min(
            abs(eigenvalues[0]) for eigenvalues in simplex_eigenvalues
        ),
        "minimum_positive_eigenvalue": min(
            eigenvalues[1] for eigenvalues in simplex_eigenvalues
        ),
        "internal_curvatures": internal_curvatures,
        "maximum_internal_curvature_residual": float(
            np.max(np.abs(internal_curvatures))
        ),
        "target_hinges": TARGET_HINGES,
        "target_hinges_are_internal": set(TARGET_HINGES)
        <= topology["bulk_internal_triangles"],
        "target_plane_projector_residual": float(
            max(
                np.max(np.abs(realized - expected))
                for realized, expected in zip(
                    realized_target_riemann, expected_target_riemann
                )
            )
        ),
        "closed_boundary_not_causal_cobordism": True,
    }


def five_cycle_weyl_star_geometry_certificate() -> dict:
    """Geometry/rank certificate for the C5*C5 two-vacuum-mode candidate."""
    topology = _topology_for(FIVE_SIMPLICES)
    background = squared_intervals(FIVE_POINTS)
    relative_boundary_edges = tuple(
        (first - 1, second - 1) for first, second in FIVE_BOUNDARY_EDGES
    )
    boundary_tangent = length_tangent(
        FIVE_POINTS[1:], relative_boundary_edges
    )
    full_tangent = length_tangent(FIVE_POINTS, FIVE_ALL_EDGES)
    simplex_eigenvalues = [
        np.linalg.eigvalsh(simplex_length_gram(simplex, background))
        for simplex in FIVE_SIMPLICES
    ]
    internal_curvatures = np.asarray(
        [
            branch_curvature(
                triangle,
                background,
                topology["bulk_internal_triangles"],
                FIVE_SIMPLICES,
            )
            for triangle in sorted(topology["bulk_internal_triangles"])
        ]
    )
    realized_target_riemann = [
        _hinge_normal_riemann(triangle, FIVE_POINTS)
        for triangle in FIVE_TARGET_HINGES
    ]
    expected_target_riemann = [
        target_hinge_tensors(rapidity, transverse)["riemann"]
        for rapidity in TARGET_RAPIDITIES
        for transverse in TARGET_TRANSVERSE_DIRECTIONS
    ]
    return {
        **topology,
        "background_lengths": background,
        "background_action": general_complex_action(
            background, FIVE_SIMPLICES, topology["bulk_internal_triangles"]
        ),
        "boundary_length_tangent_rank": _rank(boundary_tangent),
        "boundary_nonembeddability_dimension": len(FIVE_BOUNDARY_EDGES)
        - _rank(boundary_tangent),
        "full_length_tangent_rank": _rank(full_tangent),
        "full_nonembeddability_dimension": len(FIVE_ALL_EDGES)
        - _rank(full_tangent),
        "simplex_gram_eigenvalues": simplex_eigenvalues,
        "minimum_negative_eigenvalue_magnitude": min(
            abs(eigenvalues[0]) for eigenvalues in simplex_eigenvalues
        ),
        "minimum_positive_eigenvalue": min(
            eigenvalues[1] for eigenvalues in simplex_eigenvalues
        ),
        "internal_curvatures": internal_curvatures,
        "maximum_internal_curvature_residual": float(
            np.max(np.abs(internal_curvatures))
        ),
        "target_hinges": FIVE_TARGET_HINGES,
        "target_hinges_are_internal": set(FIVE_TARGET_HINGES)
        <= topology["bulk_internal_triangles"],
        "target_plane_projector_residual": float(
            max(
                np.max(np.abs(realized - expected))
                for realized, expected in zip(
                    realized_target_riemann, expected_target_riemann
                )
            )
        ),
        "closed_boundary_not_causal_cobordism": True,
    }


def _symmetric_components(tensor: np.ndarray) -> np.ndarray:
    return np.asarray(
        [tensor[first, first] for first in range(4)]
        + [
            tensor[first, second]
            for first in range(4)
            for second in range(first + 1, 4)
        ]
    )


def five_cycle_weyl_dynamics_certificate(
    hessian_step: float = 8e-5, curvature_step: float = 1e-5
) -> dict:
    """Eliminate the six physical radial modes and screen the five boundary modes.

    The Einstein contraction is deliberately called a *unit-weight screen*.
    A physically normalized simplicial Ricci tensor also requires signed dual
    volumes.  Until that separate calculation is supplied, the two-dimensional
    kernel found here is a candidate vacuum sector rather than a final theorem.
    """
    geometry = five_cycle_weyl_star_geometry_certificate()
    topology = _topology_for(FIVE_SIMPLICES)
    internal_triangles = tuple(sorted(topology["bulk_internal_triangles"]))
    background = geometry["background_lengths"]

    relative_boundary_edges = tuple(
        (first - 1, second - 1) for first, second in FIVE_BOUNDARY_EDGES
    )
    boundary_tangent = length_tangent(
        FIVE_POINTS[1:], relative_boundary_edges
    )
    boundary_u, _, _ = np.linalg.svd(boundary_tangent, full_matrices=True)
    boundary_normal_basis = boundary_u[:, geometry["boundary_length_tangent_rank"] :]
    background_boundary = np.asarray(
        [background[edge] for edge in FIVE_BOUNDARY_EDGES]
    )
    point = np.concatenate(
        (
            np.asarray([background[edge] for edge in FIVE_INTERNAL_EDGES]),
            np.zeros(boundary_normal_basis.shape[1]),
        )
    )

    def lengths_from_coordinates(coordinates: np.ndarray) -> np.ndarray:
        lengths = varied_lengths(
            background, FIVE_INTERNAL_EDGES, coordinates[: len(FIVE_INTERNAL_EDGES)]
        )
        boundary_values = (
            background_boundary
            + boundary_normal_basis @ coordinates[len(FIVE_INTERNAL_EDGES) :]
        )
        return varied_lengths(lengths, FIVE_BOUNDARY_EDGES, boundary_values)

    def action_from_coordinates(coordinates: np.ndarray) -> complex:
        return general_complex_action(
            lengths_from_coordinates(coordinates),
            FIVE_SIMPLICES,
            topology["bulk_internal_triangles"],
        )

    hessian = richardson_complex_hessian(
        action_from_coordinates, point, hessian_step
    )
    internal_count = len(FIVE_INTERNAL_EDGES)
    central_displacement_tangent = np.asarray(
        [
            2.0
            * (FIVE_POINTS[first] - FIVE_POINTS[second])
            @ MINKOWSKI_METRIC
            for first, second in FIVE_INTERNAL_EDGES
        ]
    )
    gauge_u, _, _ = np.linalg.svd(
        central_displacement_tangent, full_matrices=True
    )
    physical_internal_basis = gauge_u[:, 4:]
    internal_hessian = hessian[:internal_count, :internal_count]
    physical_internal_hessian = (
        physical_internal_basis.T
        @ internal_hessian
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
                    FIVE_SIMPLICES,
                )
                for triangle in internal_triangles
            ]
        )

    def curvature_jacobian(step_scale: float) -> np.ndarray:
        result = np.zeros(
            (len(internal_triangles), len(point)), dtype=complex
        )
        for column in range(len(point)):
            step = step_scale * max(1.0, abs(point[column]))
            displacement = np.zeros_like(point)
            displacement[column] = step
            result[:, column] = (
                curvatures(point + displacement)
                - curvatures(point - displacement)
            ) / (2.0 * step)
        return result

    fine_curvature = curvature_jacobian(curvature_step)
    coarse_curvature = curvature_jacobian(2.0 * curvature_step)
    curvature_jacobian_ = (4.0 * fine_curvature - coarse_curvature) / 3.0
    on_shell_curvature = (
        curvature_jacobian_[:, internal_count:]
        + curvature_jacobian_[:, :internal_count]
        @ physical_internal_basis
        @ internal_response
    )

    hinge_riemann = [
        _hinge_normal_riemann(triangle, FIVE_POINTS)
        for triangle in internal_triangles
    ]
    hinge_einstein = []
    for riemann in hinge_riemann:
        ricci = np.einsum("bd,abcd->ac", MINKOWSKI_METRIC, riemann)
        scalar = np.einsum("ac,ac", MINKOWSKI_METRIC, ricci)
        hinge_einstein.append(
            ricci - 0.5 * MINKOWSKI_METRIC * scalar
        )
    unit_weight_einstein_map = np.column_stack(
        [_symmetric_components(tensor) for tensor in hinge_einstein]
    )
    unit_weight_riemann_map = np.column_stack(
        [riemann.reshape(-1) for riemann in hinge_riemann]
    )
    on_shell_einstein = unit_weight_einstein_map @ on_shell_curvature
    unsigned_dual_areas = np.asarray(
        [
            _unsigned_circumcentric_dual_area(
                triangle, FIVE_POINTS, FIVE_SIMPLICES
            )
            for triangle in internal_triangles
        ]
    )
    unsigned_dual_weighted_einstein = (
        unit_weight_einstein_map
        @ np.diag(1.0 / unsigned_dual_areas)
        @ on_shell_curvature
    )
    unsigned_dual_weighted_singular_values = np.linalg.svd(
        unsigned_dual_weighted_einstein, compute_uv=False
    )
    oriented_flag_dual_area_magnitudes = np.asarray(
        [
            _oriented_flag_dual_area_magnitude(
                triangle, FIVE_POINTS, FIVE_SIMPLICES
            )
            for triangle in internal_triangles
        ]
    )
    oriented_flag_weighted_einstein = (
        unit_weight_einstein_map
        @ np.diag(1.0 / oriented_flag_dual_area_magnitudes)
        @ on_shell_curvature
    )
    oriented_flag_weighted_singular_values = np.linalg.svd(
        oriented_flag_weighted_einstein, compute_uv=False
    )
    einstein_u, einstein_singular_values, einstein_vh = np.linalg.svd(
        on_shell_einstein, full_matrices=True
    )
    einstein_rank = int(
        np.count_nonzero(
            einstein_singular_values > einstein_singular_values[0] * 1e-5
        )
    )
    candidate_vacuum_axes = einstein_vh.T[:, einstein_rank:]
    candidate_vacuum_curvatures = on_shell_curvature @ candidate_vacuum_axes
    candidate_vacuum_riemann = [
        (
            unit_weight_riemann_map @ candidate_vacuum_curvatures[:, column]
        ).reshape(4, 4, 4, 4)
        for column in range(candidate_vacuum_curvatures.shape[1])
    ]
    candidate_electric_magnetic = [
        electric_magnetic_weyl(riemann)
        for riemann in candidate_vacuum_riemann
    ]

    invariant_one = np.zeros((len(candidate_vacuum_riemann),) * 2, dtype=complex)
    invariant_two = np.zeros_like(invariant_one)
    for first, (first_electric, first_magnetic) in enumerate(
        candidate_electric_magnetic
    ):
        for second, (second_electric, second_magnetic) in enumerate(
            candidate_electric_magnetic
        ):
            invariant_one[first, second] = np.sum(
                first_electric * second_electric
                - first_magnetic * second_magnetic
            )
            invariant_two[first, second] = 0.5 * np.sum(
                first_electric * second_magnetic
                + second_electric * first_magnetic
            )

    null_wave_roots = np.asarray((), dtype=complex)
    null_wave_electric_magnetic = []
    null_wave_right_residuals = []
    null_wave_left_residuals = []
    null_wave_invariants = []
    if candidate_vacuum_axes.shape[1] == 2:
        null_wave_roots = np.roots(
            (invariant_one[1, 1], 2.0 * invariant_one[0, 1], invariant_one[0, 0])
        )
        for root in null_wave_roots:
            riemann = candidate_vacuum_riemann[0] + root * candidate_vacuum_riemann[1]
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
            null_wave_electric_magnetic.append((electric, magnetic))
            null_wave_right_residuals.append(
                float(np.linalg.norm(magnetic - rotated_electric) / scale)
            )
            null_wave_left_residuals.append(
                float(np.linalg.norm(magnetic + rotated_electric) / scale)
            )
            null_wave_invariants.append(
                np.asarray(
                    (
                        np.sum(electric * electric - magnetic * magnetic),
                        np.sum(electric * magnetic),
                    )
                )
            )

    return {
        "geometry": geometry,
        "boundary_normal_basis": boundary_normal_basis,
        "hessian": hessian,
        "maximum_hessian_imaginary_part": float(np.max(np.abs(hessian.imag))),
        "central_displacement_tangent": central_displacement_tangent,
        "physical_internal_basis": physical_internal_basis,
        "physical_internal_hessian": physical_internal_hessian,
        "physical_internal_hessian_singular_values": np.linalg.svd(
            physical_internal_hessian, compute_uv=False
        ),
        "gauge_hessian_residual": float(
            np.max(np.abs(internal_hessian @ central_displacement_tangent))
        ),
        "physical_boundary_cross": physical_boundary_cross,
        "physical_boundary_cross_singular_values": np.linalg.svd(
            physical_boundary_cross, compute_uv=False
        ),
        "internal_response": internal_response,
        "curvature_jacobian": curvature_jacobian_,
        "on_shell_curvature": on_shell_curvature,
        "on_shell_curvature_singular_values": np.linalg.svd(
            on_shell_curvature, compute_uv=False
        ),
        "unit_weight_einstein_map": unit_weight_einstein_map,
        "on_shell_unit_weight_einstein_map": on_shell_einstein,
        "on_shell_einstein_singular_values": einstein_singular_values,
        "on_shell_einstein_rank": einstein_rank,
        "candidate_vacuum_dimension": candidate_vacuum_axes.shape[1],
        "unsigned_circumcentric_dual_areas": unsigned_dual_areas,
        "unsigned_dual_weighted_einstein_map": unsigned_dual_weighted_einstein,
        "unsigned_dual_weighted_einstein_singular_values": (
            unsigned_dual_weighted_singular_values
        ),
        "unsigned_dual_weighted_einstein_rank": int(
            np.count_nonzero(
                unsigned_dual_weighted_singular_values
                > unsigned_dual_weighted_singular_values[0] * 1e-5
            )
        ),
        "oriented_flag_dual_area_magnitudes": (
            oriented_flag_dual_area_magnitudes
        ),
        "oriented_flag_weighted_einstein_singular_values": (
            oriented_flag_weighted_singular_values
        ),
        "oriented_flag_weighted_einstein_rank": int(
            np.count_nonzero(
                oriented_flag_weighted_singular_values
                > oriented_flag_weighted_singular_values[0] * 1e-5
            )
        ),
        "candidate_vacuum_boundary_axes": candidate_vacuum_axes,
        "candidate_vacuum_curvatures": candidate_vacuum_curvatures,
        "candidate_vacuum_riemann": candidate_vacuum_riemann,
        "candidate_electric_magnetic": candidate_electric_magnetic,
        "candidate_radiation_invariant_one": invariant_one,
        "candidate_radiation_invariant_two": invariant_two,
        "null_wave_roots": null_wave_roots,
        "null_wave_electric_magnetic": null_wave_electric_magnetic,
        "null_wave_right_residuals": np.asarray(null_wave_right_residuals),
        "null_wave_left_residuals": np.asarray(null_wave_left_residuals),
        "null_wave_invariants": np.asarray(null_wave_invariants),
        "dual_volume_weighting_complete": False,
        "closed_boundary_not_causal_cobordism": True,
    }


if __name__ == "__main__":
    certificate = weyl_star_geometry_certificate()
    for key in (
        "bulk_f_vector",
        "boundary_f_vector",
        "bulk_euler_characteristic",
        "boundary_euler_characteristic",
        "bulk_internal_edges",
        "boundary_length_tangent_rank",
        "boundary_nonembeddability_dimension",
        "full_length_tangent_rank",
        "full_nonembeddability_dimension",
        "minimum_negative_eigenvalue_magnitude",
        "minimum_positive_eigenvalue",
        "maximum_internal_curvature_residual",
        "background_action",
        "target_hinges_are_internal",
        "target_plane_projector_residual",
    ):
        print(f"{key}={certificate[key]}")
    five_certificate = five_cycle_weyl_star_geometry_certificate()
    print(f"five_bulk_f_vector={five_certificate['bulk_f_vector']}")
    print(f"five_boundary_f_vector={five_certificate['boundary_f_vector']}")
    print(
        "five_boundary_nonembeddability_dimension="
        f"{five_certificate['boundary_nonembeddability_dimension']}"
    )
    print(
        "five_target_plane_projector_residual="
        f"{five_certificate['target_plane_projector_residual']:.3e}"
    )
