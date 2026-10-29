"""Direct off-flat complex Lorentzian Regge action for the PCRG 4-2 cell.

Unlike the source-mapped factorized Hessian certificate, this module starts
only from signed squared edge lengths.  Every four-simplex is reconstructed
from its length Gram matrix, its complex theta-plus angles are evaluated on
the log-minus branch, and the complete bulk plus boundary action is summed.

The initial four-simplex complex has internal edge (01) and four internal
triangles (01k).  The final two-simplex complex has no bulk triangles.  Both
have the same boundary.  Therefore their action difference is the 4-2
Pachner action used in the primary Lorentzian Hessian calculation.
"""

from itertools import combinations

import numpy as np

try:
    from scripts.lorentzian_refoliation_geometry import (
        A7_TO_HISTORY,
        LORENTZIAN_POINTS,
        MINKOWSKI_METRIC,
        complex_log_minus,
        complex_sqrt_plus,
        source_mapped_pachner_hessian_entry,
    )
except ModuleNotFoundError:  # Allow direct execution from scripts/.
    from lorentzian_refoliation_geometry import (
        A7_TO_HISTORY,
        LORENTZIAN_POINTS,
        MINKOWSKI_METRIC,
        complex_log_minus,
        complex_sqrt_plus,
        source_mapped_pachner_hessian_entry,
    )


Edge = tuple[int, int]
Triangle = tuple[int, int, int]
Simplex4 = tuple[int, int, int, int, int]

A7_INITIAL_SIMPLICES: tuple[Simplex4, ...] = tuple(
    tuple(vertex for vertex in range(6) if vertex != omitted)
    for omitted in (2, 3, 4, 5)
)
A7_FINAL_SIMPLICES: tuple[Simplex4, ...] = tuple(
    tuple(vertex for vertex in range(6) if vertex != omitted)
    for omitted in (0, 1)
)
A7_INTERNAL_TRIANGLES: tuple[Triangle, ...] = tuple(
    (0, 1, vertex) for vertex in range(2, 6)
)
PROBE_EDGES: tuple[Edge, ...] = ((0, 1), (2, 3), (4, 5))
ALL_EDGES: tuple[Edge, ...] = tuple(combinations(range(6), 2))

# A deliberately well-conditioned A7-label 4-2 witness.  It was constructed as
# a subluminal nonlinear lifting of six spatial points.  All edges and all
# tetrahedra are spacelike, every 4-simplex is Lorentzian, and its weakest
# Lorentzian Gram eigenvalue is roughly three orders of magnitude farther from
# zero (relatively) than the separately certified causal-history witness.  The
# latter remains the time-orientation certificate; this one is the appropriate
# independent differentiability/Hessian certificate.
OFF_FLAT_A7_POINTS = np.asarray(
    [
        [0.14153799984636897, 0.0704534019317674, 0.5547781057237790, -0.31253081201068444],
        [0.26204024245515173, -0.6302594477133445, -0.5745232648414881, 0.3544049846832336],
        [0.16369896144992005, 0.0343589436046405, -0.2936824523065840, -0.7319912012559121],
        [0.1754152007565262, -0.22347608745019792, 0.7968138525244398, -0.09981831244564365],
        [0.3204398028208946, 0.8804881446454673, -0.3602273302691944, 0.30818322093506545],
        [0.07499549082052236, -0.13156495501833282, -0.12315891083095226, 0.4817521200939412],
    ],
    dtype=float,
)


def normalize(vertices) -> tuple[int, ...]:
    return tuple(sorted(vertices))


def squared_intervals(points: np.ndarray) -> np.ndarray:
    """Return the symmetric matrix of signed squared Minkowski intervals."""
    points = np.asarray(points, dtype=float)
    result = np.zeros((len(points), len(points)), dtype=float)
    for first, second in combinations(range(len(points)), 2):
        difference = points[second] - points[first]
        value = float(difference @ MINKOWSKI_METRIC @ difference)
        result[first, second] = result[second, first] = value
    return result


def a7_witness_lengths() -> np.ndarray:
    """Signed lengths of the well-conditioned off-flat A7 witness."""
    return squared_intervals(OFF_FLAT_A7_POINTS)


def causal_history_a7_lengths() -> np.ndarray:
    """A7 relabeling of the separate time-oriented four-move witness."""
    points = np.asarray(
        [LORENTZIAN_POINTS[A7_TO_HISTORY[vertex]] for vertex in range(6)]
    )
    return squared_intervals(points)


def varied_lengths(
    background: np.ndarray, edges: tuple[Edge, ...], values: np.ndarray
) -> np.ndarray:
    result = np.array(background, copy=True)
    for edge, value in zip(edges, values):
        result[edge] = result[edge[::-1]] = value
    return result


def simplex_length_gram(simplex: Simplex4, lengths: np.ndarray) -> np.ndarray:
    """Affine length Gram matrix based at the first simplex vertex."""
    base = simplex[0]
    others = simplex[1:]
    return np.asarray(
        [
            [
                (lengths[base, first] + lengths[base, second] - lengths[first, second])
                / 2.0
                for second in others
            ]
            for first in others
        ],
        dtype=float,
    )


def simplex_coordinates_from_lengths(
    simplex: Simplex4, lengths: np.ndarray
) -> dict[int, np.ndarray]:
    """Construct one local Minkowski realization from a signed-length Gram matrix.

    The scalar action is invariant under the arbitrary local O(1,3) frame, so
    the spectral realization is sufficient and introduces no gluing data.
    """
    gram = simplex_length_gram(simplex, lengths)
    eigenvalues, eigenvectors = np.linalg.eigh(gram)
    negative = eigenvalues < 0.0
    if np.count_nonzero(negative) != 1 or np.count_nonzero(~negative) != 3:
        raise ValueError(f"simplex {simplex} left the nondegenerate Lorentzian sector")
    # eigh orders the one negative eigenvalue first.  Rows of coordinates obey
    # X eta X^T = U diag(lambda) U^T = G.
    coordinates = eigenvectors @ np.diag(np.sqrt(np.abs(eigenvalues)))
    result = {simplex[0]: np.zeros(4, dtype=float)}
    result.update(
        {vertex: coordinates[index] for index, vertex in enumerate(simplex[1:])}
    )
    reconstruction = coordinates @ MINKOWSKI_METRIC @ coordinates.T
    if not np.allclose(reconstruction, gram, rtol=2e-10, atol=2e-12):
        raise ArithmeticError("spectral Minkowski reconstruction failed")
    return result


def local_minkowski_product(first: np.ndarray, second: np.ndarray) -> float:
    return float(first @ MINKOWSKI_METRIC @ second)


def local_outward_unit_normal(
    face: tuple[int, int, int, int],
    opposite: int,
    coordinates: dict[int, np.ndarray],
) -> np.ndarray:
    origin = coordinates[face[0]]
    tangents = np.asarray(
        [coordinates[vertex] - origin for vertex in face[1:]], dtype=float
    )
    normal = np.linalg.svd(tangents @ MINKOWSKI_METRIC)[2][-1]
    if local_minkowski_product(normal, coordinates[opposite] - origin) > 0.0:
        normal = -normal
    norm_squared = local_minkowski_product(normal, normal)
    if norm_squared >= -1e-12:
        raise ValueError("the required tetrahedral facet is not spacelike")
    return normal / np.sqrt(-norm_squared)


def complex_dihedral_angle_from_lengths(
    simplex: Simplex4, hinge: Triangle, lengths: np.ndarray
) -> complex:
    """Evaluate theta-plus from signed lengths on the primary-source branch."""
    coordinates = simplex_coordinates_from_lengths(simplex, lengths)
    opposite = [vertex for vertex in simplex if vertex not in hinge]
    if len(opposite) != 2:
        raise ValueError("hinge must be a triangle of the four-simplex")
    first_face = tuple(vertex for vertex in simplex if vertex != opposite[0])
    second_face = tuple(vertex for vertex in simplex if vertex != opposite[1])
    first = local_outward_unit_normal(first_face, opposite[0], coordinates)
    # Reverse the second outward normal to obtain the interior angle.
    second = -local_outward_unit_normal(second_face, opposite[1], coordinates)
    first_norm = local_minkowski_product(first, first)
    second_norm = local_minkowski_product(second, second)
    product = local_minkowski_product(first, second)
    discriminant = first_norm * second_norm - product**2
    # The facet normals are Hodge-dual to the projected edge vectors used in
    # the primary-source theta-plus formula.  In Lorentzian signature that
    # duality reverses the square-root term.  The plus sign below is therefore
    # the theta-plus angle of the projected wedge (using a minus sign here
    # produces theta-minus).  The opposite-edge derivative provides an
    # independent check: d theta-plus/ds_opposite =
    # -sqrt(V_h)/(24 sqrt(V_sigma)).
    argument = (
        product + 1j * complex_sqrt_plus(discriminant)
    ) / (complex_sqrt_plus(first_norm) * complex_sqrt_plus(second_norm))
    return complex(-1j * complex_log_minus(argument))


def triangle_area(triangle: Triangle, lengths: np.ndarray) -> complex:
    """Complex principal area from the signed squared side lengths."""
    first, second, third = triangle
    a = lengths[first, second]
    b = lengths[first, third]
    c = lengths[second, third]
    area_squared = (2.0 * (a * b + a * c + b * c) - a**2 - b**2 - c**2) / 16.0
    return complex_sqrt_plus(area_squared)


def triangles_of(simplex: Simplex4) -> set[Triangle]:
    return {normalize(triangle) for triangle in combinations(simplex, 3)}


def incident_angle_sum(
    triangle: Triangle, simplices: tuple[Simplex4, ...], lengths: np.ndarray
) -> complex:
    return sum(
        (
            complex_dihedral_angle_from_lengths(simplex, triangle, lengths)
            for simplex in simplices
            if set(triangle) <= set(simplex)
        ),
        0j,
    )


def complex_regge_action(
    simplices: tuple[Simplex4, ...],
    internal_triangles: tuple[Triangle, ...],
    lengths: np.ndarray,
) -> complex:
    """Return S from iS=sum_h sqrt(V_h) epsilon_h, including the boundary."""
    all_triangles = set().union(*(triangles_of(simplex) for simplex in simplices))
    internal = set(internal_triangles)
    if not internal <= all_triangles:
        raise ValueError("declared bulk triangle is absent from the complex")
    i_action = 0j
    for triangle in sorted(all_triangles):
        angle_sum = incident_angle_sum(triangle, simplices, lengths)
        if triangle in internal:
            constant = 2.0 * np.pi
        else:
            # The boundary integer k counts the branch sectors meeting at this
            # corner and is not uniformly one.  Within a light-cone-regular
            # component it is locally constant and is fixed by cancellation of
            # the discrete real (-pi) crossings.  This makes the Lorentzian
            # boundary action real rather than injecting a spurious i*pi*A.
            boundary_k = int(np.rint(-angle_sum.real / np.pi))
            if boundary_k not in (0, 1, 2):
                raise ValueError("unexpected boundary light-cone sector")
            constant = boundary_k * np.pi
        curvature = constant + angle_sum
        i_action += triangle_area(triangle, lengths) * curvature
    return -1j * i_action


def initial_action(lengths: np.ndarray) -> complex:
    return complex_regge_action(
        A7_INITIAL_SIMPLICES, A7_INTERNAL_TRIANGLES, lengths
    )


def final_action(lengths: np.ndarray) -> complex:
    return complex_regge_action(A7_FINAL_SIMPLICES, (), lengths)


def pachner_action(lengths: np.ndarray) -> complex:
    return initial_action(lengths) - final_action(lengths)


def internal_lorentzian_curvatures(lengths: np.ndarray) -> np.ndarray:
    """Real Lorentzian curvature boosts -i epsilon at the four bulk hinges."""
    deficits = np.asarray(
        [
            2.0 * np.pi
            + incident_angle_sum(triangle, A7_INITIAL_SIMPLICES, lengths)
            for triangle in A7_INTERNAL_TRIANGLES
        ],
        dtype=complex,
    )
    curvatures = -1j * deficits
    if np.max(np.abs(curvatures.imag)) > 2e-9:
        raise ValueError("curvature left the light-cone-regular Lorentzian sheet")
    return curvatures.real


def numerical_complex_hessian(function, point: np.ndarray, relative_step: float) -> np.ndarray:
    """Second-order central Hessian without discarding branch-sensitive parts."""
    point = np.asarray(point, dtype=float)
    steps = relative_step * np.maximum(1.0, np.abs(point))
    result = np.zeros((len(point), len(point)), dtype=complex)
    center = function(point)
    for first in range(len(point)):
        first_step = np.zeros_like(point)
        first_step[first] = steps[first]
        result[first, first] = (
            function(point + first_step)
            - 2.0 * center
            + function(point - first_step)
        ) / steps[first] ** 2
        for second in range(first):
            second_step = np.zeros_like(point)
            second_step[second] = steps[second]
            value = (
                function(point + first_step + second_step)
                - function(point + first_step - second_step)
                - function(point - first_step + second_step)
                + function(point - first_step - second_step)
            ) / (4.0 * steps[first] * steps[second])
            result[first, second] = result[second, first] = value
    return result


def richardson_complex_hessian(
    function, point: np.ndarray, relative_step: float
) -> np.ndarray:
    """Cancel the leading central-difference error with steps h and 2h."""
    fine = numerical_complex_hessian(function, point, relative_step)
    coarse = numerical_complex_hessian(function, point, 2.0 * relative_step)
    return (4.0 * fine - coarse) / 3.0


def source_four_two_data(points: np.ndarray, lengths: np.ndarray) -> dict:
    """Independent primary-source volume/D/Hessian evaluation in A7 labels."""
    homogeneous = np.vstack((np.ones(6), points.T))
    circuit = np.asarray(
        [
            (-1.0) ** omitted
            * np.linalg.det(np.delete(homogeneous, omitted, axis=1))
            / 24.0
            for omitted in range(6)
        ]
    )
    if circuit[0] < 0.0:
        circuit = -circuit
    if not (np.all(circuit[:2] > 0.0) and np.all(circuit[2:] < 0.0)):
        raise ValueError("witness does not realize the A7 4-2 Radon partition")
    volumes = []
    for omitted in range(6):
        simplex_points = np.delete(points, omitted, axis=0)
        volumes.append(
            abs(np.linalg.det((simplex_points[1:] - simplex_points[0]).T)) / 24.0
        )
    volumes = np.asarray(volumes)
    d_factor = float(
        -sum(
            (
                lengths[0, vertex]
                + lengths[1, vertex]
                - lengths[0, 1]
            )
            * volumes[vertex]
            for vertex in range(2, 6)
        )
        / 8.0
    )
    hessian = np.asarray(
        [
            [
                source_mapped_pachner_hessian_entry(
                    first, second, d_factor, volumes
                )
                for second in PROBE_EDGES
            ]
            for first in PROBE_EDGES
        ]
    )
    return {
        "circuit": circuit,
        "absolute_four_volumes": volumes,
        "d_factor": d_factor,
        "hessian": hessian,
    }


def flatness_curvature_certificate(relative_step: float = 1e-5) -> dict:
    """Separate the 4-2 flatness normal from embedding/gauge tangents.

    Six points in four dimensions have fourteen independent infinitesimal
    signed-length variations generated by moving their vertices (Poincare
    motions form the ten-dimensional coordinate kernel).  The remaining one
    length direction leaves the flat Cayley-Menger hypersurface.  This routine
    verifies that every bulk-curvature derivative is proportional to precisely
    that one normal and that eliminating edge (01) restores flatness for every
    boundary variation.  Thus this cell transports canonical data but carries
    no on-shell bulk curvature by itself.
    """
    points = OFF_FLAT_A7_POINTS
    background = a7_witness_lengths()
    point = np.asarray([background[edge] for edge in ALL_EDGES])
    def central_jacobian(step_scale: float) -> np.ndarray:
        steps = step_scale * np.maximum(1.0, np.abs(point))
        jacobian = np.zeros((len(A7_INTERNAL_TRIANGLES), len(ALL_EDGES)))
        for column, step in enumerate(steps):
            displacement = np.zeros_like(point)
            displacement[column] = step
            plus = internal_lorentzian_curvatures(
                varied_lengths(background, ALL_EDGES, point + displacement)
            )
            minus = internal_lorentzian_curvatures(
                varied_lengths(background, ALL_EDGES, point - displacement)
            )
            jacobian[:, column] = (plus - minus) / (2.0 * step)
        return jacobian

    curvature_jacobian = (
        4.0 * central_jacobian(relative_step)
        - central_jacobian(2.0 * relative_step)
    ) / 3.0

    vertex_tangent = np.zeros((len(ALL_EDGES), 24), dtype=float)
    for row, (first, second) in enumerate(ALL_EDGES):
        difference = points[first] - points[second]
        covector = 2.0 * difference @ MINKOWSKI_METRIC
        vertex_tangent[row, 4 * first : 4 * first + 4] = covector
        vertex_tangent[row, 4 * second : 4 * second + 4] = -covector

    tangent_u, tangent_singular_values, _ = np.linalg.svd(
        vertex_tangent, full_matrices=True
    )
    flatness_normal = tangent_u[:, -1]
    if flatness_normal[ALL_EDGES.index((0, 1))] < 0.0:
        flatness_normal = -flatness_normal
    cayley_menger = np.ones((7, 7), dtype=float)
    cayley_menger[0, 0] = 0.0
    cayley_menger[1:, 1:] = background
    cayley_gradient = []
    for first, second in ALL_EDGES:
        row = first + 1
        column = second + 1
        cofactor = (-1.0) ** (row + column) * np.linalg.det(
            np.delete(np.delete(cayley_menger, row, axis=0), column, axis=1)
        )
        # The edge variable occurs in two symmetric matrix entries.
        cayley_gradient.append(2.0 * cofactor)
    cayley_gradient = np.asarray(cayley_gradient)
    cayley_normal = cayley_gradient / np.linalg.norm(cayley_gradient)
    if cayley_normal @ flatness_normal < 0.0:
        cayley_normal = -cayley_normal
    projected_jacobian = np.outer(
        curvature_jacobian @ flatness_normal, flatness_normal
    )
    curvature_singular_values = np.linalg.svd(
        curvature_jacobian, compute_uv=False
    )

    internal_index = ALL_EDGES.index((0, 1))
    on_shell_curvature_responses = []
    for boundary_index in range(len(ALL_EDGES)):
        if boundary_index == internal_index:
            continue
        variation = np.zeros(len(ALL_EDGES))
        variation[boundary_index] = 1.0
        variation[internal_index] = -(
            flatness_normal[boundary_index] / flatness_normal[internal_index]
        )
        on_shell_curvature_responses.append(curvature_jacobian @ variation)
    return {
        "background_curvatures": internal_lorentzian_curvatures(background),
        "curvature_jacobian": curvature_jacobian,
        "curvature_singular_values": curvature_singular_values,
        "vertex_tangent_singular_values": tangent_singular_values,
        "vertex_tangent_rank": int(
            np.count_nonzero(tangent_singular_values > tangent_singular_values[0] * 1e-10)
        ),
        "flatness_normal": flatness_normal,
        "cayley_menger_determinant": float(np.linalg.det(cayley_menger)),
        "cayley_menger_normal": cayley_normal,
        "cayley_menger_normal_residual": np.max(
            np.abs(cayley_normal - flatness_normal)
        ),
        "curvature_tangent_residual": np.max(
            np.abs(curvature_jacobian @ vertex_tangent)
        ),
        "rank_one_curvature_residual": np.max(
            np.abs(curvature_jacobian - projected_jacobian)
        ),
        "maximum_on_shell_curvature_response": np.max(
            np.abs(on_shell_curvature_responses)
        ),
        "fixed_internal_boundary_curvature_response": np.max(
            np.abs(np.delete(curvature_jacobian, internal_index, axis=1))
        ),
    }


def off_flat_action_certificate(relative_step: float = 3e-5) -> dict:
    """Directly compare the differentiated complex action with the source Hessian."""
    background = a7_witness_lengths()
    point = np.asarray([background[edge] for edge in PROBE_EDGES])

    def evaluate(action, values: np.ndarray) -> complex:
        return action(varied_lengths(background, PROBE_EDGES, values))

    initial_hessian = richardson_complex_hessian(
        lambda values: evaluate(initial_action, values), point, relative_step
    )
    final_hessian = richardson_complex_hessian(
        lambda values: evaluate(final_action, values), point, relative_step
    )
    pachner_hessian = richardson_complex_hessian(
        lambda values: evaluate(pachner_action, values), point, relative_step
    )
    source = source_four_two_data(OFF_FLAT_A7_POINTS, background)
    source_hessian = source["hessian"]
    source_residuals = abs(pachner_hessian - source_hessian)
    effective_initial = (
        initial_hessian[1:, 1:]
        - np.outer(initial_hessian[1:, 0], initial_hessian[0, 1:])
        / initial_hessian[0, 0]
    )

    derivative_simplex = A7_INITIAL_SIMPLICES[0]
    derivative_hinge = tuple(derivative_simplex[:3])
    opposite_edge = tuple(
        vertex for vertex in derivative_simplex if vertex not in derivative_hinge
    )

    def angle_derivative(step_scale: float) -> complex:
        step = step_scale * max(1.0, abs(background[opposite_edge]))
        plus = np.array(background, copy=True)
        minus = np.array(background, copy=True)
        plus[opposite_edge] = plus[opposite_edge[::-1]] = (
            background[opposite_edge] + step
        )
        minus[opposite_edge] = minus[opposite_edge[::-1]] = (
            background[opposite_edge] - step
        )
        return (
            complex_dihedral_angle_from_lengths(
                derivative_simplex, derivative_hinge, plus
            )
            - complex_dihedral_angle_from_lengths(
                derivative_simplex, derivative_hinge, minus
            )
        ) / (2.0 * step)

    direct_angle_derivative = (
        4.0 * angle_derivative(1e-5) - angle_derivative(2e-5)
    ) / 3.0
    simplex_signed_volume_square = (
        np.linalg.det(simplex_length_gram(derivative_simplex, background))
        / 24.0**2
    )
    source_angle_derivative = (
        -triangle_area(derivative_hinge, background)
        / (24.0 * complex_sqrt_plus(simplex_signed_volume_square))
    )
    return {
        "background": background,
        "point": point,
        "initial_action": evaluate(initial_action, point),
        "final_action": evaluate(final_action, point),
        "pachner_action": evaluate(pachner_action, point),
        "initial_hessian": initial_hessian,
        "final_hessian": final_hessian,
        "pachner_hessian": pachner_hessian,
        "source_hessian": source_hessian,
        "source_circuit": source["circuit"],
        "source_absolute_four_volumes": source["absolute_four_volumes"],
        "source_d_factor": source["d_factor"],
        "source_residuals": source_residuals,
        "maximum_source_residual": np.max(source_residuals),
        "fixed_internal_separated_response": initial_hessian[1, 2],
        "effective_initial_boundary_hessian": effective_initial,
        "final_boundary_hessian": final_hessian[1:, 1:],
        "schur_final_residual": abs(effective_initial[0, 1] - final_hessian[1, 2]),
        "action_difference_residual": abs(
            pachner_hessian - (initial_hessian - final_hessian)
        ).max(),
        "maximum_hessian_imaginary_part": max(
            np.max(np.abs(initial_hessian.imag)),
            np.max(np.abs(final_hessian.imag)),
            np.max(np.abs(pachner_hessian.imag)),
        ),
        "direct_opposite_edge_angle_derivative": direct_angle_derivative,
        "source_opposite_edge_angle_derivative": source_angle_derivative,
        "opposite_edge_angle_derivative_residual": abs(
            direct_angle_derivative - source_angle_derivative
        ),
    }


if __name__ == "__main__":
    certificate = off_flat_action_certificate()
    print(f"initial_action={certificate['initial_action']}")
    print(f"final_action={certificate['final_action']}")
    print(f"pachner_action={certificate['pachner_action']}")
    print(f"direct_pachner_hessian=\n{certificate['pachner_hessian']}")
    print(f"source_hessian=\n{certificate['source_hessian']}")
    print(f"source_residuals={certificate['source_residuals']}")
    print(
        "fixed_internal_separated_response="
        f"{certificate['fixed_internal_separated_response']}"
    )
    print(
        "schur_final_residual="
        f"{certificate['schur_final_residual']:.3e}"
    )
    print(
        "maximum_hessian_imaginary_part="
        f"{certificate['maximum_hessian_imaginary_part']:.3e}"
    )
