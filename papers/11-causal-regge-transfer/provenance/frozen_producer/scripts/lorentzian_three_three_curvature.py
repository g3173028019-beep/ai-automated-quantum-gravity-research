"""Minimal Lorentzian curvature carrier for Phase-First Cartan-Regge Gravity.

The 4-2 cell has one internal edge, so every infinitesimal boundary variation
can be returned to the flat six-point hypersurface by solving for that edge.
The 3-3 cell is different: all fifteen edge lengths are boundary data and no
bulk edge is available to absorb the unique non-embeddable variation.  Its
initial triangle (012) and final triangle (345) therefore carry a genuine,
vertex-displacement-invariant linearized curvature response.

This module evaluates the complete complex Lorentzian action directly from
signed lengths and compares its Hessian with the primary-source factorized
3-3 formula.  It establishes a curvature vertex, not yet a continuum graviton
propagator or a refinement-complete quantum theory.
"""

from itertools import combinations

import numpy as np

try:
    from scripts.complex_lorentzian_regge import (
        ALL_EDGES,
        MINKOWSKI_METRIC,
        complex_dihedral_angle_from_lengths,
        complex_regge_action,
        incident_angle_sum,
        richardson_complex_hessian,
        squared_intervals,
        triangle_area,
        triangles_of,
        varied_lengths,
    )
    from scripts.lorentzian_refoliation_geometry import (
        induced_gram,
        outward_unit_normal,
    )
    from scripts.refoliation_graviton_history import (
        boundary_triangles as slice_boundary_triangles,
        complex_faces,
        euler_characteristic,
        replace_patch,
    )
    from scripts.first_order_angle_measure import (
        LOCAL_PAIRS,
        _DET_GRAM,
        _GRAD_DET_GRAM,
        _HESS_DET_GRAM,
    )
except ModuleNotFoundError:  # Allow direct execution from scripts/.
    from complex_lorentzian_regge import (
        ALL_EDGES,
        MINKOWSKI_METRIC,
        complex_dihedral_angle_from_lengths,
        complex_regge_action,
        incident_angle_sum,
        richardson_complex_hessian,
        squared_intervals,
        triangle_area,
        triangles_of,
        varied_lengths,
    )
    from lorentzian_refoliation_geometry import induced_gram, outward_unit_normal
    from refoliation_graviton_history import (
        boundary_triangles as slice_boundary_triangles,
        complex_faces,
        euler_characteristic,
        replace_patch,
    )
    from first_order_angle_measure import (
        LOCAL_PAIRS,
        _DET_GRAM,
        _GRAD_DET_GRAM,
        _HESS_DET_GRAM,
    )


Simplex4 = tuple[int, int, int, int, int]
Edge = tuple[int, int]

INITIAL_SIMPLICES: tuple[Simplex4, ...] = tuple(
    tuple(vertex for vertex in range(6) if vertex != omitted)
    for omitted in (3, 4, 5)
)
FINAL_SIMPLICES: tuple[Simplex4, ...] = tuple(
    tuple(vertex for vertex in range(6) if vertex != omitted)
    for omitted in (0, 1, 2)
)
INITIAL_INTERNAL_TRIANGLE = (0, 1, 2)
FINAL_INTERNAL_TRIANGLE = (3, 4, 5)
SIGN_BITS = np.asarray((0, 0, 0, 1, 1, 1), dtype=int)
PROBE_EDGES: tuple[Edge, ...] = ((0, 1), (3, 4), (2, 5))

THREE_THREE_POINTS = np.asarray(
    [
        [-0.06848396700269202, -0.4771498255239843, 0.0394530162178891, -0.20711723122589912],
        [-0.10845072326415406, 0.24713395111570013, -0.6901546536898029, 0.0025724667916458524],
        [-0.09307367772316104, 0.0498709991322696, 0.5320940458584346, -0.19278098792090692],
        [-0.02680923568740418, 0.06287851834835519, 0.06854169614319873, -0.3693729481277757],
        [-0.00851757174463905, 0.8080450777295256, -0.29117507236301454, 0.512133019430966],
        [-0.3078494287353175, -0.6907787208018662, 0.34124096783329494, 0.25456568105196986],
    ],
    dtype=float,
)

# A second 3-3 witness chosen for the complex angle-Gram Schur calculation.
# It sacrifices the particular A22 shelling labels but improves the weakest
# auxiliary relative singular value by about four orders of magnitude.  Keeping
# it independent prevents a thin auxiliary block from masquerading as a failure
# of the Lorentzian first-order algebra.
FIRST_ORDER_THREE_THREE_POINTS = np.asarray(
    [
        [0.005635860999710608, -0.19992076564728137, 0.025950285350680752, 0.3934684108294953],
        [0.08907382720436609, 0.6260513082214478, 0.25434343862185377, 0.17488684774405525],
        [-0.11645655895748219, -0.10845802590864886, -0.11260918116684002, -0.5933795487033677],
        [-0.033377923194932976, 0.040333398245549705, 0.23180408011841014, 0.44207294454825924],
        [0.12710992287687248, 0.12418489202724327, -0.5819424008630142, 0.43980760226642346],
        [-0.12216125802498172, -0.4821908069383106, 0.1824537779389094, -0.8568562566848654],
    ],
    dtype=float,
)

INITIAL_SLICE = {
    (0, 2, 3, 4),
    (0, 1, 3, 4),
    (0, 2, 4, 5),
    (0, 1, 4, 5),
}
FINAL_SLICE = {
    (1, 2, 3, 4),
    (1, 2, 4, 5),
    (1, 2, 3, 5),
    (0, 2, 3, 5),
    (0, 1, 3, 5),
}
CAUSAL_SHELLING: tuple[tuple[Simplex4, set[tuple[int, ...]], set[tuple[int, ...]]], ...] = (
    (
        (0, 1, 2, 3, 4),
        {(0, 2, 3, 4), (0, 1, 3, 4)},
        {(1, 2, 3, 4), (0, 1, 2, 4), (0, 1, 2, 3)},
    ),
    (
        (0, 1, 2, 4, 5),
        {(0, 2, 4, 5), (0, 1, 4, 5), (0, 1, 2, 4)},
        {(1, 2, 4, 5), (0, 1, 2, 5)},
    ),
    (
        (0, 1, 2, 3, 5),
        {(0, 1, 2, 5), (0, 1, 2, 3)},
        {(1, 2, 3, 5), (0, 2, 3, 5), (0, 1, 3, 5)},
    ),
)


def initial_action(lengths: np.ndarray) -> complex:
    return complex_regge_action(
        INITIAL_SIMPLICES, (INITIAL_INTERNAL_TRIANGLE,), lengths
    )


def final_action(lengths: np.ndarray) -> complex:
    return complex_regge_action(FINAL_SIMPLICES, (FINAL_INTERNAL_TRIANGLE,), lengths)


def pachner_action(lengths: np.ndarray) -> complex:
    return initial_action(lengths) - final_action(lengths)


def lorentzian_curvature(
    triangle: tuple[int, int, int],
    simplices: tuple[Simplex4, ...],
    lengths: np.ndarray,
) -> float:
    deficit = 2.0 * np.pi + incident_angle_sum(triangle, simplices, lengths)
    curvature = -1j * deficit
    if abs(curvature.imag) > 2e-9:
        raise ValueError("3-3 curvature left the regular Lorentzian sheet")
    return float(curvature.real)


def source_geometry(points: np.ndarray, lengths: np.ndarray) -> dict:
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
    if not (np.all(circuit[:3] > 0.0) and np.all(circuit[3:] < 0.0)):
        raise ValueError("witness does not realize the selected 3-3 partition")
    volumes = []
    for omitted in range(6):
        simplex_points = np.delete(points, omitted, axis=0)
        volumes.append(
            abs(np.linalg.det((simplex_points[1:] - simplex_points[0]).T)) / 24.0
        )
    volumes = np.asarray(volumes)
    # Eq. D_ij with i=0,j=1 and (-1)^(s_0+s_k).
    d_factor = float(
        sum(
            (-1.0) ** int(SIGN_BITS[vertex])
            * (
                lengths[0, vertex]
                + lengths[1, vertex]
                - lengths[0, 1]
            )
            * volumes[vertex]
            / 8.0
            for vertex in range(2, 6)
        )
    )

    def hessian_entry(first: Edge, second: Edge) -> float:
        i, j = first
        m, n = second
        sign = (-1.0) ** int(
            SIGN_BITS[i] + SIGN_BITS[j] + SIGN_BITS[m] + SIGN_BITS[n]
        )
        return float(
            d_factor
            * sign
            * volumes[i]
            * volumes[j]
            * volumes[m]
            * volumes[n]
            / (48.0 * np.prod(volumes))
        )

    hessian = np.asarray(
        [
            [hessian_entry(first, second) for second in PROBE_EDGES]
            for first in PROBE_EDGES
        ]
    )
    return {
        "circuit": circuit,
        "volumes": volumes,
        "d_factor": d_factor,
        "hessian": hessian,
    }


def lorentzian_angle_gram_saddle(
    simplex: Simplex4, lengths: np.ndarray
) -> dict:
    """Complex Barrett angle-Gram saddle for one Lorentzian simplex."""
    angles = []
    areas = []
    for first, second in LOCAL_PAIRS:
        triangle = tuple(
            vertex
            for position, vertex in enumerate(simplex)
            if position not in (first, second)
        )
        angles.append(
            complex_dihedral_angle_from_lengths(simplex, triangle, lengths)
        )
        areas.append(triangle_area(triangle, lengths))
    angles = np.asarray(angles, dtype=complex)
    areas = np.asarray(areas, dtype=complex)
    determinant = complex(_DET_GRAM(*angles))
    gradient = np.asarray(_GRAD_DET_GRAM(*angles), dtype=complex).reshape(-1)
    hessian = np.asarray(_HESS_DET_GRAM(*angles), dtype=complex)
    multiplier = np.vdot(gradient, areas) / np.vdot(gradient, gradient)
    stationarity_residual = np.linalg.norm(areas - multiplier * gradient) / np.linalg.norm(
        areas
    )
    # For I=iS, use I_aux=-mu det(G): A-mu grad(det G)=0.
    auxiliary_block = np.block(
        [
            [-multiplier * hessian, -gradient[:, None]],
            [-gradient[None, :], np.zeros((1, 1), dtype=complex)],
        ]
    )
    return {
        "simplex": simplex,
        "angles": angles,
        "areas": areas,
        "gram_determinant": determinant,
        "gradient": gradient,
        "hessian": hessian,
        "multiplier": multiplier,
        "stationarity_residual": float(stationarity_residual),
        "auxiliary_block": auxiliary_block,
        "condition_number": float(np.linalg.cond(auxiliary_block)),
    }


def lorentzian_first_order_transfer_certificate(
    relative_step: float = 3e-5,
) -> dict:
    """Eliminate complex angle/Gram variables for the causal 3-3 history."""
    background = squared_intervals(FIRST_ORDER_THREE_THREE_POINTS)
    physical_edges = ((0, 1), (0, 2), (2, 3))
    point = np.asarray([background[edge] for edge in physical_edges])
    saddles = [
        lorentzian_angle_gram_saddle(simplex, background)
        for simplex in INITIAL_SIMPLICES
    ]

    all_triangles = set().union(
        *(triangles_of(simplex) for simplex in INITIAL_SIMPLICES)
    )
    boundary_constants = {}
    for triangle in all_triangles:
        if triangle == INITIAL_INTERNAL_TRIANGLE:
            boundary_constants[triangle] = 2.0 * np.pi
        else:
            angle_sum = incident_angle_sum(
                triangle, INITIAL_SIMPLICES, background
            )
            boundary_constants[triangle] = (
                int(np.rint(-angle_sum.real / np.pi)) * np.pi
            )

    def fixed_auxiliary_i_action(values: np.ndarray) -> complex:
        lengths = varied_lengths(background, physical_edges, values)
        result = sum(
            boundary_constants[triangle] * triangle_area(triangle, lengths)
            for triangle in all_triangles
        )
        for saddle in saddles:
            simplex = saddle["simplex"]
            for angle, (first, second) in zip(saddle["angles"], LOCAL_PAIRS):
                triangle = tuple(
                    vertex
                    for position, vertex in enumerate(simplex)
                    if position not in (first, second)
                )
                result += angle * triangle_area(triangle, lengths)
            result -= saddle["multiplier"] * saddle["gram_determinant"]
        return complex(result)

    fixed_i_hessian = richardson_complex_hessian(
        fixed_auxiliary_i_action, point, relative_step
    )

    def area_derivative(
        simplex: Simplex4, edge_position: int, pair_position: int, step_scale: float
    ) -> complex:
        edge = physical_edges[edge_position]
        step = step_scale * max(1.0, abs(background[edge]))
        plus = np.array(background, copy=True)
        minus = np.array(background, copy=True)
        plus[edge] = plus[edge[::-1]] = background[edge] + step
        minus[edge] = minus[edge[::-1]] = background[edge] - step
        first, second = LOCAL_PAIRS[pair_position]
        triangle = tuple(
            vertex
            for position, vertex in enumerate(simplex)
            if position not in (first, second)
        )
        return (triangle_area(triangle, plus) - triangle_area(triangle, minus)) / (
            2.0 * step
        )

    cross_blocks = []
    auxiliary_blocks = []
    for saddle in saddles:
        simplex = saddle["simplex"]
        area_cross = np.zeros((len(physical_edges), 10), dtype=complex)
        for edge_position in range(len(physical_edges)):
            for pair_position in range(10):
                fine = area_derivative(
                    simplex, edge_position, pair_position, relative_step
                )
                coarse = area_derivative(
                    simplex, edge_position, pair_position, 2.0 * relative_step
                )
                area_cross[edge_position, pair_position] = (
                    4.0 * fine - coarse
                ) / 3.0
        cross_blocks.append(
            np.hstack((area_cross, np.zeros((len(physical_edges), 1), dtype=complex)))
        )
        auxiliary_blocks.append(saddle["auxiliary_block"])

    cross = np.hstack(cross_blocks)
    auxiliary = np.zeros((33, 33), dtype=complex)
    for index, block in enumerate(auxiliary_blocks):
        auxiliary[11 * index : 11 * (index + 1), 11 * index : 11 * (index + 1)] = block
    reduced_i_hessian = fixed_i_hessian - cross @ np.linalg.solve(auxiliary, cross.T)
    reduced_s_hessian = -1j * reduced_i_hessian

    def direct_action(values: np.ndarray) -> complex:
        return initial_action(
            varied_lengths(background, physical_edges, values)
        )

    direct_s_hessian = richardson_complex_hessian(
        direct_action, point, relative_step
    )
    return {
        "physical_edges": physical_edges,
        "saddles": saddles,
        "fixed_i_hessian": fixed_i_hessian,
        "cross": cross,
        "auxiliary": auxiliary,
        "reduced_s_hessian": reduced_s_hessian,
        "direct_s_hessian": direct_s_hessian,
        "maximum_hessian_residual": np.max(
            np.abs(reduced_s_hessian - direct_s_hessian)
        ),
        "maximum_gram_determinant": max(
            abs(saddle["gram_determinant"]) for saddle in saddles
        ),
        "maximum_stationarity_residual": max(
            saddle["stationarity_residual"] for saddle in saddles
        ),
        "maximum_reduced_imaginary_part": np.max(
            np.abs(reduced_s_hessian.imag)
        ),
        "minimum_auxiliary_relative_singular_value": min(
            np.linalg.svd(block, compute_uv=False)[-1]
            / np.linalg.svd(block, compute_uv=False)[0]
            for block in auxiliary_blocks
        ),
    }


def vertex_length_tangent(points: np.ndarray) -> np.ndarray:
    tangent = np.zeros((len(ALL_EDGES), 24), dtype=float)
    for row, (first, second) in enumerate(ALL_EDGES):
        difference = points[first] - points[second]
        covector = 2.0 * difference @ MINKOWSKI_METRIC
        tangent[row, 4 * first : 4 * first + 4] = covector
        tangent[row, 4 * second : 4 * second + 4] = -covector
    return tangent


def curvature_vertex_certificate(
    hessian_step: float = 3e-5, curvature_step: float = 1e-5
) -> dict:
    points = THREE_THREE_POINTS
    background = squared_intervals(points)
    probe_point = np.asarray([background[edge] for edge in PROBE_EDGES])

    def evaluate(values: np.ndarray) -> complex:
        return pachner_action(
            varied_lengths(background, PROBE_EDGES, values)
        )

    direct_hessian = richardson_complex_hessian(
        evaluate, probe_point, hessian_step
    )
    source = source_geometry(points, background)

    all_point = np.asarray([background[edge] for edge in ALL_EDGES])

    def curvature_jacobian(step_scale: float) -> np.ndarray:
        steps = step_scale * np.maximum(1.0, np.abs(all_point))
        result = np.zeros((2, len(ALL_EDGES)))
        for column, step in enumerate(steps):
            displacement = np.zeros_like(all_point)
            displacement[column] = step
            plus = varied_lengths(
                background, ALL_EDGES, all_point + displacement
            )
            minus = varied_lengths(
                background, ALL_EDGES, all_point - displacement
            )
            result[0, column] = (
                lorentzian_curvature(
                    INITIAL_INTERNAL_TRIANGLE, INITIAL_SIMPLICES, plus
                )
                - lorentzian_curvature(
                    INITIAL_INTERNAL_TRIANGLE, INITIAL_SIMPLICES, minus
                )
            ) / (2.0 * step)
            result[1, column] = (
                lorentzian_curvature(
                    FINAL_INTERNAL_TRIANGLE, FINAL_SIMPLICES, plus
                )
                - lorentzian_curvature(
                    FINAL_INTERNAL_TRIANGLE, FINAL_SIMPLICES, minus
                )
            ) / (2.0 * step)
        return result

    curvature_jacobian_ = (
        4.0 * curvature_jacobian(curvature_step)
        - curvature_jacobian(2.0 * curvature_step)
    ) / 3.0
    tangent = vertex_length_tangent(points)
    tangent_u, tangent_singular_values, _ = np.linalg.svd(
        tangent, full_matrices=True
    )
    flatness_normal = tangent_u[:, -1]
    if flatness_normal[0] < 0.0:
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
        cayley_gradient.append(2.0 * cofactor)
    cayley_normal = np.asarray(cayley_gradient)
    cayley_normal /= np.linalg.norm(cayley_normal)
    if cayley_normal @ flatness_normal < 0.0:
        cayley_normal = -cayley_normal
    projected = np.outer(
        curvature_jacobian_ @ flatness_normal, flatness_normal
    )
    curvature_singular_values = np.linalg.svd(
        curvature_jacobian_, compute_uv=False
    )
    return {
        "background_initial_curvature": lorentzian_curvature(
            INITIAL_INTERNAL_TRIANGLE, INITIAL_SIMPLICES, background
        ),
        "background_final_curvature": lorentzian_curvature(
            FINAL_INTERNAL_TRIANGLE, FINAL_SIMPLICES, background
        ),
        "initial_action": initial_action(background),
        "final_action": final_action(background),
        "pachner_action": pachner_action(background),
        "direct_hessian": direct_hessian,
        "source_hessian": source["hessian"],
        "maximum_source_residual": np.max(
            np.abs(direct_hessian - source["hessian"])
        ),
        "maximum_hessian_imaginary_part": np.max(
            np.abs(direct_hessian.imag)
        ),
        "source_circuit": source["circuit"],
        "source_volumes": source["volumes"],
        "source_d_factor": source["d_factor"],
        "curvature_jacobian": curvature_jacobian_,
        "curvature_singular_values": curvature_singular_values,
        "vertex_tangent_rank": int(
            np.count_nonzero(
                tangent_singular_values > tangent_singular_values[0] * 1e-10
            )
        ),
        "curvature_tangent_residual": np.max(
            np.abs(curvature_jacobian_ @ tangent)
        ),
        "rank_one_curvature_residual": np.max(
            np.abs(curvature_jacobian_ - projected)
        ),
        "physical_normal_response": curvature_jacobian_ @ flatness_normal,
        "cayley_menger_determinant": float(np.linalg.det(cayley_menger)),
        "cayley_menger_normal_residual": np.max(
            np.abs(cayley_normal - flatness_normal)
        ),
        "no_internal_edge_to_eliminate": True,
    }


def causal_curvature_history_certificate() -> dict:
    """Shell the initial 3-3 complex as a time-oriented graviton history.

    The three four-simplices enact 2->3, 3->2, and 2->3 moves.  Edge (12) is
    born in the first move, edge (04) is removed by the dynamical middle move,
    and edge (35) is born in the last move.  The bulk has an internal curvature
    triangle (012) but no internal edge, so its curvature cannot be flattened by
    a Regge edge equation inside this region.
    """
    states = [set(INITIAL_SLICE)]
    for _, past_patch, future_patch in CAUSAL_SHELLING:
        states.append(replace_patch(states[-1], past_patch, future_patch))
    if states[-1] != FINAL_SLICE:
        raise ValueError("3-3 causal shelling reached the wrong final slice")

    state_edges = [complex_faces(state, 2) for state in states]
    state_corner_triangles = [slice_boundary_triangles(state) for state in states]
    if any(
        corner != state_corner_triangles[0]
        for corner in state_corner_triangles[1:]
    ):
        raise ValueError("canonical curvature history changed its fixed corner")
    corner_edges = complex_faces(state_corner_triangles[0], 2)
    state_interior_edges = [edges - corner_edges for edges in state_edges]
    new_edges = [state_edges[index + 1] - state_edges[index] for index in range(3)]
    removed_edges = [state_edges[index] - state_edges[index + 1] for index in range(3)]

    directed_times = []
    tetrahedron_margins = []
    face_normals: dict[tuple[int, ...], list[np.ndarray]] = {}
    for simplex, past_patch, future_patch in CAUSAL_SHELLING:
        if past_patch | future_patch != {
            tuple(vertex for vertex in simplex if vertex != opposite)
            for opposite in simplex
        }:
            raise ValueError("past/future faces do not partition a shelling simplex")
        for opposite in simplex:
            face = tuple(vertex for vertex in simplex if vertex != opposite)
            eigenvalues = np.linalg.eigvalsh(induced_gram(face, THREE_THREE_POINTS))
            if eigenvalues[0] <= 0.0:
                raise ValueError("shelling face is not spacelike")
            tetrahedron_margins.append(float(eigenvalues[0]))
            normal = outward_unit_normal(face, opposite, THREE_THREE_POINTS)
            face_normals.setdefault(face, []).append(normal)
            directed_time = -normal[0] if face in past_patch else normal[0]
            if directed_time <= 0.0:
                raise ValueError("shelling face has the wrong time orientation")
            directed_times.append(float(directed_time))
    shared_normal_residuals = [
        np.linalg.norm(normals[0] + normals[1])
        for normals in face_normals.values()
        if len(normals) == 2
    ]

    bulk_simplices = set(INITIAL_SIMPLICES)
    tetrahedron_counts: dict[tuple[int, ...], int] = {}
    for simplex in bulk_simplices:
        for face in combinations(simplex, 4):
            tetrahedron_counts[face] = tetrahedron_counts.get(face, 0) + 1
    boundary_tetrahedra = {
        face for face, count in tetrahedron_counts.items() if count == 1
    }
    boundary_edges = complex_faces(boundary_tetrahedra, 2)
    bulk_edges = complex_faces(bulk_simplices, 2)
    bulk_boundary_triangles = complex_faces(boundary_tetrahedra, 3)
    bulk_triangles = complex_faces(bulk_simplices, 3)

    background = squared_intervals(THREE_THREE_POINTS)
    transfer_edges = ((0, 4), (1, 2), (3, 5))
    transfer_point = np.asarray([background[edge] for edge in transfer_edges])

    def transfer_action(values: np.ndarray) -> complex:
        return initial_action(
            varied_lengths(background, transfer_edges, values)
        )

    transfer_hessian = richardson_complex_hessian(
        transfer_action, transfer_point, 3e-5
    )
    curvature = curvature_vertex_certificate()
    curvature_columns = [ALL_EDGES.index(edge) for edge in transfer_edges]
    return {
        "states": states,
        "state_edges": state_edges,
        "state_corner_triangles": state_corner_triangles,
        "corner_edges": corner_edges,
        "state_interior_edges": state_interior_edges,
        "state_euler_characteristics": [euler_characteristic(state) for state in states],
        "new_edges": new_edges,
        "removed_edges": removed_edges,
        "minimum_directed_normal_time_component": min(directed_times),
        "minimum_tetrahedron_gram_eigenvalue": min(tetrahedron_margins),
        "maximum_shared_normal_residual": max(shared_normal_residuals),
        "boundary_tetrahedra": boundary_tetrahedra,
        "bulk_internal_edges": bulk_edges - boundary_edges,
        "bulk_internal_triangles": bulk_triangles - bulk_boundary_triangles,
        "transfer_edges": transfer_edges,
        "transfer_hessian": transfer_hessian,
        "input_to_surviving_graviton_response": transfer_hessian[0, 1],
        "input_to_late_born_edge_response": transfer_hessian[0, 2],
        "physical_initial_to_final_block": transfer_hessian[0:1, 1:3],
        "physical_transfer_rank": int(
            np.linalg.matrix_rank(transfer_hessian[0:1, 1:3], tol=1e-7)
        ),
        "curvature_response_columns": curvature["curvature_jacobian"][:, curvature_columns],
        "physical_normal_response": curvature["physical_normal_response"],
    }


def canonical_graviton_relation_certificate() -> dict:
    """Reduce the A22 mixed Hessian to one symplectic propagating relation.

    The final relative edge space is two-dimensional, while the initial one is
    one-dimensional.  Rotate the final coordinates so the first axis is the
    image of the mixed block.  The orthogonal axis is the graviton newly born in
    the last 2->3 move and carries a post-constraint rather than an incoming
    datum.  At fixed value of that new coordinate, the old-to-surviving sector
    is an exact one-degree-of-freedom symplectic map.
    """
    history = causal_curvature_history_certificate()
    hessian = history["transfer_hessian"].real
    initial_block = float(hessian[0, 0])
    mixed = np.asarray(hessian[1:, 0], dtype=float)
    final_block = np.asarray(hessian[1:, 1:], dtype=float)
    coupling = float(np.linalg.norm(mixed))
    if coupling <= 1e-10:
        raise ValueError("the curvature history has no propagating mixed block")
    propagating_axis = mixed / coupling
    new_axis = np.asarray((-propagating_axis[1], propagating_axis[0]))
    rotation = np.column_stack((propagating_axis, new_axis))
    rotated_final = rotation.T @ final_block @ rotation
    c_pp = rotated_final[0, 0]
    transfer = np.asarray(
        [
            [-initial_block / coupling, -1.0 / coupling],
            [
                coupling - c_pp * initial_block / coupling,
                -c_pp / coupling,
            ],
        ]
    )
    symplectic_form = np.asarray(((0.0, 1.0), (-1.0, 0.0)))
    return {
        "initial_hessian": initial_block,
        "mixed_block": mixed,
        "coupling": coupling,
        "propagating_axis": propagating_axis,
        "new_axis": new_axis,
        "rotated_final_hessian": rotated_final,
        "transfer_matrix": transfer,
        "transfer_determinant": float(np.linalg.det(transfer)),
        "symplectic_residual": np.max(
            np.abs(transfer.T @ symplectic_form @ transfer - symplectic_form)
        ),
        "new_mode_postconstraint": {
            "p_new_minus_c_new_prop_q_prop_minus_c_new_new_q_new": (
                1.0,
                -rotated_final[1, 0],
                -rotated_final[1, 1],
            )
        },
        "incoming_mode_count": 1,
        "outgoing_mode_count": 2,
        "propagating_rank": history["physical_transfer_rank"],
    }


if __name__ == "__main__":
    certificate = curvature_vertex_certificate()
    for key in (
        "background_initial_curvature",
        "background_final_curvature",
        "initial_action",
        "final_action",
        "pachner_action",
        "source_d_factor",
        "maximum_source_residual",
        "maximum_hessian_imaginary_part",
        "curvature_singular_values",
        "vertex_tangent_rank",
        "curvature_tangent_residual",
        "rank_one_curvature_residual",
        "physical_normal_response",
    ):
        print(f"{key}={certificate[key]}")
    print(f"direct_hessian=\n{certificate['direct_hessian']}")
    print(f"source_hessian=\n{certificate['source_hessian']}")
    history = causal_curvature_history_certificate()
    print(f"state_tetrahedron_counts={[len(state) for state in history['states']]}")
    print(f"new_edges={history['new_edges']}")
    print(f"removed_edges={history['removed_edges']}")
    print(
        "input_to_surviving_graviton_response="
        f"{history['input_to_surviving_graviton_response']}"
    )
    first_order = lorentzian_first_order_transfer_certificate()
    print(
        "first_order_hessian_residual="
        f"{first_order['maximum_hessian_residual']:.3e}"
    )
    print(
        "first_order_stationarity_residual="
        f"{first_order['maximum_stationarity_residual']:.3e}"
    )
    canonical = canonical_graviton_relation_certificate()
    print(f"canonical_transfer_matrix=\n{canonical['transfer_matrix']}")
    print(f"canonical_symplectic_residual={canonical['symplectic_residual']:.3e}")
