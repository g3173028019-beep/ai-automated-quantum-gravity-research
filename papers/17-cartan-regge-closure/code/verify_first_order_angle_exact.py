"""Exact SymPy certificate for the canonical first-order 4-2 Schur reduction.

All vertex coordinates are rational.  Dihedral sines and cosines are reconstructed from
exact barycentric face normals, so no inverse trigonometric or floating-point operations
enter the certificate.  The script proves one exact non-singular auxiliary example and the
canonical ``H_first/H_coordinate_target = -1/8`` identity.
"""

from itertools import combinations
from functools import lru_cache

import sympy as sp


PAIRS = tuple(combinations(range(5), 2))
ANGLES = sp.symbols("theta_0:10", real=True)
GRAM = sp.eye(5)
for angle, (first, second) in zip(ANGLES, PAIRS):
    GRAM[first, second] = GRAM[second, first] = -sp.cos(angle)
DET_GRAM = sp.det(GRAM)
GRADIENT = sp.Matrix([sp.diff(DET_GRAM, angle) for angle in ANGLES])
HESSIAN = sp.hessian(DET_GRAM, ANGLES)

POINTS = [
    sp.Matrix(
        [sp.Rational(7, 10), sp.Rational(1, 4), sp.Rational(1, 4), sp.Rational(1, 4)]
    ),
    sp.Matrix(
        [sp.Rational(-4, 5), sp.Rational(1, 4), sp.Rational(1, 4), sp.Rational(1, 4)]
    ),
    sp.zeros(4, 1),
    sp.Matrix([0, 1, 0, 0]),
    sp.Matrix([0, 0, 1, 0]),
    sp.Matrix([0, 0, 0, 1]),
]


def squared_length(first: int, second: int):
    return sp.expand(sum((POINTS[first][k] - POINTS[second][k]) ** 2 for k in range(4)))


def exact_area_derivative(triangle: tuple[int, int, int], edge: tuple[int, int]):
    """Derivative of a triangle area with respect to one squared edge length."""
    edge = tuple(sorted(edge))
    if not set(edge) <= set(triangle):
        return sp.S.Zero
    third = next(vertex for vertex in triangle if vertex not in edge)
    q_edge = squared_length(*edge)
    q_first = squared_length(edge[0], third)
    q_second = squared_length(edge[1], third)
    radicand = sp.factor(
        2 * (q_edge * q_first + q_edge * q_second + q_first * q_second)
        - q_edge**2
        - q_first**2
        - q_second**2
    )
    return sp.factor((q_first + q_second - q_edge) / (4 * sp.sqrt(radicand)))


@lru_cache(maxsize=None)
def exact_simplex_auxiliary_data(omitted: int):
    vertices = [vertex for vertex in range(6) if vertex != omitted]
    coordinates = [POINTS[vertex] for vertex in vertices]
    affine = sp.Matrix.hstack(
        *[sp.Matrix.vstack(sp.ones(1, 1), coordinate) for coordinate in coordinates]
    )
    inverse_affine = affine.inv()
    outward_normals = []
    for position in range(5):
        gradient = sp.Matrix(inverse_affine[position, 1:]).T
        outward_normals.append(-gradient / sp.sqrt((gradient.T * gradient)[0]))

    substitutions = {}
    areas = []
    triangles = []
    for position, (first, second) in enumerate(PAIRS):
        cosine = sp.simplify(-(outward_normals[first].T * outward_normals[second])[0])
        substitutions[sp.cos(ANGLES[position])] = cosine
        substitutions[sp.sin(ANGLES[position])] = sp.sqrt(sp.factor(1 - cosine**2))
        triangle = [vertices[k] for k in range(5) if k not in (first, second)]
        triangles.append(tuple(sorted(triangle)))
        edge_a = POINTS[triangle[1]] - POINTS[triangle[0]]
        edge_b = POINTS[triangle[2]] - POINTS[triangle[0]]
        areas.append(
            sp.sqrt(
                sp.factor(
                    (edge_a.T * edge_a)[0] * (edge_b.T * edge_b)[0]
                    - (edge_a.T * edge_b)[0] ** 2
                )
            )
            / 2
        )

    gradient = sp.simplify(GRADIENT.xreplace(substitutions))
    hessian = sp.simplify(HESSIAN.xreplace(substitutions))
    multiplier = sp.simplify(areas[0] / gradient[0])
    assert all(sp.simplify(areas[k] - multiplier * gradient[k]) == 0 for k in range(10))
    block = sp.zeros(11)
    for row in range(10):
        for column in range(10):
            block[row, column] = multiplier * hessian[row, column]
        block[row, 10] = gradient[row]
        block[10, row] = gradient[row]

    return multiplier, sp.factor(block.det()), block, tuple(triangles)


def exact_cross_vector(
    triangles: tuple[tuple[int, int, int], ...], edge: tuple[int, int]
) -> sp.Matrix:
    """Squared-length/angle mixed Hessian vector ``-dA/dq`` plus zero for mu."""
    return sp.Matrix(
        [-exact_area_derivative(triangle, edge) for triangle in triangles] + [0]
    )


def exact_simplex_block(omitted: int):
    multiplier, determinant, block, triangles = exact_simplex_auxiliary_data(omitted)
    cross = exact_cross_vector(triangles, (0, 1))
    contribution = sp.factor(-(cross.T * block.inv() * cross)[0])
    return multiplier, determinant, contribution


def exact_cluster_angle_schur_cross(
    first_edge: tuple[int, int], second_edge: tuple[int, int]
):
    """Sum the auxiliary-angle Schur cross term over the initial 4-2 cluster.

    A fixed-angle length Hessian must be added in general.  For ``(01)`` crossed
    with either ``(23)`` or ``(45)`` it is exactly zero because no triangle area
    depends on both disjoint squared edges.
    """
    result = sp.S.Zero
    for omitted in range(2, 6):
        _, _, block, triangles = exact_simplex_auxiliary_data(omitted)
        first = exact_cross_vector(triangles, first_edge)
        second = exact_cross_vector(triangles, second_edge)
        result -= (first.T * block.inv() * second)[0]
    return sp.factor(result)


def source_target():
    homogeneous = sp.Matrix.vstack(sp.ones(1, 6), sp.Matrix.hstack(*POINTS))
    circuit = [
        sp.factor(
            (-1) ** omitted
            * homogeneous[:, [column for column in range(6) if column != omitted]].det()
            / 24
        )
        for omitted in range(6)
    ]
    if circuit[0] < 0:
        circuit = [-coefficient for coefficient in circuit]
    d_factor = sp.factor(
        2
        * sum(
            circuit[vertex] * sum(POINTS[vertex][k] ** 2 for k in range(4))
            for vertex in range(6)
        )
    )
    volumes = []
    for omitted in range(6):
        vertices = [vertex for vertex in range(6) if vertex != omitted]
        edges = sp.Matrix.hstack(
            *[POINTS[vertices[k]] - POINTS[vertices[0]] for k in range(1, 5)]
        )
        volumes.append(abs(edges.det()) / 24)
    target = sp.factor(-d_factor * volumes[0] * volumes[1] / (48 * sp.prod(volumes[2:])))
    return circuit, d_factor, volumes, target


if __name__ == "__main__":
    contributions = []
    for omitted in range(2, 6):
        multiplier, determinant, contribution = exact_simplex_block(omitted)
        assert determinant != 0
        contributions.append(contribution)
        print(
            f"omitted={omitted} multiplier={multiplier} determinant={determinant} "
            f"schur_contribution={contribution}"
        )
    schur_hessian = sp.factor(sum(contributions))
    circuit, d_factor, volumes, target = source_target()
    ratio = sp.factor(schur_hessian / target)
    assert ratio == sp.Rational(-1, 8)
    input_cross = exact_cluster_angle_schur_cross((0, 1), (2, 3))
    output_cross = exact_cluster_angle_schur_cross((0, 1), (4, 5))
    assert input_cross == output_cross == sp.Rational(-1, 300)
    # Edges (23) and (45) never occur in the same four-simplex, so their direct
    # fixed-internal Hessian vanishes by simplex locality.  Eliminating (01)
    # therefore produces the exact separated boundary transfer below.
    transfer = sp.factor(-input_cross * output_cross / schur_hessian)
    assert transfer == sp.Rational(3, 3584)
    print(f"circuit={circuit}")
    print(f"D={d_factor}")
    print(f"volumes={volumes}")
    print(f"first_order_schur_hessian={schur_hessian}")
    print(f"source_target={target}")
    print(f"ratio={ratio}")
    print(f"input_internal_cross={input_cross}")
    print(f"internal_output_cross={output_cross}")
    print(f"separated_boundary_transfer={transfer}")
