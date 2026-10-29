"""Euclidean 4-2 Regge action and Hessian calibration for Road A.

The module reconstructs each four-simplex independently from its ten squared edge lengths,
adds the Regge boundary term, and differentiates the internal-edge equation at a flat Radon
configuration.  It supplies the geometric target for the discontinuous Cartan action; it is
not itself an independent-connection formulation.
"""

from itertools import combinations
from math import acos, factorial, pi, sqrt

import numpy as np
from scipy.linalg import logm


INITIAL_SIMPLICES = [tuple(v for v in range(6) if v != k) for k in range(2, 6)]
ALL_TRIANGLES = sorted(
    {triangle for simplex in INITIAL_SIMPLICES for triangle in combinations(simplex, 3)}
)
INTERNAL_TRIANGLES = [triangle for triangle in ALL_TRIANGLES if triangle[:2] == (0, 1)]


def squared_lengths(points: np.ndarray) -> np.ndarray:
    points = np.asarray(points, dtype=float)
    if points.shape != (6, 4):
        raise ValueError("points must have shape (6,4)")
    differences = points[:, None, :] - points[None, :, :]
    return np.einsum("...i,...i->...", differences, differences)


def simplex_coordinates(simplex: tuple[int, ...], lengths: np.ndarray) -> dict[int, np.ndarray]:
    """Reconstruct Euclidean simplex coordinates from squared edge lengths."""
    base = simplex[0]
    others = simplex[1:]
    gram = np.array(
        [
            [
                0.5 * (lengths[base, first] + lengths[base, second] - lengths[first, second])
                for second in others
            ]
            for first in others
        ]
    )
    cholesky = np.linalg.cholesky(gram)
    coordinates = np.zeros((5, 4))
    coordinates[1:] = cholesky
    return {vertex: coordinates[index] for index, vertex in enumerate(simplex)}


def outward_face_normal(
    face: tuple[int, ...] | list[int], opposite: int, coordinates: dict[int, np.ndarray]
) -> np.ndarray:
    ordered = sorted(face)
    origin = coordinates[ordered[0]]
    tangents = np.stack([coordinates[v] - origin for v in ordered[1:]])
    normal = np.linalg.svd(tangents)[2][-1]
    normal /= np.linalg.norm(normal)
    if normal @ (coordinates[opposite] - origin) > 0:
        normal = -normal
    return normal


def interior_dihedral_angle(
    simplex: tuple[int, ...], triangle: tuple[int, int, int], lengths: np.ndarray
) -> float:
    coordinates = simplex_coordinates(simplex, lengths)
    opposite = [vertex for vertex in simplex if vertex not in triangle]
    face_first = [vertex for vertex in simplex if vertex != opposite[0]]
    face_second = [vertex for vertex in simplex if vertex != opposite[1]]
    normal_first = outward_face_normal(face_first, opposite[0], coordinates)
    normal_second = outward_face_normal(face_second, opposite[1], coordinates)
    return acos(float(np.clip(-normal_first @ normal_second, -1.0, 1.0)))


def triangle_area(triangle: tuple[int, int, int], lengths: np.ndarray) -> float:
    first, second, third = triangle
    a = lengths[first, second]
    b = lengths[first, third]
    c = lengths[second, third]
    radicand = 2.0 * (a * b + a * c + b * c) - a * a - b * b - c * c
    if radicand <= 0:
        raise ValueError("triangle is degenerate")
    return 0.25 * sqrt(radicand)


def internal_deficit(
    triangle: tuple[int, int, int], lengths: np.ndarray, q_01: float
) -> float:
    varied = np.array(lengths, copy=True)
    varied[0, 1] = varied[1, 0] = q_01
    incident = [simplex for simplex in INITIAL_SIMPLICES if set(triangle) <= set(simplex)]
    return 2.0 * pi - sum(
        interior_dihedral_angle(simplex, triangle, varied) for simplex in incident
    )


def oriented_simplex_coordinates(
    omitted: int, lengths: np.ndarray, q_01: float
) -> tuple[tuple[int, ...], dict[int, np.ndarray]]:
    """Coordinates whose orientation is the 4-2 chain coefficient (-1)^omitted."""
    varied = np.array(lengths, copy=True)
    varied[0, 1] = varied[1, 0] = q_01
    simplex = tuple(vertex for vertex in range(6) if vertex != omitted)
    coordinates = simplex_coordinates(simplex, varied)
    orientation = np.diag([1.0, 1.0, 1.0, (-1.0) ** omitted])
    return simplex, {vertex: orientation @ value for vertex, value in coordinates.items()}


def face_transport(
    source_omitted: int, target_omitted: int, lengths: np.ndarray, q_01: float
) -> np.ndarray:
    """Levi-Civita SO(4) transport across the shared tetrahedron.

    Tangent vectors are identified pointwise and the source outward normal maps to minus the
    target outward normal, as required by oriented gluing.
    """
    source_simplex, source = oriented_simplex_coordinates(source_omitted, lengths, q_01)
    target_simplex, target = oriented_simplex_coordinates(target_omitted, lengths, q_01)
    shared = sorted(set(source_simplex) & set(target_simplex))
    base = shared[0]
    source_tangents = np.column_stack([source[v] - source[base] for v in shared[1:]])
    target_tangents = np.column_stack([target[v] - target[base] for v in shared[1:]])
    source_opposite = next(iter(set(source_simplex) - set(shared)))
    target_opposite = next(iter(set(target_simplex) - set(shared)))
    source_normal = outward_face_normal(shared, source_opposite, source)
    target_normal = outward_face_normal(shared, target_opposite, target)
    source_frame = np.column_stack([source_tangents, source_normal])
    target_frame = np.column_stack([target_tangents, -target_normal])
    transport = target_frame @ np.linalg.inv(source_frame)
    if not np.allclose(transport.T @ transport, np.eye(4), atol=2e-10):
        raise ValueError("face matching did not produce an orthogonal transport")
    if np.linalg.det(transport) < 1.0 - 2e-10:
        raise ValueError("face transport is not in SO(4)")
    return transport


def hodge_bivector(bivector: np.ndarray) -> np.ndarray:
    dual = np.zeros((4, 4))
    for a in range(4):
        for b in range(4):
            for c in range(4):
                for d in range(4):
                    indices = (a, b, c, d)
                    if len(set(indices)) < 4:
                        continue
                    inversions = sum(
                        indices[i] > indices[j]
                        for i in range(4)
                        for j in range(i + 1, 4)
                    )
                    dual[a, b] += 0.5 * (-1.0 if inversions % 2 else 1.0) * bivector[c, d]
    return dual


def cartan_internal_curvature_action(lengths: np.ndarray, q_01: float) -> float:
    """Sum 1/2 <star B, log Hol> on the four internal dual faces.

    The connection here is the geometric Levi-Civita transport derived from matched simplex
    coframes.  This tests the Cartan/Regge geometric-sector identity, not independent
    connection integration.
    """
    value = 0.0
    dual_nodes = set(range(2, 6))
    for dual_face in combinations(range(2, 6), 3):
        first, second, third = dual_face
        holonomy = (
            face_transport(third, first, lengths, q_01)
            @ face_transport(second, third, lengths, q_01)
            @ face_transport(first, second, lengths, q_01)
        )
        generator = np.real_if_close(logm(holonomy), tol=1000).real
        remaining = next(iter(dual_nodes - set(dual_face)))
        triangle = (0, 1, remaining)
        _, coordinates = oriented_simplex_coordinates(first, lengths, q_01)
        first_edge = coordinates[1] - coordinates[0]
        second_edge = coordinates[remaining] - coordinates[0]
        bivector = np.outer(first_edge, second_edge) - np.outer(second_edge, first_edge)
        raw_pairing = 0.25 * np.sum(hodge_bivector(bivector) * generator)
        value += (-1.0) ** remaining * raw_pairing
    return float(value)


def regge_internal_curvature_action(lengths: np.ndarray, q_01: float) -> float:
    varied = np.array(lengths, copy=True)
    varied[0, 1] = varied[1, 0] = q_01
    return float(
        sum(
            triangle_area(triangle, varied) * internal_deficit(triangle, lengths, q_01)
            for triangle in INTERNAL_TRIANGLES
        )
    )


def regge_action_with_boundary(lengths: np.ndarray, q_01: float) -> float:
    """Euclidean Regge action with the common-boundary exterior-angle term."""
    varied = np.array(lengths, copy=True)
    varied[0, 1] = varied[1, 0] = q_01
    value = 0.0
    internal = set(INTERNAL_TRIANGLES)
    for triangle in ALL_TRIANGLES:
        incident = [simplex for simplex in INITIAL_SIMPLICES if set(triangle) <= set(simplex)]
        angle_sum = sum(
            interior_dihedral_angle(simplex, triangle, varied) for simplex in incident
        )
        curvature_angle = 2.0 * pi - angle_sum if triangle in internal else pi - angle_sum
        value += triangle_area(triangle, varied) * curvature_angle
    return float(value)


def area_derivative_wrt_q(triangle: tuple[int, int, int], lengths: np.ndarray) -> float:
    """Derivative d A_(01k)/d s_01 in squared-length variables."""
    if triangle[:2] != (0, 1):
        return 0.0
    third = triangle[2]
    q = lengths[0, 1]
    a = lengths[0, third]
    b = lengths[1, third]
    radicand = 2.0 * (q * a + q * b + a * b) - q * q - a * a - b * b
    return (a + b - q) / (4.0 * sqrt(radicand))


def flat_regge_hessian(lengths: np.ndarray, relative_step: float = 1e-6) -> float:
    """Use the Schlaefli-reduced equation to compute d²S/d(s_01)² at flatness."""
    q = float(lengths[0, 1])
    step = relative_step * max(1.0, abs(q))
    hessian = 0.0
    for triangle in INTERNAL_TRIANGLES:
        deficit_derivative = (
            internal_deficit(triangle, lengths, q + step)
            - internal_deficit(triangle, lengths, q - step)
        ) / (2.0 * step)
        hessian += area_derivative_wrt_q(triangle, lengths) * deficit_derivative
    return float(hessian)


def oriented_affine_circuit(points: np.ndarray) -> np.ndarray:
    """Fix the circuit sign so the 4-2 Radon side has c_0,c_1>0."""
    homogeneous = np.vstack([np.ones(6), np.asarray(points, dtype=float).T])
    coefficients = np.array(
        [
            (-1) ** k
            * np.linalg.det(np.delete(homogeneous, k, axis=1))
            / factorial(4)
            for k in range(6)
        ]
    )
    if coefficients[0] < 0:
        coefficients = -coefficients
    if not (np.all(coefficients[:2] > 0) and np.all(coefficients[2:] < 0)):
        raise ValueError("points do not realize the {0,1}|{2,3,4,5} Radon partition")
    return coefficients


def source_mapped_regge_target(points: np.ndarray) -> float:
    """A8 affine target before source-D normalization and Euclidean sign matching.

    A20 proves that this function's coordinate ``D`` is ``8*D_source``.  Hence the
    project Euclidean Hessian relation ``H_Regge/H_A8=-1/8`` combines that factor
    with the opposite Euclidean action sign; it is not an unexplained dynamical ratio.
    """
    points = np.asarray(points, dtype=float)
    coefficients = oriented_affine_circuit(points)
    d_factor = 2.0 * coefficients @ np.sum(points**2, axis=1)
    volumes = []
    for omitted in range(6):
        simplex = np.delete(points, omitted, axis=0)
        volumes.append(abs(np.linalg.det((simplex[1:] - simplex[0]).T)) / factorial(4))
    volumes = np.asarray(volumes)
    return float(-d_factor * volumes[0] * volumes[1] / (48.0 * np.prod(volumes[2:])))


def canonical_radon_configuration() -> np.ndarray:
    return np.array(
        [
            [0.7, 0.25, 0.25, 0.25],
            [-0.8, 0.25, 0.25, 0.25],
            [0.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ]
    )


def random_radon_configuration(rng: np.random.Generator) -> np.ndarray:
    """Generate a generic convex 4-2 configuration with edge (01) crossing tetrahedron 2345."""
    for _ in range(1000):
        tetrahedron = rng.normal(size=(4, 4))
        tangents = tetrahedron[1:] - tetrahedron[0]
        if np.linalg.matrix_rank(tangents, tol=1e-9) != 3:
            continue
        normal = np.linalg.svd(tangents)[2][-1]
        normal /= np.linalg.norm(normal)
        weights = rng.dirichlet(np.ones(4) * 1.5)
        intersection = weights @ tetrahedron
        tangent_shift = rng.normal(scale=0.25, size=4)
        tangent_shift -= (tangent_shift @ normal) * normal
        crossing_direction = normal + tangent_shift
        positive, negative = rng.uniform(0.45, 1.4, size=2)
        points = np.vstack(
            [
                intersection + positive * crossing_direction,
                intersection - negative * crossing_direction,
                tetrahedron,
            ]
        )
        volumes = []
        for omitted in range(6):
            simplex = np.delete(points, omitted, axis=0)
            volumes.append(abs(np.linalg.det((simplex[1:] - simplex[0]).T)))
        if min(volumes) < 2e-3:
            continue
        try:
            oriented_affine_circuit(points)
            lengths = squared_lengths(points)
            q = lengths[0, 1]
            for sign in (-1.0, 1.0):
                varied = np.array(lengths, copy=True)
                varied[0, 1] = varied[1, 0] = q * (1.0 + sign * 2.5e-3)
                for simplex in INITIAL_SIMPLICES:
                    simplex_coordinates(simplex, varied)
        except (ValueError, np.linalg.LinAlgError):
            continue
        return points
    raise RuntimeError("failed to generate a stable Radon 4-2 configuration")


if __name__ == "__main__":
    configuration = canonical_radon_configuration()
    edge_lengths = squared_lengths(configuration)
    q_background = edge_lengths[0, 1]
    print(f"flat_first_derivative_check_step=1e-5")
    step = 1e-5 * q_background
    derivative = (
        regge_action_with_boundary(edge_lengths, q_background + step)
        - regge_action_with_boundary(edge_lengths, q_background - step)
    ) / (2.0 * step)
    hessian = flat_regge_hessian(edge_lengths)
    target = source_mapped_regge_target(configuration)
    print(f"flat_first_derivative={derivative:.16g}")
    print(f"regge_hessian={hessian:.16g}")
    print(f"source_target={target:.16g}")
    print(f"ratio={hessian / target:.16g}")
