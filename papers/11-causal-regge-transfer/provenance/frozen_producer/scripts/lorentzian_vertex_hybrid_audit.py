"""Vertex-local hybrid-curvature audit of the C5*C6 Lorentzian star.

A27 contracted all hinge tensors into one Einstein screen using three diagnostic
weights.  This module implements the distinct vertex-hybrid prescription forced by
the circumcentric Regge volume decomposition.  For a hinge h incident on the sole
interior vertex v,

    R_h Delta V_hv  proportional to  epsilon_h A_hv,

because the dual area in the sectional curvature cancels the dual area in the
hybrid volume.  A_hv is computed by signed circumcentric subdivision inside the
hinge.  The normal-plane Riemann generator is also multiplied by its plane-signature
so that its scalar contraction has the same sign for Euclidean and Lorentzian
normal planes.

This is a finite vertex average, not a continuum pointwise curvature theorem and
not a causal propagation certificate.
"""

import numpy as np

try:
    from scripts.complex_lorentzian_regge import (
        MINKOWSKI_METRIC,
        squared_intervals,
        triangle_area,
        varied_lengths,
    )
    from scripts.lorentzian_weighted_weyl_star import (
        INTERNAL_EDGES,
        POINTS,
        transverse_angle_deformed_points,
        weighted_weyl_star_certificate,
    )
    from scripts.lorentzian_weyl_star import _hinge_normal_riemann, _symmetric_components
except ModuleNotFoundError:
    from complex_lorentzian_regge import (
        MINKOWSKI_METRIC,
        squared_intervals,
        triangle_area,
        varied_lengths,
    )
    from lorentzian_weighted_weyl_star import (
        INTERNAL_EDGES,
        POINTS,
        transverse_angle_deformed_points,
        weighted_weyl_star_certificate,
    )
    from lorentzian_weyl_star import _hinge_normal_riemann, _symmetric_components


def hinge_vertex_area_share(
    triangle: tuple[int, int, int], points: np.ndarray, vertex: int = 0
) -> complex:
    """Signed circumcentric portion of a triangle assigned to one vertex.

    In affine coordinates based at ``vertex``, the circumcenter coefficients c
    obey G c = diag(G)/2.  The two circumcentric flag triangles incident on the
    vertex occupy (c_1+c_2)/2 of the oriented hinge area.  The principal complex
    square root distinguishes spacelike from timelike hinge area.
    """
    others = [index for index in triangle if index != vertex]
    if len(others) != 2:
        raise ValueError("the requested vertex must occur exactly once in the hinge")
    tangents = points[others] - points[vertex]
    gram = tangents @ MINKOWSKI_METRIC @ tangents.T
    determinant = float(np.linalg.det(gram))
    if abs(determinant) < 1e-13:
        raise ValueError("null hinge has no nondegenerate hybrid area")
    circumcenter_coefficients = np.linalg.solve(gram, 0.5 * np.diag(gram))
    total_area = 0.5 * np.sqrt(complex(determinant))
    return complex(0.5 * np.sum(circumcenter_coefficients) * total_area)


def normal_plane_signature(
    triangle: tuple[int, int, int], points: np.ndarray
) -> int:
    origin = points[triangle[0]]
    tangents = points[list(triangle[1:])] - origin
    normal_plane = np.linalg.svd(tangents @ MINKOWSKI_METRIC)[2][-2:].T
    determinant = float(
        np.linalg.det(normal_plane.T @ MINKOWSKI_METRIC @ normal_plane)
    )
    if abs(determinant) < 1e-13:
        raise ValueError("null normal plane is outside the audit")
    return 1 if determinant > 0.0 else -1


def signature_corrected_hinge_riemann(
    triangle: tuple[int, int, int], points: np.ndarray
) -> np.ndarray:
    return normal_plane_signature(triangle, points) * _hinge_normal_riemann(
        triangle, points
    )


def _einstein_tensor(riemann: np.ndarray) -> np.ndarray:
    ricci = np.einsum("bd,abcd->ac", MINKOWSKI_METRIC, riemann)
    scalar = np.einsum("ac,ac", MINKOWSKI_METRIC, ricci)
    return ricci - 0.5 * MINKOWSKI_METRIC * scalar


def regge_edge_moment_matrix(
    triangles: tuple[tuple[int, int, int], ...],
    points: np.ndarray,
    step_scale: float = 1e-6,
) -> np.ndarray:
    """Linear map from hinge deficits to internal squared-length equations.

    At a flat background Schlaefli cancellation reduces the linearized internal
    equation to sum_h (partial A_h/partial s_e) delta epsilon_h = 0.
    """
    background = squared_intervals(points)
    result = np.zeros((len(INTERNAL_EDGES), len(triangles)), dtype=complex)
    for row, edge in enumerate(INTERNAL_EDGES):
        step = step_scale * max(1.0, abs(background[edge]))
        plus = varied_lengths(background, (edge,), (background[edge] + step,))
        minus = varied_lengths(background, (edge,), (background[edge] - step,))
        for column, triangle in enumerate(triangles):
            if set(edge) <= set(triangle):
                result[row, column] = (
                    triangle_area(triangle, plus) - triangle_area(triangle, minus)
                ) / (2.0 * step)
    return result


def vertex_hybrid_audit(points: np.ndarray | None = None) -> dict:
    if points is None:
        points = POINTS
    points = np.asarray(points, dtype=float)
    base = weighted_weyl_star_certificate(points=points)
    triangles = base["internal_triangles"]
    area_shares = np.asarray(
        [hinge_vertex_area_share(triangle, points) for triangle in triangles]
    )
    signatures = np.asarray(
        [normal_plane_signature(triangle, points) for triangle in triangles]
    )
    riemann_tensors = [
        signature_corrected_hinge_riemann(triangle, points)
        for triangle in triangles
    ]
    riemann_map = np.column_stack(
        [riemann.reshape(-1) for riemann in riemann_tensors]
    )
    einstein_map = np.column_stack(
        [_symmetric_components(_einstein_tensor(riemann)) for riemann in riemann_tensors]
    )
    screen = (
        einstein_map
        @ np.diag(area_shares)
        @ base["on_shell_curvature"]
    )
    _, singular_values, right = np.linalg.svd(screen, full_matrices=True)
    threshold = singular_values[0] * 1e-5
    rank = int(np.count_nonzero(singular_values > threshold))
    null_axes = right.T[:, rank:]
    null_riemann = riemann_map @ np.diag(area_shares) @ base["on_shell_curvature"] @ null_axes
    edge_moment_map = regge_edge_moment_matrix(triangles, points)
    edge_equation_residual = edge_moment_map @ base["on_shell_curvature"]
    return {
        "base": base,
        "area_shares": area_shares,
        "normal_plane_signatures": signatures,
        "screen": screen,
        "singular_values": singular_values,
        "rank": rank,
        "nullity": null_axes.shape[1],
        "null_axes": null_axes,
        "null_riemann_singular_values": np.linalg.svd(
            null_riemann, compute_uv=False
        ),
        "edge_moment_singular_values": np.linalg.svd(
            edge_moment_map, compute_uv=False
        ),
        "edge_equation_residual_singular_values": np.linalg.svd(
            edge_equation_residual, compute_uv=False
        ),
        "maximum_edge_equation_residual": float(
            np.max(np.abs(edge_equation_residual))
        ),
        "spacelike_hinge_count": int(np.count_nonzero(signatures < 0)),
        "timelike_hinge_count": int(np.count_nonzero(signatures > 0)),
    }


if __name__ == "__main__":
    for label, points in (
        ("base", POINTS),
        ("deformed", transverse_angle_deformed_points()),
    ):
        result = vertex_hybrid_audit(points)
        print(f"{label}_area_shares={result['area_shares']}")
        print(f"{label}_normal_plane_signatures={result['normal_plane_signatures']}")
        print(f"{label}_singular_values={result['singular_values']}")
        print(f"{label}_rank={result['rank']}")
        print(f"{label}_nullity={result['nullity']}")
        print(
            f"{label}_null_riemann_singular_values="
            f"{result['null_riemann_singular_values']}"
        )
        print(
            f"{label}_edge_moment_singular_values="
            f"{result['edge_moment_singular_values']}"
        )
        print(
            f"{label}_edge_equation_residual_singular_values="
            f"{result['edge_equation_residual_singular_values']}"
        )
        print(
            f"{label}_maximum_edge_equation_residual="
            f"{result['maximum_edge_equation_residual']}"
        )
