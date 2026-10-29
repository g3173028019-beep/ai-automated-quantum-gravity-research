"""Gauge-reduced first-order Regge angle block on a Euclidean four-simplex.

For ten four-dimensional dihedral angles ``theta_ij`` define the angle Gram matrix

    G_ii = 1,             G_ij = -cos(theta_ij).

The Barrett first-order simplex term uses a multiplier ``mu det(G)``.  At a geometric
simplex its angle equation is ``A_ij = mu * d det(G)/d theta_ij``, where ``A_ij`` is the
area of the triangle complementary to vertices i and j.  The bordered Hessian

    C = [[mu Hess(det G), grad(det G)], [grad(det G)^T, 0]]

is the local auxiliary Gaussian block.  This module evaluates that block and its exact
scale law.  It is a gauge-reduced first-order Regge calculation, not a derivation of the
Haar/BFV measure of an independently integrated Cartan connection.
"""

from itertools import combinations

import numpy as np
import sympy as sp

try:
    from scripts.regge_four_two_action import (
        ALL_TRIANGLES,
        INITIAL_SIMPLICES,
        INTERNAL_TRIANGLES,
        area_derivative_wrt_q,
        interior_dihedral_angle,
        triangle_area,
    )
except ModuleNotFoundError:  # Allow ``python scripts/first_order_angle_measure.py``.
    from regge_four_two_action import (
        ALL_TRIANGLES,
        INITIAL_SIMPLICES,
        INTERNAL_TRIANGLES,
        area_derivative_wrt_q,
        interior_dihedral_angle,
        triangle_area,
    )


FINAL_SIMPLICES = [
    tuple(vertex for vertex in range(6) if vertex != omitted) for omitted in (0, 1)
]
LOCAL_PAIRS = tuple(combinations(range(5), 2))


def _build_angle_gram_derivatives():
    angles = sp.symbols("theta_0:10", real=True)
    gram = sp.eye(5)
    for angle, (first, second) in zip(angles, LOCAL_PAIRS):
        gram[first, second] = gram[second, first] = -sp.cos(angle)
    determinant = sp.det(gram)
    gradient = [sp.diff(determinant, angle) for angle in angles]
    hessian = sp.hessian(determinant, angles)
    return (
        sp.lambdify(angles, determinant, "numpy"),
        sp.lambdify(angles, gradient, "numpy"),
        sp.lambdify(angles, hessian, "numpy"),
    )


_DET_GRAM, _GRAD_DET_GRAM, _HESS_DET_GRAM = _build_angle_gram_derivatives()


def simplex_dihedral_angles_and_areas(
    simplex: tuple[int, ...], lengths: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Return local-pair ordered dihedral angles and complementary triangle areas."""
    if len(simplex) != 5 or len(set(simplex)) != 5:
        raise ValueError("simplex must contain five distinct vertices")
    angles = []
    areas = []
    for first, second in LOCAL_PAIRS:
        triangle = tuple(
            vertex for position, vertex in enumerate(simplex) if position not in (first, second)
        )
        angles.append(interior_dihedral_angle(simplex, triangle, lengths))
        areas.append(triangle_area(triangle, lengths))
    return np.asarray(angles), np.asarray(areas)


def angle_gram_matrix(angles: np.ndarray) -> np.ndarray:
    """Construct G_ii=1 and G_ij=-cos(theta_ij)."""
    angles = np.asarray(angles, dtype=float)
    if angles.shape != (10,):
        raise ValueError("angles must have shape (10,)")
    gram = np.eye(5)
    for angle, (first, second) in zip(angles, LOCAL_PAIRS):
        gram[first, second] = gram[second, first] = -np.cos(angle)
    return gram


def barrett_auxiliary_block(simplex: tuple[int, ...], lengths: np.ndarray) -> dict:
    """Evaluate the geometric saddle and bordered angle/multiplier Hessian.

    ``stationarity_residual`` is relative.  The determinant is best handled through
    ``slogdet`` because the block can be ill-conditioned near area-length singularities.
    """
    angles, areas = simplex_dihedral_angles_and_areas(simplex, lengths)
    determinant = float(_DET_GRAM(*angles))
    gradient = np.asarray(_GRAD_DET_GRAM(*angles), dtype=float)
    hessian = np.asarray(_HESS_DET_GRAM(*angles), dtype=float)
    multiplier = float(gradient @ areas / (gradient @ gradient))
    residual = float(np.linalg.norm(areas - multiplier * gradient) / np.linalg.norm(areas))
    block = np.block(
        [
            [multiplier * hessian, gradient[:, None]],
            [gradient[None, :], np.zeros((1, 1))],
        ]
    )
    sign, logabsdet = np.linalg.slogdet(block)
    singular_values = np.linalg.svd(block, compute_uv=False)
    return {
        "simplex": simplex,
        "angles": angles,
        "areas": areas,
        "gram": angle_gram_matrix(angles),
        "gram_determinant": determinant,
        "gradient": gradient,
        "hessian": hessian,
        "multiplier": multiplier,
        "stationarity_residual": residual,
        "block": block,
        "determinant_sign": float(sign),
        "logabsdet": float(logabsdet),
        "singular_values": singular_values,
        "relative_min_singular_value": float(singular_values[-1] / singular_values[0]),
    }


def cluster_auxiliary_determinant(lengths: np.ndarray, side: str) -> dict:
    """Return the product determinant data for the block-local 4-2 clusters."""
    if side == "initial":
        simplices = INITIAL_SIMPLICES
    elif side == "final":
        simplices = FINAL_SIMPLICES
    else:
        raise ValueError("side must be 'initial' or 'final'")
    blocks = [barrett_auxiliary_block(simplex, lengths) for simplex in simplices]
    sign = float(np.prod([entry["determinant_sign"] for entry in blocks]))
    logabsdet = float(sum(entry["logabsdet"] for entry in blocks))
    return {
        "side": side,
        "simplices": simplices,
        "blocks": blocks,
        "determinant_sign": sign,
        "logabsdet": logabsdet,
    }


def simplex_radial_cross_block(simplex: tuple[int, ...], lengths: np.ndarray) -> np.ndarray:
    """Mixed q=s_01 versus (theta,mu) Hessian row for one simplex.

    At fixed angles the first-order term is ``-sum A_ij theta_ij``.  Hence the first ten
    entries are ``-d A_ij/dq`` and the multiplier entry vanishes.
    """
    cross = np.zeros(11)
    for position, (first, second) in enumerate(LOCAL_PAIRS):
        triangle = tuple(
            sorted(
                vertex
                for local_position, vertex in enumerate(simplex)
                if local_position not in (first, second)
            )
        )
        cross[position] = -area_derivative_wrt_q(triangle, lengths)
    return cross


def first_order_radial_schur_hessian(lengths: np.ndarray) -> dict:
    """Eliminate all local angle/multiplier blocks in the flat initial 4-2 cluster.

    The fixed-angle q-q entry vanishes at the flat saddle because all four q-dependent
    triangle coefficients are internal deficits.  The reduced physical Hessian is therefore
    ``-sum B_sigma C_sigma^{-1} B_sigma^T``.  This is a direct finite first-order realization
    of the Regge radial Hessian, not merely the determinant identity for an abstract block.
    """
    contributions = []
    for simplex in INITIAL_SIMPLICES:
        data = barrett_auxiliary_block(simplex, lengths)
        cross = simplex_radial_cross_block(simplex, lengths)
        contribution = -float(cross @ np.linalg.solve(data["block"], cross))
        contributions.append(
            {
                "simplex": simplex,
                "cross": cross,
                "contribution": contribution,
                "condition_number": float(np.linalg.cond(data["block"])),
            }
        )
    return {
        "fixed_angle_radial_hessian": 0.0,
        "contributions": contributions,
        "schur_hessian": float(sum(entry["contribution"] for entry in contributions)),
    }


def initial_first_order_action_at_geometric_saddle(lengths: np.ndarray) -> float:
    """Evaluate the full boundary first-order action without calling the Regge action."""
    internal = set(INTERNAL_TRIANGLES)
    value = sum(
        triangle_area(triangle, lengths) * (2.0 * np.pi if triangle in internal else np.pi)
        for triangle in ALL_TRIANGLES
    )
    for simplex in INITIAL_SIMPLICES:
        data = barrett_auxiliary_block(simplex, lengths)
        value -= float(data["areas"] @ data["angles"])
        value += data["multiplier"] * data["gram_determinant"]
    return float(value)


def scaled_squared_lengths(lengths: np.ndarray, linear_scale: float) -> np.ndarray:
    if linear_scale <= 0:
        raise ValueError("linear_scale must be positive")
    return np.asarray(lengths, dtype=float) * linear_scale**2


if __name__ == "__main__":
    try:
        from scripts.regge_four_two_action import (
            canonical_radon_configuration,
            flat_regge_hessian,
            regge_action_with_boundary,
            source_mapped_regge_target,
            squared_lengths,
        )
    except ModuleNotFoundError:
        from regge_four_two_action import (
            canonical_radon_configuration,
            flat_regge_hessian,
            regge_action_with_boundary,
            source_mapped_regge_target,
            squared_lengths,
        )

    edge_data = squared_lengths(canonical_radon_configuration())
    first_order = first_order_radial_schur_hessian(edge_data)
    regge = flat_regge_hessian(edge_data)
    target = source_mapped_regge_target(canonical_radon_configuration())
    print(f"first_order_action={initial_first_order_action_at_geometric_saddle(edge_data):.16g}")
    print(f"regge_action={regge_action_with_boundary(edge_data, edge_data[0, 1]):.16g}")
    print(f"first_order_schur_hessian={first_order['schur_hessian']:.16g}")
    print(f"regge_hessian={regge:.16g}")
    print(f"source_target_ratio={first_order['schur_hessian'] / target:.16g}")
