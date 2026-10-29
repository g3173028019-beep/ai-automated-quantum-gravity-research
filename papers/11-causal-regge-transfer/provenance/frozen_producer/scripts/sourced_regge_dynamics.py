"""Nonlinear sourced Einstein-Regge dynamics for the fixed-boundary 4-2 cell.

The sole bulk metric coordinate is q=s_01.  With the full Regge boundary term included,
Schlaefli's identity reduces the vacuum equation to

    E_R(q) = sum_{h superset (01)} (d A_h/dq) epsilon_h(q) = 0.

A signed worldline/source term ``j sqrt(q)`` gives ``E_R + j/(2 sqrt(q)) = 0``.  In units
with c=hbar=1, ``j`` can be read as ``8 pi G eta m`` for a particle segment on edge (01),
where ``eta`` records the Euclidean/Lorentzian action-sign convention.  The Euclidean cell
does not by itself fix that continuation.
"""

from math import sqrt

import numpy as np
from scipy.optimize import brentq, minimize_scalar

try:
    from scripts.regge_four_two_action import (
        INITIAL_SIMPLICES,
        INTERNAL_TRIANGLES,
        cartan_internal_curvature_action,
        canonical_radon_configuration,
        flat_regge_hessian,
        internal_deficit,
        regge_action_with_boundary,
        regge_internal_curvature_action,
        simplex_coordinates,
        squared_lengths,
    )
except ModuleNotFoundError:  # Allow direct execution from scripts/.
    from regge_four_two_action import (
        INITIAL_SIMPLICES,
        INTERNAL_TRIANGLES,
        cartan_internal_curvature_action,
        canonical_radon_configuration,
        flat_regge_hessian,
        internal_deficit,
        regge_action_with_boundary,
        regge_internal_curvature_action,
        simplex_coordinates,
        squared_lengths,
    )


def varied_lengths(lengths: np.ndarray, q_01: float) -> np.ndarray:
    varied = np.asarray(lengths, dtype=float).copy()
    varied[0, 1] = varied[1, 0] = q_01
    return varied


def area_derivative_at_q(
    triangle: tuple[int, int, int], lengths: np.ndarray, q_01: float
) -> float:
    """Return d A_(01k)/d q_01 and zero for the fixed-area triangles."""
    if triangle[:2] != (0, 1):
        return 0.0
    third = triangle[2]
    q_0k = lengths[0, third]
    q_1k = lengths[1, third]
    radicand = (
        2.0 * (q_01 * q_0k + q_01 * q_1k + q_0k * q_1k)
        - q_01**2
        - q_0k**2
        - q_1k**2
    )
    if radicand <= 0:
        raise ValueError("the varied internal triangle is degenerate")
    return float((q_0k + q_1k - q_01) / (4.0 * sqrt(radicand)))


def vacuum_regge_equation(lengths: np.ndarray, q_01: float) -> float:
    """The exact fixed-boundary Regge equation after Schlaefli cancellation."""
    return float(
        sum(
            area_derivative_at_q(triangle, lengths, q_01)
            * internal_deficit(triangle, lengths, q_01)
            for triangle in INTERNAL_TRIANGLES
        )
    )


def sourced_regge_equation(lengths: np.ndarray, q_01: float, source: float) -> float:
    if q_01 <= 0:
        raise ValueError("q_01 must be positive")
    return float(vacuum_regge_equation(lengths, q_01) + source / (2.0 * sqrt(q_01)))


def sourced_action(lengths: np.ndarray, q_01: float, source: float) -> float:
    """The 8 pi G-rescaled action S_R(q)+j sqrt(q)."""
    return float(regge_action_with_boundary(lengths, q_01) + source * sqrt(q_01))


def source_response_curve(lengths: np.ndarray, q_01: float) -> float:
    """Signed source j(q) for which q is a saddle."""
    if q_01 <= 0:
        raise ValueError("q_01 must be positive")
    return float(-2.0 * sqrt(q_01) * vacuum_regge_equation(lengths, q_01))


def is_valid_q(lengths: np.ndarray, q_01: float) -> bool:
    if q_01 <= 0:
        return False
    varied = varied_lengths(lengths, q_01)
    try:
        for simplex in INITIAL_SIMPLICES:
            simplex_coordinates(simplex, varied)
    except (ValueError, np.linalg.LinAlgError):
        return False
    return True


def valid_q_interval(lengths: np.ndarray) -> tuple[float, float]:
    """Numerically bracket the connected nondegenerate interval containing the flat q."""
    q_background = float(lengths[0, 1])
    if not is_valid_q(lengths, q_background):
        raise ValueError("background q is not a valid initial 4-2 geometry")

    lower_invalid = max(np.finfo(float).eps, q_background * 0.5)
    while is_valid_q(lengths, lower_invalid) and lower_invalid > 1e-12 * q_background:
        lower_invalid *= 0.5
    lower_valid = q_background
    for _ in range(80):
        midpoint = 0.5 * (lower_invalid + lower_valid)
        if is_valid_q(lengths, midpoint):
            lower_valid = midpoint
        else:
            lower_invalid = midpoint

    upper_valid = q_background
    upper_invalid = q_background * 1.2
    for _ in range(80):
        if not is_valid_q(lengths, upper_invalid):
            break
        upper_valid = upper_invalid
        upper_invalid *= 1.4
    else:
        raise RuntimeError("failed to bracket the upper degeneracy")
    for _ in range(80):
        midpoint = 0.5 * (upper_valid + upper_invalid)
        if is_valid_q(lengths, midpoint):
            upper_valid = midpoint
        else:
            upper_invalid = midpoint

    width = upper_valid - lower_valid
    margin = max(1e-12 * q_background, 1e-10 * width)
    return lower_valid + margin, upper_valid - margin


def equation_hessian(lengths: np.ndarray, q_01: float, source: float) -> float:
    """Second derivative of the sourced action, evaluated by a centered difference."""
    step = 2e-6 * max(1.0, abs(q_01))
    return float(
        (
            sourced_regge_equation(lengths, q_01 + step, source)
            - sourced_regge_equation(lengths, q_01 - step, source)
        )
        / (2.0 * step)
    )


def find_sourced_saddles(
    lengths: np.ndarray,
    source: float,
    samples: int = 600,
    q_bounds: tuple[float, float] | None = None,
) -> list[dict]:
    """Find every simple saddle in the valid q interval by sign-change bracketing."""
    valid_lower, valid_upper = valid_q_interval(lengths)
    if q_bounds is None:
        lower, upper = valid_lower, valid_upper
    else:
        lower = max(valid_lower, q_bounds[0])
        upper = min(valid_upper, q_bounds[1])
        if lower >= upper:
            raise ValueError("q_bounds do not overlap the valid geometry interval")
    grid = np.linspace(lower, upper, samples)
    values = [sourced_regge_equation(lengths, q_01, source) for q_01 in grid]
    roots = []
    for q_left, q_right, f_left, f_right in zip(grid[:-1], grid[1:], values[:-1], values[1:]):
        if not np.isfinite(f_left) or not np.isfinite(f_right) or f_left * f_right > 0:
            continue
        root = brentq(
            lambda q_01: sourced_regge_equation(lengths, q_01, source),
            q_left,
            q_right,
            xtol=5e-14,
            rtol=5e-14,
        )
        if roots and abs(root - roots[-1]["q_01"]) < 1e-9 * max(1.0, abs(root)):
            continue
        deficits = {
            triangle: internal_deficit(triangle, lengths, root)
            for triangle in INTERNAL_TRIANGLES
        }
        roots.append(
            {
                "q_01": float(root),
                "length_01": float(sqrt(root)),
                "equation_residual": sourced_regge_equation(lengths, root, source),
                "action": sourced_action(lengths, root, source),
                "hessian": equation_hessian(lengths, root, source),
                "deficits": deficits,
                "regge_internal_curvature": regge_internal_curvature_action(lengths, root),
                "cartan_internal_curvature": cartan_internal_curvature_action(lengths, root),
            }
        )
    return roots


def positive_source_fold(lengths: np.ndarray) -> dict:
    """Locate the first positive-j fold on the q>q_flat Euclidean branch."""
    q_background = float(lengths[0, 1])
    _, upper = valid_q_interval(lengths)
    result = minimize_scalar(
        lambda q_01: -source_response_curve(lengths, q_01),
        bounds=(q_background + 1e-10 * max(1.0, q_background), upper),
        method="bounded",
        options={"xatol": 1e-14},
    )
    if not result.success:
        raise RuntimeError("failed to locate the positive source fold")
    q_critical = float(result.x)
    source_critical = source_response_curve(lengths, q_critical)
    return {
        "q_critical": q_critical,
        "source_critical": source_critical,
        "equation_hessian": equation_hessian(lengths, q_critical, source_critical),
    }


def linear_response(lengths: np.ndarray, source: float) -> float:
    """Return delta q from the flat physical Hessian."""
    q_background = float(lengths[0, 1])
    return float(-source / (2.0 * sqrt(q_background) * flat_regge_hessian(lengths)))


if __name__ == "__main__":
    edge_data = squared_lengths(canonical_radon_configuration())
    q_flat = edge_data[0, 1]
    probe_source = 1e-5
    print(f"formal_name=Phase-First Cartan-Regge Gravity (PCRG)")
    print(f"q_flat={q_flat:.16g}")
    print(f"vacuum_equation={vacuum_regge_equation(edge_data, q_flat):.16g}")
    print(f"flat_hessian={flat_regge_hessian(edge_data):.16g}")
    print(f"linear_delta_q(source={probe_source})={linear_response(edge_data, probe_source):.16g}")
    _, q_upper = valid_q_interval(edge_data)
    for index, saddle in enumerate(
        find_sourced_saddles(edge_data, probe_source, q_bounds=(q_flat, q_upper))
    ):
        print(
            f"saddle_{index}: q={saddle['q_01']:.16g} "
            f"hessian={saddle['hessian']:.16g} "
            f"residual={saddle['equation_residual']:.3e}"
        )
    fold = positive_source_fold(edge_data)
    print(f"fold_q={fold['q_critical']:.16g}")
    print(f"fold_source={fold['source_critical']:.16g}")
    print(f"fold_hessian={fold['equation_hessian']:.3e}")
