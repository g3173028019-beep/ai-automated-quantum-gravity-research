"""Algebraic nine-hinge target for a vacuum spin-two Regge wave sector.

A24 shows that two independent hinge deficits are insufficient: their Einstein
tensors have no nonzero kernel.  This module reverses the problem.  It constructs
nine nondegenerate Lorentzian normal bivectors whose signed curvature combination
contains an explicit two-polarization, Ricci-flat, null wave sector.

This is a local algebraic design target.  It does not assert that the nine hinges
have already been realized by one causal simplicial four-ball or that their
coefficients solve internal Regge edge equations.
"""

import numpy as np


ETA = np.diag((-1.0, 1.0, 1.0, 1.0))
E0, EX, EY, EZ = np.eye(4)
EP = (EX + EY) / np.sqrt(2.0)
RAPIDITIES = (0.5, 1.0, 0.0)
TRANSVERSE_DIRECTIONS = (EX, EY, EP)
HINGE_LABELS = tuple(
    f"eta={rapidity}:{direction}"
    for rapidity in RAPIDITIES
    for direction in ("x", "y", "p")
)


def normal_bivector(rapidity: float, transverse: np.ndarray) -> np.ndarray:
    timelike = np.cosh(rapidity) * E0 + np.sinh(rapidity) * EZ
    if not np.isclose(timelike @ ETA @ timelike, -1.0, atol=2e-14):
        raise ArithmeticError("boosted normal is not unit timelike")
    if not np.isclose(transverse @ ETA @ transverse, 1.0, atol=2e-14):
        raise ArithmeticError("transverse normal is not unit spacelike")
    if not np.isclose(timelike @ ETA @ transverse, 0.0, atol=2e-14):
        raise ArithmeticError("normal-plane frame is not orthogonal")
    return np.outer(timelike, transverse) - np.outer(transverse, timelike)


def hinge_tensors(rapidity: float, transverse: np.ndarray) -> dict:
    bivector = normal_bivector(rapidity, transverse)
    covariant_bivector = ETA @ bivector @ ETA
    riemann = np.einsum("ab,cd->abcd", covariant_bivector, covariant_bivector)
    ricci = np.einsum("bd,abcd->ac", ETA, riemann)
    scalar = float(np.einsum("ac,ac", ETA, ricci))
    einstein = ricci - 0.5 * ETA * scalar
    return {
        "bivector": bivector,
        "riemann": riemann,
        "ricci": ricci,
        "scalar": scalar,
        "einstein": einstein,
    }


def symmetric_components(tensor: np.ndarray) -> np.ndarray:
    return np.asarray(
        [tensor[first, first] for first in range(4)]
        + [
            tensor[first, second]
            for first in range(4)
            for second in range(first + 1, 4)
        ]
    )


def electric_magnetic_weyl(riemann: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Electric/magnetic parts for observer e_0; input is Ricci-flat in use."""
    electric = np.asarray(
        [[riemann[i, 0, j, 0] for j in range(1, 4)] for i in range(1, 4)]
    )
    epsilon = np.zeros((3, 3, 3), dtype=float)
    epsilon[0, 1, 2] = epsilon[1, 2, 0] = epsilon[2, 0, 1] = 1.0
    epsilon[0, 2, 1] = epsilon[2, 1, 0] = epsilon[1, 0, 2] = -1.0
    magnetic = np.einsum(
        "ikl,klj->ij", epsilon, riemann[1:, 1:, 1:, 0]
    ) / 2.0
    return electric, magnetic


def _rapidity_difference(rapidity: float) -> float:
    # sinh(eta)(cosh(eta)-sinh(eta)) is the mismatch B-E of one
    # finite-boost vacuum polarization block.
    return float(
        np.sinh(rapidity) * (np.cosh(rapidity) - np.sinh(rapidity))
    )


def wave_coefficients() -> tuple[np.ndarray, np.ndarray]:
    """Return exact-by-construction plus/cross coefficients in the nine hinges."""
    first_weight = _rapidity_difference(RAPIDITIES[1])
    second_weight = -_rapidity_difference(RAPIDITIES[0])
    zero_weight = first_weight + second_weight
    plus = np.asarray(
        (
            first_weight,
            -first_weight,
            0.0,
            second_weight,
            -second_weight,
            0.0,
            -zero_weight,
            zero_weight,
            0.0,
        )
    )
    cross = np.asarray(
        (
            -first_weight,
            -first_weight,
            2.0 * first_weight,
            -second_weight,
            -second_weight,
            2.0 * second_weight,
            zero_weight,
            zero_weight,
            -2.0 * zero_weight,
        )
    )
    return plus, cross


def algebraic_wave_target_certificate(rotation_angle: float = 0.37) -> dict:
    hinges = [
        hinge_tensors(rapidity, transverse)
        for rapidity in RAPIDITIES
        for transverse in TRANSVERSE_DIRECTIONS
    ]
    einstein_map = np.column_stack(
        [symmetric_components(hinge["einstein"]) for hinge in hinges]
    )
    riemann_map = np.column_stack(
        [hinge["riemann"].reshape(-1) for hinge in hinges]
    )
    einstein_rank = int(np.linalg.matrix_rank(einstein_map, tol=1e-11))
    _, _, right = np.linalg.svd(einstein_map, full_matrices=True)
    vacuum_basis = right.T[:, einstein_rank:]
    vacuum_riemann_singular_values = np.linalg.svd(
        riemann_map @ vacuum_basis, compute_uv=False
    )

    plus_coefficients, cross_coefficients = wave_coefficients()
    plus_riemann = (riemann_map @ plus_coefficients).reshape(4, 4, 4, 4)
    cross_riemann = (riemann_map @ cross_coefficients).reshape(4, 4, 4, 4)
    plus_electric, plus_magnetic = electric_magnetic_weyl(plus_riemann)
    cross_electric, cross_magnetic = electric_magnetic_weyl(cross_riemann)

    wave_amplitude = abs(plus_electric[0, 0])
    plus_template = np.diag((1.0, -1.0, 0.0))
    cross_template = np.asarray(
        ((0.0, 1.0, 0.0), (1.0, 0.0, 0.0), (0.0, 0.0, 0.0))
    )
    normalized_plus = plus_electric / plus_electric[0, 0]
    normalized_cross = cross_electric / cross_electric[0, 1]

    spatial_rotation = np.asarray(
        (
            (np.cos(rotation_angle), -np.sin(rotation_angle), 0.0),
            (np.sin(rotation_angle), np.cos(rotation_angle), 0.0),
            (0.0, 0.0, 1.0),
        )
    )
    rotated_plus = spatial_rotation @ plus_template @ spatial_rotation.T
    expected_rotated_plus = (
        np.cos(2.0 * rotation_angle) * plus_template
        + np.sin(2.0 * rotation_angle) * cross_template
    )
    rotated_cross = spatial_rotation @ cross_template @ spatial_rotation.T
    expected_rotated_cross = (
        -np.sin(2.0 * rotation_angle) * plus_template
        + np.cos(2.0 * rotation_angle) * cross_template
    )

    def null_invariants(electric: np.ndarray, magnetic: np.ndarray) -> np.ndarray:
        return np.asarray(
            (
                np.sum(electric * electric) - np.sum(magnetic * magnetic),
                np.sum(electric * magnetic),
            )
        )

    return {
        "hinge_labels": HINGE_LABELS,
        "hinges": hinges,
        "einstein_map": einstein_map,
        "riemann_map": riemann_map,
        "einstein_singular_values": np.linalg.svd(
            einstein_map, compute_uv=False
        ),
        "einstein_rank": einstein_rank,
        "vacuum_coefficient_nullity": len(hinges) - einstein_rank,
        "vacuum_riemann_singular_values": vacuum_riemann_singular_values,
        "plus_coefficients": plus_coefficients,
        "cross_coefficients": cross_coefficients,
        "plus_einstein_residual": float(
            np.max(np.abs(einstein_map @ plus_coefficients))
        ),
        "cross_einstein_residual": float(
            np.max(np.abs(einstein_map @ cross_coefficients))
        ),
        "plus_riemann": plus_riemann,
        "cross_riemann": cross_riemann,
        "plus_electric": plus_electric,
        "plus_magnetic": plus_magnetic,
        "cross_electric": cross_electric,
        "cross_magnetic": cross_magnetic,
        "wave_amplitude": wave_amplitude,
        "plus_template_residual": float(
            np.max(np.abs(normalized_plus - plus_template))
        ),
        "cross_template_residual": float(
            np.max(np.abs(normalized_cross - cross_template))
        ),
        "right_moving_plus_residual": float(
            np.max(np.abs(plus_magnetic - cross_electric))
        ),
        "right_moving_cross_residual": float(
            np.max(np.abs(cross_magnetic + plus_electric))
        ),
        "plus_null_invariants": null_invariants(plus_electric, plus_magnetic),
        "cross_null_invariants": null_invariants(cross_electric, cross_magnetic),
        "spin_two_plus_rotation_residual": float(
            np.max(np.abs(rotated_plus - expected_rotated_plus))
        ),
        "spin_two_cross_rotation_residual": float(
            np.max(np.abs(rotated_cross - expected_rotated_cross))
        ),
        "realized_by_causal_regge_complex": False,
        "internal_edge_equations_imposed": False,
    }


if __name__ == "__main__":
    certificate = algebraic_wave_target_certificate()
    for key in (
        "einstein_singular_values",
        "einstein_rank",
        "vacuum_coefficient_nullity",
        "vacuum_riemann_singular_values",
        "plus_einstein_residual",
        "cross_einstein_residual",
        "wave_amplitude",
        "right_moving_plus_residual",
        "right_moving_cross_residual",
        "plus_null_invariants",
        "cross_null_invariants",
        "spin_two_plus_rotation_residual",
        "spin_two_cross_rotation_residual",
    ):
        print(f"{key}={certificate[key]}")
    print(f"plus_electric=\n{certificate['plus_electric']}")
    print(f"plus_magnetic=\n{certificate['plus_magnetic']}")
    print(f"cross_electric=\n{certificate['cross_electric']}")
    print(f"cross_magnetic=\n{certificate['cross_magnetic']}")
