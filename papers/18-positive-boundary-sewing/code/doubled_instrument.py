"""Exact/numerical audit of the boundary programme doubled-boundary two-packet calibration datum.

This verifies finite matrix identities only.  It does not evaluate a CEZ/PT kernel and
does not establish a continuum or physical quantum-gravity sector.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np


TOL = 1.0e-12


def matrix_unit(n: int, i: int, j: int) -> np.ndarray:
    out = np.zeros((n, n), dtype=complex)
    out[i, j] = 1.0
    return out


def swap_matrix(d: int) -> np.ndarray:
    out = np.zeros((d * d, d * d), dtype=complex)
    for i in range(d):
        for j in range(d):
            out[j * d + i, i * d + j] = 1.0
    return out


def theta(a: np.ndarray, swap: np.ndarray) -> np.ndarray:
    """Anti-linear multiplicative swap reflection J a J^{-1}."""

    return swap @ a.conj() @ swap.conj().T


def heisenberg(a: np.ndarray, kraus: np.ndarray) -> np.ndarray:
    return kraus.conj().T @ a @ kraus


def choi_heisenberg(kraus_ops: list[np.ndarray]) -> np.ndarray:
    """Choi matrix for E(a)=sum K^* a K, with K:H_in->H_out."""

    n_out, n_in = kraus_ops[0].shape
    choi = np.zeros((n_out * n_in, n_out * n_in), dtype=complex)
    for i in range(n_out):
        for j in range(n_out):
            eij = matrix_unit(n_out, i, j)
            block = sum(heisenberg(eij, k) for k in kraus_ops)
            choi[
                i * n_in : (i + 1) * n_in,
                j * n_in : (j + 1) * n_in,
            ] = block
    return choi


def trace_distance(rho: np.ndarray, sigma: np.ndarray) -> float:
    delta = 0.5 * ((rho - sigma) + (rho - sigma).conj().T)
    return 0.5 * float(np.sum(np.abs(np.linalg.eigvalsh(delta))))


def os_gram(rho: np.ndarray, swap: np.ndarray) -> tuple[np.ndarray, list[np.ndarray]]:
    ident = np.eye(2, dtype=complex)
    half_basis = [
        np.kron(ident, matrix_unit(2, i, j)) for i in range(2) for j in range(2)
    ]
    gram = np.empty((len(half_basis), len(half_basis)), dtype=complex)
    for i, a in enumerate(half_basis):
        for j, b in enumerate(half_basis):
            gram[i, j] = np.trace(rho @ theta(a, swap) @ b)
    return gram, half_basis


def build_report() -> str:
    lam = 0.7
    delta = 1.3
    ident2 = np.eye(2, dtype=complex)
    ident4 = np.eye(4, dtype=complex)
    x = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex)
    u = np.kron(x, x)
    k0 = ident4 / np.sqrt(2.0)
    k1 = u / np.sqrt(2.0)
    v = np.vstack([k0, k1])
    swap = swap_matrix(2)

    omega = np.array([np.sqrt(lam), 0.0, 0.0, np.sqrt(1.0 - lam)], dtype=complex)
    rho = np.outer(omega, omega.conj())

    isometry_residual = np.linalg.norm(v.conj().T @ v - ident4, ord=2)
    effect0 = k0.conj().T @ k0
    effect1 = k1.conj().T @ k1
    normalization_residual = np.linalg.norm(effect0 + effect1 - ident4, ord=2)

    choi0 = choi_heisenberg([k0])
    choi1 = choi_heisenberg([k1])
    min_choi = min(
        float(np.min(np.linalg.eigvalsh(0.5 * (choi0 + choi0.conj().T)))),
        float(np.min(np.linalg.eigvalsh(0.5 * (choi1 + choi1.conj().T)))),
    )

    z0, z1 = 0.37, -0.23
    center_image = heisenberg(z0 * ident4, k0) + heisenberg(z1 * ident4, k1)
    sample_center_residual = np.linalg.norm(center_image - 0.5 * (z0 + z1) * ident4)

    sigma0 = k0 @ rho @ k0.conj().T
    sigma1 = k1 @ rho @ k1.conj().T
    p0 = float(np.trace(sigma0).real)
    p1 = float(np.trace(sigma1).real)
    post0 = sigma0 / p0
    post1 = sigma1 / p1
    posterior_fiber_trace_residual = max(
        abs(np.trace(post0) - 1.0), abs(np.trace(post1) - 1.0)
    )
    posterior_fiber_distance = trace_distance(post0, post1)
    posterior_fiber_distance_formula = abs(2.0 * lam - 1.0)
    recorded_post0 = np.zeros((8, 8), dtype=complex)
    recorded_post1 = np.zeros((8, 8), dtype=complex)
    recorded_post0[:4, :4] = post0
    recorded_post1[4:, 4:] = post1
    posterior_recorded_outcome_distance = trace_distance(recorded_post0, recorded_post1)

    rng = np.random.default_rng(20260710)
    a0 = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4))
    a1 = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4))
    b = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4))
    reflected_cut = heisenberg(theta(a0, swap), k0) + heisenberg(theta(a1, swap), k1)
    cut_then_reflect = theta(heisenberg(a0, k0) + heisenberg(a1, k1), swap)
    sample_reflection_covariance_residual = np.linalg.norm(reflected_cut - cut_then_reflect)
    sample_theta_multiplicativity_residual = np.linalg.norm(
        theta(a0 @ b, swap) - theta(a0, swap) @ theta(b, swap)
    )
    sample_theta_star_residual = np.linalg.norm(
        theta(a0.conj().T, swap) - theta(a0, swap).conj().T
    )

    gram, half_basis = os_gram(rho, swap)
    gram_hermitian_residual = np.linalg.norm(gram - gram.conj().T)
    gram_eigs = np.linalg.eigvalsh(0.5 * (gram + gram.conj().T))
    expected_os = np.sort(
        np.array(
            [
                lam,
                np.sqrt(lam * (1.0 - lam)),
                np.sqrt(lam * (1.0 - lam)),
                1.0 - lam,
            ]
        )
    )
    os_spectrum_residual = np.linalg.norm(np.sort(gram_eigs) - expected_os)

    h = np.diag([0.0, delta]).astype(complex)
    generator_eigs = []
    for i in range(2):
        for j in range(2):
            eij = matrix_unit(2, i, j)
            image = h @ eij + eij @ h
            generator_eigs.append(float(image[i, j].real))
    generator_spectrum_residual = np.linalg.norm(
        np.sort(generator_eigs) - np.array([0.0, delta, delta, 2.0 * delta])
    )

    # Invalid controls: each is required to trigger a positive failure signal.
    k1_oversized = 1.1 * u / np.sqrt(2.0)
    oversized_effect = k0.conj().T @ k0 + k1_oversized.conj().T @ k1_oversized
    oversized_superunital_excess = float(
        np.max(np.linalg.eigvalsh(oversized_effect - ident4))
    )
    dropped_outcome_normalization_loss = np.linalg.norm(effect0 - ident4, ord=2)

    omega_singular = np.array([1.0, 0.0, 0.0, 0.0], dtype=complex)
    rho_singular = np.outer(omega_singular, omega_singular.conj())
    singular_gram, _ = os_gram(rho_singular, swap)
    singular_os_min = float(
        np.min(np.linalg.eigvalsh(0.5 * (singular_gram + singular_gram.conj().T)))
    )

    alpha = 1.0 + 1.0j
    # The obsolete form is linear, rather than conjugate-linear, in its first slot.
    # Use the half-algebra identity so the base value is nonzero and the variance defect
    # cannot be accidentally hidden by an orthogonal matrix unit.
    test_a = half_basis[0] + half_basis[3]
    test_b = half_basis[0] + half_basis[3]

    def obsolete_form(a: np.ndarray, bb: np.ndarray) -> complex:
        return np.trace(rho @ theta(a, swap).conj().T @ bb)

    obsolete_conjugate_linearity_defect = abs(
        obsolete_form(alpha * test_a, test_b)
        - np.conj(alpha) * obsolete_form(test_a, test_b)
    )

    checks = {
        "isometry": isometry_residual <= TOL,
        "instrument_normalization": normalization_residual <= TOL,
        "choi_nonnegative": min_choi >= -TOL,
        "sample_center_map": sample_center_residual <= TOL,
        "posterior_fiber_normalization": posterior_fiber_trace_residual <= TOL,
        "posterior_fiber_distance_formula": abs(
            posterior_fiber_distance - posterior_fiber_distance_formula
        )
        <= TOL,
        "posterior_recorded_outcome_distance": abs(
            posterior_recorded_outcome_distance - 1.0
        )
        <= TOL,
        "sample_reflection_covariance": sample_reflection_covariance_residual <= TOL,
        "sample_theta_multiplicative": sample_theta_multiplicativity_residual <= TOL,
        "sample_theta_star": sample_theta_star_residual <= TOL,
        "os_hermitian": gram_hermitian_residual <= TOL,
        "os_positive_definite": float(np.min(gram_eigs)) > TOL,
        "os_spectrum_formula": os_spectrum_residual <= TOL,
        "generator_spectrum": generator_spectrum_residual <= TOL,
        "invalid_oversized_detected": oversized_superunital_excess > TOL,
        "invalid_dropped_outcome_detected": dropped_outcome_normalization_loss > TOL,
        "invalid_singular_os_detected": (
            singular_os_min >= -TOL and abs(singular_os_min) <= TOL
        ),
        "obsolete_reflection_detected": obsolete_conjugate_linearity_defect > TOL,
    }

    lines = [
        "the boundary programme doubled-boundary two-packet finite operational check",
        f"lambda={lam:.16g}",
        f"Delta={delta:.16g}",
        f"tolerance={TOL:.3e}",
        "scope=exact finite calibration; no CEZ/PT kernel and no maturity promotion",
        f"isometry_residual={isometry_residual:.16e}",
        f"normalization_residual={normalization_residual:.16e}",
        f"minimum_outcome_choi_eigenvalue={min_choi:.16e}",
        f"sample_center_map_residual={sample_center_residual:.16e}",
        f"outcome_probability_0={p0:.16e}",
        f"outcome_probability_1={p1:.16e}",
        f"posterior_fiber_trace_residual={posterior_fiber_trace_residual:.16e}",
        f"posterior_fiber_trace_distance={posterior_fiber_distance:.16e}",
        f"posterior_fiber_trace_distance_formula={posterior_fiber_distance_formula:.16e}",
        f"posterior_recorded_outcome_trace_distance={posterior_recorded_outcome_distance:.16e}",
        f"sample_reflection_covariance_residual={sample_reflection_covariance_residual:.16e}",
        f"sample_theta_multiplicativity_residual={sample_theta_multiplicativity_residual:.16e}",
        f"sample_theta_star_residual={sample_theta_star_residual:.16e}",
        f"os_gram_hermitian_residual={gram_hermitian_residual:.16e}",
        "os_gram_eigenvalues=" + ",".join(f"{x:.16e}" for x in gram_eigs),
        f"os_spectrum_formula_residual={os_spectrum_residual:.16e}",
        "os_generator_eigenvalues="
        + ",".join(f"{x:.16e}" for x in np.sort(generator_eigs)),
        f"os_generator_spectrum_residual={generator_spectrum_residual:.16e}",
        f"invalid_oversized_superunital_excess={oversized_superunital_excess:.16e}",
        f"invalid_dropped_outcome_normalization_loss={dropped_outcome_normalization_loss:.16e}",
        f"invalid_singular_os_minimum={singular_os_min:.16e}",
        f"obsolete_reflection_conjugate_linearity_defect={obsolete_conjugate_linearity_defect:.16e}",
    ]
    lines.extend(f"check_{name}={'PASS' if ok else 'FAIL'}" for name, ok in checks.items())
    lines.append(f"overall={'PASS' if all(checks.values()) else 'FAIL'}")
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true", help="write the checked output file")
    args = parser.parse_args()
    report = build_report()
    print(report, end="")
    if args.write:
        Path(__file__).with_name("doubled_instrument_output.txt").write_text(
            report, encoding="utf-8"
        )


if __name__ == "__main__":
    main()
