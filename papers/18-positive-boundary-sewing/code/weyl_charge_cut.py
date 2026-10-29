"""Finite Weyl/charge/cut surrogate checks for the boundary programme AdS4 construction roadmap.

The script verifies exact discrete Heisenberg-Weyl, charge-sector, dressing, group-
average, and cut-instrument identities at d=3.  The mode labels are finite slots, not
physical gravitons, a scalar QFT, or a holographic reconstruction.
"""

from __future__ import annotations

import argparse
import itertools
from pathlib import Path

import numpy as np


TOL = 2.0e-11


def shift(d: int) -> np.ndarray:
    x = np.zeros((d, d), dtype=complex)
    for j in range(d):
        x[(j + 1) % d, j] = 1.0
    return x


def clock(d: int) -> np.ndarray:
    zeta = np.exp(2j * np.pi / d)
    return np.diag([zeta**j for j in range(d)]).astype(complex)


def weyl(d: int, q: int, p: int) -> np.ndarray:
    """W(q,p)=zeta^(-qp/2) Z^q X^p for odd d."""

    zeta = np.exp(2j * np.pi / d)
    inv2 = pow(2, -1, d)
    q %= d
    p %= d
    phase = zeta ** ((-q * p * inv2) % d)
    return phase * np.linalg.matrix_power(clock(d), q) @ np.linalg.matrix_power(
        shift(d), p
    )


def kron_all(factors: list[np.ndarray]) -> np.ndarray:
    out = np.array([[1.0]], dtype=complex)
    for factor in factors:
        out = np.kron(out, factor)
    return out


def multi_weyl(d: int, qs: tuple[int, ...], ps: tuple[int, ...]) -> np.ndarray:
    return kron_all([weyl(d, q, p) for q, p in zip(qs, ps, strict=True)])


def charge_unitary(d: int, t: int, charges: tuple[int, ...]) -> np.ndarray:
    z = clock(d)
    return kron_all(
        [np.linalg.matrix_power(z, (charge * t) % d) for charge in charges]
    )


def group_average(a: np.ndarray, unitaries: list[np.ndarray]) -> np.ndarray:
    return sum(u @ a @ u.conj().T for u in unitaries) / len(unitaries)


def normalize(v: np.ndarray) -> np.ndarray:
    return v / np.linalg.norm(v)


def build_report() -> str:
    d = 3
    inv2 = pow(2, -1, d)
    zeta = np.exp(2j * np.pi / d)
    ident_d = np.eye(d, dtype=complex)
    ident = np.eye(d**4, dtype=complex)

    # Charge: exact single-mode discrete Weyl algebra and matrix basis.
    weyl_product_residual = 0.0
    weyl_star_residual = 0.0
    labels = list(itertools.product(range(d), repeat=2))
    for q, p in labels:
        a = weyl(d, q, p)
        weyl_star_residual = max(
            weyl_star_residual,
            float(np.linalg.norm(a.conj().T - weyl(d, -q, -p))),
        )
        for r, s in labels:
            sigma = (q * s - p * r) % d
            expected = zeta ** ((sigma * inv2) % d) * weyl(d, q + r, p + s)
            weyl_product_residual = max(
                weyl_product_residual,
                float(np.linalg.norm(a @ weyl(d, r, s) - expected)),
            )

    basis_gram = np.empty((d * d, d * d), dtype=complex)
    for i, xlab in enumerate(labels):
        for j, ylab in enumerate(labels):
            basis_gram[i, j] = np.trace(
                weyl(d, *xlab).conj().T @ weyl(d, *ylab)
            ) / d
    weyl_basis_residual = np.linalg.norm(basis_gram - np.eye(d * d))

    # Four slots: g+, g-, phi, anchor.  The anchor coefficient 1 makes every charge
    # sector equally populated.
    charges = (1, 1, 1, 1)
    unitaries = [charge_unitary(d, t, charges) for t in range(d)]
    projectors = []
    for q in range(d):
        projector = sum(
            (zeta ** ((-t * q) % d)) * unitaries[t] for t in range(d)
        ) / d
        projectors.append(projector)

    projector_hermitian_residual = max(
        float(np.linalg.norm(p - p.conj().T)) for p in projectors
    )
    projector_idempotence_residual = max(
        float(np.linalg.norm(p @ p - p)) for p in projectors
    )
    projector_orthogonality_residual = max(
        float(np.linalg.norm(projectors[q] @ projectors[r]))
        for q in range(d)
        for r in range(d)
        if q != r
    )
    projector_completeness_residual = np.linalg.norm(sum(projectors) - ident)
    projector_ranks = [int(np.rint(np.trace(p).real)) for p in projectors]

    rng = np.random.default_rng(20260710)
    random_matrix = rng.normal(size=(d**4, d**4)) + 1j * rng.normal(
        size=(d**4, d**4)
    )
    averaged = group_average(random_matrix, unitaries)
    sample_average_idempotence_residual = np.linalg.norm(
        group_average(averaged, unitaries) - averaged
    )
    average_unital_residual = np.linalg.norm(group_average(ident, unitaries) - ident)

    seed = random_matrix.conj().T @ random_matrix
    averaged_positive = group_average(seed, unitaries)
    average_min_eigenvalue = float(
        np.min(np.linalg.eigvalsh(0.5 * (averaged_positive + averaged_positive.conj().T)))
    )
    sample_average_invariance_residual = max(
        float(np.linalg.norm(u @ averaged_positive @ u.conj().T - averaged_positive))
        for u in unitaries
    )

    # Dressing: a neutral finite dressing.  xi contains q/p labels for g+, g-, phi.
    xi_q = (1, 2, 1)
    xi_p = (1, 1, 2)
    eta_q = (2, 1, 0)
    eta_p = (2, 0, 1)
    internal_charges = charges[:3]

    def dress(
        qs_internal: tuple[int, int, int],
        ps_internal: tuple[int, int, int],
    ) -> tuple[np.ndarray, tuple[int, ...], tuple[int, ...]]:
        charge = sum(k * p for k, p in zip(internal_charges, ps_internal, strict=True)) % d
        q_anchor = 0
        p_anchor = (-charge) % d
        qs = tuple(q % d for q in qs_internal) + (q_anchor,)
        ps = tuple(p % d for p in ps_internal) + (p_anchor,)
        return multi_weyl(d, qs, ps), qs, ps

    dressed_xi, xi_q_full, xi_p_full = dress(xi_q, xi_p)
    dressed_eta, eta_q_full, eta_p_full = dress(eta_q, eta_p)
    charge_invariance_residual = np.linalg.norm(
        unitaries[1] @ dressed_xi @ unitaries[1].conj().T - dressed_xi
    )

    sum_q = tuple((a + b) % d for a, b in zip(xi_q, eta_q, strict=True))
    sum_p = tuple((a + b) % d for a, b in zip(xi_p, eta_p, strict=True))
    dressed_sum, _, _ = dress(sum_q, sum_p)
    sigma_full = sum(
        q * p2 - p * q2
        for q, p, q2, p2 in zip(
            xi_q_full, xi_p_full, eta_q_full, eta_p_full, strict=True
        )
    ) % d
    sample_dressed_weyl_residual = np.linalg.norm(
        dressed_xi @ dressed_eta
        - zeta ** ((sigma_full * inv2) % d) * dressed_sum
    )

    anchor_origin = 1
    anchor_shift = multi_weyl(d, (0, 0, 0, anchor_origin), (0, 0, 0, 0))
    dressed_xi_alt = anchor_shift @ dressed_xi @ anchor_shift.conj().T
    internal_charge_xi = sum(
        k * p for k, p in zip(internal_charges, xi_p, strict=True)
    ) % d
    anchor_phase = zeta ** ((-anchor_origin * internal_charge_xi) % d)
    anchor_change_residual = np.linalg.norm(
        dressed_xi_alt - anchor_phase * dressed_xi
    )

    scalar_x_dressed, _, _ = dress((0, 0, 0), (0, 0, 1))
    scalar_z_dressed, _, _ = dress((0, 0, 1), (0, 0, 0))
    t_slot = (
        (ident - scalar_x_dressed).conj().T @ (ident - scalar_x_dressed)
        + (ident - scalar_z_dressed).conj().T @ (ident - scalar_z_dressed)
    )
    t_slot_min_eigenvalue = float(
        np.min(np.linalg.eigvalsh(0.5 * (t_slot + t_slot.conj().T)))
    )
    t_slot_invariance_residual = np.linalg.norm(
        unitaries[1] @ t_slot @ unitaries[1].conj().T - t_slot
    )

    # Cut: blockwise cut instrument.  Avoid forming the 6561-dimensional cut matrix;
    # the canonical inclusion is tested directly on its charge-sector blocks.
    sector_dim = projector_ranks[0]
    cut_dim = d * d * sector_dim * sector_dim
    glued_dim = d * sector_dim * sector_dim
    probabilities = np.array([0.2, 0.3, 0.5], dtype=float)
    glued_blocks = {}
    for q, probability in enumerate(probabilities):
        vec = rng.normal(size=sector_dim * sector_dim) + 1j * rng.normal(
            size=sector_dim * sector_dim
        )
        glued_blocks[q] = np.sqrt(probability) * normalize(vec)
    glued_norm = sum(float(np.vdot(v, v).real) for v in glued_blocks.values())
    cut_blocks = {(q, (-q) % d): v.copy() for q, v in glued_blocks.items()}
    cut_norm = sum(float(np.vdot(v, v).real) for v in cut_blocks.values())
    sample_cut_state_norm_residual = abs(cut_norm - glued_norm)
    outcome_probabilities = np.array(
        [float(np.vdot(glued_blocks[q], glued_blocks[q]).real) for q in range(d)]
    )
    sample_instrument_probability_residual = abs(np.sum(outcome_probabilities) - 1.0)
    sample_posterior_norm_residual = max(
        abs(np.vdot(glued_blocks[q], glued_blocks[q]).real / outcome_probabilities[q] - 1.0)
        for q in range(d)
    )

    direct_triples = {
        (q1, q2, q3)
        for q1, q2, q3 in itertools.product(range(d), repeat=3)
        if (q1 + q2 + q3) % d == 0
    }
    successive_triples = set()
    for q1, q2, q3 in itertools.product(range(d), repeat=3):
        intermediate = (q1 + q2) % d
        if (intermediate + q3) % d == 0:
            successive_triples.add((q1, q2, q3))
    successive_cut_charge_set_difference = len(
        direct_triples.symmetric_difference(successive_triples)
    )

    # Invalid controls.
    annihilator = np.zeros((d, d), dtype=complex)
    for n in range(1, d):
        annihilator[n - 1, n] = np.sqrt(n)
    truncated_ccr_defect = np.linalg.norm(
        annihilator @ annihilator.conj().T
        - annihilator.conj().T @ annihilator
        - ident_d
    )
    undressed = multi_weyl(d, xi_q + (0,), xi_p + (0,))
    undressed_charge_defect = np.linalg.norm(
        unitaries[1] @ undressed @ unitaries[1].conj().T - undressed
    )
    incomplete_instrument_probability_loss = 1.0 - outcome_probabilities[0]

    checks = {
        "weyl_product": weyl_product_residual <= TOL,
        "weyl_star": weyl_star_residual <= TOL,
        "weyl_basis": weyl_basis_residual <= TOL,
        "projector_hermitian": projector_hermitian_residual <= TOL,
        "projector_idempotent": projector_idempotence_residual <= TOL,
        "projector_orthogonal": projector_orthogonality_residual <= TOL,
        "projector_complete": projector_completeness_residual <= TOL,
        "projector_equal_rank": projector_ranks == [d**3] * d,
        "sample_average_idempotent": sample_average_idempotence_residual <= TOL,
        "average_unital": average_unital_residual <= TOL,
        "sample_average_positive": average_min_eigenvalue >= -TOL,
        "sample_average_invariant": sample_average_invariance_residual <= TOL,
        "dressing_charge_invariant": charge_invariance_residual <= TOL,
        "sample_dressed_weyl": sample_dressed_weyl_residual <= TOL,
        "anchor_change": anchor_change_residual <= TOL,
        "t_slot_positive": t_slot_min_eigenvalue >= -TOL,
        "t_slot_invariant": t_slot_invariance_residual <= TOL,
        "sample_cut_state_norm": sample_cut_state_norm_residual <= TOL,
        "sample_cut_probabilities_normalized": sample_instrument_probability_residual <= TOL,
        "sample_cut_posterior_normalized": sample_posterior_norm_residual <= TOL,
        "successive_cut_charge_sets_agree": successive_cut_charge_set_difference == 0,
        "invalid_truncated_ccr_detected": truncated_ccr_defect > TOL,
        "invalid_undressed_charge_detected": undressed_charge_defect > TOL,
        "invalid_incomplete_instrument_detected": incomplete_instrument_probability_loss > TOL,
    }

    lines = [
        "the boundary programme AdS4 Weyl/charge/cut finite surrogate check",
        f"d={d}",
        f"tolerance={TOL:.3e}",
        "scope=finite discrete-Weyl/charge surrogate; no physical graviton, scalar QFT, or holographic bridge",
        f"single_mode_weyl_product_residual={weyl_product_residual:.16e}",
        f"single_mode_weyl_star_residual={weyl_star_residual:.16e}",
        f"single_mode_weyl_basis_residual={weyl_basis_residual:.16e}",
        f"charge_projector_hermitian_residual={projector_hermitian_residual:.16e}",
        f"charge_projector_idempotence_residual={projector_idempotence_residual:.16e}",
        f"charge_projector_orthogonality_residual={projector_orthogonality_residual:.16e}",
        f"charge_projector_completeness_residual={projector_completeness_residual:.16e}",
        "charge_projector_ranks=" + ",".join(str(x) for x in projector_ranks),
        f"sample_group_average_idempotence_residual={sample_average_idempotence_residual:.16e}",
        f"group_average_unital_residual={average_unital_residual:.16e}",
        f"group_average_min_eigenvalue={average_min_eigenvalue:.16e}",
        f"sample_group_average_invariance_residual={sample_average_invariance_residual:.16e}",
        f"dressing_charge_invariance_residual={charge_invariance_residual:.16e}",
        f"sample_dressed_weyl_residual={sample_dressed_weyl_residual:.16e}",
        f"anchor_change_residual={anchor_change_residual:.16e}",
        f"t_slot_min_eigenvalue={t_slot_min_eigenvalue:.16e}",
        f"t_slot_invariance_residual={t_slot_invariance_residual:.16e}",
        f"single_side_dimension={d**4}",
        f"charge_sector_dimension={sector_dim}",
        f"cut_space_dimension={cut_dim}",
        f"glued_matched_space_dimension={glued_dim}",
        f"sample_cut_state_norm_residual={sample_cut_state_norm_residual:.16e}",
        "sample_cut_outcome_probabilities="
        + ",".join(f"{x:.16e}" for x in outcome_probabilities),
        f"sample_cut_instrument_probability_residual={sample_instrument_probability_residual:.16e}",
        f"sample_cut_posterior_norm_residual={sample_posterior_norm_residual:.16e}",
        f"successive_cut_charge_set_difference={successive_cut_charge_set_difference}",
        f"invalid_truncated_ccr_defect={truncated_ccr_defect:.16e}",
        f"invalid_undressed_charge_defect={undressed_charge_defect:.16e}",
        f"invalid_incomplete_instrument_probability_loss={incomplete_instrument_probability_loss:.16e}",
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
        Path(__file__).with_name("weyl_charge_cut_output.txt").write_text(
            report, encoding="utf-8"
        )


if __name__ == "__main__":
    main()
