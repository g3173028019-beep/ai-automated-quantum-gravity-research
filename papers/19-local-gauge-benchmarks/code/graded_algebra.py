"""Finite algebra audit for the local benchmark graded algebra equicausal BRST star algebra.

This checks only the coefficient identities behind the exact free Ward condition,
the graded-dual field sign, the scalar direct sum, the equation-ideal compatibility,
the first-order graded commutator, the ordinary kernel adjoint/star anti-involution,
and the single hbar normalization.  The
equicausal closure theorem and all wavefront/equicontinuity claims are analytic
source-specializations and are not proved by this script.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from sympy import I, Matrix, Rational, diag, expand, simplify, symbols


def reduce_source_relations(expr, p, k, d, r):
    """Use only P K=0, D P=0, and D K=R."""

    return simplify(expand(expr).subs(d * k, r).subs(p * k, 0).subs(d * p, 0))


def graded_difference(plus: Matrix, minus: Matrix, ghost_numbers) -> Matrix:
    return Matrix(
        plus.rows,
        plus.cols,
        lambda i, j: simplify(
            plus[i, j]
            - (-1 if (ghost_numbers[i] * ghost_numbers[j]) % 2 else 1)
            * minus[i, j]
        ),
    )


def graded_transpose(matrix: Matrix, parities) -> Matrix:
    """Form-level graded argument exchange for a homogeneous component matrix."""

    return Matrix(
        matrix.cols,
        matrix.rows,
        lambda i, j: (-1 if (parities[i] * parities[j]) % 2 else 1)
        * matrix[j, i],
    )


def main() -> str:
    p, k, d, r = symbols("P K D R", commutative=False)
    l2p, l2m, l1p, l1m = symbols(
        "lambda_2_plus lambda_2_minus lambda_1_plus lambda_1_minus",
        commutative=False,
    )
    rho_p, rho_m = 2 * l1p, 2 * l1m

    # Grouped Green-matrix source order: (Theta,beta,eta,bar_eta).
    lambda_plus = Matrix(
        [
            [-l2p, k * rho_p, 0, 0],
            [rho_p * d, 0, 0, 0],
            [0, 0, 0, rho_p],
            [0, 0, -rho_p, 0],
        ]
    )
    lambda_minus = Matrix(
        [
            [-l2m, k * rho_m, 0, 0],
            [rho_m * d, 0, 0, 0],
            [0, 0, 0, -rho_m],
            [0, 0, rho_m, 0],
        ]
    )
    gamma = Matrix(
        [
            [0, 0, k, 0],
            [0, 0, 0, 0],
            [0, 0, 0, 0],
            [0, 1, 0, 0],
        ]
    )
    gamma_star = Matrix(
        [
            [0, 0, 0, 0],
            [0, 0, 0, 1],
            [-d, 0, 0, 0],
            [0, 0, 0, 0],
        ]
    )

    ward_left = gamma * lambda_plus
    ward_right = lambda_plus * gamma_star
    expected_ward = Matrix(
        [
            [0, 0, 0, k * rho_p],
            [0, 0, 0, 0],
            [0, 0, 0, 0],
            [rho_p * d, 0, 0, 0],
        ]
    )
    check_exact_plus_covariance_ward = ward_left == ward_right == expected_ward

    # The bare equality above is not itself Rejzner's two-argument functional
    # Ward identity.  The graded-dual parity matrix supplies the needed sign.
    parity_rows = diag(1, 1, -1, -1)
    plus_ward_defect = (
        parity_rows * gamma * lambda_plus
        + lambda_plus * gamma_star * parity_rows
    )
    check_rejzner_two_argument_ward_with_Pi = plus_ward_defect == Matrix.zeros(4)

    # Hadamard-pair's Schur-signed minus covariance has the same plus-pattern and obeys
    # the same two-argument Ward identity.  The unschur-signed minus kernel does not.
    ghost_numbers = (0, 0, 1, -1)
    schur_minus = Matrix(
        4,
        4,
        lambda i, j: (-1 if (ghost_numbers[i] * ghost_numbers[j]) % 2 else 1)
        * lambda_minus[i, j],
    )
    schur_minus_defect = (
        parity_rows * gamma * schur_minus
        + schur_minus * gamma_star * parity_rows
    )
    bare_minus_defect = (
        parity_rows * gamma * lambda_minus
        + lambda_minus * gamma_star * parity_rows
    )
    check_schur_minus_ward = schur_minus_defect == Matrix.zeros(4)
    check_bare_minus_rejected = bare_minus_defect != Matrix.zeros(4)

    # Hadamard-pair Hermiticity is the ordinary sesquilinear adjoint, not graded exchange.
    # The bundle pairings are normalized to +1 on eta and -1 on bar_eta.
    b_odd = diag(1, -1)
    odd_indices = [2, 3]
    w_plus_odd = b_odd * lambda_plus.extract(odd_indices, odd_indices)
    w_minus_odd = b_odd * schur_minus.extract(odd_indices, odd_indices)
    odd_parities = (1, 1)
    check_q2_odd_form_ordinary_hermitian = (
        w_plus_odd.T == w_plus_odd and w_minus_odd.T == w_minus_odd
    )
    check_invalid_graded_kernel_dagger_rejected = (
        graded_transpose(w_plus_odd, odd_parities) != w_plus_odd
    )
    w_q3_odd = (
        w_plus_odd + graded_transpose(w_minus_odd, odd_parities)
    ) / 2
    check_q3_odd_form_ordinary_hermitian = w_q3_odd.T == w_q3_odd

    # Two self-adjoint odd linear generators obey psi*theta=-theta*psi.  Reversal
    # under dagger and ordinary Hermiticity of W give the ungraded anti-involution.
    theta_psi, hbar_test = symbols("theta_psi hbar_test")
    dagger_theta_star_psi = -theta_psi + hbar_test * w_q3_odd[0, 1]
    psi_star_theta = -theta_psi + hbar_test * w_q3_odd[1, 0]
    check_two_odd_generator_star_anti_involution = (
        simplify(dagger_theta_star_psi - psi_star_theta) == 0
    )

    # An illicit antighost/ghost sign flip destroys the exact Ward identity.
    invalid_plus = lambda_plus.copy()
    invalid_plus[3, 2] = rho_p
    check_invalid_odd_sign_rejected = gamma * invalid_plus != invalid_plus * gamma_star

    # The functional differential is the graded transpose of gamma, not gamma
    # componentwise.  In this source order s(bar_eta)=-beta.
    s_on_fields = parity_rows * gamma
    check_graded_dual_field_sign = (
        s_on_fields[0, 2] == k
        and s_on_fields[3, 1] == -1
        and s_on_fields * s_on_fields == Matrix.zeros(4)
    )

    # A BRST-inert scalar spectator preserves the exact Ward identity.
    lambda_phi = symbols("lambda_phi", commutative=False)
    lambda_plus_5 = diag(lambda_plus, lambda_phi)
    gamma_5 = diag(gamma, 0)
    gamma_star_5 = diag(gamma_star, 0)
    check_scalar_direct_sum_ward = gamma_5 * lambda_plus_5 == lambda_plus_5 * gamma_star_5

    # Equation-ideal stability: gamma^star P = P gamma at lambda=1/2.
    p_q = Matrix(
        [
            [-p, k, 0, 0],
            [d, Rational(1, 2), 0, 0],
            [0, 0, 0, -r],
            [0, 0, r, 0],
        ]
    )
    p_defect = gamma_star * p_q - p_q * gamma
    p_defect = p_defect.applyfunc(lambda x: reduce_source_relations(x, p, k, d, r))
    check_equation_ideal_brst_stability = p_defect == Matrix.zeros(4)

    # The antisymmetric part of the unit covariance gives the exact Green-matrix causal
    # kernel, hence the linear graded star commutator is i*hbar times Peierls.
    g2, g1 = symbols("G_2 G_1", commutative=False)
    causal = Matrix(
        [
            [-g2, 2 * k * g1, 0, 0],
            [2 * g1 * d, 0, 0, 0],
            [0, 0, 0, 2 * g1],
            [0, 0, -2 * g1, 0],
        ]
    )
    gd = graded_difference(lambda_plus, lambda_minus, ghost_numbers)
    gd = gd.subs({l2p: l2m + I * g2, l1p: l1m + I * g1}).applyfunc(simplify)
    check_linear_star_commutator_kernel = gd == I * causal

    # Abstract finite crossing symmetrization.  If W+ - W- = iG and G is
    # graded antisymmetric, W=(W+ + tau_gr W-)/2 has W-tau_gr(W)=iG.
    parities = (0, 0, 1, 1)
    x = Matrix(4, 4, symbols("x0:16"))
    y = Matrix(4, 4, symbols("y0:16"))
    g_form = (y - graded_transpose(y, parities)) / 2
    minus_form = x - I * g_form
    w_q3 = (x + graded_transpose(minus_form, parities)) / 2
    crossing_defect = (
        w_q3 - graded_transpose(w_q3, parities) - I * g_form
    ).applyfunc(simplify)
    check_q2_pair_crossing_symmetrization = crossing_defect == Matrix.zeros(4)

    # Hadamard-pair already put c0 into W_unit.  Multiplying once by hbar_LB gives hbar_H;
    # multiplying by c0 again is the rejected double count.
    c0, hbar_lb, hbar_h = symbols("c_0 hbar_LB hbar_H", nonzero=True)
    first_contraction = (hbar_lb * c0).subs(hbar_lb, hbar_h / c0)
    double_count = (hbar_lb * c0**2).subs(hbar_lb, hbar_h / c0)
    check_single_hbar_transport = simplify(first_contraction - hbar_h) == 0
    check_invalid_double_c0_rejected = simplify(double_count - hbar_h) != 0

    checks = {
        "check_exact_plus_covariance_ward": check_exact_plus_covariance_ward,
        "check_rejzner_two_argument_ward_with_Pi": check_rejzner_two_argument_ward_with_Pi,
        "check_schur_minus_ward": check_schur_minus_ward,
        "check_bare_minus_rejected": check_bare_minus_rejected,
        "check_q2_odd_form_ordinary_hermitian": (
            check_q2_odd_form_ordinary_hermitian
        ),
        "check_invalid_graded_kernel_dagger_rejected": (
            check_invalid_graded_kernel_dagger_rejected
        ),
        "check_q3_odd_form_ordinary_hermitian": (
            check_q3_odd_form_ordinary_hermitian
        ),
        "check_two_odd_generator_star_anti_involution": (
            check_two_odd_generator_star_anti_involution
        ),
        "check_invalid_odd_sign_rejected": check_invalid_odd_sign_rejected,
        "check_graded_dual_field_sign_and_nilpotence": check_graded_dual_field_sign,
        "check_scalar_direct_sum_ward": check_scalar_direct_sum_ward,
        "check_equation_ideal_brst_stability": check_equation_ideal_brst_stability,
        "check_linear_star_commutator_kernel": check_linear_star_commutator_kernel,
        "check_q2_pair_crossing_symmetrization": check_q2_pair_crossing_symmetrization,
        "check_single_hbar_transport": check_single_hbar_transport,
        "check_invalid_double_c0_rejected": check_invalid_double_c0_rejected,
    }
    overall = all(checks.values())

    lines = [
        "the boundary programme local benchmark graded algebra equicausal BRST-star symbolic audit",
        "scope=finite graded covariance/Ward/sign/normalization identities only",
        "analytic_input=equicausal closure and Hadamard wavefront estimates are source-specialized",
        "nonclaims=no unrestricted microcausal closure; no physical faithfulness; no time-slice; no interactions/QME",
        "source_order=(Theta,beta,eta,bar_eta)",
        "bare_ward=gamma*lambda_plus=lambda_plus*gamma_star=M",
        "two_argument_ward=Pi*gamma*lambda_plus+lambda_plus*gamma_star*Pi=0",
        "crossing_kernel=W_graded-algebra=(W_plus+tau_gr(W_minus_signed))/2",
        "kernel_adjoint=conjugate_tau_0_not_conjugate_tau_gr",
        "peierls_derivative_order=right derivative on first factor; left derivative on second",
        "star_derivative_order=left derivative on first factor; right derivative on second",
        "functional_sign=s_0_ell=(-1)^|ell| ell*gamma",
        "normalization=hbar_LB*(c0*W_H)=hbar_H*W_H exactly once",
        "operator_placeholders=noncommuting",
        "causal_difference=lambda_plus-(-1)^gh*lambda_minus=i*G",
    ]
    lines.extend(f"{name}={'PASS' if value else 'FAIL'}" for name, value in checks.items())
    lines.append(f"overall={'PASS_ALGEBRA_ONLY' if overall else 'FAIL'}")
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    output = main()
    print(output, end="")
    if args.write:
        target = Path(__file__).with_name(
            "graded_algebra_output.txt"
        )
        target.write_text(output, encoding="utf-8")
