"""Exact coefficient audit for the local benchmark Hadamard-pair BRST-lift convention map.

This checks the finite 4x4 row/field rescalings, the graded pairing duality,
the transported covariance blocks and causal Green matrix, and the free BRST matrix
identity.  It does not construct the host Hadamard covariances, prove the Fredholm
pair theorem, prove BRST-cohomology faithfulness, or establish a microlocal/Wick
algebra.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from sympy import I, Matrix, Rational, diag, expand, simplify, symbols


def reduce_source_relations(expr, p, k, d, r):
    """Apply only P K=0, D P=0, and D K=R to a scalar block entry."""

    return simplify(
        expand(expr).subs(d * k, r).subs(p * k, 0).subs(d * p, 0)
    )


def main() -> str:
    # Operator blocks are deliberately noncommuting.  Only the constant gauge
    # parameter is central, so an accidental reordering cannot make a check pass.
    p, k, d, r = symbols("P K D R", commutative=False)
    lambda_aux = symbols("lambda_aux")

    # Wrochna-Zahn order and matched alpha=1: (u,b,c,cbar).
    l_wz = Matrix(
        [
            [p, 2 * k, 0, 0],
            [d, -1, 0, 0],
            [0, 0, 0, 2 * r],
            [0, 0, 2 * r, 0],
        ]
    )
    s_cfg = diag(2, -1, 1, -1)
    s_row = diag(Rational(-1, 2), Rational(1, 2), Rational(1, 2), Rational(1, 2))

    # Reordered source fields: (Theta,beta,eta,bar_eta).
    source_half = Matrix(
        [
            [-p, k, 0, 0],
            [d, Rational(1, 2), 0, 0],
            [0, 0, 0, -r],
            [0, 0, r, 0],
        ]
    )
    check_operator_map = s_row * l_wz * s_cfg == source_half

    # Pairing audit in grouped order (Theta,beta,eta,bar_eta).  Normalize the
    # Hack/vector forms to one.  T_g=2 and
    # B_G(T_g x,T_g y)=2 B_H(x,y) imply B_G=B_H/2 in the common coordinates.
    b_wz = diag(Rational(1, 2), 1, 1, 1)
    b_q1 = diag(-1, -1, 1, -1)
    half = Rational(1, 2)
    check_pairing_duality = b_q1 * s_row == half * s_cfg.T * b_wz
    check_pairing_duality_invalid_no_half = b_q1 * s_row != s_cfg.T * b_wz

    lambda_2_plus, lambda_2_minus = symbols(
        "lambda_2_plus lambda_2_minus", commutative=False
    )
    lambda_1_plus, lambda_1_minus = symbols(
        "lambda_1_plus lambda_1_minus", commutative=False
    )
    lambda_wz_plus = Matrix(
        [
            [lambda_2_plus, 2 * k * lambda_1_plus, 0, 0],
            [lambda_1_plus * d, 0, 0, 0],
            [0, 0, 0, lambda_1_plus],
            [0, 0, lambda_1_plus, 0],
        ]
    )
    lambda_wz_minus = Matrix(
        [
            [lambda_2_minus, 2 * k * lambda_1_minus, 0, 0],
            [lambda_1_minus * d, 0, 0, 0],
            [0, 0, 0, -lambda_1_minus],
            [0, 0, -lambda_1_minus, 0],
        ]
    )
    lambda_q1_plus = s_cfg.inv() * lambda_wz_plus * s_row.inv()
    lambda_q1_minus = s_cfg.inv() * lambda_wz_minus * s_row.inv()
    rho_plus = 2 * lambda_1_plus
    rho_minus = 2 * lambda_1_minus
    expected_q1_plus = Matrix(
        [
            [-lambda_2_plus, k * rho_plus, 0, 0],
            [rho_plus * d, 0, 0, 0],
            [0, 0, 0, rho_plus],
            [0, 0, -rho_plus, 0],
        ]
    )
    expected_q1_minus = Matrix(
        [
            [-lambda_2_minus, k * rho_minus, 0, 0],
            [rho_minus * d, 0, 0, 0],
            [0, 0, 0, -rho_minus],
            [0, 0, rho_minus, 0],
        ]
    )
    check_covariance_plus_transport = lambda_q1_plus == expected_q1_plus
    check_covariance_minus_transport = lambda_q1_minus == expected_q1_minus

    # At form level the same operator transport is exactly one half of the WZ
    # covariance pulled back by the inverse row/source map.  This is the already
    # required Gerard-to-Hack factor, not an additional rescaling of lambda_WZ.
    check_covariance_form_half_pullback = (
        b_q1 * lambda_q1_plus
        == half * s_row.inv().T * b_wz * lambda_wz_plus * s_row.inv()
    )
    wrong_same_side_transport = s_cfg.inv() * lambda_wz_plus * s_cfg.inv()
    check_invalid_same_side_covariance_transport_rejected = (
        wrong_same_side_transport != expected_q1_plus
    )

    # The WZ vector operator is D_1=2R, hence G_(D_1)=G_R/2.  Write g_1 for
    # G_(D_1); the transported Green-matrix vector causal propagator is c_R=2g_1.
    g_2, g_1 = symbols("G_2 G_1", commutative=False)
    green_wz = Matrix(
        [
            [g_2, 2 * k * g_1, 0, 0],
            [g_1 * d, 0, 0, 0],
            [0, 0, 0, g_1],
            [0, 0, g_1, 0],
        ]
    )
    green_q1 = s_cfg.inv() * green_wz * s_row.inv()
    c_r = 2 * g_1
    expected_green_q1 = Matrix(
        [
            [-g_2, k * c_r, 0, 0],
            [c_r * d, 0, 0, 0],
            [0, 0, 0, c_r],
            [0, 0, -c_r, 0],
        ]
    )
    check_exact_q1_half_green_matrix = green_q1 == expected_green_q1

    ghost_numbers = (0, 0, 1, -1)

    def graded_difference(plus: Matrix, minus: Matrix) -> Matrix:
        return Matrix(
            4,
            4,
            lambda i, j: simplify(
                plus[i, j]
                - (
                    -1
                    if (ghost_numbers[i] * ghost_numbers[j]) % 2
                    else 1
                )
                * minus[i, j]
            ),
        )

    graded_q1_difference = graded_difference(lambda_q1_plus, lambda_q1_minus)
    reduced_graded_q1_difference = graded_q1_difference.subs(
        {
            lambda_2_plus: lambda_2_minus + I * g_2,
            lambda_1_plus: lambda_1_minus + I * g_1,
        }
    ).applyfunc(simplify)
    check_graded_causal_difference = (
        reduced_graded_q1_difference == I * expected_green_q1
    )

    gamma_wz = Matrix(
        [
            [0, 0, 2 * k, 0],
            [0, 0, 0, 0],
            [0, 0, 0, 0],
            [0, 1, 0, 0],
        ]
    )
    gamma_source = Matrix(
        [
            [0, 0, k, 0],
            [0, 0, 0, 0],
            [0, 0, 0, 0],
            [0, 1, 0, 0],
        ]
    )
    check_gamma_map = s_cfg.inv() * gamma_wz * s_cfg == gamma_source
    check_gamma_nilpotent = gamma_source * gamma_source == Matrix.zeros(4)
    check_compact_minus_one_contractible = (
        gamma_source[:, 3] == Matrix.zeros(4, 1)
        and gamma_source[3, 1] == 1
    )

    source_lambda = Matrix(
        [
            [-p, k, 0, 0],
            [d, lambda_aux, 0, 0],
            [0, 0, 0, -r],
            [0, 0, r, 0],
        ]
    )
    gamma_star_source = Matrix(
        [
            [0, 0, 0, 0],
            [0, 0, 0, 1],
            [-d, 0, 0, 0],
            [0, 0, 0, 0],
        ]
    )
    gamma_star_wz = Matrix(
        [
            [0, 0, 0, 0],
            [0, 0, 0, 1],
            [d, 0, 0, 0],
            [0, 0, 0, 0],
        ]
    )
    check_gamma_star_row_transport = (
        s_row * gamma_star_wz * s_row.inv() == gamma_star_source
    )
    brst_defect = gamma_star_source * source_lambda - source_lambda * gamma_source
    reduced_defect = brst_defect.applyfunc(
        lambda entry: reduce_source_relations(entry, p, k, d, r)
    )
    check_brst_intertwiner = reduced_defect == Matrix.zeros(4)

    checks = {
        "check_operator_row_field_map": check_operator_map,
        "check_pairing_duality_with_single_half": check_pairing_duality,
        "check_pairing_duality_invalid_no_half_rejected": (
            check_pairing_duality_invalid_no_half
        ),
        "check_covariance_plus_transport": check_covariance_plus_transport,
        "check_covariance_minus_transport": check_covariance_minus_transport,
        "check_covariance_form_half_pullback": check_covariance_form_half_pullback,
        "check_invalid_same_side_covariance_transport_rejected": (
            check_invalid_same_side_covariance_transport_rejected
        ),
        "check_exact_q1_half_green_matrix": check_exact_q1_half_green_matrix,
        "check_graded_causal_difference": check_graded_causal_difference,
        "check_gamma_configuration_map": check_gamma_map,
        "check_gamma_star_row_transport": check_gamma_star_row_transport,
        "check_gamma_nilpotent": check_gamma_nilpotent,
        "check_gamma_star_L_equals_L_gamma_all_lambda": check_brst_intertwiner,
        "check_compact_ghost_minus_one_pair_contractible": (
            check_compact_minus_one_contractible
        ),
    }
    overall = all(checks.values())

    lines = [
        "the boundary programme local benchmark Ward transport BRST-lift symbolic audit",
        "scope=exact finite pairing/normalization maps, covariance blocks, causal matrix, and free BRST identities",
        "nonclaims=no host Hadamard construction; no Fredholm proof; no BRST faithfulness; no Wick algebra",
        "wz_field_order=(u,b,c,cbar)",
        "source_reordered_field_order=(Theta,beta,eta,bar_eta)",
        "S_cfg=diag(2,-1,1,-1)",
        "S_row=diag(-1/2,1/2,1/2,1/2)",
        "pairing_identity=B_Q*S_row=(1/2)*S_cfg^T*B_WZ",
        "covariance_transport=lambda_Q=S_cfg^(-1)*lambda_WZ*S_row^(-1)",
        "form_transport=Lambda_Q=(1/2)*(S_row^(-1))^*Lambda_WZ*S_row^(-1)",
        "rho_plusminus=2*lambda_1_plusminus",
        "source_relations_used=P*K=0,D*P=0,D*K=R",
        "operator_placeholders=noncommuting",
        "causal_difference=lambda_plus-(-1)^gh*lambda_minus=i*G",
    ]
    lines.extend(
        f"{name}={'PASS' if value else 'FAIL'}" for name, value in checks.items()
    )
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
            "ward_transport_output.txt"
        )
        target.write_text(output, encoding="utf-8")
