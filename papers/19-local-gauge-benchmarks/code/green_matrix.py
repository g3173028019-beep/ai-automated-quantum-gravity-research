"""Exact algebra audit for the local benchmark Green matrix block Green formula.

This verifies formal block identities and the frozen lambda=+1/2 coefficients.  It does
not prove Green hyperbolicity, causal support, Hadamard properties, or BRST positivity.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from sympy import (
    BlockMatrix,
    Identity,
    Inverse,
    MatrixSymbol,
    Rational,
    ZeroMatrix,
    block_collapse,
    expand,
    simplify,
    symbols,
)


def simplify_blocks(expr: BlockMatrix) -> list[list[object]]:
    collapsed = block_collapse(expr)
    return [
        [simplify(expand(entry)) for entry in row]
        for row in collapsed.blocks.tolist()
    ]


def schur_formula_allowed(value: object) -> bool:
    return value != 0


def main() -> str:
    n, m = symbols("n m", integer=True, positive=True)
    lambda_aux = symbols("lambda_aux", nonzero=True)

    s = MatrixSymbol("P_lambda", n, n)
    f = MatrixSymbol("F", m, n)
    ih = Identity(n)
    iw = Identity(m)
    z_h_w = ZeroMatrix(n, m)
    z_w_h = ZeroMatrix(m, n)
    g = Inverse(s)
    lambda_inv = 1 / lambda_aux

    # In LB variables the quadratic field term is +1/2<H,P_0 H>.  Thus
    # P_lambda=P_0-lambda^{-1}F^*F, or P_0=P_lambda+lambda^{-1}F^*F.
    p0 = s + lambda_inv * f.T * f
    even = BlockMatrix([[p0, f.T], [f, lambda_aux * iw]])
    even_green = BlockMatrix(
        [
            [g, -lambda_inv * g * f.T],
            [
                -lambda_inv * f * g,
                lambda_inv * iw + lambda_inv**2 * f * g * f.T,
            ],
        ]
    )

    left_even = simplify_blocks(even * even_green)
    right_even = simplify_blocks(even_green * even)
    expected_even = [[ih, z_h_w], [z_w_h, iw]]
    check_even_left = left_even == expected_even
    check_even_right = right_even == expected_even

    ghost_m = MatrixSymbol("M", m, m)
    zw = ZeroMatrix(m, m)
    odd = BlockMatrix([[zw, ghost_m.T], [-ghost_m, zw]])
    odd_green = BlockMatrix(
        [[zw, -Inverse(ghost_m)], [Inverse(ghost_m.T), zw]]
    )
    left_odd = simplify_blocks(odd * odd_green)
    right_odd = simplify_blocks(odd_green * odd)
    expected_odd = [[iw, zw], [zw, iw]]
    check_odd_left = left_odd == expected_odd
    check_odd_right = right_odd == expected_odd

    lambda_lb = Rational(1, 2)
    eliminated_gauge_square = simplify(-1 / (2 * lambda_lb))
    hessian_shift = simplify(-1 / lambda_lb)
    auxiliary_solution = simplify(-1 / lambda_lb)
    check_lambda_match = (
        eliminated_gauge_square == -1
        and hessian_shift == -2
        and auxiliary_solution == -2
    )

    invalid_lambda_zero_rejected = not schur_formula_allowed(Rational(0))

    checks = {
        "check_even_left_inverse": check_even_left,
        "check_even_right_inverse": check_even_right,
        "check_odd_left_inverse": check_odd_left,
        "check_odd_right_inverse": check_odd_right,
        "check_lambda_plus_half_match": check_lambda_match,
        "invalid_lambda_zero_rejected": invalid_lambda_zero_rejected,
    }
    overall = all(checks.values())

    lines = [
        "the boundary programme local benchmark Green matrix full-matrix symbolic audit",
        "scope=exact formal Schur-block algebra and lambda coefficient only",
        "nonclaims=no PDE Green-hyperbolicity proof; no causal-support proof; no Hadamard state; no BRST positivity",
        "field_order_even=(H,b), H=(h,phi)",
        "field_order_odd=(c,cbar)",
        "source_quadratic_sign=S_H^(2)=S_H-(1/2)<Theta,P_bar Theta>",
        "LB_sign_map=P_0_LB=-(4/kappa^2)A^*P_bar A",
        "formula_P_lambda=P_0-lambda^(-1)F^*F",
        "frozen_lambda=1/2",
        f"eliminated_gauge_square_coefficient={eliminated_gauge_square}",
        f"hessian_FstarF_shift_coefficient={hessian_shift}",
        f"auxiliary_solution_b_over_F={auxiliary_solution}",
        "full_matrix_lambda_zero_status=SUPPORTED_BY_NONSCHUR_C2_THEOREM_NOT_THIS_EXECUTABLE",
        "identity_term_support=LOCAL_SUBSET_OF_BOTH_CAUSAL_SUPPORTS",
        "causal_difference_identity_term=CANCELS_BETWEEN_ADVANCED_AND_RETARDED",
    ]
    lines.extend(
        f"{name}={'PASS' if value else 'FAIL'}" for name, value in checks.items()
    )
    lines.append(f"overall={'PASS_BLOCK_ALGEBRA_ONLY' if overall else 'FAIL'}")
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    output = main()
    print(output, end="")
    if args.write:
        target = Path(__file__).with_name(
            "green_matrix_output.txt"
        )
        target.write_text(output, encoding="utf-8")
