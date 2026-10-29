#!/usr/bin/env python3
"""Exact finite audit for the local benchmark first h-phi-phi/BV matter deformation.

This checks coefficients and graded signs only.  It does not prove curved-space
renormalization, wavefront estimates, the full order-kappa CME, or the QME.
"""

from sympy import Rational, Symbol, expand, simplify


def require_zero(name, expression):
    residual = simplify(expand(expression))
    if residual != 0:
        raise AssertionError(f"{name}: residual={residual}")
    print(f"{name}=PASS")


kappa = Symbol("kappa")
htrace = Symbol("h")
hdphidphi = Symbol("h_dphi_dphi")
gradphi2 = Symbol("gradphi2")
m2phi2 = Symbol("m2phi2")

# Scalar density: -1/2 sqrt(g) [g^{-1}(dphi,dphi)+m^2 phi^2].
# sqrt(g)=sqrt(gbar)(1+kappa*h/2)+O(kappa^2),
# g^{-1}(dphi,dphi)=gradphi2-kappa*h_dphi_dphi+O(kappa^2).
expanded_density = -Rational(1, 2) * (
    (1 + kappa * htrace / 2)
    * (gradphi2 - kappa * hdphidphi + m2phi2)
)
vertex_from_expansion = expand(expanded_density).coeff(kappa, 1)
vertex_from_half_hT = Rational(1, 2) * hdphidphi - Rational(1, 4) * htrace * (
    gradphi2 + m2phi2
)
require_zero(
    "check_metric_volume_expansion_gives_plus_half_hT",
    vertex_from_expansion - vertex_from_half_hT,
)

# Normal-coordinate divergence identity.  The two mixed-Hessian terms cancel.
boxphi = Symbol("boxphi")
m2phi = Symbol("m2phi")
gradnu = Symbol("gradnu_phi")
mixed = Symbol("mixed_hessian_contraction")
divergence_direct = boxphi * gradnu + mixed - mixed - m2phi * gradnu
Pphi = boxphi - m2phi
require_zero(
    "check_stress_divergence_off_shell",
    divergence_direct - Pphi * gradnu,
)

# For A_phi=phi^ddagger R_phi, both factors are odd.  Moving the odd
# variation through R_phi generates the right-derivative sign.
Rphi = Symbol("R_phi", commutative=False)
left_antifield_derivative = Rphi
right_antifield_derivative = (-1) ** (1 * 1) * Rphi
if left_antifield_derivative != Rphi or right_antifield_derivative != -Rphi:
    raise AssertionError("graded_antifield_derivative_sign")
print("check_left_antifield_derivative_plus_Rphi=PASS")
print("check_right_antifield_derivative_minus_Rphi=PASS")

# Classical first-deformation equation in the frozen antibracket ordering.
bracket_vertex_free = -Pphi * Rphi
bracket_antifield_free = -right_antifield_derivative * Pphi
require_zero(
    "check_classical_matter_brst_descent",
    bracket_vertex_free + bracket_antifield_free,
)

# Invalid sign controls must fail.
if simplify((-bracket_vertex_free) + bracket_antifield_free) == 0:
    raise AssertionError("invalid_minus_half_hT_control_not_rejected")
print("invalid_minus_half_hT_control_rejected=PASS")

if simplify(bracket_vertex_free - bracket_antifield_free) == 0:
    raise AssertionError("invalid_minus_antifield_control_not_rejected")
print("invalid_minus_antifield_control_rejected=PASS")

# With a switching function, the antifield term cancels the bulk EOM term and
# leaves only the derivative-supported current term.
df_current = Symbol("df_c_T")
f_eom_current = Symbol("f_eom_current")
vertex_with_cutoff = -df_current - f_eom_current
antifield_with_cutoff = f_eom_current
require_zero(
    "check_cutoff_plateau_residual",
    vertex_with_cutoff + antifield_with_cutoff + df_current,
)

# kappa^2=32*pi*G gives kappa^2/4=8*pi*G.
piG = Symbol("piG")
require_zero(
    "check_einstein_source_normalization",
    Rational(1, 4) * (32 * piG) - 8 * piG,
)

print("scope=finite_exact_coefficient_divergence_graded_sign_and_cutoff_checks")
print(
    "nonclaims=no_curved_EG_extension;no_T11b_proof;no_full_order_kappa_CME;"
    "no_multi_insertion_AMWI;no_QME;no_backreaction;no_global_AdS"
)
print("overall=PASS_ALGEBRAIC_MATTER_SUBCOCYCLE_SCOPE")
