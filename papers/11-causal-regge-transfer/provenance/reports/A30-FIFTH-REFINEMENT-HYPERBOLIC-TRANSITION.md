# Fifth refinement verdict: the carrier becomes hyperbolic

Date: 2026-07-15 (Beijing / UTC+8)

Status: numerical finite-cell result with an exact reduced classification; no continuum or
quantum-gravity claim.

## Abstract

A29 followed the PCRG transverse spin-two sector through longitudinal subdivision levels
`m=1,2,3,4`.  The mixed Hessian retained rank two, but its canonical phase drifted strongly.
The follow-up is the complete `m=5` Hessian on the same nested Lorentzian family,
now containing `C13*C6` and 78 four-simplices.  Two independent Richardson base steps give
the same qualitative result.  The transmitted spin-two singular pair remains present and
degenerate, while the canonical eigenvalues leave the unit circle and become a reciprocal
real pair near `1.29025` and `0.77504`, duplicated for the two polarizations.  Symplecticity
is preserved to about `2e-16`.

Thus the bare Regge kernel does not merely fail to approach a stable phase.  On this tested
family it crosses from an elliptic finite oscillator to a hyperbolic transfer between
`m=4` and `m=5`.  The finite spin-two representation survives; stable graviton dynamics
does not on this trajectory.

## 1. Scope and method

The geometry, symmetry reduction and complete Regge-action Hessian are those defined in
`scripts/multilevel_causal_spin2.py` and A29.  The exact transverse `D6` symmetry permits
restriction to the real spin-two pair.  Four internal vertex-displacement tangents are
removed, the resolved physical internal block is eliminated by a Schur complement, and
the surviving time-reflection-symmetric principal function is

```text
S_m(q_i,q_f) = 1/2 a_m (q_i^2+q_f^2) + b_m q_i.q_f.          (1.1)
```

For either polarization its canonical characteristic equation is

```text
lambda^2 - 2 chi_m lambda + 1 = 0,
chi_m = -a_m/b_m.                                             (1.2)
```

Therefore `|chi_m|<1` is elliptic, `|chi_m|=1` parabolic and `|chi_m|>1`
hyperbolic.  This classification is invariant under scalar field normalization and
canonical conjugation.

## 2. Independent finite-difference runs

The complete calculation was repeated at Richardson base steps `1.8e-4` and `1.2e-4`.
The reported values are:

| base step | leading mixed singulars | largest detail singular | coarse overlap | unit-circle residual | symplectic residual |
|---:|---:|---:|---:|---:|---:|
| `1.8e-4` | `0.479963944, 0.479963804` | `1.01e-7` | `0.8919179` | `0.2902555` | `2.22e-16` |
| `1.2e-4` | `0.479964367, 0.479964196` | `1.85e-7` | `0.8919179` | `0.2902481` | `1.76e-16` |

At the finer step the four canonical eigenvalues are

```text
1.29024266 - 1.16e-6 i,
1.29024809 - 0.84e-6 i,
0.77504801 + 0.70e-6 i,
0.77504474 + 0.51e-6 i.                                      (2.1)
```

The small imaginary parts set the numerical differentiation scale; the departure of the
moduli from one is about `0.29025`, six orders of magnitude larger.  The physical internal
Hessian is resolved: at the finer step its smallest retained singular value is
`0.31048966`.  The thinnest simplex is also nondegenerate, with minimum negative Gram
eigenvalue magnitude `0.00327578` and minimum positive one `0.413006`.  The transition is
not attributable to an unresolved Schur denominator or a null simplex at the displayed
precision.

There is an important geometric limitation.  This refinement subdivides selected cap
edges without co-refining the central bulk vertex, so the minimum negative Gram magnitude
falls from `0.132328` at `m=1` to `0.00327578` at `m=5`.  A five-point descriptive fit gives
approximately `m^-2.30`.  The sequence is therefore becoming anisotropic and tending
toward sliver simplices even though every tested member is nondegenerate.  The present data
do not separate a genuine bare-Regge renormalization instability from a pathology of this
non-shape-regular refinement family.  A shape-regular bulk-and-boundary refinement is an
independent test, not a reinterpretation of these numbers.

## 3. What is exact and what is numerical

Once `a_m,b_m` are given, equation (1.2), reciprocal eigenvalue pairing and the
elliptic/hyperbolic classification are exact algebra.  The `m=5` coefficients themselves
are numerical derivatives of the complete finite Regge action.  Agreement of the two step
sizes is evidence of robustness, not an analytic theorem about all refinements.

`scripts/verify_fifth_refinement_transition.wl` independently checks the reduced
classification from rounded imported coefficients.  It is deliberately labelled a
consistency audit; it does not masquerade as a second Hessian computation.

## 4. Relation to prior Regge results

Improved and perfect actions were introduced to restore continuum dynamics and gauge
symmetry under coarse graining:

https://arxiv.org/abs/0907.4323

Four-dimensional linearized Regge calculus is known to retain triangulation-dependent
one-loop structure:

https://arxiv.org/abs/1110.6866

The present result is narrower.  It supplies one explicit causal Lorentzian,
spin-resolved family in which the bare transfer changes stability class while the carrier
rank and symplectic relation persist.  No claim of general Regge instability or literature
priority beyond this concrete computation is made.

## 5. Verdict

The following finite statement is numerically certified on the declared family:

> At `m=5`, the twofold-degenerate PCRG spin-two transfer remains rank two and symplectic,
> but is hyperbolic rather than elliptic.

This strengthens A29's rejection of the bare Regge kernel on the declared nested family as
a cylindrically consistent macroscopic dynamics.  It does not reject all shape-regular
Regge refinement schemes.  It also narrows the surviving asset: PCRG has an exact spin-two
representation and a finite causal carrier through several cells, but not a stable
refinement trajectory.  Calling that carrier a physical continuum graviton is no longer
defensible without a new improved/perfect action that independently reaches a fixed point.

## 6. Reproduction boundary

The two complete Hessian runs were executed with the existing multilevel code and the
stated base steps.  The lightweight classification audit is:

```powershell
& 'D:\Wolfram Research\Mathematica\14.0\wolframscript.exe' `
  -file scripts\verify_fifth_refinement_transition.wl
```

Because the full `m=5` calculation is expensive, the rounded numerical certificate is
recorded here rather than silently recomputed during every documentation check.
