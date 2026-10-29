# Nested refinement flow verdict: the carrier survives, the bare continuum claim does not

## Abstract

A28 found a causal two-polarization Regge oscillator and preserved it through one
longitudinal edge subdivision.  That was an existence result, not a continuum result.  Here
the same Lorentzian region is refined by uniformly subdividing a time-reflected pair of cap
edges into `m=1,2,3,4` segments, producing `C_(2m+3)*C6` cone complexes with 30, 42, 54 and
66 four-simplices.  Exact transverse dihedral symmetry reduces the calculation to physical
central-radial and spin-two boundary directions before differentiating the action.

At every level the cap-to-cap mixed spin-two Hessian has numerical rank two, its two
polarizations remain degenerate, and the reduced canonical spectrum stays on the unit
circle.  The lattice spin-two carrier therefore survives two nested refinements.  However,
its invariant oscillation phase flows

```text
0.815035, 0.832011, 0.634608, 0.388072,
```

while the overlap with the injected coarse mode falls to `0.93110`.  The last result is
stable under a 50 percent finite-difference step change.  A field rescaling or endpoint
canonical redefinition cannot remove this phase flow.  Thus the bare Regge kernel fails its
first multi-level cylindrical-consistency test.  The finite lattice graviton remains a real
result; a macroscopic or quantum-gravity limit now requires a genuinely improved/perfect
quadratic action, not another reinterpretation of the same bare kernel.

## 1. Exact reduction before numerics

Let the two real spin-two samples at a longitudinal site be

```text
u_2(phi) = (cos(2 phi), sin(2 phi)).                            (1.1)
```

The regular transverse hexagon has `D6` symmetry.  The exact Mathematica commutant audit in
`scripts/verify_causal_spin2.wl` gives

```text
Comm_C6(u_2) = {alpha I + beta J},
Comm_D6(u_2) = {alpha I}.                                      (1.2)
```

Consequently spin two cannot mix with the other transverse irreducible sectors, and every
reflection-invariant plus/cross block is scalar.  It is therefore legitimate to restrict
the Regge Hessian to spin-two boundary variations before numerical differentiation.  This
removes no spin-two response; it removes exactly decoupled scalar, vector and lattice
harmonic sectors.

The `m`-segment complex has

```text
q-cycle vertices                    2m+3,
four-simplices                      6(2m+3),
relative squared lengths per cap    7m+1,
embedded cap tangent rank           4m,
physical cap modes                  3m+1,
spin-two samples per cap            2m,
physical central-radial modes       2m+5.                       (1.3)
```

The four central vertex-displacement directions are removed by the left null quotient of
the internal radial length tangent.  The physical internal Hessian is then eliminated by a
Schur complement.  This leaves spin-two self blocks `A_m,C_m` and mixed block `B_m`.

For `m=2`, the symmetry-reduced implementation independently reproduces the unreduced A28
calculation:

```text
rank(B_2) = 2,
coarse overlap = 0.9999028,
theta_2 = 0.8320106.                                           (1.4)
```

The agreement establishes that the reduction is not manufacturing the later flow.

## 2. Nested causal geometries

The inserted vertices lie exactly on the original time-reflected cap-edge pair.  Thus
`m=4` contains the `m=2` midpoint vertices and is a genuine nested refinement.  Every
simplex remains Lorentzian and the boundary normals retain the decomposition

```text
past / future / side tetrahedra = 6(m+1) / 6(m+1) / 6.         (2.1)
```

The positive simplex Gram margin stays `0.5`.  The negative Gram magnitude decreases as
the subdivided simplices become thinner:

| `m` | q cycle | simplices | minimum negative Gram magnitude |
|---:|---:|---:|---:|
| 1 | 5 | 30 | 0.132328 |
| 2 | 7 | 42 | 0.0252817 |
| 3 | 9 | 54 | 0.0100139 |
| 4 | 11 | 66 | 0.00530682 |

No value is zero.  Conditioning is therefore monitored explicitly rather than hidden.

## 3. Transfer and phase flow

The numerical certificate is:

| `m` | leading `B_m` singular value | largest detail singular value | coarse overlap | `theta_m` |
|---:|---:|---:|---:|---:|
| 1 | 0.840003 | -- | 1.000000 | 0.815035 |
| 2 | 0.672133 | 3.4e-8 | 0.999903 | 0.832011 |
| 3 | 0.573840 | 6.6e-8 | 0.970068 | 0.634608 |
| 4 | 0.516755 | 8.6e-8 | 0.931103 | 0.388070 |

The two leading singular values agree within the finite-difference error at each level.
All additional longitudinal detail directions are pre/post constrained: they do not form
new propagating bands.  The surviving canonical maps have symplectic residuals below
`8e-16`; their unit-circle defects are below `8e-8`.

The `m=4` physical internal Hessian has thirteen resolved singular values, the smallest
being `0.32235`.  Repeating the complete Richardson Hessian at steps `1.2e-4` and `1.8e-4`
gives

```text
theta_4 = (0.38807015, 0.38807475)  at 1.2e-4,
theta_4 = (0.38806967, 0.38806974)  at 1.8e-4,                 (3.1)
```

with identical rank, coarse overlap `0.9311026`, and internal spectrum to the displayed
precision.  The phase flow is not a differencing artifact caused by the thinnest simplex.

## 4. Why normalization cannot rescue the flow

Time reflection and `D6` reduce the surviving quadratic principal function to

```text
S_m(q_i,q_f) = 1/2 a_m (q_i^2+q_f^2) + b_m q_i.q_f.           (4.1)
```

For either polarization the canonical characteristic polynomial is

```text
lambda^2 + 2(a_m/b_m) lambda + 1,
cos(theta_m) = -a_m/b_m.                                      (4.2)
```

A boundary wave-function renormalization `q -> Z_m q` multiplies both `a_m` and `b_m` by
`Z_m^2`; it cannot change (4.2).  A canonical endpoint redefinition conjugates the transfer
matrix and likewise cannot change its eigenvalues.  Numerically,

```text
cos(theta_m) = 0.685843, 0.673391, 0.805304, 0.925640.         (4.3)
```

Therefore the observed flow is a genuine change of the effective dynamics.  It cannot be
renamed away as a basis choice, polarization convention, or field normalization.

## 5. Comparison with established Regge results

Bahr and Dittrich construct improved and perfect Regge actions precisely to recover
continuum dynamics and gauge symmetry on a fixed triangulation:

https://arxiv.org/abs/0907.4323

Their analysis of broken and restored diffeomorphism symmetry explains why a naive discrete
action should not be expected to be cylindrically consistent under refinement:

https://arxiv.org/abs/0909.5688

Dittrich and Steinhaus show that four-dimensional linearized Regge calculus has genuine
triangulation dependence and that a triangulation-independent measure is substantially more
obstructed than in three dimensions:

https://arxiv.org/abs/1110.6866

The present calculation does not supersede those general results.  Its contribution is a
specific Lorentzian, causal, spin-resolved refinement flow in which the physical carrier,
pre/post constraints, canonical phase and conditioning are measured together.

## 6. Verdict

### Proved within the declared algebra

- `D6` symmetry permits exact restriction to the spin-two Hessian sector.
- The canonical phase is invariant under scalar boundary field renormalization and
  canonical conjugation.

### Numerically certified

- A rank-two, twofold-degenerate, stable spin-two carrier survives through `m=4`.
- All added longitudinal detail samples remain pre/post constrained rather than becoming
  propagating modes.
- The phase and coarse embedding fail to stabilize over `m=1,2,3,4`.
- The `m=4` flow is stable under the stated step-size variation.

### Death sentence, with exact scope

The **bare Regge realization of PCRG as a cylindrically consistent macroscopic dynamics is
rejected**.  Its phase-valued quantum refinement framework cannot be built on an unmodified
kernel whose classical spin-two eigenphase changes this strongly under exact nested
subdivision.

This does **not** erase the finite causal lattice graviton.  It changes its status from a
candidate continuum carrier to input data for a renormalization problem.  Road A remains
alive only if its next dynamics is an improved/perfect action whose coarse-grained spin-two
principal function has a refinement fixed point.  That is a new construction, not an
unfinished proof of the bare model.

## 7. Next admissible move

Derive a quadratic perfect kernel by integrating the fine spin-two and physical radial
variables while allowing the minimal symmetry-compatible counterterm basis.  Fix the
counterterms by demanding equality of the invariant transfer trace, not merely equality of
determinants or ranks.  Reject the road if no local or exponentially decaying counterterm
basis reaches a fixed point without destroying the Regge edge equations or causal
Lorentzian domain.  Only after that classical fixed point exists is a measure, ghost,
Cartan-connection or quantum refinement calculation meaningful.

## 8. Reproduction

```powershell
python scripts\multilevel_causal_spin2.py
& 'D:\Wolfram Research\Mathematica\14.0\wolframscript.exe' -file scripts\verify_causal_spin2.wl
```

## 9. Dated follow-up: `m=5` changes the stability class

The statements above remain the certificate for `m=1,2,3,4`, but they are not the end of
the refinement trajectory.  On 2026-07-15 the complete calculation was extended to
`m=5`.  The mixed spin-two block remains rank two and the reduced canonical relation
remains symplectic, yet its reciprocal eigenvalues move off the unit circle to
approximately `1.29025` and `0.77504`, duplicated for plus/cross.  Two Richardson base
steps, `1.8e-4` and `1.2e-4`, agree on this classification.

Consequently “stable carrier” must be read as a finite statement only through `m=4`.
The new result strengthens, rather than weakens, this document's rejection of the bare
continuum dynamics: between `m=4` and `m=5` the tested transfer becomes hyperbolic.  Full
data, conditioning checks and the exact reduced classification are recorded in
`A30-FIFTH-REFINEMENT-HYPERBOLIC-TRANSITION.md`.
