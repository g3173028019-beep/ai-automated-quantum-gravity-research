# A28. Causal spin-two verdict: the Ricci-nullity claim dies and the lattice graviton survives

## Abstract

A27 found a two-dimensional kernel after contracting the on-shell hinge curvatures of a
Lorentzian `C5*C6` join-star with two inverse-dual-area Einstein diagnostics.  That result
was numerically robust but physically misidentified.  The circumcentric hybrid-volume
derivation cancels the dual area and weights a hinge by its signed area portion at the
interior vertex.  With the corresponding Lorentzian normal-plane signature, the vertex
Einstein screen has rank one and nullity six, not nullity two.  Moreover all seven boundary
curvature modes already solve the true internal Regge equations.  The dimension of a finite
Dirichlet solution space was therefore being confused with the number of local helicities.

The same correction exposes a different and stronger result.  The star boundary has a
nondegenerate causal decomposition into past spacelike, future spacelike and timelike-side
tetrahedra.  Exact transverse `C6` symmetry decomposes its four cap-to-cap modes into a
two-dimensional spin-two representation and two nonradiative representations.  The
Hamilton principal Hessian is block diagonal on this decomposition.  On a time-reflection
symmetric witness the spin-two block is a stable symplectic oscillator with two exactly
degenerate polarizations.  A longitudinal `C7*C6` subdivision preserves the same two-mode
transmitted subspace with more than `0.9999` overlap; the added spin-two detail modes become
pre/post constrained.  The old Ricci-nullity argument is dead, while a more physical finite
lattice-graviton carrier survives.  Its finite-spacing electric/magnetic impedance is not
yet lightlike, so this is not a continuum graviton or a quantum-gravity completion.

## 1. Fatal correction to A27

For a hinge `h` incident on the sole interior vertex `v`, circumcentric Regge curvature and
hybrid volume have the schematic forms

```text
K_h proportional to epsilon_h / A_h^*,
Delta V_hv proportional to A_hv A_h^*.
```

Consequently the integrated vertex curvature obeys

```text
R_v Delta V_v proportional to
    sum_h epsilon_h A_hv sigma(N_h) U_h tensor U_h.          (1.1)
```

Here `A_hv` is the signed circumcentric portion of the hinge assigned to `v`, `U_h` is a
unit normal bivector, and `sigma(N_h)` is the signature of its normal plane.  The factor
`sigma(N_h)` is required because `U_ab U_cd` alone has opposite scalar contraction on
Euclidean and Lorentzian normal planes.  Equation (1.1) follows independently from the
hybrid-volume cancellation and from demanding a signature-independent scalar-curvature
normalization.

On the A27 witness the corrected screen singular values are

```text
(5.01164218e-1,
 1.66804230e-7, 1.19853805e-7, 8.40598283e-8,
 2.33174475e-8, 7.77329602e-9, 3.47464095e-24).             (1.2)
```

Thus its resolved rank is one and its nullity is six.  The 50-degree transverse deformation
again gives rank one and nullity six.  This is incompatible with A27's claim that two
inverse-dual screens isolated a physically preferred two-dimensional vacuum sector.

The actual flat-background Regge equation is sharper and does not require a pointwise
tensor reconstruction:

```text
delta E_e = sum_(h contains e) (partial A_h/partial s_e) delta epsilon_h = 0.  (1.3)
```

The eleven-by-forty-one area-moment matrix in (1.3) has eleven resolved singular values,
the smallest being `0.286114`.  After the seven physical internal radial variables are
eliminated, all seven boundary curvature modes satisfy all eleven equations with maximum
residual `1.16e-7`.  They are finite-region vacuum Regge solutions.  Counting the dimension
of this Dirichlet solution space cannot count continuum helicities.

The scientific verdict on A27 is therefore negative and specific:

- the numerical nullities reported there are reproducible;
- the two inverse-dual diagnostics were not the integrated vertex Ricci prescription;
- their common nullity two is not evidence for two helicities;
- the underlying Regge dynamics is not falsified by this correction.

## 2. The boundary was causal all along

Each boundary tetrahedron is the join of one `q`-cycle edge with one `r`-cycle edge.  Its
outward normal gives an invariant classification.  For the A27 embedding the thirty
tetrahedra split exactly as

```text
past-spacelike q strips       12,
future-spacelike q strips     12,
timelike-side strip            6.                              (2.1)
```

No normal is close to null.  Fixing the endpoint cones and timelike strip leaves eight
relative edge lengths on each spacelike cap.  The one-vertex embedded tangent map has rank
four on either cap, hence four relative configuration modes remain.  Eliminating the seven
physical central radial modes gives cap-to-cap singular values

```text
(1.90456128, 1.29995110, 0.89152243, 0.59434833),            (2.2)
```

so the complete transfer has rank four.  The earlier statement that the closed `S^3`
boundary precluded a causal interpretation was too strong: a total cobordism boundary is
closed while still admitting past, future and side pieces.

## 3. The physical two-dimensional sector comes from spin, not Ricci rank

Restore a regular transverse hexagon.  The relative cap data then carries the real
irreducible representations `m=0,1,2,3` of `C6`.  The pair

```text
cos(2 phi), sin(2 phi)                                          (3.1)
```

is the two-dimensional spin-two representation.  The Regge action is `C6` invariant, so
its Hessian commutes with this representation.  A real symmetric self-block on this
irreducible pair is proportional to the identity.  A mixed block commuting with `C6` can
more generally be `alpha I + beta J`, where `J` is the quarter-turn complex structure; the
transverse reflection of the regular witness sets `beta=0`.  Therefore plus/cross
degeneracy and decoupling from the scalar, vector and highest-lattice-harmonic sectors
follow from a finite-group `C6`/`D6` commutant calculation; they are not inferred from a
singular-value coincidence.

On the exact `C6` witness the raw harmonics project into the cap normal quotient with
residual below `4.7e-16`.  Spin-to-complement transfer leakage is below `1.0e-7`, and

```text
H_if^(2) = 0.68629428 I_2 + O(1e-7).                           (3.2)
```

The original asymmetric longitudinal geometry gives a symplectic squeeze.  This is not a
stationary wave test because the two caps are inequivalent.  A second witness was therefore
chosen before evaluating the action, using only causal and simplex-conditioning margins.
It is invariant under time reflection `q2 <-> q5`, `q3 <-> q4` and under transverse `C6`.
All thirty four-simplices are nondegenerate Lorentzian simplices; the smallest magnitudes
of the negative and positive simplex Gram eigenvalues are `0.1323` and `0.5`.

For this witness,

```text
H_ii^(2) = -0.57611043 I_2 + O(7e-7),
H_if^(2) =  0.84000285 I_2 + O(2e-7),
H_ff^(2) = -0.57611038 I_2 + O(5e-7).                          (3.3)
```

The Hamilton principal function gives a four-dimensional canonical map with determinant
one and symplectic residual below `6e-17`.  Its eigenvalues are

```text
lambda = exp(+/- i theta),
theta = 0.81503467,
```

each twice degenerate.  The maximum unit-circle defect is `9.4e-8`.  Repeating the complete
Hessian at steps `8e-5` and `1.6e-4` changes `theta` by only `7.1e-7`.  A separate symmetric
geometry deformation retains unit-circle evolution and spin-two decoupling while changing
the phase, as physical background dependence requires.

An independent exact Mathematica audit solves the commutant equations, verifies the
`C6` and `D6` statements above, and proves for scalar two-polarization blocks

```text
det(lambda I - M) = [lambda^2 + ((a+c)/b) lambda + 1]^2,       (3.4)
M^T Omega M = Omega.
```

Substitution of (3.3) gives `theta=0.81503435` using only the displayed finite-difference
digits.  The small difference from the full-precision value is therefore a rounding effect,
not an independent discrepancy.

This closes a finite-cell statement: PCRG has a causal, stable, two-polarization lattice
spin-two oscillator selected by spacetime dynamics and transverse rotations.

## 4. Curvature and the remaining wave defect

The positive-frequency multiplier is

```text
lambda_+ = 0.68584313 + 0.72774934 i.                          (4.1)
```

Combining initial and final spin-two boundary data with (4.1), solving the internal radial
equations, and reconstructing (1.1) gives a relative vertex Einstein residual
`2.85e-7`.  The electric tensors are transverse and traceless to the same accuracy and
have the plus/cross templates.  Their magnetic partners obey a single fitted discrete
impedance relation

```text
B_plus  = Z E_cross,
B_cross = -Z E_plus,
Z = (-3.7e-9 + 0.80024446 i),                                 (4.2)
```

with normalized fitting residual `1.06e-6`.

Equation (4.2) is strong evidence that the modes form one spin-two dynamical pair, but it is
also the remaining obstruction.  The impedance has a quarter-phase convention and modulus
different from one.  Without a controlled refinement law mapping it to the continuum
`B=J E` relation, the mode is not yet a lightlike Regge gravitational wave.  The previous
21--29 percent A27 residual was not the right observable, but the corrected observable has
not magically become continuum exact.

## 5. First longitudinal refinement

Subdivide a time-reflected pair of cap edges.  The result is a `C7*C6` cone with forty-two
four-simplices and boundary classification

```text
past spacelike / future spacelike / timelike side = 18 / 18 / 6.  (5.1)
```

Each cap has fifteen relative lengths.  Removing the rank-eight two-vertex embedded tangent
space leaves seven modes, four of which are spin-two samples at two longitudinal sites.  The
nine physical internal radial Hessian singular values range from `10.0093` to `0.375508`.

The four-by-four spin-two mixed block has singular values

```text
(0.672132955, 0.672132849, 1.05e-7, 1.70e-9).                 (5.2)
```

Thus only one longitudinal combination for each polarization transfers.  The two added
detail modes are pre/post constrained, precisely the type of varying-phase-space behavior
expected in canonical Pachner evolution.  The surviving two-dimensional subspace overlaps
the exact midpoint injection of the coarse mode with singular values above `0.9999028` on
both caps.  Its reduced canonical eigenvalues are

```text
exp(+/- i 0.8320111),                                         (5.3)
```

twice degenerate, with symplectic residual `2.7e-15` and unit-circle defect `6.7e-9`.
The `0.017` phase shift from (3.3) is a finite-spacing/refinement effect, not disappearance
of the carrier.  Because survivor and constrained self-blocks still mix at order `0.40`,
this calculation is a pre/post constrained refinement, not a full two-band propagator.

## 6. Comparison with primary literature

Höhn's canonical linearized Regge analysis identifies vertex-displacement gauge, lattice
gravitons and the distinct creation/removal roles of Pachner moves.  The rank-two survivor
and null refinement-detail transfer in (5.2) fit that pre/post-constraint picture rather
than a fixed-dimensional canonical map:

https://arxiv.org/abs/1411.5672

The general varying-phase-space formalism of Dittrich and Höhn explains why singular mixed
Hessians do not license inversion on every boundary variable:

https://arxiv.org/abs/1108.1974

The hybrid-volume cancellation in (1.1) is the tensorial extension of the Voronoi--Delaunay
Regge construction, in which hinge curvature is deficit per dual area and curvature is
assigned through hybrid cells:

https://arxiv.org/abs/0804.0279

The edge-local, volume-weighted character of simplicial Ricci curvature is developed in:

https://arxiv.org/abs/1107.2458

PCRG's new contribution is not those established ingredients separately.  It is the
explicit Lorentzian join-star witness in which causal boundary classification, internal
Regge elimination, finite-group spin-two projection, stable canonical evolution and a
controlled subdivision certificate are calculated together.

## 7. Verdict ledger

### Proved within the declared finite algebra

- The Regge edge equation is the area-moment equation (1.3), not an arbitrary global
  inverse-dual Einstein screen.
- Exact `C6` invariance isolates the two-dimensional `m=2` representation and forces
  plus/cross degeneracy of symmetric Hessian blocks.
- A nonsingular mixed spin-two block defines a symplectic canonical map.

### Numerically certified

- The A27 inverse-dual two-mode interpretation fails under the derived vertex-hybrid screen.
- The `C5*C6` boundary has a robust past/future/side causal decomposition.
- The time-symmetric spin-two map is stable, symplectic and twofold degenerate.
- Its integrated vertex curvature is Einstein-null at differentiation accuracy.
- The `C7*C6` subdivision preserves the coarse two-mode survivor above `0.9999` overlap.

### Conditional

- The finite oscillator is called a **candidate lattice graviton** because its spin,
  causality and Regge dynamics are controlled, while its continuum dispersion is not.

### Open

- convergence of the impedance modulus and phase to the continuum lightlike value;
- a multi-layer refinement dispersion curve rather than one survivor phase;
- positive-frequency normalization and a retarded Green function;
- nonlinear stability and interaction of the spin-two sector;
- the independent Cartan connection/ghost/measure lift and quantum contour.

## 8. Reproduction

```powershell
python scripts\lorentzian_vertex_hybrid_audit.py
python scripts\causal_join_star_transfer.py
python scripts\refined_causal_spin2.py
& 'D:\Wolfram Research\Mathematica\14.0\wolframscript.exe' -file scripts\verify_causal_spin2.wl
python tests\test_two_cell_amplitude.py
```
