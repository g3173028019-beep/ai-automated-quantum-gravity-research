# Verification and evidence

[Home](README.md) · [Paper index](docs/PAPER-INDEX.md) · [Reproduction guide](docs/REPRODUCTION.md)

The original archive contains 20 full research packets. Later analytical manuscripts are documented separately below. The manuscripts state their arguments, calculation files expose supporting algebra or numerics, and local verification notes distinguish executed checks, mathematical inputs and retained evidence. The [historical verification record](docs/VERIFICATION-HISTORY.md) preserves the full delivery narrative.

## Evidence used in this archive

| Evidence | What to inspect | Scope |
|---|---|---|
| Analytical proof | Definitions, hypotheses, domains, proof steps and cited theorems | Establishes the stated conclusion under its assumptions |
| Exact algebra | Rational identities, symbolic operators, integer relations and finite field equations | Checks the encoded identities, without replacing an analytic domain or PDE argument |
| Certified estimate | Interval, residual and tail bounds, their inputs and error propagation | Distinguish newly generated enclosures from recombination of saved bounds |
| Numerical calculation | Model, discretization, parameters, convergence checks and observable | Supports the specified finite or approximated calculation |
| Retained certificate or input | Provenance and role in the downstream argument | Replaying a certificate differs from regenerating its producer's data |
| Packaging and layout | Source/PDF consistency, citations, bylines, metadata and rendered pages | Checks the artifact, without certifying its physical interpretation |

A successful finite check applies to the object it encodes. It does not establish an unencoded infinite-dimensional limit, Fredholm or wavefront theorem, interacting quantum master equation, or physical reconstruction.

## Evidence across the 20 papers

The [reproduction guide](docs/REPRODUCTION.md) gives every default command and its output location. Each linked packet supplies its local requirements and verification record.

| Paper | Principal evidence and boundary |
|---|---|
| [01 · Projection sign rigidity](papers/01-projection-sign-rigidity/README.md) | Analytic genericity and fixed-data holonomy criterion; portable exact revision check. Historical Wolfram subset enumeration is optional. |
| [02 · Half-wormholes and Wiener chaos](papers/02-half-wormholes-wiener-chaos/README.md) | Probability/combinatorial derivations and two portable exact checks; prior variance formulas are credited. Historical Monte Carlo is separate. |
| [03 · Energy-weighted thermal recovery](papers/03-jt-energy-weighted-recovery/README.md) | Infinite-input analytic recovery bound, supported by charge/Gaussian algebra, finite coherent examples and interval evaluations. |
| [04 · Heralded controlled-phase gate](papers/04-dbi-heralded-controlled-phase/README.md) | Complete-instrument and uniform conditional-fidelity arguments; the default replay recombines retained continuum residual partitions rather than reintegrating them. |
| [05 · Integer flux families](papers/05-integer-flux-families/README.md) | Exact assignments, F-term rank, canonical masses and pointwise lattice calculations; global geometric capacities remain conditional inputs. |
| [06 · Positive-energy local clocks](papers/06-maxwell-local-clocks/README.md) | Full-state clock comparison and resource/error estimates using the packaged accepted electromagnetic reference data. |
| [07 · Clock-foliated effective reference](papers/07-clock-foliated-reference/README.md) | Coupled harmonic constraints, positivity and finite-band Gaussian/evolution calculations for the specified modified effective action. |
| [08 · Einstein–scalar mixed quartic](papers/08-einstein-scalar-mixed-quartic/README.md) | Off-shell quartic/source identities and a separate exact finite-backreaction family, with the split-signature limitation explicit. |
| [09 · Canonical Einstein–scalar Hadamard states](papers/09-canonical-einstein-scalar-hadamard/README.md) | Analytic physical-state existence proof on the stated compact-Cauchy, critical-potential backgrounds, including singular Fredholm repair; the algebra script supports normalization and examples. |
| [10 · Controlled quantum-weight networks](papers/10-quantum-weight-networks/README.md) | Finite cell/network proofs and exact binary-symplectic checks; arbitrary-graph and collective-limit statements are not extrapolated from a finite table. |
| [11 · Finite causal transfer](papers/11-causal-regge-transfer/README.md) | Finite quotient/endpoint analysis and replay from frozen complete-action Hessian blocks; optional full-action producers are separate. |
| [12 · Complexity-weighted moments and corona states](papers/12-complexity-moment-corona/README.md) | Analytic determinacy, critical escape and quotient/reachability arguments; finite geometric-series values illustrate proved bounds. |
| [13 · Thermal point-process statistics](papers/13-thermal-statistics-reconstruction/README.md) | Exponential bounds, conditional reconstruction and explicit counterexamples; finite diagnostics do not close the separate Airy density problem. |
| [14 · Operational instrument geometry](papers/14-operational-instrument-geometry/README.md) | Endpoint normalization, phase/reference obstructions and complete-tree arguments, with exact algebra and finite covariance checks. |
| [15 · Positive Cartan carriers](papers/15-positive-cartan-carriers/README.md) | Explicit Clifford carrier, positive-selector obstruction and restricted invariant/action analysis, supported by exact finite certificates. |
| [16 · Fuzzy half-turn covariance](papers/16-fuzzy-half-turn-limits/README.md) | Finite spectrum/covariance and sector-resolved limit proofs; heat-trace sequences are diagnostics rather than proofs of the limits. |
| [17 · Finite Cartan–Regge closure](papers/17-cartan-regge-closure/README.md) | Local metric quotient, exact Schur blocks and nonlinear source calculations; numerical fold search has its own limited scope. |
| [18 · Positive boundary sewing](papers/18-positive-boundary-sewing/README.md) | Typed framework, operator-valued CP criterion and two finite operational families; exhaustive small tests and seeded samples are distinguished. |
| [19 · Local gauge benchmark](papers/19-local-gauge-benchmarks/README.md) | Source-specialized Green/Ward/BRST composition and finite sign checks; physical-state existence is imported, ordinary-dagger descent is obstructed and the interacting anomaly remains uncomputed. |
| [20 · AI-assisted research case study](papers/20-ai-research-case-study/README.md) | Preserved corpus reaggregation and keyword sensitivity, with historical PDF counts retained as evidence; no causal audit effect or model ranking is inferred. |

Paper 07 modifies a clock-foliated effective action and constructs a finite-band Gaussian. Paper 09 proves the scoped canonical coupled physical-state theorem. They are distinct constructions; neither finite algebra command proves the other's analytic result. Paper 19 uses the physical-state theorem as an input and addresses additional local benchmark interfaces.

## Analytical evidence in papers 21–24

- [21 · Weyl relative entropy](papers/21-weyl-diamond-relative-entropy/manuscript.md). Free linearized Weyl gravity in a round four-dimensional Minkowski diamond. The exclusion concerns the stated class of bulk densities and local edge jets, not all descriptions of entropy.
- [22 · Lorentzian Regge Airy asymptotics](papers/22-lorentzian-regge-airy/manuscript.md). An explicit six-vertex circuit and a localized real contour. It does not establish a global contour prescription or a refinement limit.
- [23 · Gravitational memory and boosts](papers/23-gravitational-memory-boost/manuscript.md). Shell-flat, rotation-equivariant quasifree spin-two representations. No universal Lorentz-symmetry claim and no result for general smooth overlapping radial bases.
- [24 · Calibrated de Sitter tidal diffusion](papers/24-calibrated-tidal-diffusion/manuscript.md). Free linearized Einstein gravity in the Bunch–Davies state, with the specified tracer calibration. Neither detector backreaction nor a state-construction theorem is proved.

The statements and derivations are in the Markdown files. The earlier twenty-packet calculation replays do not verify these four proofs, and no new external review is claimed by this integration.

## Analytical evidence in papers 25–31 and companions 29B and 29B+31

Papers [25 and 26](docs/PAPER-INDEX.md#25-characteristic-data-and-radiative-symplectic-rank) give complete scoped Markdown derivations. The typeset papers below provide TeX sources and paper-specific reviews. Internal cross-checks and clean compilation do not constitute external expert review.

- [27 · Parity-odd Weyl-cubic ADM gravity](papers/27-parity-odd-weyl-cubic-adm/VERIFICATION.md): Check the normal-acceleration elimination and boundary term, the unrestricted bracket using evolution but not multiplier equations, and the reconstructed metric map. The radiative and field-redefinition-invariant vertex tests have different roles.
- [28 · The Einstein tensor from Lorentzian Regge variation](papers/28-lorentzian-regge-einstein-tensor/VERIFICATION.md): Check the sign at timelike and spacelike joints, Hessian symmetry at the smooth basepoint, the separate bulk/face/hinge Taylor bounds, and the affine interpolation example.
- [29 · Physical polarization in the vacuum high-frequency limit](papers/29-high-frequency-graviton-polarization/VERIFICATION.md): Check the screen quotient and Hermitian test, the covariant matrix symbol and compensated source terms, and the two opposite-helicity exact vacuum families.
- [29B · Physical polarization of high-frequency vacuum gravity](papers/29b-high-frequency-graviton-polarization/VERIFICATION.md): Check the exact null constraints, compensated shear transport in the constructed overlap, and calibrated frequency-filtered curvature extraction. These are analytic checks, independent of the forward-limit verification for Paper 29.

- [30 · Nariai rigidity under null completeness](papers/30-nariai-null-completeness/VERIFICATION.md): Check the covering/seed from Galloway–Ling, the normal-connection identity and Killing-field integral, the equal pole-period obstruction, and Lafuente López’s global identification. Cross-checked analytically; no numerical dependency.
- [31 · Classically matched vacuum waves with opposite quantum transfer](papers/31-graviton-metaplectic-transfer/VERIFICATION.md): Check the smooth trace-free retuning, the even-repeat metaplectic sign, nonzero-band gravitational constraint return, fixed Weyl dual smears, and exact same-collar composition. Cross-checked analytically; no numerical dependency.
- [29B+31 · Polarization defects and opposite graviton transfers](papers/29b-31-polarization-defects-graviton-transfer/README.md): See the [manuscript](papers/29b-31-polarization-defects-graviton-transfer/manuscript.pdf) for the exact vacuum shear normal form, differentiable endpoint averaging, simple-angle retuning, finite-wavelength matching, local optical normalization and physical band bracket. No numerical dependency.

## Reusing evidence and reading the history

Unchanged completed evidence can be retained without rerunning it for a documentation update. A changed model, input class, state, time interval or analytic hypothesis needs the corresponding argument or calculation for that setting. Packet-specific `VERIFICATION.md` and `PROVENANCE.json` files identify these boundaries.

The [historical record](docs/VERIFICATION-HISTORY.md) reports the terminal/desktop builds, rendered-page inspections, portability runs and archive checks performed for the recorded deliveries. Earlier paper 01/02 editions remain preserved with their edition history. Saved archive-check results describe their recorded packaging state; they are not an instruction to rerun the checker whenever prose changes.

This documentation-only revision read packet instructions, requirements and relevant output-writing code. It did **not** execute scientific calculations, optional producers, manuscript builds, new hash inventories or the archive checker. Guide commands are instructions for readers, not new execution receipts.

The recorded review consists of automated calculations and AI review. Human expert review has not occurred; no comprehensive novelty certificate or publication-priority determination is supplied.
