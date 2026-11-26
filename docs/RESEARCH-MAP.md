# Research map

[Home](../README.md) · [Paper index](PAPER-INDEX.md) · [Reproduction](REPRODUCTION.md)

The archive connects explicit information tasks, coupled-field state construction, finite geometry/control, and limits of reconstruction. A channel error bound, a Hadamard state-existence theorem, an exact finite algebra and a descriptive process study answer different scientific questions. The routes below organize the collection by those questions while retaining each action, state space, input class and observable.

## Geometry–matter information tasks

**Question:** what information task can a specified dynamics perform, and what preparation or apparatus does it require?

| Paper | Result to read for | Essential scope |
|---|---|---|
| [03. JT energy-weighted recovery](../papers/03-jt-energy-weighted-recovery/README.md) | A channel approximation uniform on a mean-energy-bounded class, including coherent tails and external references | The specified fixed-topology Hamiltonian, initial Gaussian and time interval |
| [04. DBI heralded gate](../papers/04-dbi-heralded-controlled-phase/README.md) | A controlled-phase instrument with explicit uniform success and conditional-fidelity bounds | The independent-field half-line model, charge code, radial herald and local corrections |
| [06. Maxwell local clocks](../papers/06-maxwell-local-clocks/README.md) | Autonomous local pulse controls with positive apparatus energy and a certified entangled output | The reduced geometry–Maxwell model, joint initial spectral preparation and stated large resource scales |

Suggested route: **04 → 06 → 03** for gate, apparatus and channel estimates. These are different models. Shared techniques include complete instruments, reference-system tests, full-output error estimates and resource accounting; a theorem transfers only after its target assumptions are checked.

## Coupled fields: states, stability and interactions

[34 · Normalizable quadrupole states off the vacuum stress orbit](../papers/34-ads3-vacuum-stress/README.md) An explicit normalizable quantum quadrupole family realizes stress one-point profiles outside the AdS₃ vacuum stress orbit. A Hill spectral criterion characterizes classical vacuum-orbit representatives of prescribed stress shapes.

[36 · Full counting statistics of neutral emission in a quantum charged throat](../papers/36-charged-throat-emission/README.md) A single-channel charged-throat emission kernel yields an exact scaled cumulant generator and a large-deviation rate for the emission count.

**Question:** which compatible physical state, quadratic dynamics or finite interaction can actually be constructed?

| Paper | Positive construction | Distinction to retain |
|---|---|---|
| [09. Canonical Einstein–scalar Hadamard states](../papers/09-canonical-einstein-scalar-hadamard/README.md) | Positive gauge-invariant quasifree Hadamard covariances on the physical quotient, with a spacetime lift | A finite canonical real scalar multiplet with fixed positive-definite target metric; smooth on-shell globally hyperbolic four-dimensional background, compact Cauchy surface and potential gradient zero along the scalar image. Nonconstant fields and a nonzero, not necessarily positive, Hessian are included. |
| [09+21, Edition 2. Local curvature and modular dynamics](../papers/09-21-coupled-modular-kernel/README.md) | Positive physical states and convergent modular reconstruction on the full local algebra | Any smooth globally hyperbolic solution of the canonical finite-scalar theory with a fixed positive target kinetic form and arbitrary smooth potential. The original round four-field curvature-completion and entropy example is retained. |
| [07. Preferred-clock reference model](../papers/07-clock-foliated-reference/README.md) | Positive physical inhomogeneous quadratic sectors and an explicit 141-oscillator finite-band Gaussian, preserving the homogeneous action | A modified effective action with altered characteristics. This is quadratic stability and a finite-band state, not paper 09’s canonical metric-Hadamard theorem or nonlinear stability. |
| [08. Native mixed quartic and Walker geometry](../papers/08-einstein-scalar-mixed-quartic/README.md) | An off-shell mixed contact including auxiliary feedback, plus a separate exact matter-supported finite solution | The finite solution has real split signature `(2,2)`. A Lorentzian real lift and full quartic chart descent are not supplied. |
| [19. Local Einstein–scalar gauge benchmark](../papers/19-local-gauge-benchmarks/README.md) | A nonminimal Green matrix, compatible free BRST/equicausal interfaces and explicit dagger/interaction limits | An imported action and local AdS-interior/compact-host construction. State existence is credited as an input; a full interacting QME is not proved. |

Suggested routes: **09 → 19** for physical state construction and its use in a local gauge benchmark; **07** for the distinct effective-model harmonic repair; **08** for native interactions and the split-signature finite solution. None substitutes for the hypotheses or conclusion of another.

## Gravity observables: from states to readouts

[33 · Time-resolved scalar returns from a marked timelike tube](../papers/33-vacuum-wave-returns/README.md) Complete temporal kernels in finitely many angular channels recover a marked homogeneous-orbit metric and its twist. Vacuum eliminates the conformal ambiguity at a nonumbilic tube, but exact inhomogeneous vacuum ambiguity survives finite angular projection.

### Physical degrees of freedom and state comparisons

- [25 · Characteristic data and symplectic rank](../papers/25-characteristic-data-symplectic-rank/manuscript.md): in a solved linear characteristic sector the terminal sphere record is governed by two moments, while its kernel contains infinitely many radiative symplectic pairs. Read beside [19](../papers/19-local-gauge-benchmarks/) for the separate issue of free gauge reduction; neither proof assumes the other.
- [26 · Neutral de Sitter mass susceptibility](../papers/26-neutral-de-sitter-mass-susceptibility/manuscript.md): removing a homogeneous scalar position pole gives the stated neutral comparison, with local normality but divergent mass-squared entropy susceptibility. [21](../papers/21-weyl-diamond-relative-entropy/manuscript.md) compares coherent Weyl states on a different algebra; [24](../papers/24-calibrated-tidal-diffusion/manuscript.md) concerns a spin-two tidal record, not the scalar mass problem.


| Comparison | What is shared | What does not follow |
|---|---|---|
| [09](../papers/09-canonical-einstein-scalar-hadamard/README.md) → [21 · Weyl relative entropy](../papers/21-weyl-diamond-relative-entropy/manuscript.md) | Both work with physical, gauge-invariant degrees of freedom. | Constructing a coupled Hadamard state does not give the free Minkowski diamond’s modular form. |
| [09](../papers/09-canonical-einstein-scalar-hadamard/README.md) → [18](../papers/18-positive-boundary-sewing/README.md) → [23 · Gravitational memory and boosts](../papers/23-gravitational-memory-boost/manuscript.md) | Separate the existence of a positive representation from the operations it supports. | Positivity or an abstract sewing test does not imply implementation of Lorentz boosts on a finite-energy memory representation. |
| [06](../papers/06-maxwell-local-clocks/README.md) → [24 · Calibrated de Sitter tidal diffusion](../papers/24-calibrated-tidal-diffusion/manuscript.md) | Both make the physical readout or its calibration explicit. | The Maxwell apparatus is not the de Sitter tracer protocol. |
| [21 · Weyl relative entropy](../papers/21-weyl-diamond-relative-entropy/manuscript.md) ↔ [24 · Calibrated de Sitter tidal diffusion](../papers/24-calibrated-tidal-diffusion/manuscript.md) | Both use linearized Weyl curvature rather than gauge-fixed metric components. | A local-algebra relative entropy in Minkowski does not determine a long-time de Sitter diffusion coefficient. |

## Moments, low-order data and reconstruction

[35 · Positive string poles and indefinite finite-transfer graviton amplitudes](../papers/35-typeii-full-polarization/README.md) For the full 35-dimensional ten-dimensional physical projectile polarization space and a fixed real probe, positive Regge-subtracted massive poles coexist with an indefinite finite-transfer coefficient.

[Paper 32](../papers/32-typeii-scale-recovery/README.md) applies sparse recovery to fixed-transfer graviton scattering: the nonforward coefficients expose scale moments absent at forward transfer, and an explicit ambiguity pair proves the optimal precision.

**Question:** when do partial data identify a matrix, a scalar law, or a random spectrum?

| Paper | Main distinction |
|---|---|
| [01. Real projection recovery](../papers/01-projection-sign-rigidity/README.md) | Idempotency restores generic signs from order-one/two data; uniqueness at an observed feasible fiber is weaker than uniform injectivity. |
| [02. Gaussian replicas and Wiener chaos](../papers/02-half-wormholes-wiener-chaos/README.md) | Replica/sector formulas and finite-model moment indeterminacy; canonical discrete extensions do not exhaust all representing measures. |
| [12. Complexity-weighted moments and corona states](../papers/12-complexity-moment-corona/README.md) | Fixed-regulator determinacy, critical mass escape and actual regulator-reachable states are separate properties. |
| [13. Thermal point-process statistics](../papers/13-thermal-statistics-reconstruction/README.md) | Determining one heat statistic, reconstructing a point-process law from joint clocks, and obtaining a pure state are different tasks. |
| [18. Positive boundary systems and sewing](../papers/18-positive-boundary-sewing/README.md) | Positive scalar data do not automatically provide noncommutative positivity, reflection, complete positivity, normalization or quotient compatibility. |

Suggested route: **01 → 02 → 13**, then **12 → 18** for regulator limits and a typed boundary framework. The Gaussian-polynomial law in 02 is distinct from the point-process heat statistics in 13. The latter retains a conditional semigroup method while leaving the stochastic-Airy tail-density target unresolved. Neither moment growth nor failure of a sufficient Carleman condition alone proves indeterminacy.

The current revisions of 01/02 apply local characterization, attribution, proof-exposition and interpretation corrections while preserving their valid main results. Earlier source versions and copies are lineage, not extra independent papers.

## Finite control, operational geometry and positive carriers

**Question:** what algebra or geometry is generated by specified finite operations, and which extra structure fails to follow?

| Paper | Result to read for | Boundary |
|---|---|---|
| [10. Quantum-weight networks](../papers/10-quantum-weight-networks/README.md) | Spin(5) cells, a chiral Spin(10) pair, and full `su(4^n)` control on every stated connected network with n>=3 | Independent controls and bonds are hypotheses; iid rms convergence is not operator-norm or interacting-Gibbs convergence. |
| [14. Resolved instruments](../papers/14-operational-instrument-geometry/README.md) | Covariant endpoint normalization, phase non-identifiability and a spectral-reference obstruction with calibrated completion | Finite operational statements do not derive a spacetime coframe or gravitational action. |
| [15. Positive Cartan carriers](../papers/15-positive-cartan-carriers/README.md) | A positive carrier together with a signed Lorentz selector obstruction, including the semidefinite boundary | Classical Clifford/action ingredients are credited; a physical selector or coupling scale is not selected. |
| [16. Fuzzy half-turn limits](../papers/16-fuzzy-half-turn-limits/README.md) | Exact finite covariance, a uniform polar-window defect, four-copy heat-trace limit and persistent cross-pole modes | A selected-sector limit is not a full operator-topology contraction or interacting continuum theory. |

Suggested routes: **14 → 15** for operational versus carrier/selector geometry, and **10 → 16** for finite algebra versus collective or spectral limits. These are mathematical reading routes, not an assertion that the papers define the same gravity model.

## Three distinct finite Regge questions

| Paper | Question answered | What remains separate |
|---|---|---|
| [17. Fixed-boundary Cartan–Regge closure](../papers/17-cartan-regge-closure/README.md) | Exact independent-angle elimination, metric quotient matching, determinant placement and nonlinear sourced curvature on one 4–2 complex | A branch-scoped Euclidean finite-complex result is not a Lorentzian propagation theorem or a connection path-integral construction. |
| [11. Causal transfer and partial refinement](../papers/11-causal-regge-transfer/README.md) | A finite angular-weight-two transfer, architectural rank bound and subdivision dependence with corrected endpoint identification | The tested family is not a vanishing-mesh continuum limit; the fifth level’s hyperbolic transfer is not a universal continuum no-go. |
| [22 · Lorentzian Regge Airy asymptotics](../papers/22-lorentzian-regge-airy/manuscript.md) | Cubic saddle coalescence and a uniform Airy approximation on the actual localized real contour | One Lorentzian 4–2 circuit. Do not identify its contour with 17’s Euclidean quotient or treat it as a refinement theorem for 11. |

Suggested route: **17 → 11 → 22** to compare exact local elimination, finite causal transfer and the resolution of a saddle merger. Papers 17 and 11 use Schur reduction; the issue in 22 is uniform semiclassical asymptotics. Their geometries, quotients and observables differ. Shared terminology does not merge these papers.

## Integer compactification data

[05. Integer STU flux families](../papers/05-integer-flux-families/README.md) constructs exact F-term solutions and positive masses in the complex-structure/axiodilaton sectors, then relates a weak-coupling branch to integral charge capacity. Read the divisor construction and mass spectrum before the pointwise lattice minima. The inherited global geometry and untreated Kähler directions delimit the claim; the within-family optimum is not a global landscape optimum.

## Research-process evidence

[20. Claim correction and artifact growth](../papers/20-ai-research-case-study/README.md) studies preserved artifact states, actual correction cases and sensitivity of a filename classification. It is a descriptive methods paper. File growth, a mutable timestamp cutoff or a changed classification fraction does not establish model superiority, a causal audit effect or scientific value by itself.

This provides a route into how the archive’s arguments were corrected, not a ranking or validation of the other 19 papers. [Declarations](../DECLARATIONS.md) records operator-reported model chronology and unresolved credit; [History](../HISTORY.md) separates source dates and editions.

## The Einstein dynamics and correspondence limits

[37 · Winding rank and shape confinement in homogeneous Einstein–scalar tori](../papers/37-winding-rank/README.md) Full winding rank locks the classical shape and yields a compact fixed-clock quantum resolvent. Nonzero deficient rank leaves an escape direction and excludes point spectrum on the marked space.

- [27 · Parity-odd Weyl-cubic ADM gravity](../papers/27-parity-odd-weyl-cubic-adm/README.md) starts with a continuum curvature interaction and asks what remains of the geometric constraint algebra after order reduction. Compare [19](../papers/19-local-gauge-benchmarks/) for a distinct free gauge and causality question.
- [28 · The Einstein tensor from Lorentzian Regge variation](../papers/28-lorentzian-regge-einstein-tensor/README.md) asks whether the complete signed variation converges along genuinely refining Lorentzian meshes. It extends the archive's finite Regge questions [11](../papers/11-causal-regge-transfer/), [17](../papers/17-cartan-regge-closure/) and [22](../papers/22-lorentzian-regge-airy/) in a different direction; those papers establish neither its hypotheses nor its limit.
- [29 · Physical polarization in the vacuum high-frequency limit](../papers/29-high-frequency-graviton-polarization/README.md) asks what physical tensor information survives a high-frequency weak limit. [25](../papers/25-characteristic-data-symplectic-rank/) addresses a related loss of radiative information at a characteristic endpoint, using a different construction and no high-frequency limit.
- [29B · Physical polarization of high-frequency vacuum gravity](../papers/29b-high-frequency-graviton-polarization/README.md) reorganizes the generalized-wave transport theorem and adds an independent characteristic construction: two independent screen covariances at fixed resolved scalar geometry, with frequency-selected tidal readouts in the marked overlap. The two theorems have different hypotheses.

## Exact vacuum rigidity and quantum transfer

- [30 · Nariai rigidity under null completeness](../papers/30-nariai-null-completeness/README.md): Under the stated four-dimensional Einstein vacuum, positive cosmological constant, two-sided null-completeness, compact nonspherical semidefinite Cauchy hypotheses, the Lorentzian universal cover is the Nariai product.
- [31 · Classically matched vacuum waves with opposite quantum transfer](../papers/31-graviton-metaplectic-transfer/README.md): Two smooth exact vacuum plane waves match full classical Jacobi and longitudinal returns but have opposite calibrated physical two-helicity band Weyl brackets. For (N) concatenated intervals, the adjacent brackets distinguish (2^N) relative lift histories while the exterior map retains their product.
- [29B+31+Q19 · Causal counts for polarization and signed transfer](../papers/29b-31-polarization-defects-graviton-transfer/README.md) links these complementary observations on one exact plane-wave class: the lifted physical transfer follows from short marked intervals, while fixed-length chords resolve the background tidal polarization.

The two manuscripts use different stocks. The Lorentzian universal-cover rigidity result is classical; it does not furnish a quantum state for the sandwich-wave comparison.

## Tensor spectral geometry

- [38 · From Euclidean tensor prices to EPRL amplitudes: a rank obstruction at the four-to-two move](../papers/38-tensor-prices-eprl-boundary/README.md): Coherent cut-price operators reproduce finite corrected Euclidean EPRL amplitudes and their gluing. A marked 4↔2 boundary block has rank five against a rank-one coarse seam, with exact cutoff-limit values at γ=1/3. Corrected Euclidean EPRL at the stated rational Immirzi parameters, weights and boundary marks. The spin-cutoff limit keeps the two complexes fixed; it is not a graph continuum limit.
- [39 · Spectral Steiner costs on real trees and a Poisson field on the p-adic boundary](../papers/39-real-tree-spectral-boundary-fields/README.md): A compatible spectral assignment recovers real-tree distances and finite Steiner-hull lengths. A Poisson field has explicit non-Gaussian connected boundary kernels below threshold and a Gaussian white-noise limit at critical scaling. A classical tree-indexed probability field and a tensor spectral encoding, with the existing finite Poisson model explicitly credited. No spacetime or gravity construction is asserted.
- [40 · Null-volume shared-label rigidity](../papers/40-null-volume-shared-label-rigidity/README.md): Classical CGHS constraints and a static-kernel classification fix when the chosen shared-label sampler can reproduce physical null volumes. In four dimensions, spherical flat tubes need exactly three labels, whereas every regular Schwarzschild exterior box has unbounded exact cost. The joint physical volume sampler is assumed, and the bound is on exact conditional-independence labels. It neither constrains arbitrary tensor encodings nor detects a horizon from rank alone.
- [41 · Turaev–Viro fusion, recoupling and Wilson moduli](../papers/41-tensor-prices-tv-fusion-and-wilson/README.md): Physical cylinder tensors recover quantum dimensions and signed recoupling data. Joint Wilson projector prices also reconstruct the known marked SU(2) torus character variety in the large-level limit. Exact finite physical fusion and coherent operator tests, plus a separate marked torus limit. The native recovery is not a new character variety or a derivation of dynamics.
- [42 · Open Turaev–Viro Dehn fillings](../papers/42-tensor-prices-tv-dehn-filling/README.md): Open fillings with the same nonzero closed three-sphere amplitude have different exact tensor prices. Even caps with the same separate prices can have different sewing behavior. Finite doubled SU(2) state sums with marked whole-torus boundary components. This is a boundary sewing result, not a claim about distinct closed manifolds or the physical zero-cosmological-constant limit.

## Boundary sewing, matter cones and metric probes

- [44 · Gravitational region sewing](../papers/44-gravitational-region-sewing/README.md): The minimum independent register for a bracket-preserving linear dressing is the cross-pairing rank. Linked linear graviton charges require ten pairs per independent link for a single electric–magnetic pairing. A propagating polarization instead has continuous spectral load, yielding quantitative finite-reservoir obstructions on planar and Schwarzschild cuts.
- [45 · Vacuum stress and its massless pole](../papers/45-vacuum-stress-massless-pole/README.md): The large-window vacuum variance is an explicit positive functional times the massless-shell atom. Only a helicity-zero stress tensor structure survives on that shell. A separate scalar–tensor model distinguishes Ricci focusing and the bare source from the complete linear gravitational charge.
- [46 · The common Einstein light cone](../papers/46-einstein-scalar-common-cone/README.md): For positive strongly hyperbolic scalar modes with independent gradients on an open attained family, equality of the scalar and Einstein characteristic sets forces a positive nonlinear sigma model with a potential. Two exact constrained data families distinguish this property from hyperbolicity alone.
- [47 · Spacelike area tomography](../papers/47-spacelike-area-tomography/README.md): Twelve spacelike area types recover the marked metric with two independent consistency checks; ten are locally sufficient and minimal among continuous pointwise scalar readings. On a regular intrinsic patch, the differential of the inverse gives area profiles for every compact conserved linearized Einstein test.

These are separate model-scoped results. No map edge identifies a finite canonical register with a radiative bath, a stress-pole readout with a graviton pole, or a marked area inverse with nonlinear quantization.

## Using the routes

The [paper index](PAPER-INDEX.md) links all manuscripts and distinguishes the original PDF/source packets, later analytical Markdown manuscripts and typeset analytical articles. Use packet notes for commands and dependencies, and the [verification overview](../VERIFICATION.md) for the distinction between proof, exact finite checks, numerical diagnostics and retained inputs. No additional scientific rerun is implied by this navigation layer.

The 53 archive items are selected entry points into a wider corpus. Related source drafts are consolidated where they answer the same question; independent scoped questions remain distinct. Human authors and affiliations remain blank, and human expert review has not occurred. Neither a completed packet nor its private source date certifies global novelty, public priority, journal acceptance or full quantum gravity.
