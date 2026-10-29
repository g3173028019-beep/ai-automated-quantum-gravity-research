# Research map

[Home](../README.md) · [Paper index](PAPER-INDEX.md) · [Reproduction](REPRODUCTION.md)

The archive connects explicit information tasks, coupled-field state construction, finite geometry/control, and limits of reconstruction. A channel error bound, a Hadamard state-existence theorem, an exact finite algebra and a descriptive process study answer different scientific questions. The routes below group all 20 papers by those questions while retaining each action, state space, input class and observable.

## Geometry–matter information tasks

**Question:** what information task can a specified dynamics perform, and what preparation or apparatus does it require?

| Paper | Result to read for | Essential scope |
|---|---|---|
| [03. JT energy-weighted recovery](../papers/03-jt-energy-weighted-recovery/README.md) | A channel approximation uniform on a mean-energy-bounded class, including coherent tails and external references | The specified fixed-topology Hamiltonian, initial Gaussian and time interval |
| [04. DBI heralded gate](../papers/04-dbi-heralded-controlled-phase/README.md) | A controlled-phase instrument with explicit uniform success and conditional-fidelity bounds | The independent-field half-line model, charge code, radial herald and local corrections |
| [06. Maxwell local clocks](../papers/06-maxwell-local-clocks/README.md) | Autonomous local pulse controls with positive apparatus energy and a certified entangled output | The reduced geometry–Maxwell model, joint initial spectral preparation and stated large resource scales |

Suggested route: **04 → 06 → 03** for gate, apparatus and channel estimates. These are different models. Shared techniques include complete instruments, reference-system tests, full-output error estimates and resource accounting; a theorem transfers only after its target assumptions are checked.

## Coupled fields: states, stability and interactions

**Question:** which compatible physical state, quadratic dynamics or finite interaction can actually be constructed?

| Paper | Positive construction | Distinction to retain |
|---|---|---|
| [09. Canonical Einstein–scalar Hadamard states](../papers/09-canonical-einstein-scalar-hadamard/README.md) | Positive gauge-invariant quasifree Hadamard covariances on the physical quotient, with a spacetime lift | A finite canonical real scalar multiplet with fixed positive-definite target metric; smooth on-shell globally hyperbolic four-dimensional background, compact Cauchy surface and potential gradient zero along the scalar image. Nonconstant fields and a nonzero, not necessarily positive, Hessian are included. |
| [07. Preferred-clock reference model](../papers/07-clock-foliated-reference/README.md) | Positive physical inhomogeneous quadratic sectors and an explicit 141-oscillator finite-band Gaussian, preserving the homogeneous action | A modified effective action with altered characteristics. This is quadratic stability and a finite-band state, not paper 09’s canonical metric-Hadamard theorem or nonlinear stability. |
| [08. Native mixed quartic and Walker geometry](../papers/08-einstein-scalar-mixed-quartic/README.md) | An off-shell mixed contact including auxiliary feedback, plus a separate exact matter-supported finite solution | The finite solution has real split signature `(2,2)`. A Lorentzian real lift and full quartic chart descent are not supplied. |
| [19. Local Einstein–scalar gauge benchmark](../papers/19-local-gauge-benchmarks/README.md) | A nonminimal Green matrix, compatible free BRST/equicausal interfaces and explicit dagger/interaction limits | An imported action and local AdS-interior/compact-host construction. State existence is credited as an input; a full interacting QME is not proved. |

Suggested routes: **09 → 19** for physical state construction and its use in a local gauge benchmark; **07** for the distinct effective-model harmonic repair; **08** for native interactions and the split-signature finite solution. None substitutes for the hypotheses or conclusion of another.

## Moments, low-order data and reconstruction

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

## Two distinct finite Regge questions

| Paper | Question answered | What remains separate |
|---|---|---|
| [17. Fixed-boundary Cartan–Regge closure](../papers/17-cartan-regge-closure/README.md) | Exact independent-angle elimination, metric quotient matching, determinant placement and nonlinear sourced curvature on one 4–2 complex | A branch-scoped Euclidean finite-complex result is not a Lorentzian propagation theorem or a connection path-integral construction. |
| [11. Causal transfer and partial refinement](../papers/11-causal-regge-transfer/README.md) | A finite angular-weight-two transfer, architectural rank bound and subdivision dependence with corrected endpoint identification | The tested family is not a vanishing-mesh continuum limit; the fifth level’s hyperbolic transfer is not a universal continuum no-go. |

Suggested route: **17 → 11** to compare exact local elimination with causal transfer. Both use Schur reduction, but their geometries, quotients and observables differ. Shared terminology does not merge these papers.

## Integer compactification data

[05. Integer STU flux families](../papers/05-integer-flux-families/README.md) constructs exact F-term solutions and positive masses in the complex-structure/axiodilaton sectors, then relates a weak-coupling branch to integral charge capacity. Read the divisor construction and mass spectrum before the pointwise lattice minima. The inherited global geometry and untreated Kähler directions delimit the claim; the within-family optimum is not a global landscape optimum.

## Research-process evidence

[20. Claim correction and artifact growth](../papers/20-ai-research-case-study/README.md) studies preserved artifact states, actual correction cases and sensitivity of a filename classification. It is a descriptive methods paper. File growth, a mutable timestamp cutoff or a changed classification fraction does not establish model superiority, a causal audit effect or scientific value by itself.

This provides a route into how the archive’s arguments were corrected, not a ranking or validation of the other 19 papers. [Declarations](../DECLARATIONS.md) records operator-reported model chronology and unresolved credit; [History](../HISTORY.md) separates source dates and editions.

## Using the routes

The [paper index](PAPER-INDEX.md) links every current PDF, editable main source, packet README and verification record. Use packet notes for commands and dependencies, and the [verification overview](../VERIFICATION.md) for the distinction between proof, exact finite checks, numerical diagnostics and retained inputs. No additional scientific rerun is implied by this navigation layer.

The 20 packets are selected entry points into a wider corpus. Related source drafts are consolidated where they answer the same question; independent scoped questions remain distinct. Human authors and affiliations remain blank, and human expert review has not occurred. Neither a completed packet nor its private source date certifies global novelty, public priority, journal acceptance or full quantum gravity.
