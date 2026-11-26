# Reproduce a calculation

[Home](../README.md) · [Paper index](PAPER-INDEX.md) · [Verification](../VERIFICATION.md) · [Delivery history](VERIFICATION-HISTORY.md)

All 20 original packets have a default portable entry point. Reading an analytical proof, checking finite algebra, replaying a saved certificate and rebuilding its original producer are different operations. Choose a row below, then use the packet's `README.md` and `VERIFICATION.md` for its scientific scope.

## Setup

Run each listed command from the corresponding paper directory. Follow that packet's Python minimum and `requirements.txt`: documented minimum versions range from Python 3.8 to 3.11, and dependencies differ. A separate environment per packet keeps version pins or ranges independent. For a packet with third-party dependencies, a PowerShell setup example is:

```powershell
cd papers/03-jt-energy-weighted-recovery
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
python code/reproduce.py
```

Dependency installation may use the network. The default scientific commands use packaged definitions or inputs locally after setup. Papers 01, 02, 12 and 20 use only Python's standard library and need no third-party package installation. The local requirements' exact pins and version ranges remain authoritative.

## Papers 01–10

Output paths are relative to the paper directory. A packaged recorded run is labeled separately when a command prints rather than writes it.

| Paper directory | Default command | Dependencies | Output and meaning |
|---|---|---|---|
| [01-projection-sign-rigidity](../papers/01-projection-sign-rigidity/README.md) | `python verification/check_revision.py` | Standard library; Python ≥3.10 | `verification/revision-checks.json`: exact fixed-data revision check, including 1024 star-gauge assignments. The generic theorem rests on its analytic proof. |
| [02-half-wormholes-wiener-chaos](../papers/02-half-wormholes-wiener-chaos/README.md) | `python verification/check_revision.py`<br>`python verification/check_chaos_extras.py` | Standard library; Python ≥3.10 | `verification/revision-checks.json`, `verification/chaos-extras-results.json`: exact variance reductions, convolution/sign controls and quadrature moments. No Monte Carlo is required. |
| [03-jt-energy-weighted-recovery](../papers/03-jt-energy-weighted-recovery/README.md) | `python code/reproduce.py` | mpmath, NumPy, SymPy | `results/reproduction.json`: charge/Gaussian identities, finite coherent/reference checks and interval bound evaluations. |
| [04-dbi-heralded-controlled-phase](../papers/04-dbi-heralded-controlled-phase/README.md) | `python replay.py` | mpmath, SymPy | `results/gate.json`: saved-input/certificate checks, residual recombination and conditional gate bounds. It does not reintegrate the continuum residual. |
| [05-integer-flux-families](../papers/05-integer-flux-families/README.md) | `python reproduce.py` | SymPy | `results/construction_certificate.json`, `results/divisor_assignments.json`: exact assignments, rank, masses and charge calculations. |
| [06-maxwell-local-clocks](../papers/06-maxwell-local-clocks/README.md) | `python scripts/construct_local_clocks.py` | NumPy, SciPy | `results/local-clocks.json`: clock/resource/error estimates using the bundled reference extraction. |
| [07-clock-foliated-reference](../papers/07-clock-foliated-reference/README.md) | `python scripts/construct_inhomogeneous.py` | NumPy, SciPy, SymPy | `results/construction.json`, `results/controlled_ell1_evolution.csv`, `results/quadratic_action.txt`: effective-model harmonic reduction, positivity and finite-band evolution. |
| [08-einstein-scalar-mixed-quartic](../papers/08-einstein-scalar-mixed-quartic/README.md) | `python code/native_mixed_quartic.py`<br>`python code/finite_backreaction.py` | SymPy | `results/native-mixed-quartic.json`, `results/finite-backreaction.json`: quartic/source checks and separate finite-metric field identities. |
| [09-canonical-einstein-scalar-hadamard](../papers/09-canonical-einstein-scalar-hadamard/README.md) | `python code/canonical_algebra_check.py` | SymPy | **Prints JSON to stdout.** `results/algebra-results.json` is the packaged recorded run, not a file automatically rewritten by this command. It supports normalization, charged adjoints, finite-rank defects and the example background. |
| [10-quantum-weight-networks](../papers/10-quantum-weight-networks/README.md) | `python code/verify.py` | SymPy | `outputs/verification.json`: exact cell/network and binary-symplectic checks using `inputs/model.json`. `outputs/established_verification.json` is retained reference output. |

Paper 07 modifies a specified clock-foliated effective action and constructs a finite-band Gaussian. Paper 09 proves physical Hadamard-state existence for the canonical coupled system under its compact-Cauchy and critical-potential assumptions. These are different constructions; neither finite algebra command replaces the other's analytic hypotheses or proof.

## Papers 11–20

| Paper directory | Default command | Dependencies | Output and meaning |
|---|---|---|---|
| [11-causal-regge-transfer](../papers/11-causal-regge-transfer/README.md) | `python reproduce.py` | NumPy, SymPy, SciPy | `results/replay.json`: frozen-file checks, quotient/endpoint algebra, geometry and spectra reconstructed from packaged Hessian blocks. Full Hessian regeneration is optional below. |
| [12-complexity-moment-corona](../papers/12-complexity-moment-corona/README.md) | `python code/verify_constructive.py` | Standard library; Python ≥3.8 | `reproduced-output/calculation-output.json` and stdout: finite geometric-formula illustrations. `saved-output/calculation-output.json` preserves the earlier output. Asymptotic and quotient results are analytic. |
| [13-thermal-statistics-reconstruction](../papers/13-thermal-statistics-reconstruction/README.md) | `python verification/check_bounds.py`<br>`python verification/check_channels.py` | SymPy, SciPy | `verification/checks.json`, `verification/channels.json`: symbolic identities, finite determinantal examples, quadrature diagnostics and channel controls. |
| [14-operational-instrument-geometry](../papers/14-operational-instrument-geometry/README.md) | `python code/verify.py` | NumPy, SymPy | `data/results.json`: normalization, phase/tree/reference checks and the finite covariance example. |
| [15-positive-cartan-carriers](../papers/15-positive-cartan-carriers/README.md) | `python code/verify.py` | NumPy, SymPy | `data/results.json`: carrier/selector, invariant and restricted-action certificates constructed from definitions. |
| [16-fuzzy-half-turn-limits](../papers/16-fuzzy-half-turn-limits/README.md) | `python verification/check_fuzzy.py` | SymPy | `verification/checks.json`: four exact spin representations, symbolic oscillator parameters and finite heat-trace diagnostics. |
| [17-cartan-regge-closure](../papers/17-cartan-regge-closure/README.md) | `python code/replay.py` | NumPy, SciPy, SymPy | `outputs/verification.json`: exact geometry/Schur blocks and recomputed nonlinear saddles/fold search from `inputs/calibration.json`; compares decisive quantities with retained `outputs/established_verification.json`. |
| [18-positive-boundary-sewing](../papers/18-positive-boundary-sewing/README.md) | `python reproduce.py` | NumPy | `outputs/doubled_instrument.txt`, `outputs/weyl_charge_cut.txt`, `outputs/diagnostic.txt`, `outputs/receipt.json`: finite operational families and independent ancillary/restriction/compression controls. |
| [19-local-gauge-benchmarks](../papers/19-local-gauge-benchmarks/README.md) | `python reproduce.py` | SymPy | `outputs/green_matrix.txt`, `outputs/ward_transport.txt`, `outputs/graded_algebra.txt`, `outputs/matter_vertex.txt`, `outputs/receipt.json`: 49 printed coefficient/sign rows, not a new PDE or interacting-QME proof. |
| [20-ai-research-case-study](../papers/20-ai-research-case-study/README.md) | `python reproduce.py` | Standard library; Python ≥3.10 | `outputs/results.json`: saved hash-relation recount, corpus totals and filename-rule sensitivity. Historical PDF parsing is retained evidence, not repeated against unbundled PDFs. |

## Optional heavy producers and historical calculations

These are separate from the defaults. Their inclusion does not mean they were rerun for the latest documentation or packaging change.

**Continuum residual generation for paper 04.** The retained `reference/certify_trial_residual.py` documents interval Taylor-cell integration and an analytic tail. The default `replay.py` consumes saved `inputs/` partitions. The reference producer is not a drop-in replay command: its code expects `continuum-trial.json` and resumable `residual-partition-<branch>.json` beside the script. Regeneration therefore requires a separately prepared working copy with those source-defined inputs and branch/slice settings. A bounded time slice may save an incomplete resumable partition; that is not a completed residual certificate. The packet's README and verification note explain this retained-input boundary.

**Complete-action Hessians for paper 11.** From that packet directory, the optional local producer is:

```powershell
python run_full.py --level 1
python run_full.py --level 5
```

These commands write `results/full-level-1.json` and `results/full-level-5.json`. Level 5 recomputes the substantially heavier complete-action Hessian. The default replay uses frozen blocks and does not call this producer. The optional command uses the preserved original implementation and local imports, not an independent second Regge implementation. Frozen inputs are not overwritten.

**Historical Wolfram Language supplements for papers 01 and 02.** The unchanged scripts remain optional after the current portable exact Python checks. With Wolfram Language available, from the chosen packet directory:

```powershell
wolframscript -file verification/legacy_checks.wls
```

The projection supplement contains exponential subset-sum enumeration; the chaos supplement includes 400,000-sample Monte Carlo comparisons. Inspect the script and workload before choosing to run it. These historical checks were not repeated for the current revisions. Their finite or stochastic calculations do not replace the manuscripts' analytic genericity, moment or limit arguments.

## Papers 21–26: reading and checking the derivations

- [21 · Weyl relative entropy](../papers/21-weyl-diamond-relative-entropy/manuscript.md): An explicit coherent-state quadratic form in electric and magnetic radial curvature, and an obstruction to representing it by a regular finite-jet local curvature density.
- [22 · Lorentzian Regge Airy asymptotics](../papers/22-lorentzian-regge-airy/manuscript.md): A nonzero cubic on a light-cone-regular branch leads to a uniform localized Airy approximation and opposite Fresnel signs for the separated saddles.
- [23 · Gravitational memory and boosts](../papers/23-gravitational-memory-boost/manuscript.md): Memory is implementable in an explicit finite-energy example, but Doppler leakage into unsqueezed radial modes rules out joint implementation with a nontrivial standard boost in the stated class.
- [24 · Calibrated de Sitter tidal diffusion](../papers/24-calibrated-tidal-diffusion/manuscript.md): Repeatedly recalibrating a freely falling finite-aperture tracefree tidal record yields a controlled nonzero diffusion scale.

These Markdown manuscripts are the editable derivations. They do not include a separate calculation runner or a LaTeX/PDF build; the packet build instructions below apply to papers 01–20.

- [25 · Characteristic data](../papers/25-characteristic-data-symplectic-rank/manuscript.md): Recompute the endpoint-moment rank and check the symplectic Gram pairing on disjoint radial intervals. The Raychaudhuri check is a separate second-order statement.
- [26 · Neutral de Sitter states](../papers/26-neutral-de-sitter-mass-susceptibility/manuscript.md): Check the neutral-mode comparison and the susceptibility lower bound on bounded local approximants. This is an analytic derivation, not a simulation.

## Papers 27–31 and companion 29B: rebuild the analytical articles

With a standard TeX Live installation, open the chosen paper directory and run `latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex`; the generated file is `main.pdf`. The repository copy is named `manuscript.pdf`.

- [27 · Parity-odd Weyl-cubic ADM gravity](../papers/27-parity-odd-weyl-cubic-adm/README.md): Check the normal-acceleration elimination and boundary term, the unrestricted bracket using evolution but not multiplier equations, and the reconstructed metric map. The radiative and field-redefinition-invariant vertex tests have different roles.
- [28 · The Einstein tensor from Lorentzian Regge variation](../papers/28-lorentzian-regge-einstein-tensor/README.md): Check the sign at timelike and spacelike joints, Hessian symmetry at the smooth basepoint, the separate bulk/face/hinge Taylor bounds, and the affine interpolation example.
- [29 · Physical polarization in the vacuum high-frequency limit](../papers/29-high-frequency-graviton-polarization/README.md): Check the screen quotient and Hermitian test, the covariant matrix symbol and compensated source terms, and the two opposite-helicity exact vacuum families.
- [29B · Physical polarization of high-frequency vacuum gravity](../papers/29b-high-frequency-graviton-polarization/README.md): Check the exact characteristic data, transverse transport, and normalized frequency-filtered Weyl readouts on the specified overlap. See its [verification note](../papers/29b-high-frequency-graviton-polarization/VERIFICATION.md) for the proof and rebuild scope; no numerical replay is required.

These are analytical derivations, not numerical runs. Rebuilding confirms the source and layout, not the mathematical argument.

- [30 · Nariai rigidity under null completeness](../papers/30-nariai-null-completeness/README.md): `cd papers/30-nariai-null-completeness` and run `pdflatex -halt-on-error main.tex` twice. Check the covering/seed from Galloway–Ling, the normal-connection identity and Killing-field integral, the equal pole-period obstruction, and Lafuente López’s global identification. Cross-checked analytically; no numerical dependency.
- [31 · Classically matched vacuum waves with opposite quantum transfer](../papers/31-graviton-metaplectic-transfer/README.md): `cd papers/31-graviton-metaplectic-transfer` and run `pdflatex -halt-on-error main.tex` twice. Check the smooth trace-free retuning, the even-repeat metaplectic sign, nonzero-band gravitational constraint return, fixed Weyl dual smears, and exact same-collar composition. Cross-checked analytically; no numerical dependency.

- [32 · Type-II scale recovery](../papers/32-typeii-scale-recovery/README.md): Rebuild the standalone `main.tex`. Check the Laurent selection rule, the three-scale signed rational witness in the PDF, and the proof of the matching full-record lower bound. No numerical data or extra TeX sources are required.

## Rebuild a manuscript

With an existing compatible TeX distribution providing `latexmk`, `pdflatex` and BibTeX, run from the repository root:

```powershell
./tools/build.ps1 -Paper 03-jt-energy-weighted-recovery
```

The helper builds the selected paper and copies `build/main.pdf` to its `manuscript.pdf`. Omitting `-Paper` builds every manuscript directory selected by the helper. To build one paper directly from its directory:

```powershell
latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=build main.tex
Copy-Item -LiteralPath build/main.pdf -Destination manuscript.pdf -Force
```

Some papers include additional local TeX parts; use `main.tex` as the entry point and keep the packet files together. Papers 01 and 02 also provide a local `build.ps1` implementing the source/BibTeX/PDF-copy workflow. Building changes build artifacts and the delivery copy; it is separate from scientific reproduction.

## Find the evidence you need

| File or location | Purpose |
|---|---|
| `manuscript.pdf` | Reading copy |
| `main.tex`, local included TeX files, `references.bib` | Editable scientific argument and bibliography |
| `README.md`, `requirements.txt` | Packet-specific setup, command and dependencies |
| `VERIFICATION.md` | Executed coverage, inherited evidence and unresolved limits |
| `PROVENANCE.json`, retained input/output records | Source identity and the role of packaged evidence |
| The output path in the tables above | Result produced or printed by the selected entry point |

Saved `archive-checks.json` and the [verification history](VERIFICATION-HISTORY.md) record prior packaging work. `tools/check_archive.py` and `tools/render_review.py` are packaging/layout helpers, separate from scientific entry points; they need not be rerun for an unchanged calculation or a prose-only documentation update.

This guide was updated by reading local packet instructions and relevant code. None of the scientific, heavy-producer, build or archive-check commands above was executed for this documentation-only revision.

## Building Papers 33–37

The TeX for these five articles is self-contained. In each paper directory run `latexmk -pdf -halt-on-error manuscript.tex` or run `pdflatex -halt-on-error manuscript.tex` twice. The build produces a local `manuscript.pdf`.

- [33 · Time-resolved scalar returns from a marked timelike tube](../papers/33-vacuum-wave-returns/README.md)
- [34 · Normalizable quadrupole states off the vacuum stress orbit](../papers/34-ads3-vacuum-stress/README.md)
- [35 · Positive string poles and indefinite finite-transfer graviton amplitudes](../papers/35-typeii-full-polarization/README.md)
- [36 · Full counting statistics of neutral emission in a quantum charged throat](../papers/36-charged-throat-emission/README.md)
- [37 · Winding rank and shape confinement in homogeneous Einstein–scalar tori](../papers/37-winding-rank/README.md)

## Building the tensor spectral geometry papers

In each linked packet, run `latexmk -pdf -halt-on-error manuscript.tex` or run `pdflatex -halt-on-error manuscript.tex` twice. The peer-review status and any standalone exact replay are described in the packet verification notes.

- [38 · From Euclidean tensor prices to EPRL amplitudes: a rank obstruction at the four-to-two move](../papers/38-tensor-prices-eprl-boundary/README.md)
- [39 · Spectral Steiner costs on real trees and a Poisson field on the p-adic boundary](../papers/39-real-tree-spectral-boundary-fields/README.md)
- [40 · Null-volume shared-label rigidity](../papers/40-null-volume-shared-label-rigidity/README.md)
- [41 · Turaev–Viro fusion, recoupling and Wilson moduli](../papers/41-tensor-prices-tv-fusion-and-wilson/README.md)
- [42 · Open Turaev–Viro Dehn fillings](../papers/42-tensor-prices-tv-dehn-filling/README.md)

## Building Papers 44–47

Each paper is a self-contained analytical article. Run `latexmk -pdf -halt-on-error -jobname=manuscript main.tex` in the packet directory, or run `pdflatex -halt-on-error -jobname=manuscript main.tex` three times. This is a document rebuild, not a numerical reproduction of the theorem.

- [44 · Gravitational region sewing](../papers/44-gravitational-region-sewing/README.md)
- [45 · Vacuum stress and its massless pole](../papers/45-vacuum-stress-massless-pole/README.md)
- [46 · The common Einstein light cone](../papers/46-einstein-scalar-common-cone/README.md)
- [47 · Spacelike area tomography](../papers/47-spacelike-area-tomography/README.md)
