# Verification of the local-clock construction

Run on 29 October 2026, Beijing:

```powershell
python 'scripts/construct_local_clocks.py'
```

Exit code 0, 23 checks passed, using Python 3.11.15, NumPy 2.4.3 and SciPy 1.17.1. The saved machine-readable output is `results/local-clocks.json`, including both the shorter reduced-kinetic point and the explicit rest-energy-inclusive point.

New checks cover the raised-cosine pulse area and half-areas, local unitary endpoint, clock energy positivity and available negative reserve, focused minimum uncertainty and exact focus-width formula, independent kinetic-energy variance, the analytic Gaussian-tail regime, nonoverlapping pulses, applicability of the retained nonlinear comparison, total entanglement margin, and the charged reflected/negative-momentum bound. Hermite quadrature at orders 64 and 128 followed by adaptive time integration gives inner error integrals 0.00102610224930 and 0.00102610224942, respectively. This independently checks the scale and conservatism of the analytic clock error; it does not replace the analytic proof.

The new proof is a full-Hamiltonian Duhamel comparison, not a free-clock trajectory substituted into the interacting solution. It needs no grid approximation of the full nonlinear geometry Hilbert space. Self-adjointness follows from tensor sums and bounded perturbations; positive apparatus energy follows from the operator inequality h_r>=0. Exact positive-dust support follows from a single specified spectral projection commuting with the autonomous total Hamiltonian. These operator arguments are stated explicitly in DERIVATION.md; the script checks their finite formulas, not a formal proof assistant transcription.

The calculation reads only the bundled extraction of the accepted reference JSON for its phase constant and switching times, and reconstructs the new parameter point and old accepted bound coefficients algebraically. It does not import or rerun the old producer. The JSON records the accepted input's SHA-256. The retained nonlinear Duhamel estimate, compact-shape norm and source graph estimate are reused rather than independently re-proven. Their physical regime remains nu<=1e-12, the same completed finite-code H, and the specified cutoff radial/compact-shape seed.

The primary, rest-energy-inclusive theorem is rounded conservatively to D<.04014, N>.459 and dust gap>.459 E0. It also checks the total rest cost .08 E0 and the corresponding spectral/clock-tail regime. The shorter reduced-kinetic point gives D<.04855, N>.451 and dust gap>.399 E0. The Gaussian outer-tail logarithms are below -31211 and -31220, respectively. Their machine exponentials underflow, but the analytic positive remainders are included in the proof and are vastly smaller than the displayed rounding margins. This run uses double precision with a wide strict margin, not a new outward interval or formal proof certificate.

Material assumptions remain: the one initial joint total-H spectral preparation, the specified chirped Gaussian packets, ideal closed-system local control potentials and positive offsets, the accepted Q-dependent completion, and the new reduced apparatus Hamiltonian. Rest masses are explicitly charged in the primary point, but quadratic nonrelativistic dispersion is assumed exactly; small velocity does not establish uniform agreement with relativistic propagation over the very long T. No covariant apparatus stress tensor, unknown-input quantum channel, hardware preparation error, concurrent environment or human review has been verified. Clock backreaction and negative-momentum amplitude are bounded in norm; no unbounded final clock-energy expectation has been inferred from trace closeness.

## Complete manuscript compilation

The complete English manuscript is `main.tex`, with three verified primary records in `references.bib` and build/use instructions in `README.md`. It contains the trusted/completed model and seed, weighted echo and sign derivation, explicit nonlinear comparison and initial graph constants with their analytic derivation, the four actual clocks, full interacting comparison, initial spectral resource, rest-energy-inclusive theorem and limitations. No inaccessible internal document is the sole support for a mathematical estimate.

Ran from `paper` on 29 October 2026:

```powershell
latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=build main.tex
Copy-Item -LiteralPath build/main.pdf -Destination manuscript.pdf -Force
```

Latexmk 4.86a, pdfTeX 1.40.27 / TeX Live 2025 and BibTeX 0.99d completed with exit code 0. The final PDF has 11 pages. The converged `build/main.log` has no warnings, unresolved citations/cross-references, multiply defined labels, or overfull/underfull boxes. The initial pass's unresolved-reference notices were resolved by latexmk's normal BibTeX/recompile sequence. A Perl host-locale fallback notice did not affect compilation.

PyMuPDF inspection confirmed 11 pages and an empty PDF Author field. Human bylines and affiliations are blank. The exact operator uncertainty quotation and distinctions between the assembly date, historical work and public disclosure are included. The parent task is performing integrated visual QA; these compile/text checks do not claim a completed independent visual review. No PDF-creation marker was rerun, and no publication checkout was modified.

Bibliographic verification used official APS pages for [Sørensen and Mølmer (1999)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.82.1971), [Salecker and Wigner (1958)](https://journals.aps.org/pr/abstract/10.1103/PhysRev.109.571), and [Vidal and Werner (2002)](https://link.aps.org/doi/10.1103/PhysRevA.65.032314). Citation purpose is limited to the established closed-loop principle, general quantum-clock context and negativity definition, respectively.

## Integrated final manuscript inspection

The desktop editor compiler also confirmed success. The parent inspected all 11 rendered pages in contact sheets and pages 5, 9 at readable resolution. No clipping, overlap, missing glyphs or unresolved references were observed. Human byline and PDF author metadata are blank; the corrected chronology and full uncertainty statement are present. Manuscript-facing text contains no internal project codename. No unchanged scientific calculation was rerun for this layout pass.

## Private repository delivery replay

The default commands in README.md were executed successfully from this repository packet on 29 October 2026 Beijing time. They use no external project directory or network input. This reuses completed scientific coverage and verifies portability; it is not an additional independent proof. The manuscript PDF, LaTeX and bibliography match the previously completed edition byte for byte.
