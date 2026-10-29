# Verification of the finished manuscript and package

Date: 29 October 2026, Beijing. Human expert review: none. Human authors and affiliations remain blank. No original research source, completed recovery or publication checkout was edited.

## Analytic content

The manuscript supplies explicit definitions and proofs for the single-internal-angular-channel rank bound, canonical scalar transfer identities, endpoint-identification repair and exact independent-endpoint counterexample, positive conserved quadratic form, and scalar counterterm/common-dilation matching theorem. The angular commutants are derived by finite matrix equations. Their scope is finite real or complex linear algebra with the stated symmetry, support and nondegeneracy hypotheses.

The rank bound uses the unique internal radial angular-weight-two pair and the absence of a simplex containing moving vertices from both caps. Adding internal angular copies or nonlocal terms changes the hypotheses. A positive invariant does not select a Hamiltonian logarithm branch, quantum energy positivity or a physical vacuum.

## Executed for this package

- `python -B reproduce.py`: exit code 0. Verified all 20 frozen-file hashes, exact algebra and counterexample, independently rebuilt the explicit geometry and causal-normal counts for levels 1–5, reconstructed saved transfer matrices, checked their spectra and symplectic residuals, and verified the scalar repairs. The receipt is `results/replay.json`.
- `python -B run_full.py --help`: exit code 0. This loaded the packaged complete-action producer/import dependencies and displayed the explicit level selector. It verifies import portability and entry-point availability, **not** a new full numerical action run.
- `latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=build main.tex`: completed successfully using pdfLaTeX and BibTeX. The finished PDF has nine pages. Final-log inspection found no undefined references/citations, missing characters, overfull boxes or underfull boxes. The Windows wrapper emitted a nonfatal Perl locale fallback warning.

The final PDF is copied to `manuscript.pdf`. PDF Author is explicitly empty in the source. Integrated rendered-page visual QA is assigned to the parent delivery; it is not represented here as already completed by this worker.

The tested runtime had NumPy 2.4.3, SymPy 1.14.0 and SciPy 1.17.1. `pdfinfo` confirms nine A4 pages and an empty Author field. A focused AI synthesis check of the structural locality/rank proof, detail momentum signs, endpoint proposition and scalar matching theorem found no further bearing error; this is an AI check, not human peer review. The final text explicitly restricts “kernel” to the quadratic principal function.

## Retained full-producer evidence

`data/level-1-reproduction.json` and `data/level-5-reproduction.json` are unchanged complete-action outputs from the completed recovery on the same day. Their SHA-256 values are recorded in `PROVENANCE.json`; optional producer dependencies are copied byte-for-byte. The first run took 12.70 seconds and gives mixed singular values about 0.840003 with elliptic arguments about 0.815035. The fifth took 423.37 seconds and gives mixed singulars 0.47996436720 and 0.47996419596, reciprocal-real eigenvalues near 1.29025 and 0.77504 with imaginary residuals at most 1.17e-6, overlap about 0.8919179, symplectic residual 1.76e-16 and smallest physical internal singular 0.31048966.

These full calculations were **not rerun for this manuscript package**. The default replay verifies their saved matrices, rather than independently differentiating the action. The generic positive-phase output field at the fifth level is deliberately not treated as a physical oscillation phase. Its tiny values are arguments of numerical imaginary residuals on a reciprocal-real spectrum.

Numerical levels 2–4, the original second finite-difference-step runs, curvature-impedance reconstruction and earlier branch diagnostics are inherited evidence. The connected original reports are frozen under provenance. The article does not claim to recertify all of those numerical calculations independently.

## Remaining limits

There is no second independent complete Regge action implementation, interval-certified floating-point error analysis, physical transport of the endpoint identification through every survivor quotient, full phase-space refinement intertwiner, local fine-lattice realization of the scalar counterterm, or full gluing/constraint preservation result. Maximum mesh diameter does not shrink in the tested partial refinement. No asymptotic continuum rejection is asserted.

Physical positive frequency, retarded support, a massless pole, unique continuum helicity, nonlinear stability, a ghost-free quantum measure and a full continuum limit remain open. The bibliography preserves credit to established canonical Regge, Lorentzian action and perfect-action work; the novelty comparison is targeted, with no comprehensive global priority claim.
