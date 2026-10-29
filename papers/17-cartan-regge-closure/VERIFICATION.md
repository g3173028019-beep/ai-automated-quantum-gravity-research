# Verification record

Date: 29 October 2026, Beijing. The portable replay completed successfully with bundled inputs and local modules. It was invoked by absolute script path from an unrelated working directory, with no historical project import or network dependency.

## Fresh portability results

Python dependencies: NumPy 2.4.3, SciPy 1.17.1, SymPy 1.14.0. Exit code 0.

- Metric matching: initial matching rank 25, boundary-combined rank 39 and nullity 1; final matching rank 6, combined rank 20 and nullity 0.
- Rebuilt all four exact bordered auxiliary blocks. Their nonzero determinants and Schur contributions reproduce the manuscript table.
- Exact metric-channel Hessian `-224/16875`, independent coordinate target `1792/16875`, ratio `-1/8`.
- At source `1e-6`, roots `2.2500252370920344` and `2.255001222506141`; sourced equation residuals below `1.2e-16`.
- Their Hessians have opposite signs: approximately `-0.0131419932222` and `+0.0132542946421`.
- Internal Regge/Cartan curvature sums agree within `1.6e-15` on the stated near-flat logarithm branch.
- Fold witness `q=2.25251922313574`, source `5.028351131077464e-5`; finite-difference Hessian residual `6.86e-8`.

The fresh output is `outputs/verification.json`. The prior established output was copied without changing bytes to `outputs/established_verification.json`; its provenance hash is recorded. Exact blocks, hessian and ratio agree. Root positions agree within `2e-10`, Hessians within `1e-7`, and fold source within `1e-10`. These tolerances compare a portable numerical replay with saved evidence, not a new physical error model.

## Analytic content and reused evidence

The manuscript explains the metric-to-edge quotient, the local geometric angle–Gram branch, Schur reduction, exact scale law and principal-log branch restriction. The numerical source equation and flat response follow from the stated Regge action with the boundary term.

The replay reconstructs the inspected original producer rather than a second independently implemented first-order geometric model. The historical thirteen-geometry first-order and twenty-five-geometry Regge-target ensembles were not repeated because their producer logic and inputs were unaffected. Their role is labeled as inherited numerical evidence. A previous monolithic legacy test timed out and is not counted as passed.

## Build and LaTeX cleanup

`latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=build main.tex` completed with exit code 0 under TeX Live 2025. Final PDF: 9 pages. The final log has no undefined/multiply defined references, citation warnings or overfull/underfull boxes. PDF Author is explicitly blank. A shell-level Perl locale fallback did not affect compilation.

The full source was checked for source-path dependencies, draft placeholders, label/citation errors and obsolete internal terminology. Integrated PDF visual QA and upload are the parent delivery process's responsibility; this worker does not count them as already completed.

## Material limits

The principal-log identity is exact only in its declared branch or with an explicit continuous lift. The scalar optimizer and sign-change root search provide finite numerical witnesses, not a global proof that every stationary branch has been found. The work does not derive an independent connection measure, simplicity/ghost/BFV Jacobian, global Lorentzian contour, two-helicity causal propagator or continuum limit. Human expert review and exhaustive global novelty verification have not occurred.
