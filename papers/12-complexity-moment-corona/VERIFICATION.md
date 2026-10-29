# Verification of the complete article and package

Completed: 29 October 2026 (Beijing). Human authors and affiliations: blank. No human expert review is claimed.

## Article coverage

The complete English article gives the central definitions and arguments for finite-exponential-moment determinacy; weighted/translated heat clocks and the exact heat normalizer; sufficient sector counting/intensity bounds; the true weighted abscissa; the finite-critical total-variation branch and divergent-critical escape branch; positive-moment divergence under escape; two-sided order estimates; convergent-profile constants and exponential scaling; scalar and bounded-sector corona descent; ambient ultrafilter separation; the word/total-ideal counterexample; the further word quotient and isometric central embedding; and actual long-block Abel reachability.

Classical imported results are explicitly identified: Stieltjes--Carleman determinacy, Banach--Alaoglu compactness, standard positive state extension and basic C*-quotient facts. The manuscript does not claim a new general moment or Banach-limit theorem. The extra regular-profile scaling corollary is proved directly by Riemann sums; it is presented as a standard consequence under the additional convergence hypothesis, not a general consequence of comparability.

The translated heat-clock support gap is repaired by the explicit exponential-moment ratio. The false general inference that mass escape may leave a positive complexity moment bounded is absent. The order estimates are not promoted to asymptotic constants without a residual limit. Fixed-word microlocal rays are not silently killed by a word-c0 pullback. The further J_W quotient is defined and proved directly; preservation of this ideal by physical operations is not asserted.

The code and table reuse the already verified exact-formula construction. Finite numeric rows illustrate the analytic bound and do not prove a limit by extrapolation. Explicit subnets, rather than unjustified convergent sequences, are used for compactness on the full nonseparable bounded-sequence algebra.

## Bibliographic verification

Six entries are cited; all six are used and all citation keys resolve.

- Simon (1998), *Advances in Mathematics* 137, 82–203, DOI 10.1006/aima.1998.1728: primary author/arXiv paper checked at Proposition 1.5, Corollary 4.5, Section 3 and Appendix C. The arXiv posting is 1999; the publication date remains 1998.
- Patterson (1976), *Acta Mathematica* 136, 241–273, DOI 10.1007/BF02392046: primary paper checked at Section 3, especially normalized positive measures and boundary concentration. Its geometric conclusions are not transplanted to the present arbitrary graded ensemble.
- Lorentz (1948), *Acta Mathematica* 80, 167–190, DOI 10.1007/BF02393648: primary scanned article checked at Section 1 and Theorem 1, including the explicit parity/Banach-limit example.
- Alaoglu (1940), *Annals of Mathematics* 41(1), 252–267, DOI 10.2307/1968829: publisher/JSTOR issue record verifies the entry; compactness is a standard imported result, not freshly reproved here.
- Murphy (1990), *C*-Algebras and Operator Theory*, Academic Press, ISBN 9780125113601: official Elsevier product metadata and publisher preview verify the book. The preview identifies the author as Gerard J. Murphy; the shop page's current “Gerald” spelling was not copied. The verification page is https://shop.elsevier.com/books/c-algebras-and-operator-theory/unknown/978-0-08-092496-0 . No third-party book or paper is bundled.
- Collier, Eberhardt and Zhang (2023), *SciPost Physics* 15, 151, DOI 10.21468/SciPostPhys.15.4.151: primary arXiv version 2304.13650v3 checked at equation (2.23), Appendix B and the fixed-topology scope. Only the density convention is an input to the heat example.

The DOI resolver calls for Lorentz/Alaoglu did not return pages in the browser tool; the primary scan and publisher issue record supplied the needed evidence. The ScienceDirect book endpoint returned 403; the official Elsevier product record and publisher preview supplied metadata. This is a targeted reference check, not a global priority certificate.

## Executed portability check

`python code/verify_constructive.py` completed with Python 3.11.15, using only standard-library modules. The script was then copied into an isolated package-local directory under `build/portable-check/code/` and run from that copied root. It required no old project path, network, provider or external input.

The fresh copied-location output, the normal package reproduction output and the saved recovery output were byte-identical, SHA-256 `ff80506f7610d22a8749d9c74cd4342d1f60d7ed096900f7398be48017382758`. The copied saved input/output bytes are preserved. `saved-output/package-check.json` records the successful comparison, interpreter version, labels and bibliography checks.

The script evaluates geometric block masses and the long-block observable using exact sum formulas with floating-point exponentials. It checks the analytic deviation bound with a small rounding allowance. The independently evaluated complement uses `expm1`; late witness values round to endpoints, while positive analytic error bounds remain. No Monte Carlo, arbitrary precision loop or physical-model simulation was run.

## Compilation and cleanup

Command: `latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=build main.tex`.

Toolchain: latexmk 4.86a, pdfTeX 1.40.27, TeX Live 2025, BibTeX. The completed document is 11 A4 pages. `build/main.pdf` was copied to `manuscript.pdf`. `pdfinfo` confirms an empty Author field.

The first compile exposed a missing xcolor package used by hyperlink color expressions; this was fixed directly. A long publisher URL caused bibliography underfull boxes; the verified book is now identified by ISBN in the bibliography, with the full verification URL retained here. Final `build/main.log` has no LaTeX warnings, undefined citations/references, overfull boxes or underfull boxes. The Perl wrapper emits a locale-fallback warning on this Windows installation; the compiler and bibliography build completed successfully.

Cleanup checked 48 unique labels, resolved all referenced labels and all six citation keys, and found no duplicate label, unused bibliography entry or draft placeholder. Scientific names and public article text contain no internal project codename. The manuscript supplies the requested chronology and exact authorship/credit quote. The separate completed batch marker was not rerun.

Root performs integrated PDF visual QA and repository upload. This worker's verification is source/proof, compilation, metadata and code portability; it does not claim a fresh page-by-page visual review.

## Material limits

The article is a mathematical methods synthesis with explicit repairs and examples, not a global originality certificate. Actual physical weighted critical asymptotics, a selected common operation quotient, unbounded-domain control, a physical witness/pair and semiclassical survival remain unconstructed. No journal acceptance, experimental realization or human specialist validation is claimed. Original manuscripts, recovery evidence, audits and the publication checkout were not edited.
