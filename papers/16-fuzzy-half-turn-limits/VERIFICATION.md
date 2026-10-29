# Verification

Edition: 29 October 2026, Beijing. Human expert review has not occurred.

## Mathematical coverage

The paper includes analytic proofs of the exact finite spectrum and half-turn covariance, the fixed-polar-window bound, the dominated-convergence four-sector heat-trace limit, the two persistent cross-pole states, and the Weyl-symbol parity formula. The generator-level comparison with the prior fuzzy harmonic regulator is explicit. The Fourier formula is proved for Schwartz kernels and extended by Hilbert--Schmidt/L2 density; no unstated continuum embedding is inferred.

`python verification/check_fuzzy.py` passed under Python 3.11.15 and SymPy 1.14.0. It independently constructs spin generators with exact radicals and checks complete matrices at j=1/2,1,3/2,2, with theta=1 and symbolic Omega:

- every entry of the self-dual spectral formula;
- exact covariance and projective square;
- both cross-pole eigenvectors and their energy;
- the axial difference from the generator form of the prior regulator.

The output `verification/checks.json` also contains finite heat-trace diagnostics. Those decimals do not prove the limit. The Weyl integral identity is verified analytically in the text rather than by a sampled transform.

## References and scope

Primary arXiv/publisher metadata were checked for Langmann--Szabo, Grosse--Wulkenhaar, Chu--Madore--Steinacker, and Ydri. The Langmann--Szabo Proposition 1 order-two convention and Ydri section 4.2 regulator formulas were inspected at formula level in the recovery. The relevant generator comparison was replayed exactly here. This is a targeted comparison, not a comprehensive novelty certificate.

The old large-spin plotting producers were not rerun; their output is not needed by the new proofs. There is no Monte Carlo result, interacting continuum proof, global fuzzy Fourier theorem, or operator-topology limit in this packet.

## Compilation and cleanup

`latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=build main.tex` completed successfully with TeX Live 2025, pdfTeX 1.40.27, latexmk 4.86a and BibTeX. The delivered manuscript has **6 pages**.

The first build exposed a missing xcolor load, which was fixed. The final LaTeX log contains no undefined references/citations, LaTeX warnings, overfull boxes, or underfull boxes. The Perl launcher emitted a locale fallback warning; this did not affect the successful PDF build and is distinct from the clean final TeX log. The PDF author field is blank.

The source was checked for placeholders, duplicate labels, citation-key consistency, old-path dependencies and outward internal codenames. The original scientific archives and publication checkout were not modified. Parent integration owns visual PDF QA and upload; this packet does not claim those steps were performed by the worker.
