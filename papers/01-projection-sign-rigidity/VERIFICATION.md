# Verification of the revised projection paper

## Fresh and analytic coverage

The feasible-fiber equivalence follows directly from `H_R={h in H:Ah=R}` and the complete-support sign/holonomy gauge correspondence. The manuscript includes the explicit rank-one example; `python verification/check_revision.py` passed exact rational idempotency, the collision of the two displayed holonomy fields and exhaustive enumeration of all 1024 gauge-fixed fields, with exactly one target solution.

The generic interior-rank theorem is unchanged. Its all-size validity rests on the manuscript's analytic nonvanishing and finite-zero-locus proof, already assessed in the completed recovery; a finite example does not prove genericity. Historical supplementary Wolfram checks and expensive original experiments were not rerun.

## Build and static checks

TeX Live 2025 latexmk/pdflatex/BibTeX compiled successfully to `build/main.pdf`; the copied `manuscript.pdf` has **6 pages**. The final log has no LaTeX warnings, undefined references/citations or overfull/underfull boxes. PDF Author is blank. Source cleanup preserves the exact uncertainty statement and contains no TODO/draft placeholders or duplicate labels. No new external project paths are needed by the paper or default check. Integrated visual QA is performed by the repository integrator.

## Limits

No complete current novelty search, human expert review, noisy/efficient decoder, exceptional-stratum classification, Hermitian phase theorem or continuum limit is certified. The rank-one explanatory example does not test or refute the generic interior-rank theorem.
