# 38. From Euclidean tensor prices to EPRL amplitudes: a rank obstruction at the four-to-two move

[Home](../../README.md) · [Paper index](../../docs/PAPER-INDEX.md) · [PDF](manuscript.pdf) · [TeX](manuscript.tex) · [Verification](VERIFICATION.md)

**Result.** Coherent cut-price operators reproduce finite corrected Euclidean EPRL amplitudes and their gluing. A marked 4↔2 boundary block has rank five against a rank-one coarse seam, with exact cutoff-limit values at γ=1/3.

**Scope.** Corrected Euclidean EPRL at the stated rational Immirzi parameters, weights and boundary marks. The spin-cutoff limit keeps the two complexes fixed; it is not a graph continuum limit.

## Build and replay

From this directory run `latexmk -pdf -halt-on-error manuscript.tex` or run `pdflatex -halt-on-error manuscript.tex` twice.

- `python replay_rank_two.py` checks selected exact formulae.
- `python replay_rank_bounds.py` checks 2,050 finite rational/cutoff rank inequalities and 1,000 minimum-pair bounds.

The manuscript is an AI-led research artifact with internal checks. Human author and affiliation fields are blank; no external expert review or peer review has been claimed. See [repository declarations](../../DECLARATIONS.md).
