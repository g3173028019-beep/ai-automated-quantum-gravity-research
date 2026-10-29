# Finite Cartan–Regge Schur closure

The complete English paper is `manuscript.pdf` (9 pages), with editable `main.tex` and `references.bib`.

The paper presents the fixed-boundary metric quotient on a four-dimensional 4–2 cell, its geometric holonomy identity on a stated logarithm branch, the exact independent-angle Schur elimination, the convention-mapped Regge determinant factor, and nonlinear sourced curvature response. It includes the rational calibration coordinates and enough local geometric construction to define the exact certificate. It is separate from any causal spin-two/refinement paper.

## Reproduce

Use Python 3.10 or later and the dependencies in `requirements.txt`. From this directory:

```text
python code/replay.py
```

The replay reads `inputs/calibration.json`, checks it against the bundled exact geometry, and uses only the local modules in `code/`. It rebuilds the matching ranks and four bordered angle–Gram blocks, checks the exact Hessian and target ratio, recomputes both nonlinear source saddles at `j=1e-6`, and locates the numerical fold. It writes `outputs/verification.json` and compares its decisive quantities against the unchanged established output in `outputs/established_verification.json`.

The script works from another working directory when called by its absolute path. There is no old project path, network, raw archive or provider dependency in the default replay. Historical source identifiers and hashes occur only in `PROVENANCE.json`.

## Build

```text
latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=build main.tex
```

Copy `build/main.pdf` to `manuscript.pdf` after rebuilding a changed source. `VERIFICATION.md` states actual execution and inherited coverage, including the limits of the numerical fold search. `PROVENANCE.json` records code adaptation and source/output hashes.

Human authors, affiliations and PDF Author are blank. Human expert review has not occurred. Neither the finite one-loop quotient nor its source response is presented as a completed quantum theory, a propagation theorem or a derived connection measure.
