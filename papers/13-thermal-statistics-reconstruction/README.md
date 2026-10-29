# Thermal Point-Process Statistics: Exponential Bounds, Reconstruction, and Their Limits

The ten-page `manuscript.pdf` combines bounded-statistic determinacy, explicit limits of scalar and joint thermal identification, determinantal channel asymptotics, and a conditional semigroup-density method. Its complete TeX sources are `main.tex` and `references.bib`.

The central positive result is the exponential bound E exp(tZ) <= exp[(m/b)(exp(tb)-1)] under bounded nonnegative weights and factorial-correlation domination, with its concentration consequence. A positive-spectrum lognormal Cox family has identical complete mixed-clock moments but distinct laws when the correlation condition is omitted. A two-level example separates a determinate scalar law from spectral reconstruction. The paper proves a sufficient superlinear-growth channel criterion and disproves its converse with an exactly pure projection DPP. Full joint clocks reconstruct a point-process law within the specified correlation/integrability class; this does not imply a pure state or one-dimensional GNS.

The semigroup lemma is conditional on explicit pathwise hypotheses. It provides an absolute-continuity/Jacobian reduction and does not close the literal stochastic-Airy tail-density target. Generic finite projection sign rigidity is outside this paper. The methods are credited classical ingredients organized around explicit assumptions, not a claimed new general theory or a black-hole dynamics theorem.

## Reproduce

With Python 3.11 or later and the packages in `requirements.txt`, run from this directory:

```text
python verification/check_bounds.py
python verification/check_channels.py
```

The scripts need no network or external project files. They write `verification/checks.json` and `verification/channels.json`. `verification/reference-bounds.json` preserves the earlier bounded-statistic output bytes as provenance/reference evidence; the reproduction does not need to load it.

The checks cover the exact signed lognormal moment identity, two-level example, finite determinantal MGF examples, Cox Laplace quadrature diagnostics, cumulant coefficients, the unrestricted-clock counterexample, symbolic pure projection blocks, and the Bernoulli GNS distinction. Channel-ratio decimals illustrate an analytically proved limit and are not an infinite-dimensional certificate.

Build with:

```text
latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=build main.tex
```

`manuscript.pdf` is copied from the successful final `build/main.pdf`. `VERIFICATION.md` records actual coverage and build issues. `PROVENANCE.json` identifies source lineage; those historical paths are not runtime dependencies.

Human author names and affiliations are blank, and no human expert review has occurred. The paper records the operator-reported chronology and unresolved attribution statement without awarding authorship, ownership or publication priority.
