# Verification

Edition: 29 October 2026, Beijing. Additional AI review; no human expert review.

## Analytic results and focused joining check

The text supplies the factorial-generating-functional proof of the Poisson-type MGF bound, entire scalar MGF and concentration inequality. The Cox construction and its distinct-law Laplace proof are explicit; its extension to all mixed-temperature moments uses the same finite conditional-polynomial argument. A two-level deterministic example separates scalar-law determinacy from reconstruction.

The trace-channel section includes the ordered-clock gap identity, a Schatten-norm proof of the p>1 sufficient condition, and a pure rank-one-block projection counterexample at p=1. The synthesis agent independently checked this dilation and the full-joint-clock reconstruction proof. The latter explicitly includes finite mixed moments, factorial partition inversion, finite exponential tilting, and absolutely convergent local Laplace functionals. It does not assume positive exponential moments for unbounded Airy heat weights.

The Bernoulli example proves unique law does not imply pure state. The semigroup lemma is an analytic conditional reduction, with Gaussian inverse-coordinate/slice coverage and literal Airy density bounds left unproved.

## Portable computations actually run

Both commands exited successfully:

```text
python verification/check_bounds.py
python verification/check_channels.py
```

Runtime: Python 3.11.15, SymPy 1.14.0, SciPy 1.17.1.

`checks.json` records exact signed integer moments for the lognormal Stieltjes family, exact same-clock/two-clock rational weights, finite independent-Bernoulli determinantal MGF examples including a negative but uniformly bounded lower energy, and Cox Laplace quadratures at two tolerances. The Laplace differences are numerical illustrations; analytic uniqueness proves distinct laws. Quadrature diagnostics are not interval bounds.

`channels.json` records cumulants 1–4 from the Bernoulli log-mgf, exact 121/4 versus 10 for the unrestricted-clock counterexample, symbolic idempotency/trace of a pure rank-one block, the positive Bernoulli GNS Gram determinant, and scaled finite sums illustrating the projection p=1 channel limit. No finite decimal computation is used to prove full process reconstruction or a channel limit.

The bounds script derives from the already verified recovery producer. Only its report's prose pointer was changed to refer to the accompanying manuscript; mathematical inputs and calculations were preserved. `reference-bounds.json` retains the original saved output bytes and their source hash in provenance. No old source path is imported at runtime.

## References, compilation and remaining scope

Primary metadata were checked for Soshnikov, Lyons, Heyde, Janssen--Mirbabayi--Zograf, and Ramirez--Rider--Virag. Heyde's title, author, volume25/issue2/pages392–393 and DOI were confirmed directly through the Crossref primary metadata API after the DOI page was unavailable in the browser tool. The classical tools are credited without a global novelty claim.

`latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=build main.tex` completed with TeX Live 2025, pdfTeX 1.40.27, latexmk 4.86a and BibTeX. The PDF has **10 pages**. The final TeX log has no undefined references/citations, LaTeX warnings, overfull boxes or underfull boxes. The Perl launcher gave a locale fallback warning outside the TeX log; the build succeeded. PDF Author is blank.

No stochastic-Airy simulation, quantitative coarea theorem, physical gravity-completion identification, evaporation dynamics, or generic projection sign-rigidity proof is part of this paper. Parent integration owns final visual QA and upload; worker compilation does not claim those steps. Original evidence and publication checkout were preserved.
