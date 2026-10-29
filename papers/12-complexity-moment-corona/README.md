# Complexity-Weighted Moments, Critical Mass Escape, and Selected Corona States

This is a complete mathematical-methods article with a portable reproduction package. It treats fixed-regulator determinacy, the critical-denominator alternative, justified moment rates, corona states, the word/total ideal distinction and actual Abel-regulator reachability. Classical ingredients are credited; no new general moment theorem or physical gravitational completion is asserted.

Human authors and affiliations are blank. The compiled PDF also has a blank Author field.

## Contents

- `main.tex`, `references.bib`: complete article source and six references.
- `manuscript.pdf`: compiled article; `build/` contains the compilation log and auxiliary files.
- `code/verify_constructive.py`: exact geometric-formula numerical illustrations, using only the Python standard library.
- `saved-output/calculation-output.json`: preserved recovery output, copied without changing its bytes.
- `reproduced-output/calculation-output.json`: fresh output from the portable script.
- `VERIFICATION.md`: actual proof, reference, build and portability coverage and limitations.
- `PROVENANCE.json`: source identifiers, hashes, version merge, chronology and artifact hashes.

## Reproduce

From this directory:

```text
python code/verify_constructive.py
latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=build main.tex
```

The script prints the calculation and writes `reproduced-output/calculation-output.json`. It needs no input data, network or historical path. Its finite values illustrate proved bounds; they do not establish asymptotics by extrapolation. To update the top-level PDF after a rebuild, copy `build/main.pdf` to `manuscript.pdf`.

The article has no figures or external TeX includes. It builds with ordinary TeX packages listed in the preamble. Python dependencies are listed in `requirements.txt`; no third-party Python package is required.

## Scientific scope

The text gives the central definitions and proofs directly. An exact heat normalizer supplies a safe translated-clock argument. The critical series distinguishes an ordinary probability limit from escape; two-sided residual bounds give order estimates, while a residual limit gives a constant and an exponential scaling law. The further word ideal corrects a pullback failure caused by fixed-word microlocal escape. Parity and long blocks distinguish ambient from genuinely reachable state separation.

The principal-series density is imported from the cited Virasoro-TQFT source. No physical critical asymptotics, common sewing/dynamics quotient, physical witness or selected gravitational state pair is supplied. Thermal point-process statistics are a separate article/question.

Earlier work was Fable 5-led; July 2026 onward results were predominantly produced by GPT models according to the operator. DeepSeek and multiple GPT models also contributed review/preparation. Human expert review has not occurred.

"I genuinely do not know how authorship, intellectual credit, or ownership should be attributed for these results."
