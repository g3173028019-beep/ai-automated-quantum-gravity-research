# Positive boundary systems and completely positive sewing

Complete English mathematical framework/method paper, combining a typed boundary constitution with doubled-matrix and finite Weyl/charge/cut calibrations. It contains definitions, GNS and quotient scope, observer restriction, a corrected scalar-versus-amplified positivity criterion, exact finite constructions/proofs, counterexamples and numerical certificates. It makes no continuum or gravitational-source identification.

- `manuscript.pdf`: compiled paper.
- `main.tex`, `references.bib`: complete local LaTeX source.
- `code/`, `reproduce.py`, `outputs/`: portable finite evidence.
- `PROVENANCE.json`, `code-provenance.json`: source identities and derivative-code transformations.
- `VERIFICATION.md`: actual coverage and limits.

Run from this directory:

```text
python reproduce.py
latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=build main.tex
```

Python requires NumPy; `requirements.txt` lists the requirement. Default reproduction uses no original research path, network or model provider. It writes finite stdout/receipts to `outputs/`. Two historical producers were retained with descriptive labels/output names; their formulas, parameters and negative controls were preserved. The independent diagnostic adds the transpose Choi, restriction, compression and Weyl-dictionary controls.

The paper credits standard operator-algebra/quantum-operation results and claims the bounded compatibility/calibration package. Scalar state Grams do not certify complete positivity. The finite matched-charge quotient uses its declared block algebra. Human byline, affiliation and PDF Author remain blank. No human expert review or global novelty certificate is claimed. Earlier work was Fable 5-led; July onward predominantly GPT-produced according to the operator, with DeepSeek and multiple GPT review/preparation contributions. The full attribution statement and exact operator quotation appear in the manuscript and provenance.
