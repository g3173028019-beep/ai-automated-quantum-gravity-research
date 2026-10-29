# A local Einstein–scalar BV/BRST benchmark

Complete English source-specialized benchmark paper. It gives the local action and convention dictionary, full nonminimal Green matrix for every real auxiliary parameter including zero, Ward-repair/normalization transport, free equicausal algebra, an ordinary-dagger obstruction, the first matter BV subcocycle and the explicit interaction stop line.

The reduced physical Hadamard state is an imported input. The separate canonical Einstein–scalar physical-state theorem is cited, not reissued as a new state-existence theorem. The benchmark's nonminimal, adjoint and interaction questions remain separate.

Run from this directory:

```text
python reproduce.py
latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=build main.tex
```

Python requires SymPy; see `requirements.txt`. The four programs in `code/` run offline and write their algebraic evidence to `outputs/`. Default reproduction needs no historical directory, source papers or model provider. `manuscript.pdf`, local TeX/bibliography, provenance and verification records complete the packet.

The finite programs check 49 printed coefficient/sign rows. They do not prove analytic support, Fredholm or wavefront inputs, the full nonlinear interaction, an anomaly cancellation or QME. Those distinctions are stated in the paper. Human byline, affiliation and PDF Author remain blank. No human expert review or comprehensive novelty certification occurred. Earlier work was Fable 5-led; July onward predominantly GPT-produced according to the operator, with DeepSeek and multiple GPT review/preparation contributions. Full chronology and unresolved attribution appear in the manuscript and provenance.
