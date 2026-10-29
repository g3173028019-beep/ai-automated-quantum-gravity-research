# Positive physical Hadamard states for critical Einstein-scalar multiplets

The independent fourteen-page manuscript proves the scoped physical state-existence theorem for linearized canonical Einstein gravity coupled to a finite fixed positive scalar multiplet on smooth on-shell four-dimensional globally hyperbolic backgrounds with compact Cauchy surface and critical potential along the background image. Nonconstant flat-direction backgrounds and a nonzero compatible potential Hessian are included.

This is distinct from the clock-foliated effective-model paper. The proof explicitly derives the common-normalized coupled operators, compatible projections, physical positivity, regular and singular Fredholm correction, physical quotient and spacetime Hadamard lift. The pseudodifferential, elliptic and normally hyperbolic analytic inputs are cited precisely.

## Read and build

manuscript.pdf is the completed English paper. main.tex includes operator.tex, fredholm.tex and example.tex, with references.bib.

~~~powershell
latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=build main.tex
Copy-Item -LiteralPath build/main.pdf -Destination manuscript.pdf -Force
~~~

## Reproduce the supporting algebra

Install requirements.txt, then run:

~~~powershell
python code/canonical_algebra_check.py
~~~

The command prints a JSON report and requires no network or historical files. The recorded run is results/algebra-results.json. It checks normalization, charged-adjoint identities, the finite-rank charge defect and the nonconstant on-shell example. It does not replace the analytic Fredholm or microlocal proof. VERIFICATION.md and PROVENANCE.json record scope and source history.

Human names and affiliations remain blank. Earlier work was Fable 5-led; July 2026 onward results were predominantly GPT-produced, as reported by the operator. DeepSeek and multiple GPT models also participated in review and preparation. No human expert review or public-priority certification is asserted.

> I genuinely do not know how authorship, intellectual credit, or ownership should be attributed for these results.
