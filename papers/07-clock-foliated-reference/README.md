# A local clock-foliated repair of the inhomogeneous spectrum of a rotating Einstein-static scalar reference

The specified effective action preserves the exact reference and homogeneous action while making every physical inhomogeneous quadratic sector positive. It constructs a 141-oscillator finite-band Gaussian. Nonlinear stability and a continuum Hadamard state for the modified characteristic cones are not claimed.

The completed English manuscript is manuscript.pdf (11 pages), with editable main.tex and references.bib. Its bytes match the completed and visually inspected research edition. Human bylines and affiliations are blank.

## Reproduction

From this packet, after installing requirements.txt:

~~~powershell
python scripts/construct_inhomogeneous.py
~~~

The commands require no historical project directory or network. Results are written under results/. The Maxwell packet bundles its accepted switching/phase input; the other calculations reconstruct their inputs algebraically. Saved outputs distinguish exact symbolic checks from numerical evaluations and inherited analytic premises.

## Manuscript build

~~~powershell
latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=build main.tex
Copy-Item -LiteralPath build/main.pdf -Destination manuscript.pdf -Force
~~~

DERIVATION.md and VERIFICATION.md give the calculation, executed coverage and material limits. PROVENANCE.json identifies the completed source edition and its hashes. This delivery preserves the scientific scope; it is not an additional research or novelty-review round.

Earlier work was Fable 5-led; July 2026 onward results were predominantly GPT-produced, as reported by the operator. DeepSeek and multiple GPT models also participated in review and preparation. Human expert review has not occurred. Local research dates and this private repository delivery do not establish public priority.

> I genuinely do not know how authorship, intellectual credit, or ownership should be attributed for these results.
