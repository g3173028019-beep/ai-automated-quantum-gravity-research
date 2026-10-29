# Positive Cartan carriers and signed Lorentz selectors

The complete ten-page paper defines an explicit minimal Clifford carrier, proves its positive-selector obstruction including the PSD boundary, and gives the restricted action/invariant and conditional Einstein–Cartan bridge with exact portable certificates. Four-component spinors and the polynomial action family receive classical prior credit.

- Read manuscript.pdf.
- Edit main.tex, selector-and-certificates.tex and certificate-appendix.tex.
- references.bib contains the cited bibliography.
- From this directory, reproduce the compact exact certificates:

~~~powershell
python code/verify.py
~~~

The program regenerates data/results.json from definitions. No stored NPZ matrix, historical project path, network or external input is needed. Python 3.11 or later is required; requirements.txt records the tested NumPy and SymPy versions.

Rebuild with an existing compatible TeX installation:

~~~powershell
latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=build main.tex
~~~

The output is build/main.pdf; manuscript.pdf is the delivery copy. The final build has no TeX warnings or box warnings. VERIFICATION.md explains the algebraic rank interpretation, including reference-signature complexification. PROVENANCE.json retains original hashes and the repaired exactness history.

Human author and affiliation fields are blank and PDF Author metadata is empty. Human expert review has not occurred. Local July 2026 source dates and this 29 October 2026 edition do not establish public priority.

“I genuinely do not know how authorship, intellectual credit, or ownership should be attributed for these results.”
