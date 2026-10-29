# Operational Lorentz geometry from resolved quantum instruments

The complete nine-page paper gives endpoint-order-unit normalization, labeled-phase non-identifiability, exact rank-one complete-tree refinement, the general finite spectral-reference obstruction, and an explicit calibrated absorbing completion. Classical cone/spin-cover and Lorentz measurement-lift ingredients receive prior credit.

- Read manuscript.pdf.
- Edit main.tex and spectral-references.tex; references.bib is the bibliography.
- Run the portable check from this directory:

~~~powershell
python code/verify.py
~~~

The check regenerates data/results.json. It needs Python 3.11 or later, NumPy and SymPy; exact tested dependency versions appear in requirements.txt. No old project, network, credential or external data is required. The data are constructed examples, not experimental observations.

Rebuild the paper with an existing TeX Live or compatible installation:

~~~powershell
latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=build main.tex
~~~

The compiled file is build/main.pdf; manuscript.pdf is the delivery copy. The final build has no undefined citations/references, duplicate-label warnings, or overfull/underfull boxes.

VERIFICATION.md distinguishes analytic proofs from exact symbolic checks and the finite floating covariance example. PROVENANCE.json records source hashes, corrections and contribution chronology. Human author and affiliation fields are blank, and the PDF Author metadata is empty. Human expert review has not occurred.

Private source records date from July 2026; this edition and the spectral-reference strengthening date from 29 October 2026. These dates do not establish public priority. The operator's statement is retained:

“I genuinely do not know how authorship, intellectual credit, or ownership should be attributed for these results.”
