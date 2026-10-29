# Controlled quantum-weight networks

The complete English paper is `manuscript.pdf` (12 pages). Its editable source is `main.tex`, with the short included derivations `family.tex`, `network.tex` and `limits.tex`, and bibliography `references.bib`.

The paper combines the finite cell and chiral pair with the decisive network question. It proves the positive-weight family under the specified transverse update, gives an exact counterexample for a different noncommuting update, classifies every connected network with the stated independent controls and bond types, and separates the iid collective limit from operator-norm convergence and interacting-Gibbs fluctuations. The mean-field model is identified as a transverse Curie–Weiss model with spectator degrees of freedom.

## Reproduce

Use Python 3.10 or later with the packages in `requirements.txt`. From the packet directory:

```text
python code/verify.py
```

This reads `inputs/model.json`, performs exact SymPy and integer binary-symplectic calculations, and writes `outputs/verification.json`. It needs no historical project path, network, account, raw conversation or third-party paper. The same script also works when invoked by absolute path from a different working directory.

`outputs/established_verification.json` preserves the previously established output bytes. The default replay checks exact dimensions 10/45/4095/65535 for vector-bond chains and 10/255/4095/65535 for clock-bond chains. The input JSON describes the local controls and Pauli realization; the arbitrary-graph theorem is proved analytically in the paper, not inferred from a finite table.

## Build

With LaTeX, BibTeX and latexmk installed:

```text
latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=build main.tex
```

Copy `build/main.pdf` to `manuscript.pdf` after a changed build. `VERIFICATION.md` distinguishes executed checks, analytic results and material limits. `PROVENANCE.json` records source/output hashes and the contribution history; its historical paths are provenance identifiers, not runtime dependencies.

Human bylines, affiliations and PDF Author are blank. Human expert review has not occurred, and no global novelty or journal-acceptance certificate is claimed.
