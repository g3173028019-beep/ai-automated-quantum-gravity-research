# Half-Turn Covariance of a Tangential Fuzzy-Sphere Oscillator and Its Sector-Resolved Planar Limit

The six-page `manuscript.pdf` is a self-contained mathematical paper. `main.tex` and `references.bib` are its complete sources.

It proves the exact finite self-dual spectrum, right-half-turn parameter inversion, a uniform error on a fixed polar occupation square, and the four-copy full heat-trace limit. It identifies persistent cross-pole eigenvectors, derives the Weyl-symbol parity/symplectic-Fourier relation with explicit normalization, and compares the tangential kinetic family with Ydri's earlier fuzzy harmonic regulator.

The results concern a complex finite matrix Hilbert space and its stated spectral limits. They do not establish a real-field symmetry, an interacting renormalized continuum theory, a global rotation-covariant Fourier transform, or an operator-topology contraction. The prior fuzzy oscillator framework is credited; no broad first-action or publication-priority claim is made.

## Reproduce

From this packet directory, with Python 3.11 or later and SymPy installed:

```text
python verification/check_fuzzy.py
```

The script writes `verification/checks.json`, checking four exact spin representations and symbolic oscillator parameters, plus heat-trace diagnostics. It needs no network, old project directory, credentials or external input data. The analytic proofs establish the limits; finite numerical sequences only illustrate them.

Build the PDF with a standard TeX Live installation:

```text
latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=build main.tex
```

The delivered `manuscript.pdf` is copied from the successful final `build/main.pdf`. Build intermediates may be regenerated. `VERIFICATION.md` records the actual execution and limitations; `PROVENANCE.json` records source lineage and SHA-256 identities. Historical evidence paths appear only in provenance and are not runtime dependencies.

Human author names and affiliations are blank. The operator reports earlier Fable 5-led work, predominantly GPT-produced results from July 2026 onward, and DeepSeek/multiple-GPT cross-review and preparation. No human expert review has occurred. The paper preserves the operator's unresolved attribution statement.
