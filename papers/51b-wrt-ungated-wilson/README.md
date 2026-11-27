# Paper 51B. An Ungated Centered Wilson Law for Paired-Cap WRT Sewing

This is a standalone follow-up to Paper 51, which remains unchanged. The paper proves the ungated first surviving WRT coefficient for a single fixed integer-spin ribbon on the terminal pair of the four-vertex path. It also proves the literal color-ratio law for matched trees with fixed lower-center profiles. The classical normal-phase appendix is a separate comparison, not an input to the discrete coefficient.

- [PDF](manuscript.pdf)
- [LaTeX entry point](main.tex); inputs: transfer.tex, radial.tex, ratio-tail.tex, normal-phase.tex.
- Run ./build.sh with pdfLaTeX and the packages imported by main.tex for a deterministic rebuild.
- Run python3 verify_coefficients.py for the exact rational coefficient replay.

Internal review only. See the [verification record](VERIFICATION.md) and the [companion derivation supplement](../51-50b-wrt-cycle-gauge-discrimination/supplement/README.txt), including the general Hall-sector ratio law, the Gaussian constant and the finite-level replays.
