# 50. What Graph Structure Survives WRT Sewing?

**Hidden connector motifs and the stability boundary of dense cores**

This paper asks which graph changes can be distinguished by closed, oriented SU(2) Witten–Reshetikhin–Turaev invariants in one explicit class of torus sewings.

- **Fixed and growing square connectors.** An infinite family matches in its first three logarithmic coefficients and splits at the fourth. In a joint growth regime the first eligible contrast is a signed sum of positive leafless-motif weights. A diamond detects structure invisible to triangle, square and pentagon counts. A second pair reverses its actual modulus ordering as the number of twins crosses a multiple of `(log r)^(2/3)`.
- **Dense rigid cores.** After taking the level limit first, an exact vertexwise normalization removes any sublinear maximum-degree deletion from the asymptotic signed fluctuation multiplier under the stated strict template condition. At linear deletion, an explicit pair has the same local degrees, densities, site means, sine correction and exact oriented phase, but different Gaussian determinants.

The graph-size and level limits in these two parts are different. The paper does not assert arbitrary simultaneous dense-graph asymptotics or reconstruction of a graph from its WRT invariant.

- [PDF](manuscript.pdf)
- [LaTeX entry point](manuscript.tex) with [connector proofs](detection.tex) and [dense-core proofs](density.tex)
- [Exact finite certificates](replay/)
- [Proof and reproduction notes](VERIFICATION.md)

Build the PDF with `latexmk -pdf manuscript.tex`, or run `pdflatex manuscript.tex` until cross-references stabilize. Keep the three source files in the same directory.
