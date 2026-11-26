# 49 — Chern–Simons Phases and Degenerate Sewing in SU(2) WRT Invariants

This article calibrates classical action and absolute phase against exact finite-level WRT sewings. It derives the local density of the H-plumbing logarithm from surgery, computes the absolute adjoint APS phase for opposite-cap graphs, and proves an explicit central-fibre collision where Weyl grouping and the separate central charge coefficient differ. A marked tree and triangle comparison isolates the contribution of a loop holonomy.

- [Compiled manuscript](manuscript.pdf)
- [LaTeX entry point](main.tex); section sources: [cubic](cubic.tex), [APS and torsion](aps.tex), [Weyl endpoint](weyl.tex), [marked examples](marked.tex)
- [Proof and build verification](VERIFICATION.md)
- [Independent finite-sum replay](replay.py), [recorded output](replay-output.txt)

Build in standard TeX Live:

    latexmk -pdf main.tex

Replay with Python 3, NumPy and SciPy:

    OPENBLAS_NUM_THREADS=1 python replay.py

The replay supplies no analytic proof input. The triangle in the final section is a Wilson-marked amplitude, not the unmarked invariant. There is no assertion of a general path-integral measure or a complete stationary-phase expansion.
