# Verification

## Exact replay

```sh
python3 replay_tv_native_cylinder.py
python3 replay_tv_torus_symbols.py
```

The cylinder replay checks, in exact integer or rational arithmetic:

1. The vacuum rank below the truncation wall and the one-excursion first deficit (tested at levels 1–40; the unrestricted-level argument is in the manuscript).
2. The parity-matched sector closure from the diagonal fundamental seed and exact half-cut probabilities at level two.
3. The actual doubled level-two projector identity \(P_{00}Q_{00}P_{10}Q_{10}P_{00}=-P_{00}/16\) in the normalized fusion coordinates.

The torus replay checks finite-dimensional phase and odd-compression conventions for the physical charge basis and the Weyl product-to-sum rule. It does not numerically certify the infinite-level support/count theorem, which is proved by the norm-exclusion, trial-packet and trace arguments in the article.

## Mathematical audit

A separate internal pass checked the use of the normalized physical cylinder, the multiplicity ranks and Perron normalization, the doubled cut projectors and their ordered phases, and the finite simultaneous unitary separation argument. The torus operator phases, missing odd endpoints, empirical measure and absence of asymptotic lower spectral outliers were checked against the stated finite representation and proof. The candidate spectral pencils are supplied as coherent linear combinations of actual physical maps; an experiment using only independently dephased histories would not implement those pencils.

The journal and preprint references were checked by title and persistent link. The compiled PDF was checked for unresolved references, compilation errors and visual layout, including all equation pages, the signed-cycle figure and the bibliography.

## Separation from other packets

This paper concerns open physical cylinders, punctured-boundary recoupling, and torus Wilson constraints in compact Turaev–Viro theory. It does not repeat the four-dimensional EPRL construction and Pachner rank obstruction of Paper 38, the real-tree and \(p\)-adic stochastic results of Paper 39, or the separate closed-manifold surgery manuscript. It does not assert priority for the standard categorical dimensions, the quantum torus skein relations, or matrix trace separation.

All checks were internal. No external human expert review is claimed.
