# 48 · Leading logarithms in closed SU(2) Witten–Reshetikhin–Turaev invariants of capped trees

This article studies the uncolored closed invariant for fixed tree plumbings with opposite order-two caps. It gives the exact H-plumbing calculation, a framed corona family with unbounded leading log degree, the even/odd law for paths of arbitrary length, and a positive leading-coefficient theorem for arbitrary unframed paired-cap trees. An appendix gives the odd-level integer-spin counterpart for the corona family.

The claims are about the specified closed invariants and normalizations. The manuscript does not claim general graph-manifold asymptotics, a full expansion for all levels and manifolds, or priority for the occurrence of logarithms.

## Source and build

- main.tex: entry point, abstract and bibliography
- introduction.tex, conventions.tex: introduction and finite surgery conventions
- h-tree.tex, h-proof.tex: H plumbing and exact cotangent reduction
- corona.tex: fixed-framing corona theorem and logarithmic estimates
- paths.tex, trees.tex: exact path transfer and tree face law
- consequences.tex, integer-spin.tex: no-log obstruction, presentation scope and SO(3) appendix
- manuscript.pdf: compiled article

Run **latexmk -pdf main.tex** from this folder with a standard TeX Live installation. No external figure assets are required.

## Reproducibility

The analytic proofs are in the article. The following optional checks replay exact finite formulas or illustrate convergence; they are not proof inputs.

- **python3 replay_h_tree.py** (standard library): literal Kirby and compact modular H sums, both level parities, and the leading logarithm's finite constant.
- **python3 replay_corona.py** (NumPy and SciPy): direct framed corona contractions and the positive- and zero-volume regimes.
- **python3 replay_paths.py** (NumPy and SciPy): original physical +2/−2 modular cap sums against the exact positive path transfer; tests the degree-two identity across its upper wedge.
- **python3 replay_tree_sewing.py** (NumPy): exact no-alias spin multiplicities and direct finite branching-tree contractions.
- **python3 replay_so3.py** (NumPy and SciPy): separate finite integer-spin and full SU(2) contractions, with the pointed factor.
