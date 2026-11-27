# Verification

## Claims and proof checks

- **Exact model:** Ordinary unitary Lorentzian EPRL transfer, one gauge-fixed group node in each simplicial vertex, normalized compact edge identity, and exactly one factor (2s+1) for each closed face. The compact spin-zero leg is kept as a spherical principal-series transfer. The main proof integrates the exact Lorentzian coefficient, not a finite booster or a stationary-phase approximation.
- **Absolute convergence:** The centered beta-density overlap bounds the minimal magnetic boost coefficients in modulus before summing. The Bernstein identity, followed by a one-dimensional integral bound, yields the uniform trace and Hilbert–Schmidt estimates. An independent Hilbert–Schmidt/weighted-spanning-tree estimate gives the same shell exponent; the manuscript presents the direct trace-norm/hyperbolic integration proof.
- **Colour counting:** Admissibility along edges of the connected simplicial link permits only a uniformly bounded number of spin offsets per minimum spin. The high slots at a vertex contract in a single triangle. Partial traces of compact intermediate-spin projectors supply three (O((1+n)^{-1/2})) factors per vertex without a leftover unconstrained trace dimension. The local integral checks both collision and large-rapidity regions. Multiplying by the standard face weights gives the claimed summable shell estimate.
- **Refinement comparison:** The zero-tuple Gram lemma is the fixed-colour result of Paper 69; the new rank statements concern full interior sums for modified local weights. Ordinary-weight endpoint rank and the sign of the fully integrated kernel were not verified. The (A_4) test uses specified common boundary marks and does not imply a rank lower bound on its own. The observation about positive compact measures concerns a single native compressed face, not the coupled seam.
- **Literature and independence:** The analytic proof and its model, normalization, rank and endpoint scope received an independent internal technical review. This is not external peer review; the manuscript does not assert literature-wide priority or a continuum or triangulation-invariance result.

## Rebuild

With a standard TeX Live installation (pdfLaTeX, Latin Modern, AMS, geometry and hyperref):

```sh
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex
```

The packaged file is the final build of `main.tex` renamed `manuscript.pdf`. No external figures, bibliography database or numerical data are needed. Page rendering and link/reference checks supplement but do not replace the proof.

## Small exact sanity check

The following Python 3 check verifies the tetrahedral incidence, shell exponent, tetrahedral character multiplicities and the signed compact example. It does not evaluate the Lorentzian integral or certify the unperturbed scalar.

```python
from fractions import Fraction as F
from itertools import combinations

v = range(4)
edges = set(combinations(v, 2))
faces = list(combinations(v, 3))
assert all(sum(set(e) <= set(f) for f in faces) == 2 for e in edges)
assert F(len(v)) - F(3, 2)*len(faces) == -2
assert (16 + 3*4 + 8*(-2)) // 12 == 1  # apex-half invariants
assert (2 + 3*2 + 8*(-1)) // 12 == 0   # coarse-seam invariants
# Phases in units of pi/3, ordered from low to high endpoint.
a = {(0, 1): 0, (0, 2): 0, (0, 3): 0,
     (1, 2): -1, (1, 3): 1, (2, 3): -1}
angles = [a[i, j] + a[j, k] - a[i, k] for i, j, k in faces]
chi_half = [2, 1, -1, -2, -1, 1]  # exact 2 cos(k pi/3)
assert [chi_half[t % 6] for t in angles] == [1, 1, 1, -2]
```
