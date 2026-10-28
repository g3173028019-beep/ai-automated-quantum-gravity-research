# Focused closure verification

The central theorem is unchanged from the completed round-one derivation. This closure supplies its full Gaussian comparison and domain argument in the manuscript and removes old-project dependencies from the replay code. No new optimality or instrument programme was started.

The portable command `python code/reproduce.py` passed. It freshly reproduces:

- Gaussian moment-generating-function identities for `M(v)` and `F(v)`, their leading coefficients `80/3` and `1`, and the real/imaginary stability ODE and Wronskian for `z(t)`.
- Exact charge polynomial identities, monotonic weighted quotients, and the stationary cutoff/leading coefficient.
- Charge columns through input `N=100`, retaining two upper-edge raising steps; 180 coherent matter/reference vectors for both first-moment and second-moment bounds. These are reused passing check structures with a local fresh run, not new independent review.
- Finite formula evaluation with 150-decimal mpmath interval arithmetic at `C=10^6,10^12,10^18,10^24`, plus the separate second-moment examples and coherent-tail hostile controls. Large analytic cutoffs are never used as simulation dimensions.

`results/reproduction.json` contains the actual outputs and code hash. No old source or inherited numerical certificate is read by the default command. The inherited theorem inputs are now derived or specified in the article. `PROVENANCE.json` records their historical source hashes as provenance only.

The article compiled successfully with `latexmk -pdf -outdir=build main.tex`, including bibliography and resolved cross-references, into a six-page PDF. The human author and PDF author metadata are empty. Integrated visual inspection covered all six pages in a contact sheet and the central theorem and disclosure pages at readable resolution, with no clipping or overlap found. The local TeX wrapper emitted an environment locale warning and an amsmath advisory about primitive fraction notation; neither is a compilation error.

Not verified: a full interacting JT evolution simulation, experimental input preparation, ADM energy constraints, exact-JT optimality, all-time or arbitrary-Fock limits, formal proof-assistant certification, independent human expert review, or global priority/novelty. The manuscript explicitly limits these claims.
