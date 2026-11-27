# All-level WRT collisions on fixed graph-manifold backbones

[Manuscript (PDF)](manuscript.pdf) · [LaTeX source](manuscript.tex) · [Joint-limit section](joint.tex) · [Verification](VERIFICATION.md)

Varying the signs of section-fiber exchanges on a fixed bipartite graph leaves the standard oriented, uncolored SU(2) Witten-Reshetikhin-Turaev invariant unchanged at every integer shifted level \(r\ge 3\). The manifolds can still have different first Betti numbers, or different odd torsion groups at a common Betti number. The paper constructs both families, calculates the exact Maslov phase, and applies the known odd pointed invariant as an independent diagnostic. A separate theorem controls the simultaneous size and level limit of the signed finite sum on strict dense nonbipartite templates, at a sufficient exponential level.

This is a companion to [Paper 50](../50-wrt-graph-detection-and-stability/), whose motif discrimination, mesoscopic sign reversal, and linear-deletion counterexample are not replaced. The joint-limit theorem uses the earlier continuum saddle as an explicit input. It does not claim a polynomial or sharp level threshold.

To compile with pdfLaTeX and the standard AMS, microtype, geometry, hyperref and xcolor packages, run:

~~~sh
./build.sh
~~~

The script runs three passes with a fixed source epoch. To check the finite examples with Python 3.8+ and SymPy:

~~~sh
python replay.py
sha256sum -c SHA256SUMS
~~~

The finite checks supplement, rather than replace, the analytic proof.

## Companion papers

[51+50B](../../papers/51-50b-wrt-cycle-gauge-discrimination/README.md) presents even-cycle and finite-gauge separation results.
