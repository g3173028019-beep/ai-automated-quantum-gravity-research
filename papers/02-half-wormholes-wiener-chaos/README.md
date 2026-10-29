# Half-wormholes as Wiener chaos

**Correlated-replica spectroscopy, exact Pfaffian sectors and boundary moment indeterminacy**

Revised edition: **29 October 2026**, 10 pages. [Manuscript](manuscript.pdf), [LaTeX source](main.tex), [revision record](REVISION.md).

The paper studies the one-time-point, zero-dimensional Gaussian-coupling polynomial. It retains the replica identities, closed q=2 sector formula, correlated antisymmetric determinant formula, finite probabilistic truncation bounds and moment-law dichotomy at degree two versus three. The exact general half-wormhole expansion and variance sum are credited to prior work. Slud's odd-degree result and the unresolved even-degree originality comparison are stated explicitly.

Supporting results retained from the longer source include the even-q fourth-chaos formula, an elementary Gaussian Pfaffian convolution identity, and an exact four-node quadrature calibration. The odd-q calculation is a formal sign-convention control, not an odd-q Grassmann partition-function theorem. Finite quadrature exactness does not identify a measure in an indeterminate problem. Canonical N-extremal discrete measures form only a subfamily of the full Nevanlinna parametrization.

Run both short portable checks with Python >=3.10; only the standard library is needed:

```powershell
python verification/check_revision.py
python verification/check_chaos_extras.py
```

They write `verification/revision-checks.json` and `verification/chaos-extras-results.json`. Checks use exact integers, fractions and Q[sqrt(13)] arithmetic; no Monte Carlo, network or historical project path is required. The first checks 90 reductions of the prior variance sum, 12 fourth-moment sums, n=1 empty-sector behavior and direct convolution at n=1,2,3. The second checks the even/odd completion-sign controls and quadrature moments through order seven.

Build with TeX Live, latexmk and BibTeX:

```powershell
./build.ps1
```

`main.tex` and `references.bib` are self-contained. `verification/legacy_checks.wls` preserves the unchanged historical supplement; its substantial loops and Monte Carlo experiments are optional and were not rerun for this revision.

Human author/affiliation fields and PDF Author remain blank. Human expert review has not occurred. Historical source dates, this revision and demonstrable public disclosure remain distinct. No full quantum-gravity construction, comprehensive novelty certificate, settled ownership or journal acceptance is claimed. See [verification](VERIFICATION.md) and [provenance](PROVENANCE.json).
