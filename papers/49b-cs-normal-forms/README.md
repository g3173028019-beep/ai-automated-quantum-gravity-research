# Exact Chern–Simons Normal Forms for Paired-Cap Tree Plumbings

Standalone follow-up to Paper 49. It proves exact analytic normal forms for the classical Chern–Simons action on closed SU(2) graph manifolds with two opposite order-two fillings per vertex.

- **Full tree:** exact flat Jacobian scheme for nonparallel seam axes and exact action at fully orthogonal configurations.
- **Perfectly matched spine:** an all-order quartic action with arbitrary positive terminal pendant counts.
- **Connected unmatched spines:** an all-order weighted normal form for a tree of centers with arbitrarily many ray kernel variables, including all star spines.
- **Leading local integral:** a one-pendant vertex gives a logarithm; extra pendants give finite coarea. The exact scalar order-polytope volume agrees with Paper 48's independently computed closed WRT volume. No quantum reference-density identification is made.

[PDF](manuscript.pdf) · [TeX entry point](main.tex) · [Verification](VERIFICATION.md)

## Build

Run `./build.sh` with TeX Live and pdfLaTeX installed (AMS packages, geometry, mathtools, microtype, and hyperref). The script builds in a temporary directory, resolves cross-references and checks for overfull boxes; it fixes the creation timestamp. It can use an unregistered TeX Live data tree through the optional `TEXLIVE_DIST` environment variable. The proof component files are `flat.tex`, `normal.tex`, `kernel.tex`, and `oscillation.tex`. References are embedded in `main.tex`; no external bibliography or network fetch is required.

## Related manuscripts

The [project index](https://github.com/g3173028019-beep/ai-automated-quantum-gravity-research/blob/main/docs/PAPER-INDEX.md) retains Papers 48 and 49. This follow-up does not replace their finite WRT asymptotics, H-plumbing density, APS comparison, Weyl regrouping or marked loop calculation. The [companion torsion and marked cutoff analysis](https://github.com/g3173028019-beep/ai-automated-quantum-gravity-research/tree/main/papers/48-49-wrt-torsion-bridge) addresses the separate physical density comparison.
