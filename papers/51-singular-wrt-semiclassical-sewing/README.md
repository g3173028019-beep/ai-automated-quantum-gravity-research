# Singular semiclassical sewing for paired-cap Witten–Reshetikhin–Turaev invariants

[PDF](manuscript.pdf) · [LaTeX entry point](main.tex) · [Verification](VERIFICATION.md)
One English manuscript on fixed untwisted SU(2) paired-cap trees with at least two backbone vertices. The central result combines an action-resolved all-orders power–log expansion for the exact radial invariant with an independently normalized classical leading coefficient and an ungated leading joint Wilson law for the stated dormant routes.

## Build

Run `./build.sh` with pdfLaTeX, Latin Modern, geometry, AMS packages, microtype, xcolor, and hyperref. The manuscript uses inline references and has no external figures or bibliography dependencies. On successful compilation the script leaves `manuscript.pdf`.

## Verification

See [Verification](VERIFICATION.md). The independent finite-sum replay uses `python verify_h_tree.py` with `mpmath`.

## Source map

- `main.tex`, `intro.tex`, `conclusion.tex`: shared conventions, motivation, synthesis.
- `geometry.tex`: flat-character geometry, matching and relative scale charts.
- `radial.tex`, `phases.tex`: exact lattice expansion, phase filters and support.
- `classic.tex`, `angular.tex`, `hall.tex`: normalized Chern–Simons local action, finite Wilson transfers, and absolute Hall sewing.

Individual higher angular coefficients and unrestricted active-skeleton words are open.

## Companion papers

[51B](../../papers/51b-wrt-ungated-wilson/README.md) develops the ungated Wilson law. [51+50B](../../papers/51-50b-wrt-cycle-gauge-discrimination/README.md) presents even-cycle and finite-gauge separation results.
