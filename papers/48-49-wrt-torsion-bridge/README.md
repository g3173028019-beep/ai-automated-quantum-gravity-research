# Adjoint torsion and a coupled logarithmic cutoff for paired-cap WRT invariants

[Manuscript PDF](manuscript.pdf) · [Verification](VERIFICATION.md)

This manuscript compares the exact finite-level SU(2) Witten–Reshetikhin–Turaev sum with adjoint torsion on the flat-connection branches of unframed paired-cap trees. The fixed macroscopic Wilson-angle limit gives a Haar sewing measure. The true square-root adjoint torsion gives the same radial density up to a fixed normalization, including for covers with central–central seams and their SO(3) bending.

In the four-backbone-site path, independently truncating this smooth torsion measure gives twice the closed leading logarithmic weight. The finite quantum Wilson-radius marks and the reduced Chern–Simons quartic, with its physical branch marks projected to the chart, give the same triangular cutoff. The result does not construct the absolute quantum integration density or its phase by a path-integral argument.

This is a joint extension of [Paper 48](../48-wrt-leading-logarithms) and [Paper 49](../49-wrt-semiclassical-sewing). The distinct [normal-form companion](../49b-cs-normal-forms) proves more general statements about the actual reduced action; the present manuscript uses the four-site quartic and focuses on physical radii, torsion and exact quantum comparison. Neither original manuscript is superseded.

## Build

Run **./build.sh** with pdfLaTeX and the standard LaTeX packages used by main.tex (Latin Modern, geometry, AMS, mathtools, microtype, xcolor and hyperref). The build runs twice, checks references and page-box warnings, and leaves only the final manuscript.pdf. With the same TeX engine and fonts it is byte-reproducible across build directories. The inputs are main.tex, quantum.tex, torsion.tex, regulator.tex and p4.tex. No external figure files, bibliography database or shell escape are required.
