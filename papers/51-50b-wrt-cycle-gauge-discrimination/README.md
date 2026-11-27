# Cycle logarithms and gauge discrimination in paired-cap sewing

Companion article labeled **51+50B**. It preserves both original papers and is independent of companion 51B's ungated Wilson estimate. Start at main.tex; the four other TeX files are included automatically. The compiled PDF is manuscript.pdf. [PDF](manuscript.pdf) · [TeX](main.tex) · [Verification](VERIFICATION.md)

## Contents

1. For every fixed even cycle, exact signed SU(2) transfer and a positive convolution trace produce the leading coefficient in physical normalization, with an explicit $O_m(r^{2m})$ remainder at either shifted-level parity. The logarithmic relative-color law is for the reduced color sum, not a crossed-ribbon observable.
2. For a different backbone, the signed $K_{4,4}$ collision retains the original WRT and all block-interior comparison experiments but changes integral homology. Every odd pointed theory and every twist for a finite odd nilpotent gauge group agrees on this pair.
3. The order-21 semidirect group $\mathbb F_7\rtimes_2 C_3$ distinguishes the pair, including all three twists. Section 5 gives exact rank enumeration and proves minimality among odd group orders. The order-16 abelian detector is a separate even-order comparison.

## Rebuild and exact checks

- Run ./build.sh with a standard pdfLaTeX installation (geometry, lmodern, microtype, AMS, mathtools, booktabs, xcolor, hyperref). An optional TEX_ENV names a local TeX setup script. The build runs three passes, disables variable trailer IDs and uses a reproducible timestamp.
- Run python3 check_order21.py. It uses only the standard library and checks every one of the 39,366 rank pairs using both modular elimination and anchored $3\times3$ minors. Its frozen output is order21-replay.txt.
- Run python3 check_signed_homology.py if SymPy is installed. It independently compares the exact $24\times24$ relation presentation with $4A_\epsilon$ and reproduces the cyclic counts. The expected exact data are in signed_homology_check.json.

## Sources and scope

The earlier working source supplied the even-cycle and integral/odd-nilpotent derivations, preserved here in reorganized form. The published predecessor source blob IDs checked read-only are Paper 50 manuscript.tex (ef60610477c3c0461c2ac6aad00c32633c978ca9) and Paper 50B manuscript.tex (7491ac507620524fe0f905d376e62c7e6a3879fd). The order-21 calculation is an additional proof, with two exact rank reductions in the companion replay. No result requires a numerical asymptotic fit, a universal closed-WRT reconstruction claim, or the Wilson theorem in companion 51B.

## Supplement

The [complete earlier derivation note](supplement/wrt_extensions.pdf), [LaTeX source](supplement/wrt_extensions.tex), [replays and inputs](supplement/README.txt) preserve the Hall-sector ratio law, angular normal-model details, finite cycle diagnostic, and variable-cap and arbitrary-depth examples. Paper [51B](../51b-wrt-ungated-wilson/README.md) supersedes the note’s open status for the specific ungated terminal Wilson coefficient.
