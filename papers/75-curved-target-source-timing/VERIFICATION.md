# Verification: Paper 75

The manuscript and all proofs are in [main.tex](main.tex). The result is analytical: no sweeps, numerical sign checks or scripts are required.

1. **Exact background and operators.** Differentiate the Friedmann relation to obtain the stated time and spatial Ricci components; use the geodesic `y=0` to verify the wave-map equation. Expand `A²` at zero to verify the curvature and its transverse derivative. Check that the transverse quadratic action has the stated sign and that the synchronous spatial TT system is the constrained tensor sector of the same coupled stock.
2. **Complete Ward accounting.** Differentiate the one marked sourced action three times in distinct drive amplitudes, retaining the fourth action derivative and all Hessian and third-derivative contacts. Lower Euler coefficients cancel the variations of the gauge derivative. The early curvature tag has zero first and second jets. Reflection, causality and field type, rather than an unsourced on-shell identity, leave the single tagged-stress chain.
3. **Strict nonvanishing.** If every tested Laplace eigenpair gave zero stress pairing for a nonzero TT tensor, the corresponding second-order spatial differential operator would have zero principal symbol. One finite eigenpair with the correct sign and the short-time scalar Green sign suffice for strictly positive widths.
4. **Conserved cut and actual reachability.** Green's identity gives the symplectic cross-cut term, and propagating a cut through a source includes the remaining same-side pieces. The quadratic identity for `Ω(v, M_cv)` is strictly positive for every nonzero input `v`; this makes both physical Cauchy data reachable by actual source variations.
5. **Same-label inverse.** The difference of translated drives is the Euler Hessian applied to a compact variation. After the strict support and contact check, the scan derivative is convolution by the reflected fixed ramp derivative. Its complete support is inside the measured interval. The Fourier multiplier is nonzero near the origin, and compact support continues equality to all frequencies. The calibrated weight is positive on the unknown support.
6. **Reference and scope.** Paper 71 supplies the compact positive physical Hadamard covariance for the same coupled stock. It does not supply nonlinear quantum composites or interacting positivity. The leading connected tree expression has one Green leg between separated temporal bands and requires no coincident loop renormalization.

Build the standalone document with a standard pdfLaTeX installation containing Latin Modern and microtype:

```sh
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex
```

The generated `main.pdf` is distributed as [manuscript.pdf](manuscript.pdf).
