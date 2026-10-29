# Verification of the forward constructions

Executed 29 October 2026, Beijing, with Python and SymPy 1.14.0. Both scripts run without importing or editing the retained research modules. The source and JSON outputs live only in this line's directory.

```powershell
python 'code/native_mixed_quartic.py'
python 'code/finite_backreaction.py'
```

Both commands exited with code 0. The latest mixed-quartic run returned 42 checks in about 2.5 seconds; the finite-backreaction calculation returned 14 checks in about 1 second. Timings are execution observations, not performance bounds.

## Mixed quartic coverage

The generic native density follows by exact finite expansion of the minimal scalar action. The code checks its matrix identity and independently polarized four-term local expression. Four exact rational configurations conserve total four-momentum, have all four external legs off shell, and use nonzero external metric and internal longitudinal pivots. For each configuration it compares the new quadratic Einstein tensor-density current with the coefficient in the independently expanded Einstein Gamma--Gamma cubic, in all four auxiliary directions. The source comparison therefore includes 16 exact off-shell coefficients, with external wave symbols and the omitted-on-shell approximation recorded in the JSON.

The generic longitudinal-frequency \(A\)-source and its two arbitrary transverse polarization factors are separately checked symbolically. Its extension to arbitrary transverse dependence follows from the exact geometric \(R_{vv}\) identity used in the derivation. At one rational momentum, the retained auxiliary inverse is multiplied by the actual directly calculated Einstein Hessian, and physical/auxiliary quadratic orthogonality is checked. The arbitrary-source completion-of-square check fixes the absence of an extra factor \(1/2\) in the mixed feedback.

Equations (2)--(6) of DERIVATION.md are generic analytic tensor constructions. The four exact configurations independently check signs, factors and the implemented formula; they are not advertised as exhaustive symbolic verification over arbitrary external momenta. No new full Compton amplitude, factorization theorem, Lorentzian S-matrix prescription, quartic Real chart identity or complete BV lift was run or inferred. The previously completed on-shell Compton result is read-only background and is not recounted as a new verification result.

## Exact finite solution coverage

The second script keeps \(\alpha,\beta,\kappa,\epsilon,v,w\) symbolic. It derives the full inverse metric, all Christoffels and all 16 Ricci components directly from (8), without calling the saved Walker Ricci identity. It checks the full Einstein--scalar equation, vanishing scalar wave, null scalar gradient, determinant, Ricci scalar, geodesic null gradient and exact scalar-amplitude dependence. The rank-two Hessian minor is derived with the actual connection, and the cubic sign obstruction is checked independently.

The finite metric is nondegenerate on \(H\ne0\), with a real split-signature section and a complex analytic continuation. It has no demonstrated Lorentzian Real section in the retained conjugate-helicity geometry. It is a local solution, with no claimed global causal or positive-energy interpretation. Its background/zero-longitudinal data do not lie in the \(a\ne0\) Fourier carrier, and no twistor representative for these zero modes is inferred. It is an exact sector result rather than a proof of the full nonlinear carrier equivalence.

## Retained assumptions and conventions

The action uses one real massless minimally coupled scalar and \(\kappa>0\). The auxiliary coordinate normalization and curvature anchor are inherited and explicitly reproduced in DERIVATION.md; the new checks test their use. Formal rational Fourier kernels are evaluated on declared chart overlaps. No retarded inverse is identified with a cyclic inverse. No distributions at a chart pole, boundary terms at infinity, nonlinear convolution-domain theorem, ghosts, loop determinant or quantum master equation are supplied.

The historical source files, previous research rounds and publication checkout were not modified. No formation/controller epoch was resumed.

## Completed English manuscript

`main.tex` is the self-contained English manuscript, with `references.bib` and `README.md`. The native action, off-shell physical inclusion, literal curvature anchor, both source vectors, auxiliary Hessian and inverse, mixed quartic derivation, finite metric and inverse, full finite equations and rank witness are included in the paper itself. The finite split-signature/zero-longitudinal sector remains explicitly distinct from the localized carrier and a Lorentzian Real lift.

Compiled with the installed TeX Live 2025 `latexmk`/`pdflatex` and BibTeX:

```powershell
latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=build main.tex
Copy-Item -LiteralPath 'build/main.pdf' -Destination 'manuscript.pdf' -Force
```

The command exited with code 0, generating `build/main.pdf` and the delivered `manuscript.pdf`, **10 pages**. Final-log checks found no overfull/underfull boxes, undefined references or citations, duplicate-label warnings or other LaTeX warnings. Human author, title, subject and keyword metadata fields are empty; tool producer and creation-date fields remain ordinary compiler metadata. The paper has no human byline or affiliation. Bibliographic records were verified directly against the two cited primary arXiv pages; only framework claims are attributed to those sources. The required chronology and exact attribution-uncertainty quote are included.

The integrated PDF visual review is performed by the parent task. No PDF marker extraction was rerun for this manuscript.

## Integrated final manuscript inspection

The desktop editor compiler also confirmed success. The parent inspected all 10 rendered pages in contact sheets and pages 5, 8 at readable resolution. No clipping, overlap, missing glyphs or unresolved references were observed. Human byline and PDF author metadata are blank; the corrected chronology and full uncertainty statement are present. Manuscript-facing text contains no internal project codename. No unchanged scientific calculation was rerun for this layout pass.

## Private repository delivery replay

The default commands in README.md were executed successfully from this repository packet on 29 October 2026 Beijing time. They use no external project directory or network input. This reuses completed scientific coverage and verifies portability; it is not an additional independent proof. The manuscript PDF, LaTeX and bibliography match the previously completed edition byte for byte.
