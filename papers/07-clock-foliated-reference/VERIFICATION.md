# Executed verification and scope

29 October 2026, Beijing. The successful construction is a preferred-clock local effective theory with exact retained background/homogeneous action and positive quadratic energy in all inhomogeneous sectors. It is not a nonlinear stable-universe theorem or a continuum Hadamard-state theorem.

Executed command:

```powershell
python 'scripts/construct_inhomogeneous.py'
```

The run used local Python with SymPy 1.14.0 and SciPy 1.17.1. It completed with exit status zero. Saved outputs are `results/construction.json`, `results/quadratic_action.txt`, and `results/controlled_ell1_evolution.csv`. `the unbundled exploratory scalar trials` preserves the earlier unsuccessful repair trials; it is exploratory and is not needed to reproduce the final result.

| Consequential question | Actual check and result |
|---|---|
| Does the scalar representation close? | Derived radial/tangent sphere Laplacian identities and integrated their exact gradient norms; the basis contains the diagonal target/spatial mixing. |
| Are lapse and shift constraints retained? | Varied the unreduced ADM quadratic action, solved both algebraically, and retained reconstructed A,B in the coupled evolution. |
| Is clock gauge a lost physical equation? | The original and repaired omitted clock equation is identically zero on the respective reduced dynamics. |
| Do reduced matter equations match covariant equations? | Four independent direct covariant KG residuals simplify to exact zero before and after repair. |
| Is exceptional ell1 handled correctly? | Fixed the genuine residual spatial gauge psi=0; checked all five original reduced equations vanish when substituting the four physical accelerations. An independent collaborating agent confirmed the exact common gauge kernel. |
| Are original instabilities numerical noise? | Exact Sturm counts: ell1 has two positive real omega-squared roots and a complex pair; ell3 has one negative and four positive roots. Growth rates are order one and 0.2, respectively. |
| Are ell2 zeros gauge? | Exact physical internal-symmetry kernel `(0,0,1,1/2,0)`, stiffness rank four. The repaired Gram term lifts it. |
| Can radial potential tuning fix ell3? | Exact determinant factor is negative throughout the retained homogeneous elliptic window. This is an analytic obstruction. |
| Is the successful repair an all-tower result? | Exact Sylvester minors: rational positive coefficients after L=t+8, plus exact exceptional L=3 minors. Original ell>=4 positivity also checked exactly at L=t+24. |
| Is the transverse vector reduction constrained? | Shift eliminated explicitly, Killing n1 gauge direction removed, physical n1 oscillator retained; K,H positive for every higher vector harmonic by displayed formulas. |
| Did the tensor model change? | New strain correction computed from the inherited EH TT kinetic normalization: delta omega-squared=8 lambda_bar b²=.1805, with compression/Gram zero at linear TT order. The old tensor derivation was reused. |
| Are strain/compression signs independently checked? | Separate agent derived E1,Q1, their integrated squares, exceptional gauge transformation, and principal scalar speeds; saved sibling evidence `the retained structural derivation (formulas reproduced in the manuscript)`. |
| Was a state actually constructed? | Saved canonical covariance matrices for all scalar multiplets ell1..3, representing 141 oscillators with degeneracies. CCR uncertainty, purity and stationarity checked numerically. |
| Was coupled evolution actually run? | Exact linear propagator for one ell1 harmonic on tau=0..100, including reconstructed A,B. Relative energy drift 1.32e-12, below the analytic positive-energy bound. |
| Is UV propagation unchanged? | No. Calculated modified scalar speeds 1.0277823 and 2.0234863 and three unit speeds. Clock-gauge action has higher spatial derivatives, though reduced frequencies remain asymptotically proportional to k. |

The scalar band covariance stationarity residuals are `7.15e-14`, `2.44e-13`, `4.94e-13` for ell1/2/3. The ell1 purity residual is `4.32e-15`; its minimum uncertainty eigenvalue is `-2.94e-15`, consistent with a pure covariance's exact zero eigenvalues and floating-point error. The matrix-square-root construction, not a hand-tuned roundoff threshold, supplies positivity and purity analytically.

The independent check is mathematical/computational collaboration, not human peer review. No global originality or priority claim was evaluated. No original source, previous round, publication checkout, formation, controller epoch, or other research line was modified.

The self-contained paper is `manuscript.pdf` (11 pages), with editable `main.tex` and `references.bib`. `latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=build main.tex` completed successfully. The final TeX log has no overfull/underfull boxes, undefined references or LaTeX warnings. PDF inspection confirms blank title/subject/keywords/author/creator/producer metadata and no creation-date field. An initial overfull prose line and a tab-escaped fraction in the first source write were corrected before this final compilation. Parent handles visual page inspection. No PDF marker was rerun.

Material gaps and assumptions:

1. The timelike-clock patch is assumed; the covariant preferred-foliation action is not analyzed across a null or vanishing clock gradient.
2. Time derivative count is verified in the local quadratic clock gauge. No global nonlinear Hamiltonian-constraint/degree-of-freedom theorem is supplied for the higher-derivative covariant clock formulation.
3. The finite scalar state is exact for the finite quadratic spectral truncation. Its effective physical applicability assumes a suitable cutoff and suppression of omitted operators. Neither ultraviolet matching nor a quantitative EFT truncation-error bound is known.
4. Modified characteristic speeds exceed the background metric speed. Metric microcausality, Lorentz-invariant ultraviolet completion, and the original normally hyperbolic Hadamard theorem are not claimed for the repaired action.
5. The state and trajectory are physical **linearized** constructions. Nonlinear constraint-solving for those finite-amplitude inhomogeneous data and nonlinear PDE evolution remain undone; homogeneous escapes of the old action persist exactly because all new terms vanish on the homogeneous ansatz.
6. All-L positivity is a property of the specified quadratic effective model. It does not justify extending the effective description to arbitrarily large physical wavenumber.

## Integrated final manuscript inspection

The desktop editor compiler also confirmed success. The parent inspected all 11 rendered pages in contact sheets and pages 6, 10 at readable resolution. No clipping, overlap, missing glyphs or unresolved references were observed. Human byline and PDF author metadata are blank; the corrected chronology and full uncertainty statement are present. Manuscript-facing text contains no internal project codename. No unchanged scientific calculation was rerun for this layout pass.

## Private repository delivery replay

The default commands in README.md were executed successfully from this repository packet on 29 October 2026 Beijing time. They use no external project directory or network input. This reuses completed scientific coverage and verifies portability; it is not an additional independent proof. The manuscript PDF, LaTeX and bibliography match the previously completed edition byte for byte.
