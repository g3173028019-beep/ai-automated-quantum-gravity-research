# Verification of the revised Gaussian-chaos paper

## Analytic scope

The compact-density repair uses a positive-probability event 0<s_-<=S<=s_+<infinity. The Gaussian mixture is then uniformly positive on each compact interval; `(log f)_+<=f` controls the other logarithmic part. Combined with the two symmetric stretched-exponential tail lower bounds, this closes the stated Krein exposition without changing hypotheses or conclusion.

The canonical-family correction is classical moment-operator terminology. The single-draw correction fixes the probabilistic interpretation and explicitly treats n=1. These points are analytic qualifications, not numerically proved density or operator statements.

The even-q W4 proof enumerates the only two possible chaos-four monomial types and uses even-block completion-sign coherence. The convolution proof expands the Pfaffian in even minors and uses orthogonality for distinct Gaussian vertex supports. The quadrature calibration is a classical finite moment identity in the determinate n=2 example. The all-size statements rest on these proofs, not finite tests.

## Fresh portable calculations

Both Python checks passed using only the standard library:

- `verification/check_revision.py`: 90 exact central-binomial reductions for 1<=n<=12 and all k; 12 fourth-moment sum checks; n=1 empty intermediate sector; independent full Pfaffian-polynomial Gaussian expectations versus the minor-convolution formula at n=1,2,3.
- `verification/check_chaos_extras.py`: exact completion-sign enumeration for the even reference N=8,q=2 and formal odd N=12,q=3; the latter yields 17,580,480 rather than the blindly extended even formula 20,328,000. Exact Q[sqrt(13)] quadrature reproduces moments 0 through 7. Odd-q output is a formal sign control, not physical SYK evidence.

The extras script was copied byte for byte from the focused recovery calculation and rerun locally in the packet. Saved outputs record the exact values. No expensive Monte Carlo, historical Wolfram loop or source-project producer was rerun. The optional unchanged historical supplement and earlier finite checks remain inherited evidence with their original scope.

## Build and static coverage

TeX Live 2025 latexmk/pdflatex/BibTeX compiled successfully. `manuscript.pdf` has **10 pages**. The final log has no LaTeX warnings, unresolved references/citations or overfull/underfull boxes. Labels are unique, cited keys resolve and PDF Author is blank. Default reproduction and the paper need no original project directory or network. Integrated PDF visual QA is the integrator's separate check.

## Attribution and remaining gaps

Slud's publisher metadata/abstract and the author-hosted powers/products paper were inspected in the existing audit; this revision reuses those primary records rather than claiming a new comprehensive search. The q=2 reduction is prior mathematics with an independent derivation. The covariance and hypercontractive tools are established. Current global originality of the remaining consequences is unresolved. The finite model does not identify a dynamical bulk, a full microscopic ensemble from one scalar data class, or a complete quantum-gravity theory. Human expert review has not occurred.
