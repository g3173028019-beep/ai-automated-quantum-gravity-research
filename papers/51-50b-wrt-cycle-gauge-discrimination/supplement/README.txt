Extensions of paired-cap semiclassical sewing
Prepared 27 November 2026.

Build the mathematical note from this directory:
  pdflatex -interaction=nonstopmode -halt-on-error wrt_extensions.tex
  pdflatex -interaction=nonstopmode -halt-on-error wrt_extensions.tex

Reproduce the finite checks (Python: numpy and sympy):
  OPENBLAS_NUM_THREADS=1 python checks/check_even_cycles.py
  python checks/check_signed_homology.py
  python checks/check_gaussian.py

Each section is include-ready LaTeX. The source manifest records the exact
input blobs. Original repository files are not duplicated in this bundle.

Proved from stated inputs:
  - Actual ratio-resolved fixed-spin Wilson asymptotics on matched trees,
    without matched-product cutoffs; explicit P4 hard-window coefficient.
  - Closed vacuum leading asymptotic for all fixed even cycles.
  - One- and two-filter asymptotics for the exact P4 quadratic normal phase.
  - Integral H1 with variable paired cap orders; full odd-pointed scalars;
    local-ribbon blindness; odd nilpotent representation factorization;
    all twisted odd nilpotent DW collisions and odd pro-p completions.

Historical status in the original note (updated by Paper 51B):
  The fully ungated next-logarithm finite-level P4 Wilson coefficient.
  Paper 51B supplies both finite-level cancellation estimates and proves
  the ungated coefficient for the specified terminal ribbon on P4,
  including the extension to fixed integer spin. See
  ../../51b-wrt-ungated-wilson/README.md.
