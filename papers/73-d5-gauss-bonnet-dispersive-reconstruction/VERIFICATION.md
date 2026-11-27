# Verification

The manuscript starts from an ordinary fixed-spacelike-transfer hierarchy, including its twice-subtracted member, finite low-contour analyticity, an absolute real coherent-cut norm, and PSD physical SO(4) angular pairings. It does not derive those premises. For the Gauss–Bonnet application the pole-subtracted local tree amplitude is polynomial or normally convergent on the low contours in its Ward-transverse contact expansion. String comparisons are tree-level only.

## Build

With TeX Live, Latin Modern and AMS packages installed, run twice:

    pdflatex -interaction=nonstopmode -halt-on-error main.tex

The only TeX source included by `main.tex` is `references.tex`. There are no figure or network dependencies.

## Analytic and algebraic checks

- `angular-certificate.md` gives the explicit all-spin Gram lower bound and the quadratic matrix/curvature remainder estimates.
- `replay_angular.py` verifies the finite exact rational Bernstein, Bessel-series, geometric-remainder and Newton-cost inequalities used by that analytic certificate. Run with stock Python 3 and asserts enabled: `python3 replay_angular.py`. No finite-spin quadrature is part of the proof.
- `replay_symbolic.py` verifies the finite screen algebra, normalized Gauss–Bonnet trace, forward-profile and latitude identities, and sample Taylor terms of the Lagrange/resolvent calculation. Run `python3 replay_symbolic.py` with SymPy. These checks are not a machine formalization of the infinite series or a verification of an arbitrary full amplitude.
- `regularity-countermodels.md` isolates the failure of inferring physical angular compatibility or complex-disc control from just off-forward values and smooth real data.

The constant angular cap \(Q/M\le10^{-5}\) is a sufficient bound for positivity of the constructed signed operator. It is not an asserted uniform analyticity or physical bandwidth. The claimed pure-metric inequality uses the bandwidth where **both** crossed reconstruction and the angular certificate hold, with the Newton pole retained.
