# Verification record

## Proof controls

| Claim | Verification in the manuscript |
| --- | --- |
| H logarithm and exact action | The exact cap contraction is checked by the finite odd sine transform. The discarded charge and angular regions are bounded at order sqrt(r), below the closed leading term. The quaternionic harmonic chart, torsion linking square, and Jacobian of the local surgery density are calculated independently. No allocation of the finite part or formal continuum measure is used. |
| Absolute APS | The terminal-deleted oriented trace extends the adjoint system. Its ordinary and twisted intersection forms are exhibited separately, including their radicals. The global sign is fixed by the induced section/fibre convention before comparison with RT. Smooth quotient and torsion are asserted only at irreducibles with det B nonzero, not at aligned reducible characters. |
| Unmarked odd-cycle comparison | The degree-two endpoint amplitude is bounded over the entire charge range. The map from dyadic vertex logs to edge sums is invertible for odd cycles, giving a direct uniform quantum bound. The angular determinant factors over edges and has an integrable inverse square root. |
| Unmarked two-seam star | Absolute values of rho are obtained by rational seam compression and the bare-line low spectrum, not from the Gaussian signs. Exact finite leaf sine sums provide linear endpoint estimates and summable nonresonant tails. Signed Weyl pairing cancels the bulk poles; the central coefficient is computed independently by the Dirichlet half value. |
| Marked tree and triangle | The Wilson cutoff is explicit. Exact signed-lattice nonstationarity removes noncentral cap residues before passing to the compact scaled limit. Uniform dyadic tails are controlled before the angle-determinant integral. Only after evaluating the damped loop integral is its angular boundary value taken. |

An independent manuscript reader checked the proof package and flagged two corrections now incorporated: the irreducibility requirement for the smooth quotient and the integral spine contribution in the eigenline CS lift. The corrections do not change the rho or WRT coefficients. This is internal verification, not external peer review.

## Build and replay

The PDF builds to 13 letter-size pages with all cross-references resolved. All pages were rendered with Poppler and visually checked; there are no overfull boxes or missing glyphs. The replay evaluates the finite sine matrices, without using the asymptotic formula to construct the input. It checks the exact cap reciprocity formula and H cotangent sum to floating-point accuracy, and displays convergence of the three indicated limiting coefficients. The small marked-triangle coefficient follows cancellation of larger summands and shows appreciable finite-level error at the smaller test levels.
