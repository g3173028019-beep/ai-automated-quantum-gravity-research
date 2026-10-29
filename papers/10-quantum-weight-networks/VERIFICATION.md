# Verification record

Date: 29 October 2026, Beijing. Scientific outcome: the finite control statements and their stated iid limits pass the checks below. Their interpretation remains the one stated in the paper.

## Fresh portability execution

Executed the bundled `code/verify.py` by absolute path from an unrelated working directory. Exit code 0, SymPy 1.14.0; no external project import or network request is required.

Verified with exact matrices:

- One-cell Lie dimension 10, nonsingular invariant antisymmetric form, Clifford volume and polarized epsilon/metric trace identities.
- All displayed J/K/E/C bracket signs and the full six-variable symmetric covariance matrix.
- All 32 chiral weight projectors: 16 positive-parity rank-one projectors, 16 negative-parity zero projectors.
- The different noncommuting update `(X+Z) tensor I` gives dimension 15.
- The positive-weight double-commutator identity and the mean-field free-energy derivatives.

Independent integer Pauli closure returned:

| Cells | Vector bond | Clock bond |
|---|---:|---:|
| 1 | 10 | 10 |
| 2 | 45 | 255 |
| 3 | 4095 | 4095 |
| 4 | 65535 | 65535 |

## Analytic results

The manuscript contains the complete positive-weight generation argument, the block-extension lemma and spanning-tree proof for arbitrary connected networks. The all-graph proof received an additional synthesis-model algebraic cross-check during the preceding recovery. It has not received human expert review or proof-assistant formalization.

The iid CCR statement uses standard bounded-product-state moment combinatorics, with all required local moments checked. The exact state rms and full operator norm of the diagonal commutator residual are derived in the paper. The Curie–Weiss spectator reduction is an explicit representation argument, not an additional numerical hypothesis.

## Build and LaTeX cleanup

`latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=build main.tex` completed with exit code 0 using TeX Live 2025. Final output: 12 pages. The final LaTeX log has no undefined or multiply defined references, no citation warnings, and no overfull or underfull boxes. Blank PDF Author is set explicitly. The shell reported a Perl locale fallback; it did not affect the completed build.

The source was checked for duplicate labels, unresolved citations, obsolete project terminology and draft placeholders. Included theorem files are part of the final source, rather than external historical inputs. Integrated PDF visual QA and upload are handled by the parent delivery process; this worker does not claim that inspection as already performed.

## Material limits

No actual interacting finite-N Gibbs fluctuation theorem, gravitational projector, energy-gap protection, spacetime locality, Hamiltonian constraint or Lorentzian gravity is established. The product-state theorem cannot supply an operator-norm gravitational closure estimate. Historical non-Pauli bond scans were not repeated. The primary comparison is targeted; no exhaustive originality certificate or public-priority inference is made.
