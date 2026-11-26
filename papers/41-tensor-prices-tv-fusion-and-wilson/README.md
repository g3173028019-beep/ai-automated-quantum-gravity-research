# Physical Gluing and Wilson Constraints from Euclidean Tensor Prices

An English research manuscript on the compact unitary Turaev–Viro theory. It prices the input–output tensors of physical cylinders, rather than individual kinematical boundary factors.

- **Fusion:** Physical rank ratios give quantum-dimension rates and a nested-cut fusion law. For doubled SU(2)\(_k\), diagonal fundamental labels generate a parity-matched sector. The first strict deficit from the untruncated vacuum rank determines \(k\).
- **Recoupling:** Coherent combinations of actual annular projectors recover relative cut data up to channel-basis gauge. In a doubled level-two four-puncture example the signed pencil vanishes at \(-1/16\).
- **Torus Wilson constraints:** A joint penalty for three physical Wilson operators recovers the marked compact torus character variety, including its degenerate points, and the normalized Haar state count in the large-level limit.

The state-sum, fusion and Wilson operators are inputs from the cited literature. The manuscript does not derive Lorentzian gravity, assign a metric length to a defect, or assert that the native price of an isolated scalar knows its phase.

## Files and build

- [manuscript.pdf](manuscript.pdf): compiled paper.
- [manuscript.tex](manuscript.tex): main LaTeX source.
- [torus_joint_wilson.tex](torus_joint_wilson.tex), [torus_bib.tex](torus_bib.tex): included section and references.
- [replay_tv_native_cylinder.py](replay_tv_native_cylinder.py): exact integer and rational checks of the finite fusion and recoupling examples.
- [replay_tv_torus_symbols.py](replay_tv_torus_symbols.py): checks of the finite torus phase convention using exact exponents.
- [VERIFICATION.md](VERIFICATION.md): verification scope and limitations.

Run `latexmk -pdf manuscript.tex` in this directory on a standard TeX Live installation (including Latin Modern and PGF/TikZ). Alternatively run `pdflatex manuscript.tex` until the cross-references settle. The two replay scripts require Python 3 and its standard library only.

This is an AI-led research manuscript with internal mathematical and computational checks; it has not been independently peer reviewed.
