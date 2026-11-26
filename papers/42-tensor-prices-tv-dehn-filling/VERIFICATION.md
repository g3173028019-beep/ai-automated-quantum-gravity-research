# Verification

## Independent review

A separate AI reviewer checked the manuscript's boundary and modular conventions, the even- and odd-root Gauss-sum arguments, the common \(S^3\) closure, both charge-support composition examples, and the claimed scope and attribution. The review prompted closer citations to prior work on Seifert and torus-link states, an explicit identification of the normalized price with one-port order-zero Rényi entropy in this diagonal family, and the example with two different nonzero joint prices. This was an internal cross-check, not external peer review or journal acceptance.

## Exact replay

`python3 verify_support.py` requires only the Python standard library. It checks zeros of the modular matrix entries with integer polynomial arithmetic modulo the cyclotomic polynomial \(\Phi_{4r}\); it uses no floating-point zero threshold. The finite replay checks:

- the two benchmark slopes \(2\) and \(2r\) for \(4\leq r\leq21\);
- the complete vacuum even-slope support formula for odd \(r\) in that range and \(0\leq m\leq r\);
- the disjoint coloured supports at \(r=15\);
- the nonzero joint support counts at \(r=135\).

The general statements and limits follow from the analytic arguments in the manuscript, not from extrapolating this finite replay.

## PDF build

The source includes its bibliography. From this directory, with a standard TeX Live installation:

```sh
pdflatex -halt-on-error manuscript.tex
pdflatex -halt-on-error manuscript.tex
```

The compiled PDF was checked for resolved references and visually rendered page by page.

## Interpretation

The result compares actual states on disconnected, marked whole torus boundaries in the compact unitary theory. Open states that give the same closed \(S^3\) contraction are not being presented as distinct closed invariants. The root-sequence result concerns the exact support price after subtraction of an unfilled baseline at each level; it is **not** an identification of \(q\to1\) with a physical \(\Lambda\to0\) limit, and it does not preclude smoother state limits.
