# Tensor spectral prices of open Dehn fillings in Turaev–Viro theory

This paper studies the compact \(\mathrm{SU}(2)_k\) Turaev–Viro theory on Seifert cobordisms with at least three marked torus boundary components. The whole, disconnected torus boundaries are the tensor factors.

- **[PDF](manuscript.pdf)**
- **[LaTeX source](manuscript.tex)**
- **[Verification and scope](VERIFICATION.md)**

Two open Dehn fillings have different exact multipartite tensor prices even though the same completion of their remaining boundaries gives \(S^3\) in both cases. The paper computes even-slope support counts at odd root orders and gives equal-price cap experiments whose joint fillings have different nonzero prices. A second example gives a zero joint tensor when the two marked charge supports are disjoint. At a fixed slope, the exact price after subtraction of the unfilled baseline can also retain a dependence on the sequence of root orders.

The results concern a fixed ordinary spectral-price construction and the compact positive-cosmological-constant state sum. They do not recover an intrinsic metric, a dimensionful cosmological constant, or a flat-gravity continuum limit.

## Rebuild

The manuscript is self-contained and embeds its bibliography. With a standard TeX Live installation:

```sh
pdflatex -halt-on-error manuscript.tex
pdflatex -halt-on-error manuscript.tex
```


## Optional exact check

[verify_support.py](verify_support.py) checks the stated finite examples and small-root support formulas in integer cyclotomic arithmetic. Run it with `python3 verify_support.py`. The arguments in the paper do not rely on this finite check.
