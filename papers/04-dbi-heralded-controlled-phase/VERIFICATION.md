# Focused closure verification

The round-02 independent check of the model's linear P/Q instrument, reference extension, fixed phase correction, time-dependent residual theorem and robustness domain is reused. No broad review or new continuum refinement was performed. This closure independently verifies the consequential sharper fidelity bound once and replays the existing numerical certificate at its frozen time.

For `||K-D||<=eta<1`, write `v=D Psi`, `x=K Psi`, `t=||x||`. Expanding `||x-v||^2` gives `Re<v,x> >= (1+t^2-eta^2)/2`. The tangent identity `(t-sqrt(1-eta^2))^2>=0` implies conditional squared fidelity at least `1-eta^2`. This remains true with arbitrary reference systems, and purification extends it to mixed inputs. The bound is sharp even for contractions: `K=(1-eta^2)I+eta sqrt(1-eta^2)J`, with real antisymmetric `J` satisfying `J^T J=I`. The default replay checks these exact algebraic identities with SymPy, not a fitted numerical example. No generic mathematical novelty is claimed.

The default `python replay.py` was executed successfully with mpmath 1.3.0 and SymPy 1.14.0. It checks the exact copied input hashes and integrator identity; recombines saved cell upper bounds and analytic tails with rational arithmetic; checks each initial mismatch from all sixteen basis coefficients; propagates the retained residual slopes at the new time; and recalculates the target-distance and robust-box guarantees at 110 interval digits. It asserts both nominal and robust distance constants are below one before using positive success bounds. The saved cellwise integrand bounds are inherited evidence, not freshly regenerated continuum integration.

Actual newly replayed guarantees:

| Bound | Nominal | Robust box |
|---|---:|---:|
| Success lower | 0.6565527711781338 | 0.6351321846638459 |
| Conditional squared fidelity lower | 0.96400610918152 | 0.9587714318247369 |
| Conditional Bell negativity lower | 0.4640061091815201 | 0.45877143182473695 |

The input class is common radial `psi5` tensor arbitrary logical/reference state, with exact invariant-code dynamics. Source/readout errors use normalized pure-vector Hilbert norm, not an unqualified trace-distance substitution. These bounded gate statements imply no domain or unbounded-density guarantee for arbitrary perturbed sources. All failures remain in the complete flagged instrument. A successful-subnormalized diamond-distance bound is distinguished from a deterministic-channel claim.

Compiled with TeX Live 2025 `latexmk` 4.86a and pdfLaTeX, using `latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=build main.tex`. Final build exits 0 and produces a six-page PDF, copied to `manuscript.pdf`. The final log has no undefined references/citations or overfull/underfull boxes. The first compile exposed an unescaped underscore in an internal bibliography path; using `path` formatting fixed it. The preamble has `author{}` and empty PDF author metadata. Integrated visual inspection covered all six pages in a contact sheet and the fidelity proof and instrument page at readable resolution, with no clipping or overlap found. This is a GPT-led same-team nonblind verification, not human expert review.
