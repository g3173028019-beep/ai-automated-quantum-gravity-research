# Finite causal spin-two transfer and partial refinement in Lorentzian Regge calculus

The finished manuscript is `manuscript.pdf`; its source is `main.tex` with `references.bib`. Human authors, affiliations and PDF Author metadata are blank. Human expert review has not occurred.

The paper gives an explicit finite Lorentzian Regge region, its complete-action and physical quotient recipe, the angular-weight-two carrier, an architectural rank-two bound, the corrected endpoint identification statement, the retained fifth-level hyperbolic transfer and a closed scalar-kernel calibration. It does not claim a continuum no-go, uniquely identified continuum helicity, retarded propagation, a selected positive-energy vacuum or a full local perfect action.

## Portable default replay

Use Python 3.10 or newer, with the dependencies in `requirements.txt`:

```text
python -m pip install -r requirements.txt
python reproduce.py
```

The replay writes `results/replay.json`. It requires no network after dependency installation, no original project directory, no API key and no conversation data. It verifies 20 frozen-file hashes, exact commutant/canonical/endpoint algebra, the positive invariant and rank factorization, the explicit geometry at levels 1–5, the saved first- and fifth-level transfer spectra, and their real scalar calibrations.

**The default replay does not recompute the complete Regge Hessian.** It reconstructs and tests transfer matrices from packaged frozen Hessian blocks. The full action calculations were completed during the preceding recovery and are retained evidence: level 1 took 12.70 seconds, and level 5 took 423.37 seconds on that workstation. No costly fifth-level calculation was rerun for this paper delivery. The numerical table's levels 2–4 and historical step-size checks are inherited report evidence.

## Optional full calculation

The original complete-action producer and its local import dependencies are copied byte-for-byte under `provenance/frozen_producer/`. To run it explicitly:

```text
python run_full.py --level 1
python run_full.py --level 5
```

Outputs go to `results/full-level-1.json` or `results/full-level-5.json`; the frozen inputs are never overwritten. These commands use the same original numerical action implementation, not a second independent Regge implementation. Runtime depends on hardware. They are separate from the quick default so that checking the paper does not silently repeat an expensive run. The entry point's imports and help path were tested; no optional full Hessian was rerun for packaging.

## Build

With an existing TeX distribution containing `latexmk`, `pdflatex` and BibTeX:

```text
latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=build main.tex
```

The deliverable `manuscript.pdf` is the copied final `build/main.pdf`. The build produces nine pages with no final-log undefined references, missing glyphs, overfull or underfull boxes. A Windows Perl locale fallback warning from the build wrapper does not affect compilation. Integrated visual QA is handled by the parent delivery.

`PROVENANCE.json` records immutable input/output SHA-256 identifiers, original source locations and source chronology. The frozen reports live under `provenance/reports/`; internal historical identifiers are confined to provenance. No third-party full papers are bundled. `VERIFICATION.md` distinguishes what was proved, recomputed, checked from saved outputs and left unresolved.

## Attribution

Earlier work was Fable 5-led; July 2026 onward results were predominantly produced by GPT models, as reported by the operator. DeepSeek and multiple GPT models also participated in review/preparation. This is reported chronology, not an allocation of individual authorship or ownership.

“I genuinely do not know how authorship, intellectual credit, or ownership should be attributed for these results.”
