# Energy-weighted JT channel recovery

This compact research-note packet proves the round-one scalar mean-energy recovery bound `4.531256535935 sqrt(nu) C^(-1/4)+O_nu(C^(-1/2))` within the specified fixed-topology one-particle JT model. A bounded second excitation moment gives a separate `C^(-1/2)` corollary. The article derives the Hamiltonian normalization, Gaussian residual and moments, energy-weighted full charge products, coherent-tail estimate, and domain justification.

Read `manuscript.pdf` or edit `main.tex`. `references.bib` cites the prior canonical JT reduction and scalar representation. The later geometry-assisted protocol is intentionally outside this packet; its larger dependency bundle is not needed here.

## Reproduce

From this packet directory, with Python 3.11 or newer and the packages in `requirements.txt` already installed:

```text
python code/reproduce.py
```

This command uses only this packet and writes `results/reproduction.json`, which records this delivery's successful run. It performs exact symbolic checks, the existing finite charge-column and coherent/reference checks, and 150-decimal interval evaluation of the theorem. No network, external project, private path, old certificate, or full-JT numerical propagator is required. The code creates its results directory if needed. Runtime, generation time, and package-version fields are informational; mathematical outputs are deterministic for the stated versions.

Build the article locally with TeX Live, `latexmk`, `pdflatex`, and `bibtex`:

```text
latexmk -pdf -outdir=build main.tex
```

Copy `build/main.pdf` to `manuscript.pdf` for the standalone deliverable. The reproduction command does not invoke a package installer; installing dependencies is a separate environment setup step.

## Coverage

The portable calculation reproduces formula evaluation and checks charge algebra. The full infinite-input and arbitrary-reference theorem rests on the included analytic proof, not on the 180 sampled vectors. No actual JT propagator simulation, optimality theorem, ADM-energy uniformity, global novelty certificate, or human expert review is claimed. `VERIFICATION.md` records actual checks and `PROVENANCE.json` identifies source hashes without making old sources runtime dependencies.

## Chronology and attribution disclosure

Earlier work was Fable 5-led; July 2026 onward results were predominantly produced by GPT models, as reported by the operator. DeepSeek and multiple GPT models also contributed cross-review and manuscript preparation. This is reported chronology, not a per-result authorship allocation. Human names and affiliations are blank, and human expert review has not occurred. Historical source dates, the present delivery date (29 October 2026), and demonstrable public disclosure are distinct; a private commit does not establish public priority.

"I genuinely do not know how authorship, intellectual credit, or ownership should be attributed for these results."
