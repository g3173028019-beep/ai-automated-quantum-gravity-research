# Integer STU flux families

This compact research note gives an exact integral divisor family, its full real F-term rank, positive canonical masses, modular point types and exact pointwise charge minima in the stated STU period/lattice model. A capacity-28 example has S=9i, U_i=3i, quotient charge 24 and g_s=1/9. The global orientifold interpretation and capacities remain inherited conditional inputs. The nonzero superpotential and unstabilized Kähler sector are explicit.

Read manuscript.pdf or main.tex. No journal-submission or global novelty certification is claimed.

## Reproduce

With Python 3.11 or newer and the pinned local dependency installed:

```text
python -m pip install -r requirements.txt
python reproduce.py
```

The reproduction command itself uses no network, old project directory, external data, credentials or user-specific path. It rebuilds all 33 divisor assignments and nine reduced point types, symbolic rank and mass checks, primitive-pair costs and a nontrivial modular control. Results are written to results/construction_certificate.json and results/divisor_assignments.json. Dependency installation may require network if SymPy is not already installed; it is separate from offline reproduction.

No external input files are needed: the exact periods, lattice ordering and integer construction are specified in the manuscript and code. PROVENANCE.json identifies completed earlier source results by relative identifiers and SHA-256 hashes; those records are provenance, not runtime inputs. The results directory contains freshly generated portable outputs, not copies presented as newly recomputed data.

Build the note with a standard TeX installation:

```text
latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=build main.tex
```

Copy build/main.pdf to manuscript.pdf after compilation. BibTeX references are in references.bib. No third-party full papers are bundled.

## Interpretation and disclosure

Masses use the common unstabilized Kähler factor e^(K_K)=1. This is a no-scale CS/axiodilaton construction, not a fully stabilized string vacuum. Nine means reduced moduli point types under separate modular factors and torus permutations, not a complete count of flux or string-duality orbits.

Earlier work was Fable 5-led; July 2026 onward results were predominantly produced by GPT models, as reported by the operator. DeepSeek and multiple GPT models also contributed cross-review and manuscript preparation. This is reported chronology, not a per-result authorship allocation. Human names and affiliations remain blank; human expert review has not occurred. The present delivery date is 29 October 2026. Historical source dates and demonstrable public disclosure are distinct; a private commit does not establish public priority.

"I genuinely do not know how authorship, intellectual credit, or ownership should be attributed for these results."
