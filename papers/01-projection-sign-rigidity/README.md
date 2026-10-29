# Generic sign rigidity of real projection kernels from two-point data

Revised edition: **29 October 2026**, 6 pages. [Manuscript](manuscript.pdf), [LaTeX source](main.tex), [revision record](REVISION.md).

For N>=5 and 2<=r<=N-2, almost every real rank-r projection is determined by its diagonal and entrywise moduli up to diagonal-sign gauge. The generic theorem and its analytic proof are retained. The revised fixed-data criterion uses the observed feasible holonomy fiber; an exact full-support rank-one example explains why uniform injectivity is stronger. Sparse/complex/continuum, stable recovery and polynomial-time learning remain outside the theorem.

Run the portable exact revision check with Python >=3.10:

```powershell
python verification/check_revision.py
```

It uses the standard library, enumerates 1024 star-gauge assignments and writes `verification/revision-checks.json`. Build the paper with TeX Live and latexmk/BibTeX:

```powershell
./build.ps1
```

The source and bibliography are self-contained. `verification/legacy_checks.wls` is the unchanged historical supplement and is optional; its Wolfram checks were not rerun. Default reproduction needs no historical project directory or network.

Human authors and affiliations are blank; PDF Author is blank. Human expert review has not occurred. Targeted prior-source comparisons do not establish comprehensive novelty or public priority. Source history, present revision and demonstrable public disclosure are distinct; [provenance](PROVENANCE.json) pins the input edition.
