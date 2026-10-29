# Claim correction and artifact growth in AI-assisted theoretical research

Complete English descriptive research-method paper with the preserved corpus, explicit comparison/classification methods, quantitative results, rule sensitivity, mathematical correction cases, limitations and a bounded evidence-routing proposal.

Run from this directory:

```text
python reproduce.py
latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=build main.tex
```

The calculation uses Python's standard library only. It reads local `provenance/` tables and writes `outputs/results.json`. No original research path, network, raw conversation, secret or model-provider data is needed.

The byte-preserved classification and aggregate JSON are accompanied by a reduced table of path identifiers/content hashes and a record of five baseline PDF counts. Hash comparison and filename sensitivity are freshly reproduced from those tables. Historical PDF page parsing was already performed during recovery and is explicitly inherited, not rerun against absent PDFs. The old mutable-timestamp cutoff is not a contemporaneous cryptographic snapshot.

This is an observational artifact-stage study. It does not infer model superiority, a causal effect of auditing, a waste percentage or a validated general workflow optimization. Human names, affiliations and PDF Author are blank. No human expert review occurred. The operator-reported chronology and exact unresolved-attribution quote appear in the manuscript and `PROVENANCE.json`.
