# Verification record

## Delivery of 29 October 2026

The three new manuscripts compiled with pdfLaTeX, BibTeX and latexmk to six pages each. Their final logs contain no unresolved references, citations or overfull boxes. The integrated review inspected all 18 rendered pages in contact sheets and selected equation, table and disclosure pages at readable resolution. No clipping, overlap or unresolved-reference markers were found. Human bylines and PDF author metadata are empty.

Each default reproduction command was run successfully from its final repository packet, without historical project paths or network input:

- `papers/03-jt-energy-weighted-recovery`: `python code/reproduce.py` reproduces Gaussian and charge identities, the optimized leading coefficient, 180 coherent matter/reference vector checks and 150-digit interval evaluations. The infinite-dimensional result rests on the stated analytical derivation; the finite checks are implementation checks.
- `papers/04-dbi-heralded-controlled-phase`: `python replay.py` verifies copied certificate hashes, exact residual sums and initial mismatch, 110-digit gate estimates, and the sharp conditional-fidelity identity. It passes with nominal success at least 0.6565527711781338 and conditional squared fidelity at least 0.96400610918152. The robust bounds are 0.6351321846638459 and 0.9587714318247369. The saved cellwise residual bounds are inherited; the continuum producer was not rerun.
- `papers/05-integer-flux-families`: `python reproduce.py` checks the exact flux signs, charge normalization, generic family rank and mass formulas, 33 assignments, nine reduced moduli point types and 72 complementary-pair weights. The geometric source capacities and global orientifold assumptions remain inherited inputs.

Packet-specific `VERIFICATION.md` files distinguish derivation, fresh checks, inherited evidence and limits. One focused closure pass reused completed research checks; no full dynamical simulation, new broad literature review, formal proof-assistant certification or human expert review is claimed. The sharper DBI fidelity lemma and the flux branch corollary are bounded mathematical additions, not new physical protocols or global optimality results.

`archive-checks.json` records the final build, metadata, disclosure, bibliography, English-language and common credential/path-pattern checks, together with packaged-file hashes. The two earlier paper directories were compared against their pre-delivery SHA-256 manifest and remain unchanged. These archive checks do not certify physical interpretation or global novelty.

## First repository batch, 28 October 2026

The English sources were read with their stated assumptions. Edits concern blank authorship, the requested AI-workflow disclosure, provenance, bibliography portability, a related-work update, and the intermediate-sector indices described in `HISTORY.md`.

Both PDFs are rebuilt from source. Current archive checks inspect blank author fields, the disclosure, bibliography keys, unresolved PDF cross-references, placeholder bylines, accidental Chinese-language content, and common credential/private-path patterns. `archive-checks.json` records the checks and SHA-256 hashes of packaged files. Rendered pages are inspected for layout.

These are scoped manuscript and archive checks, not a fresh independent proof of every result and not human expert review.

The archive checker uses Python 3 and `pypdf`. The optional visual-review helper also uses Pillow and Poppler (`pdftoppm`). These dependencies are separate from the LaTeX-only build command.

For this edition, AI visual inspection covered the five-page projection PDF and seven-page half-wormhole PDF using rendered page contact sheets and the disclosure pages. The author lines are blank; the requested workflow statement appears in both PDFs. No clipping, overlap, or unresolved cross-reference markers were found in the inspected renders and final build checks.

## Historical automated evidence

The following coverage is inherited from retained records rather than claimed as newly executed:

- **Projection:** the analytic nonvanishing and finite-certificate argument, exact-rational subset-sum examples, and symbolic examples. Finite checks support implementations; the generic theorem rests on the analytic proof.
- **Half-wormholes:** finite Pfaffian/chaos calculations, generating-function checks, exact low-dimensional determinant identities, and identified Monte Carlo comparisons. The moment-determinacy argument is analytic, not certified by Monte Carlo.

The original scripts can be substantial. The half-wormhole script includes 400,000-sample Monte Carlo comparisons; the projection script includes exponential subset-sum checks. They require Wolfram Language and are optional for reading the papers. Their historical comments are not global-originality claims.

Human review, a comprehensive current novelty review, journal submission and acceptance remain pending.
