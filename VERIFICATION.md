# Verification record for the first repository batch

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
