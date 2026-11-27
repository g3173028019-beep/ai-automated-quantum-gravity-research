# Common Einstein cones and complete relational area probes for multifield matter

[Manuscript PDF](manuscript.pdf) · [Verification](VERIFICATION.md)

This standalone manuscript combines two precise local questions: which first-derivative scalar actions share the Einstein light cone on attained open data, and when marked metric areas capture all compact coupled linear observables.

- The classification includes the physical rank-four Gram manifold for arbitrarily many scalars, the branch-dependent parity distinction, and exact constrained-data attainment.
- Ten plane types give compact source representatives off shell for complete material labels. With hidden fields they are complete in the coupled equation quotient when a local differential reconstruction exists.
- Exact backreacting examples show a necessary and sufficient cyclic condition for arbitrary finite field number. An additional six-field example shows that a derivative four-form can turn a previously invisible scalar into an observable one without changing the metric background, direct stress response, or characteristic cones.
- A scoped curved illustration reconstructs a conformally transported scalar through a scalar-curvature zero, and distinguishes the geometric response from the full constrained combination.

The manuscript extends Papers 46 and 47 without replacing either. The curved illustration uses a specified transported theory and does not advance the flat stress-pole statement of Paper 45.

## Build

Run **./build.sh** with pdfLaTeX, Latin Modern, geometry, AMS, mathtools, microtype, xcolor and hyperref. The inputs are main.tex, cone.tex, attainment.tex, area.tex, examples.tex and curved.tex. No shell escape, external graphics or bibliography files are required. The build runs twice and leaves manuscript.pdf only; with the same engine and fonts the PDF is byte-reproducible across temporary build directories.
