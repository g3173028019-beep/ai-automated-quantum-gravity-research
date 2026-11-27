# Verification

- The release PDF is 51 pages. The internal technical pass reconciled the tree hypotheses, normalization conventions and the scope of the angular and radial statements. The abstract and Introduction now explicitly exclude the one-vertex backbone, as in the main proofs.
- The TeX sources passed a three-pass pdfLaTeX build and checks for unresolved or duplicate references and overfull boxes. Sources and required dependencies are listed in the [paper notes](README.md).
- `python verify_h_tree.py` (requires `mpmath`) compares the two-vertex H-tree finite odd-grid contraction with the closed cotangent formula at four levels, including both parities. [Replay output](verification-output.txt). This finite check is distinct from the proofs for arbitrary fixed backbones.
- The all-orders claim is for the radial/action sectors. The absolute torsion comparison and the joint Wilson law concern the specified leading Hall charts and overlap limits; higher angular coefficients and unrestricted active-skeleton probes are left open.
- The checks were performed in the project's internal research workflow. No external peer review is represented here.

The distributed PDF hash is recorded in [SHA256SUMS](SHA256SUMS).
