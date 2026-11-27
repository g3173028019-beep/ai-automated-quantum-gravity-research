# Verification

## Proof audit

- For each signed exchange the bare modular operator is the same real \(S\), since \(-I\) acts by \(S^2=I\). The cap and seam composition matrix \(T\) is invertible, has signature zero, and \(UT^{-1}U^{\mathsf T}=-A_\epsilon\). The shared radical in the signature subtraction contributes a hyperbolic summand, so no nondegeneracy of signed adjacency is assumed. For bipartite graphs the remaining signature phase vanishes.
- The odd-localized homology is computed by the cap relations, then seam Mayer-Vietoris. The free graph cycle summand comes from the \(H_0\) term. The doubled odd torsion is obtained by Smith equivalence of the signed biadjacency and its transpose. Example nonhomeomorphism follows from distinct homology, without using canonical JSJ presentation.
- The polarized surgery Gauss sum gives the odd-pointed modulus, also for singular surgery matrices. This is a known MOO/Funar formula; the paper uses it as an additional diagnostic. Andersen-Jørgensen's exact monodromy and negative-monodromy agreement is credited, including the \(I,-I\) comparison with Betti numbers 3 and 1. Kania-Bartoszynska and Funar are credited for earlier quantum-invariant collisions.
- The joint proof was checked independently: anchored odd-cycle lift count, weighted sinc envelope, absorption of reflected endpoint gains per site, odd-grid signed quadrature, \(e^{Cn}\) normalization, and the complex phase/prefactor comparison with Paper 50, Theorem 5.1. It proves a sufficient exponential level window; it does not infer a polynomial threshold from the fixed-graph argument.

An independent reader checked both the exact collision/diagnostic and simultaneous-limit sections. The only requested precision was to state the fixed thresholds on class proportions and deletion degree for the uniform constants, which is included in Theorem 5.1.

## Replay and build

- Python/SymPy replay: integral Smith form of the full cap/seam presentation versus the localized signed-adjacency reduction on the stated small families; independently computes signed Maslov contraction matrices, pointed ratios from the homologies, and exact nonvanishing at \(r=3\).
- [Replay transcript](replay-output.txt)
- PDF: pdfTeX 1.40.26, three passes, no final warnings or box diagnostics. Two clean builds in different temporary directories agree byte for byte. All pages rendered and inspected; all fonts embedded. A minimal TeX installation may use the optional environment variable `P50_PDF_MAP` to supply its local Type 1 map. The build does not require online sources.
- [Content checksums](SHA256SUMS)
