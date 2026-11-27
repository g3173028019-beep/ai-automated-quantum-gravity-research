# The causal memory of gravitational radiation from compact area measurements

- [Manuscript (PDF)](manuscript.pdf)
- [LaTeX source](manuscript.tex)
- [Verification](VERIFICATION.md)

One fixed compact normal area profile, combined with all compact temporal drive and readout tests in a finite window, determines the causal force law at a cut in a calibrated gravitational shear channel. If the profile lies at distance d from the cut, the recovered window loses d at each end. The retained shear potential can be any smooth known function of calibrated time. Exact Einstein–scalar examples include nonstationary FLRW solutions and anisotropic solutions with nonzero Weyl curvature.

The paper distinguishes finite causal memory from infinite linear memory. In the attained inverse-square family the rank is exactly n at coupling n(n+1) and otherwise infinite. These are ranks of a retarded force law, not claims about Bondi displacement memory, global spectra or lossless Hamiltonian realizations.

The manuscript follows Papers 44 and 47 and leaves both originals unchanged. It is independently readable and does not use the companion common-cone manuscript.

## Reproduction

Run the following in this directory with pdfLaTeX and the AMS, geometry, microtype, xcolor and hyperref packages:

~~~sh
./build.sh
~~~

The build runs three passes with a fixed source epoch. The optional exact arithmetic replay requires Python and SymPy:

~~~sh
python replay.py
~~~

Use the [checksums](SHA256SUMS) to verify the archived source and PDF. Byte-for-byte PDF reproduction refers to the same TeX engine and font setup.
