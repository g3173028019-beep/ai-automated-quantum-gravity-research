# Verification

[Paper entry](README.md) · [PDF](manuscript.pdf) · [TeX](manuscript.tex)

## Review and executed checks

- Read the definitions, proof, hypotheses, and boundary statements against the named source literature. The note is an internally reviewed research manuscript, not a refereed publication.
- Compiled the self-contained TeX and rendered the 11-page PDF. Checked page flow, equation and reference layout, and blank author metadata.
- Ran `python replay_rank_two.py`. Exact spin-channel and symmetry-block arithmetic at γ=1/3. The script checks selected finite entries, determinant bounds, and the limiting same-position kernel; it does not replace the analytic spin-foam proof.
- `python replay_rank_bounds.py` checks 2,050 finite rational/cutoff rank inequalities and 1,000 minimum-pair bounds.

The PDF stored with the source is 330172 bytes. The analytic statements and the finite computations are separate; successful finite replay does not certify an unencoded physical interpretation.
