# AI-Automated Quantum Gravity Research

English working manuscripts and verification material from an AI-led research project, beginning with TOE-Fable.

The first batch contains two scoped mathematical manuscripts. Author and affiliation fields are blank. The operator's workflow and attribution statement appears in [DECLARATIONS.md](DECLARATIONS.md) and in both PDFs.

## Manuscripts

1. **Generic Sign Rigidity of Real Projection Kernels from Two-Point Data**
   - [English PDF](papers/01-projection-sign-rigidity/manuscript.pdf)
   - [LaTeX source](papers/01-projection-sign-rigidity/main.tex)
   - [Scope and verification notes](papers/01-projection-sign-rigidity/README.md)
2. **Half-Wormholes as Wiener Chaos: Correlated-Replica Spectroscopy, Exact Pfaffian Sectors, and Boundary Moment Indeterminacy**
   - [English PDF](papers/02-half-wormholes-wiener-chaos/manuscript.pdf)
   - [LaTeX source](papers/02-half-wormholes-wiener-chaos/main.tex)
   - [Scope and verification notes](papers/02-half-wormholes-wiener-chaos/README.md)

These were selected because their claims are concentrated, English drafts and analytic arguments already exist, and historical automated checks are available. The broader quantum-gravity and ultraviolet-completion programme remains unfinished. Archival selection does not establish publication-level originality or human validation.

## Status

- Batch date: **28 October 2026**.
- Language: **English only**.
- Human authorship and affiliations: **pending; blank in the PDFs**.
- Human expert review: **not yet performed**.
- Journal submission or acceptance: **not claimed**.
- Source history: [HISTORY.md](HISTORY.md).
- Workflow, attribution, rights and review status: [DECLARATIONS.md](DECLARATIONS.md).

This archive is private. Raw conversations, credentials, third-party paper PDFs and unrelated project files are outside this batch.

## Reproduction

With TeX Live and `latexmk` installed, run:

```powershell
./tools/build.ps1
```

Historical Wolfram Language scripts are included under each paper's `verification/` directory. Their cost and coverage are described in the paper notes. Current build and archive checks are in [archive-checks.json](archive-checks.json), with their limits in [VERIFICATION.md](VERIFICATION.md).
