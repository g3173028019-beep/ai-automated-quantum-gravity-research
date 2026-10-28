# AI-Automated Quantum Gravity Research

English working manuscripts and verification material from an AI-led research project, beginning with TOE-Fable.

The archive contains scoped mathematical and model-physics manuscripts with minimal reproduction material. Human author and affiliation fields are blank. The current chronological contribution disclosure and unresolved attribution statement are in [DECLARATIONS.md](DECLARATIONS.md).

## Manuscripts

1. **Generic Sign Rigidity of Real Projection Kernels from Two-Point Data**
   - [English PDF](papers/01-projection-sign-rigidity/manuscript.pdf)
   - [LaTeX source](papers/01-projection-sign-rigidity/main.tex)
   - [Scope and verification notes](papers/01-projection-sign-rigidity/README.md)
2. **Half-Wormholes as Wiener Chaos: Correlated-Replica Spectroscopy, Exact Pfaffian Sectors, and Boundary Moment Indeterminacy**
   - [English PDF](papers/02-half-wormholes-wiener-chaos/manuscript.pdf)
   - [LaTeX source](papers/02-half-wormholes-wiener-chaos/main.tex)
   - [Scope and verification notes](papers/02-half-wormholes-wiener-chaos/README.md)

3. **Energy-Weighted Channel Recovery in a Fixed-Topology JT Model**
   - [English PDF](papers/03-jt-energy-weighted-recovery/manuscript.pdf)
   - [LaTeX source](papers/03-jt-energy-weighted-recovery/main.tex)
   - [Reproduction and verification](papers/03-jt-energy-weighted-recovery/README.md)
4. **A Heralded Controlled-Phase Gate from Radial Return in an Independent-Field DBI Model**
   - [English PDF](papers/04-dbi-heralded-controlled-phase/manuscript.pdf)
   - [LaTeX source](papers/04-dbi-heralded-controlled-phase/main.tex)
   - [Reproduction and verification](papers/04-dbi-heralded-controlled-phase/README.md)
5. **Integer STU Flux Families: Positive Masses and Tadpole-Controlled Coupling**
   - [English PDF](papers/05-integer-flux-families/manuscript.pdf)
   - [LaTeX source](papers/05-integer-flux-families/main.tex)
   - [Reproduction and verification](papers/05-integer-flux-families/README.md)

The 29 October packet is selected for closed, reproducible claims. JT concerns the same-mean-energy leading recovery theorem, distinct from the separately retained geometry-assisted protocol. The DBI paper includes a sharper conditional-fidelity bound for its complete heralded instrument. The flux paper gives exact integer families and positive CS/axiodilaton mass spectra in its stated model. The broader research programme is not included merely to enlarge this archive. Selection does not certify global originality or human validation.

## Status

- Latest delivery: **29 October 2026**; original batch **28 October 2026**.
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

The first two papers include historical Wolfram Language scripts under `verification/`; their cost and coverage are described in their notes. Each new packet includes a short Python reproduction command, dependencies, inputs where required, and saved outputs. See its README for the exact command. Current build and archive checks are in [archive-checks.json](archive-checks.json), with their limits in [VERIFICATION.md](VERIFICATION.md).
