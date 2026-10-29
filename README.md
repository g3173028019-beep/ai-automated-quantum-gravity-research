# AI-Automated Quantum Gravity Research

English working manuscripts and reproducible calculations from an AI-led research programme in quantum gravity, quantum information and mathematical physics.

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

6. **Positive-Energy Local Clocks for a Weighted Maxwell Geometry-Mediated Gate**
   - [English PDF](papers/06-maxwell-local-clocks/manuscript.pdf)
   - [LaTeX source](papers/06-maxwell-local-clocks/main.tex)
   - [Reproduction and verification](papers/06-maxwell-local-clocks/README.md)
7. **A Local Clock-Foliated Repair of the Inhomogeneous Spectrum of a Rotating Einstein-Static Scalar Reference**
   - [English PDF](papers/07-clock-foliated-reference/manuscript.pdf)
   - [LaTeX source](papers/07-clock-foliated-reference/main.tex)
   - [Reproduction and verification](papers/07-clock-foliated-reference/README.md)
8. **Native Mixed Quartic Interactions on an Einstein–Scalar Relative Carrier and an Exact Matter-Supported Walker Geometry**
   - [English PDF](papers/08-einstein-scalar-mixed-quartic/manuscript.pdf)
   - [LaTeX source](papers/08-einstein-scalar-mixed-quartic/main.tex)
   - [Reproduction and verification](papers/08-einstein-scalar-mixed-quartic/README.md)

The first 29 October packet contains JT mean-energy recovery, a DBI heralded gate and exact integer flux families. The subsequent priority packet adds local autonomous Maxwell control, a clock-foliated effective reference model with positive inhomogeneous quadratic energy, and native Einstein–scalar quartic/backreaction constructions. Each paper states its model and assumptions. Delivery order is not a ranking of the broader research, a complete corpus inventory, or certification of global originality or human validation.

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
