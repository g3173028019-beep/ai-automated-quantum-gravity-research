# Heralded DBI controlled-phase gate

The packet constructs a fixed-time radial-return heralded CZ in a specified independent-field DBI half-line model. At dimensionless time `107.61733846478977`, every logical input has success probability at least `0.6565527711781338` and conditional squared fidelity at least `0.96400610918152`, including an arbitrary reference. The stated bounded source/readout/timing/phase box retains success at least `0.6351321846638459` and squared fidelity at least `0.9587714318247369`.

Read [manuscript.pdf](manuscript.pdf) or [main.tex](main.tex). The construction retains a complete P/Q failure flag; it is conditional, not deterministic. Exact preservation of the four-state code, a common initially uncorrelated radial source, a normalized rank-one radial readout mode and local phase controls are explicit assumptions. Their physical implementation costs are not supplied. Standard corrective `Rz` angle is minus the target phase difference, not plus it.

## Replay

With Python 3 and the packages in `requirements.txt` already installed, run:

```text
python replay.py
```

This default command is entirely local and network-free. It has no original project paths and writes `results/gate.json`. For environment setup, an optional `python -m pip install -r requirements.txt` can install the packages; installation is separate from the network-free replay.

The replay checks exact algebra of the sharp conditional fidelity bound `F_squared >= 1-eta^2`, input hashes, the saved continuum residual sums/tails, the exact initial mismatch, and all frozen-time operator and gate bounds using 110-digit outward interval arithmetic. The numerical sharpness example is not used to choose the protocol. No old continuum refinement is run.

The four JSON files under `inputs/` are exact copies of the inherited trial and complete residual partitions. The unmodified historical `reference/certify_trial_residual.py` explains the cellwise interval algorithm and its integrator hash; it is not the replay entry point and is not executed by the default command. Initial source mismatch and the analytic large-radius tail are retained. Recombining saved upper bounds does not freshly regenerate their integrand enclosures.

The JSON distinguishes exact decimal distance constants from float endpoint displays. The former define the printed fidelity guarantees; the latter are rounded outward once more for safe machine display. The robust constant is certified directly from the branch bound plus the control budget, before a redundant outward rounding of the printed nominal constant.

## Build

With an existing TeX installation:

```text
latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=build main.tex
```

The delivered `manuscript.pdf` is copied from `build/main.pdf`. The bibliography records the supplied internal source artifacts rather than inventing external provenance or a public priority claim.

## Scope and chronology

The Hamiltonians and finite nonlinear parameters are held fixed. Source/readout norm errors cover normalized pure vectors or an explicitly equivalent dilation guarantee. They do not transfer unbounded-density or RMS-domain results to arbitrary norm-close radial states. No Hamiltonian-parameter uncertainty, autonomous readout apparatus, universal state-loading interface, spatial-locality theorem, or quantum-gravity completion is claimed.

Retained scientific inputs are dated 19 September 2026. The delivery date is 29 October 2026. Earlier work was Fable 5-led; July 2026 onward results were predominantly produced by GPT models, as reported by the operator. DeepSeek and multiple GPT models also contributed cross-review and manuscript preparation. This is reported chronology, not a per-result authorship allocation. Human names and affiliations are blank, and human expert review has not occurred. Historical source dates, delivery date and demonstrable public disclosure are distinct; a private commit does not establish public priority.

"I genuinely do not know how authorship, intellectual credit, or ownership should be attributed for these results."
