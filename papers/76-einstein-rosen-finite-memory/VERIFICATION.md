# Proof and build checks

- **Exact process (§1):** a spectral projection of \(R\) commutes with the Hamiltonian and all retained instruments. Orthogonality removes off-diagonal high blocks when tracing out the discarded Fock space, including for retained–high correlated initial states. A common high label must be used throughout a history.
- **CP approximation (§2):** the two adjacent weights are positive, add to one, and reproduce \(r\); the square-root spectral instruments give the positive initialization. One backward-effect differentiation costs the actual remaining time, while the second derivative obeys the square of the same budget. Taylor interpolation then costs at most \((A g_LTh)^2/8\) in event probability. This is uniform in the number of cuts and does not average time across branches.
- **Universal lower bound (§3):** the low reset precedes the choice of fresh data. The physical coherent dwell ladder, on smooth packets in the admitted support, realizes a fixed processor with approximate unitary channel for each normal calibration preparation. Verify that all data/ref inputs satisfy the diamond guarantee, that the program is the only preparation-dependent object after reset, and that the purification is not accessed in the CPTP recycling decoder. The integer weight ladder and Fourier code give an exact finite family with orthogonal ideal outputs. A single program POVM then bounds mean success by \(d/N\).
- **One fixed preparation (§4):** the smooth periodic normalization of the high packet gives a genuinely uniform wrapped law. Identical translated low spectral packets make the lost mass exactly \(1-|\phi|^2\); it is included in the POVM, not postselected. After the middle replacement, the guessing bound uses \(d\), not \(d^2\). The \(3\eta^2+2\varepsilon\) loss includes the agreement event and the first marginal separately. The optional finite-particle extension uses the exact total-energy/simplex Jacobian, not a change of physical clock.
- **Nearest prior art (§5):** the generic barycentric and recycling arguments are credited; the sharp bounds here concern the stated physical process and control model. The proof does not identify an optimal precision bound for one fixed initialization.

Run the following for exact arithmetic tests of the Fourier-code packing, wrapped-law Parseval sum, and interpolation weights and constants. These finite checks do not replace the analytic continuity or control assumptions.

    python3 check_finite_identities.py

The self-contained manuscript needs standard pdfLaTeX (amsmath, amssymb, amsthm, mathtools, geometry, lmodern, microtype, hyperref). For example:

    pdflatex -interaction=nonstopmode -halt-on-error -jobname=manuscript main.tex
    pdflatex -interaction=nonstopmode -halt-on-error -jobname=manuscript main.tex
    pdflatex -interaction=nonstopmode -halt-on-error -jobname=manuscript main.tex
