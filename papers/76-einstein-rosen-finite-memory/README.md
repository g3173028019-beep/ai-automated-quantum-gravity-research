# Paper 76: Sharp finite-memory laws for quantum Einstein–Rosen waves

- [Manuscript PDF](manuscript.pdf)
- [Complete LaTeX source](main.tex)
- [Verification guide](VERIFICATION.md)
- [Optional exact arithmetic replay](check_finite_identities.py)

For the exact polarized Einstein–Rosen radial-mode process, a single simulator that works for all admitted preparations needs persistent memory dimension of order \(1+A g_LTW/\sqrt{\varepsilon}\). \(T\) is the total permitted physical waiting time, \(g_L\) is the span of the retained \(e^{-bH_{0,L}}\) spectrum, and \(W\) is the admitted width of the discarded spectral variable. The positive classical construction keeps one randomly rounded label throughout the run. The converse covers arbitrary finite quantum memories, including mixed programs, with the same engine and the whole adaptive tester class. It is not a bound on a simulator chosen afresh for each preparation.

A separate two-window experiment uses one fixed smooth discarded-mode state and no tester quantum ancilla or feedback. It forces linear time growth at fixed accuracy, but does not settle the small-error exponent for a one-preparation simulator. The split is between radial modes, not spatial regions. The bounds charge free physical time and persistent memory, not the resolution or gate cost of the measuring apparatus.

[Papers 18](../18-positive-boundary-sewing), [53](../53-coherent-quantum-shear-sewing), and [68](../68-einstein-rosen-spectral-memory) are preserved. Paper 68 supplies the exact spectral process and distinguishes exact infinite-rank memory from memoryless resets. Paper 76 supplies a uniform finite-budget approximation law and a two-window quantitative converse.

The manuscript compares the programming mechanism with [Hashimoto–Yoshida–Murao, arXiv:2611.02441](https://arxiv.org/abs/2611.02441) and earlier mixing/programming results. The distinctive claims are the continuous physical process upper bound for correlated preparations, the joint time/error lower bound for actual finite-energy mode probes, and the separate fixed-state witness.
