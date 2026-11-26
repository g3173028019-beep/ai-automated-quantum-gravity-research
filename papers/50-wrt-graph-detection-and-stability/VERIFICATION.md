# Proof and reproduction notes

## Claim-to-proof map

| Claim | Where proved | Key checks |
|---|---|---|
| Exact oriented WRT charge sum; finite square reduction | §2 | Actual cap and seam cobordisms, exact Maslov correction for singular adjacency as well, odd Gauss cap identity, even/odd reflected signs and all-level min transfer. |
| Fixed-connector fourth logarithm | Theorem 3.1 | Entire finite-cap/lattice correction through this order is retained; the tied deletion invariants cancel. Ordered-gap prototype certificate uses full signed monomials, not just alternating-cycle counts. |
| Leading localization and first eligible motif comparison | Proposition 4.1 and Theorem 4.2 | Exact common reflected measure and conditional transfer eliminate lower matching sectors before a limit. Neumann centering removes leaves, with strict positivity proved also for leafless motifs containing bridges. |
| Joint modulus sign reversal | Theorem 4.3 | Degree-five finite contrast, hard and pin cancellations, adaptive tail at arbitrarily slowly diverging graft count, unmarked hard-sixth and marked soft-fifth constants. No sign at the exact critical prefactor. |
| Signed rigorous density saddle | Theorem 5.1 | Sequential limits. Full exterior is bounded in absolute value; neighborhood Gamma gate avoids an invalid unrestricted product exponential moment. Refined saddle compared after exact sitewise, not pair-density, normalization. |
| Linear-degree counterexample | Proposition 5.2 | Refined quotient treated as a bare template. The `(2,2)` and `(4)` masks have quotient signature −2 and identical dimensional prefactors. The signed leading coefficient ratio is `99/(70 sqrt(2))`. |
| Finite-level caveat | End of §5 | The two-pinned-class tripartite sector is a signed-sum diagnostic, not a bound on the total or a simultaneous limit theorem. |

## Exact finite certificates

From the paper directory (Python 3.10+, and a C++17 compiler):
```sh
cd replay
python3 verify_wrt_universal_twin_family.py
python3 verify_wrt_motif_diamond_pair.py
python3 verify_wrt_motif_cancellation_sixth.py
python3 verify_wrt_motif_soft_sign_transition.py
```

The first replay compiles the bundled C++17 enumerators in a temporary folder. It checks the common low-order counting data and the exact signed third-gap signature `g! g(g+1)` times the cycle-minus-twice-path-plus-matching contrast, including the explicit small graft limits. Other replays use rational arithmetic and verify the cycle-matched diamond weight `1/467775`, its modulus coefficient `7072/4455`, the sixth hard coefficient `2449558508/945`, and marked contraction `1/1890` with multiplier `27010/63` of the negative soft integral. Run transcripts are included for convenience.

All listed replays succeeded. The analytic interpretation and limiting claims are proved in the manuscript, not inferred from numerical fits. The integrated source was cross-read for normalization, parity, the finite signed subtraction, the threshold scope and the density bounds. The PDF was compiled through stable references and checked visually at all 19 pages. This documents internal verification, not journal peer review.
