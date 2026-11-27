#!/usr/bin/env python3
"""Exact finite checks only; not a replacement for the proof in main.tex."""
from fractions import Fraction


def ceil(x):
    return -(-x.numerator // x.denominator)


for n in range(2, 100):
    exact = Fraction(sum((n - abs(k)) ** 2 for k in range(1 - n, n)), n**3)
    assert exact == Fraction(2 * n**2 + 1, 3 * n**2)

# x = Omega T W / pi, so no transcendental evaluations are needed.
# Each j represents an exact DFT mode of the length-L ladder.
for x_num in range(1, 63):
    for x_den in range(1, 17):
        x = Fraction(x_num, x_den)
        K = max(1, ceil(x))
        assert 0 < x / K <= 1
        for m in range(1, 11):
            L = m * K + 1
            packing = x * L / (2 * K)
            N = 1 + int(packing)
            assert N <= L
            assert N >= m * x / 2
            assert Fraction((N - 1) * 2 * K, L) <= x
            assert len({j % L for j in range(N)}) == N

for lo, hi in ((-10, -1), (-2, 8), (0, 1), (7, 19)):
    h = hi - lo
    for numerator in range(101):
        f = Fraction(numerator, 100)
        r = lo + f * h
        assert (1-f) * lo + f * hi == r
        assert (1-f) * (lo-r)**2 + f * (hi-r)**2 == (r-lo) * (hi-r)
        assert (r-lo) * (hi-r) / 2 <= Fraction(h*h, 8)

print("PASS: exact Parseval sums, DFT packing/residues, affine moment and variance.")
