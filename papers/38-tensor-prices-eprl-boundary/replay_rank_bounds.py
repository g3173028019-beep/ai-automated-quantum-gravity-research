"""Exact-arithmetic checks of the rank bounds in manuscript.tex.

This checks the closed formulas, not the oriented tensor contractions or
their representation-theory proof. No external Python packages are used.
"""
from fractions import Fraction as F
from math import gcd


def bound_twospin(alpha: int, beta: int, cutoff: int):
    """alpha=2*j_s^+, beta=2*j_s^- for a minimal positive spin s."""
    assert alpha > beta >= 1 and gcd(alpha, beta) == 1 and cutoff >= 1
    dims = [(alpha * m + 1) * (beta * m + 1) for m in range(cutoff + 1)]
    volume = sum(d * d for d in dims)
    first = sum(x * y for x, y in zip(dims, dims[1:]))
    second = sum(x * y for x, y in zip(dims, dims[2:]))
    neighbors = [((dims[m - 1] if m else 0) +
                  (dims[m + 1] if m < cutoff else 0))
                 for m in range(cutoff + 1)]
    bnum = 2 * volume - 1 - dims[-1] ** 2 + 2 * second
    assert bnum == sum(d * d for d in neighbors)
    assert 2 * first == sum(d * z for d, z in zip(dims, neighbors))
    assert 4 * first ** 2 <= bnum * volume

    r0 = alpha + beta + 1
    dim = dims[1]
    delta = F(1, 5) if (alpha, beta) == (2, 1) else F(1)
    anum = r0 * volume - (r0 - 1) + 2 * second
    a = F(anum, dim * dim * volume)
    c = F(2 * first, dim * volume)
    b_lower = (bnum - delta * (volume - 1)) / (dim * dim * volume)
    b_upper = F(bnum, dim * dim * volume)

    diff_lower = F((r0 - 2) * (volume - 1) + dims[-1] ** 2,
                   dim * dim * volume)
    det_lower = ((r0 - 2 - 3 * delta) * (volume - 1) + dims[-1] ** 2
                 ) / (dim * dim * volume)
    assert a - b_upper == diff_lower > 0
    assert a + 3 * b_lower - 4 * c * c >= det_lower > 0
    return a, (b_lower, b_upper), c


def min_pair_bound(m: int):
    t = F((3*m + 4) * (2*m + 1) * (m + 1),
          (3*m + 1) * (2*m + 3) * (m + 2))
    assert F(2, 3) <= t <= 1
    return t * t / (t * t + (1 - t) ** 2)


if __name__ == "__main__":
    count = 0
    for a in range(2, 12):
        for b in range(1, a):
            if gcd(a, b) != 1:
                continue
            for cap in range(1, 51):
                bound_twospin(a, b, cap)
                count += 1
    for m in range(1000):
        assert min_pair_bound(m) >= F(4, 5)
    print(f"PASS: {count} finite rank-bound checks; 1000 min-pair bounds.")
    for cap in (1, 2, 10, 100):
        a, b, c = bound_twospin(2, 1, cap)
        print(f"gamma=1/3, N={cap}: a={a}, b in [{b[0]}, {b[1]}], c={c}")

