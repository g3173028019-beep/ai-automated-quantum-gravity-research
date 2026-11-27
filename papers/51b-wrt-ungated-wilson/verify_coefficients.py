#!/usr/bin/env python3
"""Replay the exact rational integrals in the 51B coefficient."""
from fractions import Fraction as Q
from itertools import zip_longest


def chebyshev_second(count):
    prev, cur = [0], [1]
    out = [cur]
    for _ in range(count):
        nxt = [a-b for a,b in zip_longest([0]+[2*c for c in cur],
                                          prev, fillvalue=0)]
        prev, cur = cur, nxt
        out.append(cur)
    return out


U = chebyshev_second(48)
for J in range(25):
    P = U[2*J][:]
    P[0] -= (-1)**J
    assert P[0] == 0 and not any(P[1::2])
    # Fubini in the split rho-integral makes each z^p contribute c/p.
    direct = sum((Q(c,p) for p,c in enumerate(P) if p), Q(0))
    odd = 2*(-1)**(J+1)*sum((Q(1,j) for j in range(1,J+1,2)), Q(0))
    assert direct == odd, (J,direct,odd)
    assert sum((Q(c,p+1) for p,c in enumerate(U[2*J])), Q(0)) == Q(1,2*J+1)
    if J == 1:
        assert P == [0,0,4] and direct == 2
        assert Q(4,3)*(1-Q(1,2)) == Q(2,3)
        assert 1-Q(1,2) == Q(1,2)
print("Verified exact integrals for integer spins 0 through 24.")
print("Spin one: full 2; window [1,2] centered 2/3 and vacuum 1/2.")
