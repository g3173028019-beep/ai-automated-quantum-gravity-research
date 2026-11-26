"""Exact phase-exponent sanity checks for the torus Wilson formulas.

No floating-point values of roots of unity are used.  A monomial key (h,a,b)
represents eta**h U**a V**b, with eta = exp(pi*i/(2*r)) and UV=eta**2 VU.
These checks do not replace the analytic support/count proof.
"""

from collections import Counter
from itertools import product


def term(r, h, a, b):
    return Counter({(h % (4*r), a % (2*r), b % (2*r)): 1})


def linear(*pairs):
    out = Counter()
    for c, x in pairs:
        for key, val in x.items():
            out[key] += c*val
    return Counter({key: val for key, val in out.items() if val})


def mul(r, x, y):
    out = Counter()
    for (h, a, b), c in x.items():
        for (i, d, e), f in y.items():
            key = ((h+i-2*b*d) % (4*r), (a+d) % (2*r), (b+e) % (2*r))
            out[key] += c*f
    return Counter({key: val for key, val in out.items() if val})


def phase(r, x, h):
    return mul(r, term(r, h, 0, 0), x)


def L(r, a, b):
    return linear((1, term(r, -a*b, a, b)), (1, term(r, -a*b, -a, -b)))


def adjoint(r, x):
    return Counter({((-h-2*a*b) % (4*r), -a % (2*r), -b % (2*r)): c
                    for (h, a, b), c in x.items()})


def parity(r, x):
    return Counter({(h, -a % (2*r), -b % (2*r)): c for (h, a, b), c in x.items()})


for r in (3, 4, 5, 9, 12, 37):
    for a, b in product(range(-2, 3), repeat=2):
        x = L(r, a, b)
        assert adjoint(r, x) == x
        assert parity(r, x) == x
        for c, d in product(range(-2, 3), repeat=2):
            det = a*d-b*c
            expected = linear((1, phase(r, L(r, a+c, b+d), det)),
                              (1, phase(r, L(r, a-c, b-d), -det)))
            assert mul(r, x, L(r, c, d)) == expected

    for n in range(1, r):
        for step in (-1, 1):
            m = n+step
            if not 1 <= m < r:  # missing odd basis entries
                continue
            dehn_exponent = (m*m-1) - (n*n-1)
            weyl_exponent = 2*step*n+1
            assert (dehn_exponent-weyl_exponent) % (4*r) == 0

    for b in range(2*r):
        assert sum((2*n+b) % (2*r) == 0 for n in range(2*r)) <= 2

print("PASS: exact Weyl product-to-sum, adjoint/parity, Dehn-twist phases, and odd trace support")
