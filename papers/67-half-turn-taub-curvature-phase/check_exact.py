"""Exact half-turn detector coefficient and centered fluctuation checks."""
from collections import Counter
from fractions import Fraction as F
from itertools import product

# (monomial, integer frequency pair, coefficient):
# remove the common factors i*d1, i*d2,-i*d1 from c1,c2,bar(c1).
c1 = (("zm", (-1, 0), 1), ("Zp", (1, 0), -1))
c2 = (("ym", (0, -1), 1), ("Yp", (0, 1), -1))
C1 = (("Zm", (1, 0), 1), ("zp", (-1, 0), -1))
line = Counter()
for legs in product(c1, c2, C1, C1, C1):
    freq = tuple(sum(x[1][j] for x in legs) for j in range(2))
    if freq == (-2, 1):
        coeff = 1
        for x in legs:
            coeff *= x[2]
        line[tuple(sorted(x[0] for x in legs))] += coeff
assert line == {
    tuple(sorted(("zm", "Zm", "Yp", "zp", "zp"))): -3,
    tuple(sorted(("Zp", "zp", "Yp", "zp", "zp"))): -1,
}
D = (1 + F(9, 10)) / 2
assert D == F(19, 20)
assert F(-1, 2) + 3 * F(1, 2) == 1  # Weyl intensity constant

w = (1, -1, 2, -2)
n = (5, 1, 0, 2)
v = (2, 0, -1, 0)
dst = tuple(x - y for x, y in zip(n, v))
assert dst == (3, 1, 1, 2)
dot = lambda x, y: sum(a*b for a, b in zip(x, y))
assert dot(w, n) == dot(w, dst) == dot(w, v) == 0
assert 3*n[1]+n[0]+1 == 9 and n[0]*(n[0]-1)*(n[2]+1) == 20

I, d = (3, 1, 1, 2), (1, -1, 0, 0)
Dw = [x*y for x, y in zip(I, w)]
Dd = [x*y for x, y in zip(I, d)]
assert dot(w, I) == 0 and dot(w, Dw) == 16
sigmad = [x - y*F(dot(w, Dd), dot(w, Dw)) for x, y in zip(Dd, Dw)]
assert dot(d, sigmad) == 3
assert (sum(sigmad[:2]), sum(sigmad[2:])) == (F(3, 2), F(1, 2))
print("Curvature projection D=19/20, Weyl intensity +1, neutral line weight 9*sqrt(20).")
print("Centered band covariance coefficients: 3, (3/2,1/2) in the two frequencies.")
