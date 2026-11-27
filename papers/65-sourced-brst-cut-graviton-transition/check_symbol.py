"""Exact rational check of the six coefficients in the spatial Einstein cubic."""
from fractions import Fraction as F
from itertools import permutations

K = ((1, 0, 0), (0, 1, 0), (-1, -1, 0))
P = (
    ((0, 0, 0), (0, 1, 0), (0, 0, -1)),
    ((1, 0, 0), (0, 0, 0), (0, 0, -1)),
    ((F(1, 2), F(-1, 2), 0), (F(-1, 2), F(1, 2), 0), (0, 0, -1)),
)

def matvec(m, v):
    return [sum(x * y for x, y in zip(row, v)) for row in m]

def dot(a, b):
    return sum(x * y for x, y in zip(a, b))

def trace_product(a, b):
    return sum(a[i][j] * b[j][i] for i in range(3) for j in range(3))

for p, k in zip(P, K):
    assert sum(p[i][i] for i in range(3)) == 0
    assert matvec(p, k) == [0, 0, 0]
assert [sum(k[i] for k in K) for i in range(3)] == [0, 0, 0]

checks = []
for a, b, c in permutations(range(3)):
    k = K[c]
    val = -2 * dot(matvec(P[a], k), matvec(P[c], matvec(P[b], k)))
    val += trace_product(P[a], P[c]) * dot(k, matvec(P[b], k))
    checks.append(F(val))

assert checks == [F(5, 2), F(1, 2), F(5, 2), F(1, 2), F(3, 2), F(3, 2)]
assert sum(checks) == 9
print("Ordered coefficients:", ", ".join(map(str, checks)))
print("Exact symmetrized coefficient:", sum(checks))
