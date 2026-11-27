#!/usr/bin/env python3
"""Exact finite checks for the move and magnetic zero-weight argument, not quadrature."""
from collections import Counter
from itertools import combinations, product
from fractions import Fraction

A, B = frozenset((0, 1)), frozenset(range(2, 6))
four = [A | (B - {b}) for b in B]
two = [{a} | B for a in A]
def simplices_at(side, width):
    return Counter(frozenset(t) for t4 in side for t in combinations(t4, width))
def boundary_tets(side):
    f = simplices_at(side, 4)
    return {face for face, k in f.items() if k == 1}
def triangles(tets):
    return {frozenset(t) for tet in tets for t in combinations(tet, 3)}

ext = boundary_tets(four)
assert ext == boundary_tets(two) and len(ext) == 8
external_faces = triangles(ext)
assert len(external_faces) == 16
assert {k for k in simplices_at(two, 3)} - external_faces == set()
assert set(simplices_at(four, 3)) - external_faces == {A | {b} for b in B}
two_int = {t for t, n in simplices_at(two, 4).items() if n > 1}
four_int = {t for t, n in simplices_at(four, 4).items() if n > 1}
assert two_int == {B} and four_int == {A | set(p) for p in combinations(B, 2)}
assert all(sum(t in external_faces for t in triangles([T])) == 2 for T in four_int)

def magnetic_weights(two_j):
    # SU(2): multiplicity of a scalar = number of weight-0 tuples
    # less number of weight-1 tuples at any fixed integer total spin parity.
    m = range(-two_j, two_j + 1, 2)
    tuples = list(product(m, repeat=4))
    w0 = [t for t in tuples if sum(t) == 0]
    return w0, len(w0) - sum(sum(t) == 2 for t in tuples)

for two_j, rank in [(1, 16), (2, 81)]:
    w0, d = magnetic_weights(two_j)
    assert d == two_j + 1 and d**4 == rank
    for t in w0:
        c = Counter(t)
        assert all(c[m] == c[-m] for m in c) and c[0] % 2 == 0
    assert sum(Fraction(1, 2) for _ in external_faces) == 8
    print(f"j={Fraction(two_j,2)}: exact exterior faces {len(external_faces)}; "
          f"weight-zero tuples {len(w0)}; coarse cap {d}; refined rank {rank}")
w0, _ = magnetic_weights(3)
assert (3, -1, -1, -1) in w0
print("Higher spins need the analytic anchor; no magnetic-positivity extrapolation.")
