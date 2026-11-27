#!/usr/bin/env python3
"""Replay finite algebraic constants in the TT kernel; not a dynamics simulation."""
from itertools import permutations

import sympy as s

q = [s.Matrix((1, 0, 0)),
     s.Matrix((-s.S.Half, s.sqrt(3)/2, 0)),
     s.Matrix((-s.S.Half, -s.sqrt(3)/2, 0))]
z = s.Matrix((0, 0, 1))


def E(k):
    n = z.cross(k)
    assert k.dot(k) != 0
    return z*z.T - n*n.T/k.dot(k)


def contractions(k, e):
    return [tuple(map(s.simplify, (
        -(e[a]*k[c]).dot(e[c]*e[b]*k[c]),
        s.trace(e[a]*e[c])*k[c].dot(e[b]*k[c])/2,
    ))) for a, b, c in permutations(range(3))]


m = [-E(x) for x in q]
assert contractions(q, m) == [(s.Rational(3,16),s.Rational(15,32))]*6
assert sum(s.trace(m[i]*m[j]*m[k])
           for i, j, k in permutations(range(3))) == -s.Rational(27,4)

single_terms = []
for perm in permutations(range(3)):
    k = [(q[j]+q[perm[j]])/2 for j in range(3)]
    single_terms.append(sum(sum(t) for t in contractions(k, [E(x) for x in k])))
assert single_terms == [s.Rational(-63,16),0,0,
                        s.Rational(-63,64),s.Rational(-63,64),0]
assert -s.Rational(1,4)*2*sum(single_terms)/(10*6**3) == s.Rational(7,5120)
assert s.Rational(63,32)/(10*6**3) == s.Rational(7,7680)
assert (s.Rational(2,3)*(-s.Rational(27,4)))/(2*s.Rational(63,16)) == -s.Rational(4,7)
print("Exact TT contractions, all opposite-sign pairings, short-time factors and weak-ruler ratio passed.")
