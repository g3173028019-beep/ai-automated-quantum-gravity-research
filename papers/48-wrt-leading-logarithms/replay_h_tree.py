#!/usr/bin/env python3
"""Numerical checks of the exact surgery formulas; no external dependencies.

The literal six-label Kirby sum is evaluated at small levels.  At larger
levels use the finite modular graph directly, not the derived Gauss formula.
"""
from cmath import exp
from itertools import product
from math import fsum, isclose, log, pi, sin, sqrt, tan

EDGES = ((0, 1), (0, 2), (0, 3), (1, 4), (1, 5))
WEIGHTS = (0, 0, 2, -2, 2, -2)
DEG = (3, 3, 1, 1, 1, 1)


def close(z, expected, rel_tol=6e-10, abs_tol=5e-10):
    assert abs(z - expected) <= max(abs_tol, rel_tol * abs(expected)), (z, expected)


def data(r):
    idx = tuple(range(1, r))
    sq = sqrt(2 / r)
    S = [[sq * sin(pi * a * b / r) for b in range(r)] for a in range(r)]
    T = [exp(pi * 1j * (a * a - 1) / (2 * r)) for a in range(r)]
    return idx, sq, S, T


def kirby_literal(r):
    idx, sq, S, T = data(r)
    s = [0] + [sin(pi * a / r) for a in idx]
    ss = [[sin(pi * a * b / r) for b in range(r)] for a in range(r)]
    weight = [[T[a] ** w * s[a] ** (2 - d) if a else 0 for a in range(r)]
              for w, d in zip(WEIGHTS, DEG)]
    val = 0j
    for labels in product(idx, repeat=6):
        term = 1 + 0j
        for j, a in enumerate(labels):
            term *= weight[j][a]
        for j, k in EDGES:
            term *= ss[labels[j]][labels[k]]
        val += term
    f_graph = val / s[1] ** 7
    f_pos = sum(T[a] * (s[a] / s[1]) ** 2 for a in idx)
    f_neg = sum(T[a].conjugate() * (s[a] / s[1]) ** 2 for a in idx)
    close(abs(f_pos), 1 / S[1][1])
    return f_graph / (f_pos ** 3 * f_neg ** 3)


def graph_direct(r):
    idx, sq, S, T = data(r)
    leaf_p = [sum(S[a][k] * T[k] ** 2 * S[k][1] for k in idx) for a in idx]
    leaf_m = [sum(S[a][k] * T[k] ** -2 * S[k][1] for k in idx) for a in idx]
    for a in idx:
        expected = 4 * ((a + r) % 2) * sin(pi * a / (2 * r)) ** 2 / r
        close(leaf_p[a - 1] * leaf_m[a - 1], expected)
    u = [0] + [leaf_p[a - 1] * leaf_m[a - 1] / S[1][a] for a in idx]
    s_u = [0] + [sum(S[b][a] * u[a] for a in idx) for b in idx]
    for b in idx:
        close(s_u[b], (-1) ** (b + 1))
    return sum(u[a] * s_u[a] for a in idx)


def cot_sum(r):
    return fsum(1 / tan(pi * n / (2 * r)) for n in range(1, r, 2))


if __name__ == "__main__":
    for r in range(3, 9):
        exact_tau = (-1) ** r * cot_sum(r) / sin(pi / r)
        close(kirby_literal(r), exact_tau)
    for r in range(3, 100):
        exact_z = (-1) ** r * sqrt(2 / r) * cot_sum(r)
        close(graph_direct(r), exact_z)
    euler = 0.57721566490153286061
    for r in (10000, 10001, 100000, 100001):
        finite_part = pi * cot_sum(r) / r - log(r)
        assert isclose(finite_part, euler + log(4 / pi), abs_tol=1e-7)
    print("PASS: literal surgery (levels r=3..8); modular graph (3..99); both-parity logarithmic finite part")
