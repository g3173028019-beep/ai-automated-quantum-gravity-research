#!/usr/bin/env python3
"""Finite DST replay of the corona plumbing. NumPy/SciPy required.

This is only a replay of exact modular sums, not a basis for the proof.
The value printed is (2/log r)^m times the factored, phase-removed graph
contraction, whose theorem limit is the volume of the log-charge polytope.
"""
from math import log, pi, sqrt
from itertools import product
import numpy as np
from scipy.fft import dst


def phase_removed_literal(r, parents, w, h):
    """Literal contraction, after +/-2 Gauss leaves and overall factors."""
    m = len(parents)
    n = np.arange(1, r, dtype=float)
    x = n / r
    odd = (n % 2 == 1)
    ortho = sqrt(2 / r)
    b = (pi / (2 * r)) / np.tan(pi * x / 2)

    def sine(v):
        return dst(v, type=1, norm="ortho") / ortho

    children = [[] for _ in range(m)]
    for i in range(1, m):
        children[parents[i]].append(i)
    down = [None] * m
    for i in range(m - 1, -1, -1):
        d = len(children[i]) + (i != 0)
        f = np.cos(pi * x / 2) ** 2 * (pi * x / np.sin(pi * x)) ** (d + 1)
        k = (4 / pi) * sine(odd * b * np.exp(0.5j * pi * h[i] * n * n / r))
        down[i] = odd * f / n ** (d + 1) * np.exp(0.5j * pi * w[i] * n * n / r) * k
        for j in children[i]:
            down[i] *= (r / pi) * sine(down[j])
    return np.sum(down[0]) * (2 / log(r)) ** m


def small_literal(r, parents, w, h):
    """Independently evaluate the actual six/several-vertex Kirby graph directly."""
    m = len(parents)
    nback = 2 * m
    edges = [(i, parents[i]) for i in range(1, m)] + [(i, m + i) for i in range(m)]
    degrees = [0] * nback
    for a, b in edges:
        degrees[a] += 1
        degrees[b] += 1
    labels = np.arange(1, r, dtype=int)
    S = sqrt(2 / r) * np.sin(pi * labels[:, None] * labels[None, :] / r)
    th = np.exp(0.5j * pi * (labels**2 - 1) / r)
    plus = S @ (th**2 * S[:, 0])
    minus = S @ (th**-2 * S[:, 0])
    weights = [plus * minus * th**v / S[:, 0]**d for v, d in zip(w + h, degrees)]
    full = 0j
    for cs in product(range(r - 1), repeat=nback):
        elem = np.prod([weights[i][a] for i, a in enumerate(cs)])
        elem *= np.prod([S[cs[i], cs[j]] for i, j in edges])
        full += elem
    W = sum(w + h)
    phase = (-1) ** (r + W) * np.exp(0.5j * pi * W * (r * r - 1) / r)
    den = phase * 2 ** (m + 0.5) * pi ** (1 - 2*m) * r ** (2*m - 1.5)
    return full / den * (2 / log(r)) ** m


if __name__ == "__main__":
    cases = [
        ("H untwisted", (-1,), (0,), (0,), 1),
        ("H one spine twist", (-1,), (2,), (0,), 1/2),
        ("H one pendant twist", (-1,), (0,), (2,), 1/2),
        ("H both twisted", (-1,), (1,), (-1,), 0),
        ("m=2 chain", (-1, 0), (0, 0), (0, 0), 1/2),
        ("m=2 pendant clipped", (-1, 0), (0, 0), (1, 0), 1/8),
        ("m=2 spine clipped", (-1, 0), (1, 0), (0, 0), 3/8),
        ("m=2 adjacent pendants", (-1, 0), (0, 0), (1, -1), 0),
        ("m=4 path", (-1, 0, 1, 2), (0,)*4, (0,)*4, 5/24),
        ("m=4 star", (-1, 0, 0, 0), (0,)*4, (0,)*4, 1/4),
    ]
    for tag, p, w, h, volume in cases:
        if len(p) <= 2:
            for r in range(3, 11):
                assert abs(phase_removed_literal(r, p, w, h)-small_literal(r, p, w, h)) < 1e-9
        a = phase_removed_literal(4096, p, w, h)
        b = phase_removed_literal(65536, p, w, h)
        print(f"{tag:25s} r=4096 {a:+.6f}, r=65536 {b:+.6f}, limit {volume:.6f}")
