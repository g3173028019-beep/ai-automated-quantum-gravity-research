#!/usr/bin/env python3
"""Literal +2/-2 capped zero-framed SU(2) paths versus exact finite elimination.

NumPy and SciPy; relative floating-point checks, not proof. Z(S^2 x S^1)=1.
The graph's signature is zero, including the odd paths with nonzero nullity.
"""
from math import pi, sqrt

import numpy as np
from scipy.fft import dst


def sine(vec):
    return dst(vec, type=1, norm="ortho")


def literal_paths(r, limit):
    a = np.arange(1, r, dtype=np.int64)
    port = sqrt(2 / r) * np.sin(pi * a / r)
    # Actual physical cap twists theta_a^(+2), theta_a^(-2).
    plus_twist = np.exp(1j * pi * ((a*a - 1) % (2*r)) / r)
    plus = sine(port * plus_twist)
    minus = sine(port * plus_twist.conjugate())
    leaf_pair = plus * minus
    end, interior = leaf_pair / port, leaf_pair / port**2
    answers = [np.sum(leaf_pair)]
    propagated = end
    for length in range(2, limit + 1):
        against = sine(propagated)
        answers.append(np.sum(end * against))
        propagated = interior * against
    return np.array(answers)


def eliminated_paths(r, limit):
    n = np.arange(1, r, 2)
    defects = np.minimum(n, r - n)
    U = 1 / (2 * np.sin(pi * n / (2*r))**2)
    W = 1 / np.tan(pi * n / (2*r))
    M = np.minimum.outer(defects, defects)
    half = np.ones_like(U)
    answers = [1.0]
    for length in range(2, limit + 1):
        if length % 2 == 0:
            answers.append((-1)**r * sqrt(2/r) * np.dot(W, half))
        else:
            answers.append(np.dot(U, half))
            half = M @ (U * half)
    return np.array(answers)


def kernel_relative_error(r):
    l = np.arange(1, r, 2)
    a = np.arange(1, r)
    U = 1 / (2 * np.sin(pi * l / (2*r))**2)
    sine_matrix = np.sin(pi * np.outer(a, l) / r)
    defects = np.minimum(a, r-a)
    exact = r/2 * np.minimum.outer(defects, defects)
    return np.max(np.abs((sine_matrix * U) @ sine_matrix.T - exact) / exact)


if __name__ == "__main__":
    for r in (3, 4, 17, 18, 41, 42, 127, 128, 257, 258):
        a = literal_paths(r, 10)
        b = eliminated_paths(r, 10)
        rel = np.max(np.abs(a - b) / np.maximum(1, np.abs(b)))
        kernel = kernel_relative_error(r)
        assert rel < 1e-11 and kernel < 1e-11
        print(f"r={r:3}: literal +2/-2 paths 1..10 max relative {rel:.3g}; "
              f"whole-grid kernel max relative {kernel:.3g}")

