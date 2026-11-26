"""Finite unframed tree sewing checks (not used to prove the asymptotics)."""
from itertools import product
from math import fsum, pi, prod, sin, sqrt
import numpy as np


def prob_zero(labels):
    shift = [(n - 1) // 2 for n in labels]
    return sum(sum(k) == 0 for k in product(*(range(-h, h + 1) for h in shift))) / prod(labels)


def finite_bracket(r, labels):
    return (4 / r) * fsum(
        float(np.cos(pi*l/(2*r))**2) *
        prod(sin(pi*n*l/r)/(n*sin(pi*l/r)) for n in labels)
        for l in range(1, r, 2)
    )


def tree_contract(adj, weight, S):
    def visit(v, p):
        z = weight(len(adj[v])).copy()
        for a in adj[v]:
            if a != p:
                z *= S @ visit(a, v)
        return z
    return np.sum(visit(0, -1))


def finite_tree(r, adj):
    a = np.arange(1, r)
    S = sqrt(2/r)*np.sin(pi*np.outer(a, a)/r)
    cap = S @ (S[:, 0]*np.exp(1j*pi*(a**2-1)/r))
    actual = tree_contract(adj, lambda degree: abs(cap)**2 / S[:, 0]**degree, S)
    n = np.arange(1, r, 2)
    R = (-1)**r*sqrt(2/r)*np.sin(pi*np.outer(n,n)/r)
    z = tree_contract(
        adj, lambda degree: (4/r)*np.cos(pi*n/(2*r))**2 /
                           (sqrt(2/r)*np.sin(pi*n/r))**degree, R
    )
    return actual, z


if __name__ == "__main__":
    checks = 0
    for r in (7, 8, 15, 16, 31, 32):
        for d in (1, 2, 3, 4):
            for labels in product((1, 3, 5, 7), repeat=d):
                if 1 + sum(n-1 for n in labels) >= r:
                    continue
                assert abs(finite_bracket(r, labels)-prob_zero(labels)) < 2e-12
                checks += 1

    # Three-leaf star and a branched five-node tree.
    trees = [
        ((1,2,3), (0,), (0,), (0,)),
        ((1,2,3), (0,4), (0,), (0,), (1,)),
    ]
    for adj in trees:
        for r in (7, 8, 31, 32, 127, 128):
            x, y = finite_tree(r, adj)
            assert abs(x-y) < 3e-12*max(1, abs(y))
    print(f"PASS: {checks} no-alias spin multiplicities; branched tree S/T sewing on both parities")
