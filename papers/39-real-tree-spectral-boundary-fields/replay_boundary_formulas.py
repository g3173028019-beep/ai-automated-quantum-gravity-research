#!/usr/bin/env python3
"""Exact checks of selected tree-hull and critical shell identities."""
from fractions import Fraction
from itertools import product


def leaves(p, depth):
    return list(product(range(p), repeat=depth))


def prefix(a, b):
    return next((k for k, (x, y) in enumerate(zip(a, b)) if x != y), len(a))


def hull_length(group):
    """Unit-edge length of the minimal subtree spanning the chosen leaves."""
    if len(group) <= 1:
        return 0
    n = len(group[0])
    lca = min(prefix(group[0], a) for a in group[1:])
    return len({a[:k] for a in group for k in range(lca + 1, n+1)})


def cut_hull_length(points):
    """Sum leaf-induced split crossings in the same finite hull."""
    n = len(points[0])
    edges = {a[:k] for k in range(1, n+1) for a in points}
    return sum(0 < sum(a[:len(e)] == e for a in points) < len(points)
               for e in edges)


def critical_covariance(p, depth, cylinder_text):
    """Compute critical covariance of two indicator probes by leaf summation."""
    a, b = cylinder_text
    layer = leaves(p, depth)
    return Fraction(sum(p**prefix(x,y)
                        for x in layer if x[:len(a)] == a
                        for y in layer if y[:len(b)] == b), p**(2*depth))


for p in (2,3):
    for depth in range(1,5):
        layer = leaves(p, depth)
        for x in (layer[0], layer[-1]):
            group = [x, layer[len(layer)//2],layer[0]]
            assert hull_length(group) == cut_hull_length(group)
        same = critical_covariance(p, depth, [(0,),(0,)])
        across = critical_covariance(p, depth, [(0,),(1,)])
        assert same == Fraction(1, p) * (1 + Fraction(p-1,p)*(depth-1))
        assert across == Fraction(1, p*p)
        full = critical_covariance(p, depth, [(),()])
        assert full == 1 + Fraction(p-1,p)*depth
        print(f'p={p} depth={depth}: whole-boundary critical variance / lambda = {full}')
print('Known finite identities verified exactly. No Monte Carlo is used.')
