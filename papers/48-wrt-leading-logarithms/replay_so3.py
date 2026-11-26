#!/usr/bin/env python3
"""Numerical closed SU(2)/SO(3) checks plus exact finite parity Gauss tests."""
from itertools import product
from math import sin, tan, pi, sqrt
import numpy as np


def graph(parents, w, h):
    m = len(w)
    weights = w + h + (2, -2) * (2 * m)
    edges = [(i, parents[i]) for i in range(1, m)] + [(m+i, i) for i in range(m)]
    edges += [(2*m+2*i+j, i) for i in range(2*m) for j in range(2)]
    return weights, edges


def modular(r, so):
    labels = np.arange(1, r, 2 if so else 1, dtype=int)
    S = sqrt((4 if so else 2) / r) * np.sin(pi * labels[:, None] * labels[None, :] / r)
    th = np.exp(0.5j * pi * (labels**2 - 1) / r)
    kappa = np.dot(S[0] ** 2, th) / S[0, 0]
    return S, th, kappa, labels


def tau(weights, edges, S, th, kappa, sigma, labels, pins=()):
    parents = dict(edges)
    degrees = [0] * len(weights)
    children = [[] for _ in weights]
    for c, p in edges:
        children[p].append(c)
        degrees[c] += 1
        degrees[p] += 1
    pins = dict(pins)
    dynamic = [None] * len(weights)
    for i in range(len(weights)-1, -1, -1):
        local = S[0]**(2-degrees[i]) * th**weights[i]
        if i in pins:
            local *= (labels % 2 == 1-pins[i])
        for c in children[i]:
            local *= S @ dynamic[c]
        dynamic[i] = local
    return np.sum(dynamic[0]) * kappa**(-sigma) / S[0, 0]


def signature(weights, edges):
    Q = np.diag(weights)
    for a, b in edges:
        Q[a, b] = Q[b, a] = 1
    eigen = np.linalg.eigvalsh(Q)
    return np.sum(eigen>1e-9) - np.sum(eigen<-1e-9)


def point_sum(weights, edges, pins):
    pinned = dict(pins)
    free = [i for i in range(len(weights)) if i not in pinned]
    ii = {v:j for j,v in enumerate(free)}
    F = range(len(free))
    linear = [weights[v] for v in free]
    induced = []
    for a,b in edges:
        if a in ii and b in ii:
            induced.append((ii[a],ii[b]))
        elif a in ii:
            linear[ii[a]] += 2*pinned[b]
        elif b in ii:
            linear[ii[b]] += 2*pinned[a]
    pol = np.diag(np.array(linear,dtype=int)%2)
    for a,b in induced:
        pol[a,b] = pol[b,a] = 1
    def q(x):
        return (sum(linear[j]*x[j] for j in F) + 2*sum(x[a]*x[b] for a,b in induced))%4
    radical = [x for x in product((0,1), repeat=len(F)) if not np.any((pol @ np.array(x,dtype=int)) % 2)]
    P = (sum(weights[j]*s for j,s in pinned.items()) + 2*sum(pinned[a]*pinned[b] for a,b in edges if a in pinned and b in pinned))%4
    exponent_counts = [0]*4
    for x in product((0,1), repeat=len(F)):
        exponent_counts[(P+q(x))%4] += 1
    G = complex(exponent_counts[0]-exponent_counts[2],exponent_counts[1]-exponent_counts[3])
    assert G.real**2+G.imag**2 == (2**len(F)*len(radical) if all(q(x)==0 for x in radical) else 0)
    return G


if __name__ == "__main__":
    examples = [
        ((-1,), (0,), (0,)),
        ((-1,), (1,), (1,)),
        ((-1,), (2,), (3,)),
        ((-1,0), (0,0), (0,0)),
        ((-1,0), (1,0), (0,-2)),
        ((-1,0,1,2), (0,)*4, (0,)*4)
    ]
    for parents,w,h in examples:
        ws,ed = graph(parents,w,h)
        sig = signature(ws,ed)
        for r in range(3,80,2):
            zeta = 1j if r%4 == 3 else -1j
            ks = (1+zeta)/sqrt(2)
            sem = -(2**len(w))*ks**(-sig)*zeta**sum(w+h)
            u = tau(ws,ed,*modular(r,False)[:3],sig,modular(r,False)[3])
            o = tau(ws,ed,*modular(r,True)[:3],sig,modular(r,True)[3])
            assert np.isclose(u,sem*o,rtol=4e-10,atol=1e-8),(r,u,sem*o)
            if len(w)==1 and w+h==(0,0):
                assert np.isclose(o, sum(1/tan(pi*n/(2*r)) for n in range(1,r,2))/(2*sin(pi/r)))
            if len(w)==1:
                for mask in ((), ((0,0),), ((0,1),), ((2,0),(4,1)), tuple((j,0) for j in range(len(ws)))):
                    P = point_sum(ws,ed,mask)
                    if zeta == -1j:
                        P = P.conjugate()
                    direct = tau(ws,ed,*modular(r,False)[:3],sig,modular(r,False)[3],mask)
                    assert np.isclose(direct,2**(-len(ws)/2)*ks**(-sig)*P*o,rtol=1e-10,atol=1e-8)
    print("PASS: independent closed SO(3) and SU(2) sums for all odd r=3..79 (including singular surgery); exact H formula and parity radical-mask checks")
