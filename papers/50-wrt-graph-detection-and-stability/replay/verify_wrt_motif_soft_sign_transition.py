"""Exact multigraph hard polynomial through z^4 after marking a contraction."""
from fractions import Fraction as Q
from math import comb, factorial
from verify_wrt_motif_cancellation_sixth import DH, DG, CH, CG

def hard(edges, n=11, qmax=4):
    nbh = [[] for _ in range(n)]
    for a,b in edges:
        assert a != b
        nbh[a].append(b)
        nbh[b].append(a)
    F = [[0]*(qmax+1) for _ in range(1<<n)]
    F[0][0] = 1
    maxdeg = max(map(len, nbh))
    coef = {(m,q,d):[(-1)**l*comb(d,l)*(factorial(m+q-1)//factorial(m+q-1-l))
                        for l in range(min(d,q)+1)]
            for m in range(1,n+1) for q in range(qmax+1)
            for d in range(maxdeg+1)}
    for S in range(1,1<<n):
        m, bits = S.bit_count(), S
        while bits:
            b = bits & -bits
            bits ^= b
            i = b.bit_length()-1
            Sm = S ^ b
            d = sum((Sm>>j)&1 for j in nbh[i])  # repeats are parallel edges
            for q in range(qmax+1):
                F[S][q] += sum(c*F[Sm][q-l] for l,c in enumerate(coef[m,q,d]))
    return tuple(Q(F[-1][q], factorial(n+q)) for q in range(qmax+1))

def contract(edges, mark):
    u,v = mark
    ren = {x:i for i,x in enumerate(x for x in range(12) if x != v)}
    ren[v] = ren[u]
    return [(ren[a],ren[b]) for a,b in edges if (a,b) != mark]

def marked(edges):
    return tuple(map(sum,zip(*(hard(contract(edges,e)) for e in edges))))

D, C = [tuple(x-y for x,y in zip(marked(H),marked(G)))
        for H,G in ((DH,DG),(CH,CG))]
assert D == (0,)*4+(Q(1,4725),)
assert C == (0,)*4+(Q(-1,1890),)
total = 5*D[4]+C[4]
assert total == Q(1,1890)
factor = 2*total*Q(factorial(75),factorial(72))
assert factor == Q(27010,63)
assert Q(2449558508,945)/factor == Q(453454,75)
print("marked z^4 differences:", D[4], C[4], "composite:", total)
print("relative soft multiplier of w0:", factor)
