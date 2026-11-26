"""Exact finite counts and polynomial integral used in the motif theorem."""
from fractions import Fraction as Q
from math import comb, factorial

H = [(0,1),(0,4),(0,8),(1,4),(1,9),(2,5),(2,7),(2,11),(3,6),
     (3,7),(3,11),(4,9),(5,8),(5,11),(6,7),(6,10),(8,10),(9,10)]
G = [(0,4),(0,5),(0,9),(1,2),(1,9),(1,11),(2,6),(2,9),(3,4),
     (3,7),(3,8),(4,5),(5,7),(6,10),(6,11),(7,8),(8,10),(10,11)]

def classify(edges):
    A = [set() for _ in range(12)]
    for u, v in edges:
        assert v not in A[u] and u != v
        A[u].add(v)
        A[v].add(u)
    assert all(len(row) == 3 for row in A)
    seen, front = {0}, [0]
    while front:
        u = front.pop()
        for v in A[u] - seen:
            seen.add(v)
            front.append(v)
    assert len(seen) == 12
    cycles = [0] * 6
    for o in range(12):
        def walk(P):
            if len(P) >= 3 and P[1] < P[-1] and o in A[P[-1]]:
                cycles[len(P)] += 1
            if len(P) < 5:
                for v in A[P[-1]]:
                    if v > o and v not in P:
                        walk(P + [v])
        walk([o])
    diamonds = sum(comb(len(A[u] & A[v]), 2) for u, v in edges)
    return tuple(cycles[3:]), diamonds

def mul(A, B):
    from collections import defaultdict
    R = defaultdict(Q)
    for (i, j), u in A.items():
        for (h, k), v in B.items():
            R[i+h, j+k] += u*v
    return R

assert classify(H) == ((4, 2, 4), 1)
assert classify(G) == ((4, 2, 4), 0)
C = {(2,0):Q(1,2), (0,2):Q(1,2), (0,1):Q(-1), (0,0):Q(1,3)}
C2 = {key:Q(-value,360) for key, value in {
    (4,0):15, (2,2):90, (2,1):-180, (2,0):60,
    (0,4):15, (0,3):-60, (0,2):60, (0,0):-8}.items()}
D = mul(C, mul(C2, C2))
weight = 2 * sum(q / ((i+1)*(i+j+2)) for (i,j), q in D.items())
assert weight == Q(1,467775)
out = weight*Q(factorial(17),factorial(12))
assert out == Q(7072,4455)
print("cycle counts and diamond counts:", classify(H), classify(G))
print("exact diamond integral:", weight)
print("actual relative mesoscopic coefficient:", out)
