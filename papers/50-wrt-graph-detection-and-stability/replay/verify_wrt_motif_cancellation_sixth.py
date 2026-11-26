"""Exact graph counts and integrated cube polynomials through degree six."""
from fractions import Fraction as Q
from math import comb, factorial
from verify_wrt_motif_diamond_pair import H as DH, G as DG, classify

CH = [(0,3),(0,4),(0,6),(1,4),(1,10),(1,11),(2,8),(2,9),(2,11),
      (3,5),(3,10),(4,7),(5,9),(5,11),(6,7),(6,8),(7,9),(8,10)]
CG = [(0,3),(0,5),(0,10),(1,3),(1,4),(1,5),(2,6),(2,8),(2,10),
      (3,7),(4,9),(4,10),(5,11),(6,7),(6,11),(7,8),(8,9),(9,11)]

def cube(edges, n=12, qmax=6):
    nbh = [0]*n
    for u,v in edges:
        nbh[u] |= 1<<v
        nbh[v] |= 1<<u
    F = [[0]*(qmax+1) for _ in range(1<<n)]
    F[0][0] = 1
    coef = {(m,q,d):[(-1)**l*comb(d,l)*factorial(m+q-1)//factorial(m+q-1-l)
                       for l in range(min(d,q)+1)]
            for m in range(1,n+1) for q in range(qmax+1) for d in range(4)}
    for S in range(1,1<<n):
        m, bits = S.bit_count(), S
        while bits:
            bit = bits & -bits
            bits ^= bit
            u = bit.bit_length()-1
            Sm = S ^ bit
            d = (nbh[u] & Sm).bit_count()
            for q in range(qmax+1):
                F[S][q] += sum(a*F[Sm][q-l] for l,a in enumerate(coef[m,q,d]))
    return tuple(Q(F[-1][q],factorial(n+q)) for q in range(qmax+1))

def delta(H,G):
    return tuple(u-v for u,v in zip(H,G))

assert classify(DH) == ((4,2,4),1)
assert classify(DG) == ((4,2,4),0)
assert classify(CH) == ((0,2,5),0)
assert classify(CG) == ((0,2,6),0)
cl = [cube(G) for G in (DH,DG,CH,CG)]
assert all(z[1] == -12 for z in cl)
dd, dc = delta(*cl[:2]), delta(*cl[2:])
assert dd == (0,)*5+(Q(1,467775),Q(-3911,212837625))
assert dc == (0,)*5+(Q(-1,93555),Q(53,500500))
assert 5*dd[5]+dc[5] == 0
six = 5*dd[6]+dc[6]
assert six == Q(11933,851350500)
relative = six*Q(factorial(78),factorial(72))
assert relative == Q(2449558508,945)
print("diamond and pentagon contrasts:", dd[5:], dc[5:])
print("composite order-five exactly cancels; sixth cube weight:", six)
print("sixth relative hard coefficient:", relative)
