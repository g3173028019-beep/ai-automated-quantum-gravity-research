#!/usr/bin/env python3
"""Exact checks of finite arithmetic and formal normalization constants."""
from fractions import Fraction as Q

def mul(a,b):
    return tuple(tuple(sum(x*y for x,y in zip(row,col)) for col in zip(*b)) for row in a)

def tr(a): return tuple(zip(*a))
def eye(d): return tuple(tuple(int(i == j) for j in range(d)) for i in range(d))
def mod(a,p): return tuple(tuple(x%p for x in row) for row in a)
def powm(a,k,p=None):
    n=eye(len(a))
    for _ in range(k):
        n=mul(n,a)
        if p:n=mod(n,p)
    return n
def diag(a): return tuple(tuple(x if i==j else 0 for j in range(len(a))) for i,x in enumerate(a))

F=diag((-7,1,1,1,1))
b=((8,3,0,0),(21,8,0,0),(0,0,1,0),(0,0,0,1))
r=diag((1,1,1,-1))
f=diag((-7,1,1,1))
assert mul(mul(tr(b),f),b)==f
assert 8*8-3*21==1
assert mul(b,r)==mul(r,b)
g=mul(powm(b,10),r)
assert mod(g,3)==mod(r,3)
assert mod(g,79)==((0,60,0,0),(25,0,0,0),(0,0,1,0),(0,0,0,78))
assert powm(g,2,79)!=eye(4) and powm(g,4,79)==eye(4)
assert sum(F[i][i]*x*x for i,x in enumerate((1,2,1,1,1)))==0
assert all((a+b+c)%8!=7 for a in (0,1,4) for b in (0,1,4) for c in (0,1,4))
assert all(a+b+c != 4 for a in (0,1) for b in (0,1) for c in (0,1))

# Laurent-polynomial operations, variable r; each exponent and coefficient exact.
def scale(k,f,shift=0): return {n+shift:Q(k)*a for n,a in f.items() if a*k}
def dif(f): return {n-1:n*a for n,a in f.items() if n*a}
V={2:1,-1:-1}
assert scale(Q(-1,2),dif(dif(V)))=={0:-1,-3:1}
assert scale(-1,V,-2)=={0:-1,-3:1}
assert scale(Q(-1,2),dif(V),-1)=={0:-1,-3:Q(-1,2)}
beta_over_pi=Q(4,3)
assert 4*(1+1+4*Q(1,2)**2)==12
assert Q(12,3)==4  # integral of 12 r^{-4} from 1 to R
assert Q(4,24)*beta_over_pi**3==Q(32,81)
assert Q(3,8)*Q(32,81)==Q(4,27)
assert Q(3,8)*Q(4,24)*beta_over_pi==Q(1,12)
# Coefficients of S^2, S, 1 in -4(S+6)+(S+12)^2/6.
assert (Q(1,6),-4+Q(24,6),-24+Q(144,6))==(Q(1,6),0,0)
print('PASS: congruences, isotropy, curvature, CGB polynomial, volume and action coefficients.')
