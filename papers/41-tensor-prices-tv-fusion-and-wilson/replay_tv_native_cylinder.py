"""Exact integer replay of the compact doubled SU(2)_k cylinder tests."""
from math import comb
from fractions import Fraction

def walk(k,n):
    v=[1]+[0]*k
    for _ in range(n): v=[(v[i-1] if i else 0)+(v[i+1] if i<k else 0) for i in range(k+1)]
    return v

def catalan(n):return comb(2*n,n)//(n+1)

def two_arm_seams(k,n):
    v=walk(k,n)
    d=walk(k,2*n)[0]
    assert sum(x*x for x in v)==d
    p={(i,j):Fraction((u*w)**2,d*d) for i,u in enumerate(v) for j,w in enumerate(v) if u*w}
    assert sum(p.values())==1
    return p

for k in range(1,41):
    for n in range(1,k+1):assert walk(k,2*n)[0]==catalan(n)
    assert walk(k,2*k+2)[0]==catalan(k+1)-1
    for parity in (0,1):
        n=2*k+parity
        assert set(two_arm_seams(k,n))=={(i,j) for i in range(k+1) for j in range(k+1) if i%2==j%2==parity}
for n in range(1,20):
    assert walk(2,2*n)[0]==2**(n-1)
    d={(1,1):1} if n%2 else {(i,j):Fraction(1,4) for i in (0,2) for j in (0,2)}
    assert two_arm_seams(2,n)==d
assert two_arm_seams(2,2)[(0,2)]+two_arm_seams(2,2)[(2,0)]==Fraction(1,2)
print("PASS: exact finite-level wall, parity-matched closure, and level-2 mixed charge laws")

# k=2, doubled four-fundamental cut. Chiral F=(1/sqrt(2))*Hadamard,
# so all entries of the actual doubled recoupling are rational.
h=((1,1),(1,-1))
ix=[(0,0),(0,1),(1,0),(1,1)]
G=[[Fraction(h[a][c]*h[b][d],2) for c,d in ix] for a,b in ix]
assert all(sum(G[i][a]*G[j][a] for a in range(4))==(i==j) for i in range(4) for j in range(4))
p=lambda i:[[int(a==b==i) for b in range(4)] for a in range(4)]
q=lambda j:[[G[a][j]*G[b][j] for b in range(4)] for a in range(4)]
mul=lambda A,B:[[sum(A[i][k]*B[k][j] for k in range(4)) for j in range(4)] for i in range(4)]
X=p(0)
for Y in (q(0),p(2),q(2),p(0)): X=mul(X,Y)
assert X==[[v*Fraction(-1,16) for v in row] for row in p(0)]
print("PASS: actual signed doubled level-2 cut pencil has unique null root -1/16")
