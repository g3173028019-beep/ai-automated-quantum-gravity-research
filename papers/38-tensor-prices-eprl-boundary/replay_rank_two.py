#!/usr/bin/env python3
"""Exact integer/Rational check of the finite marked EPRL rank formula."""
from fractions import Fraction as F

def D(i): return (2*i+1)*(i+1)
def channels(a,b): return set(range(abs(a-b),a+b+1,2))  # doubled spins
def rank(m,n,r,t):
    if not channels(2*m,2*n) & channels(2*r,2*t):
        return 0
    if not channels(m,n) & channels(r,t):
        return 0
    return len(channels(3*m,3*n) & channels(3*r,3*t))

def R(r,t,N):
    return sum(D(m)*D(n)*rank(m,n,r,t)
               for m in range(N+1)
               for n in range(max(0,m-r-t),min(N,m+r+t)+1))/F(D(r)*D(t))

def K(r,t):
    return F(3*r*t+r+t+1,D(r)*D(t))

for r in range(8):
    for t in range(8):
        m=2*(r+t)+2
        assert sum(rank(m,m+h,r,t) for h in range(-r-t,r+t+1)) == 3*r*t+r+t+1

assert [[R(r,t,1) for t in (0,1)] for r in (0,1)] == [[37,2],[2,F(145,36)]]
for N in (1,2,4,16,64,256,1024):
    v=R(0,0,N); b=R(0,1,N)/v; c=R(1,1,N)/v
    det=c-b*b
    assert det>0
    print(f'cutoff index {N:4}: normalized offdiag {float(b):.8f}, diag {float(c):.8f}, det {float(det):.8f}')
assert K(1,1)-K(0,1)**2==F(1,18)
print("The exact limiting matrix entries are K00=1, K01=1/3, K11=1/6.")

# Exact finite lower bounds for the five-state block from the two symmetry sectors.
# No numerical diagonalization, and no assumed equality between the crossed
# projector expectation and its lower bound.
for N in (1,2,3,4,8,64,256):
    V=sum(D(m)**2 for m in range(N+1))
    S1=sum(D(m)*D(m+1) for m in range(N))
    S2=sum(D(m)*D(m+2) for m in range(N-1))
    B=sum(( (D(m-1) if m else 0) +(D(m+1) if m<N else 0) )**2 for m in range(N+1))
    assert B==2*V-1-D(N)**2+2*S2
    c=F(S1,3*V); a=F(4*V-3+2*S2,36*V)
    b_min=(B-F(V-1,5))/F(36*V); b_max=F(B,36*V)
    assert a-b_max==F(1,18)+F(D(N)**2-2,36*V)>0
    assert 36*(a+3*b_min-4*c*c)>=(F(7,5)*(V-1)+D(N)**2)/V>0
print("The two five-state symmetry-sector lower bounds hold at every replayed cutoff.")
