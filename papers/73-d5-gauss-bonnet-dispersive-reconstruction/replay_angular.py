"""Exact rational checks for the D5 K10 angular certificate.

This combines the archival rational lens polynomial, Bernstein split and
certificate assertions. It deliberately does not include finite-spin numerical
quadrature. Run with stock Python 3 and without the -O option.
"""
from fractions import Fraction as F
from functools import lru_cache
from math import comb, factorial, prod

@lru_cache(None)
def endpoint_poly(a,b):
    return [sum(F((-1)**i*comb(a,i)*comb(b,j-i),2**j)
                for i in range(max(0,j-b),min(a,j)+1)) for j in range(a+b+1)]

def exact_profile(K):
    # pi times unnormalized overlap is 2*pi times the polynomial below.
    conv=[F(0)]*(4*K+4)
    for l in range(K+1):
        for m in range(K+l+2):
            i=K-l+m; j=K+l+1
            beta=F(factorial(i)*factorial(j),factorial(i+j+1))
            amp=F(comb(K,l)*2**(K-l)*comb(K+l+1,m),K+l+1)*beta
            pol=endpoint_poly(2*K+2+m,K+l+1-m)
            for k,c in enumerate(pol): conv[K-l+k]+=amp*c
    norm=F(4)*sum(F((-1)**i*comb(K,i),2*i+3) for i in range(K+1))
    conv=[2*c/(norm*norm) for c in conv]  # g1*pi
    B=[F(1)]
    for k in range(1,K+1): B.append(-B[-1]/(2*k*(2*K+2*k+3)))
    B2=[sum(B[i]*B[k-i] for i in range(k+1)) for k in range(K+1)]
    P=[F(1)]
    for k in range(1,K+1): P.append(-sum(B2[i]*P[k-i] for i in range(1,k+1)))
    g=[F(0)]*len(conv)
    for k,p in enumerate(P):
        for d,c in enumerate(conv):
            dd=d
            for _ in range(k):
                c*=dd*(dd+1);dd-=2
            if c:
                assert dd>=0
                g[dd]+=(-1)**k*p*c
    for j in range(K+1):
        assert 4*sum(c*F(2**(d+2*j+3),d+2*j+3) for d,c in enumerate(g))==(j==0)
    return g

def split_half(b):
    u=list(b); l=[u[0]]; r=[u[-1]]
    while len(u)>1:
        u=[(a+b)/2 for a,b in zip(u,u[1:])]
        l.append(u[0]);r.append(u[-1])
    return l,r[::-1]

@lru_cache(None)
def J(d,a):
 g=exact_profile(10)
 for _ in range(d):g=[i*g[i] for i in range(1,len(g))]
 pol=[F(0)]*a+g;n=len(pol)-1
 b=[sum(pol[j]*2**j*F(comb(i,j),comb(n,j)) for j in range(i+1)) for i in range(n+1)]
 boxes=[b]
 for z in range(4):boxes=[u for w in boxes for u in split_half(w)]
 B=F(2,16*(n+1))*sum(abs(c) for box in boxes for c in box)
 p=exact_profile(10)
 chi=[F(1),F(32),F(1024),F(25000)]
 mx=lambda s:max(F(1,4)**s,F(1,2)**s)
 out=mx(a-d-2)*B
 for l in range(2,d):
  pL=p[:]
  for _ in range(l):pL=[i*pL[i] for i in range(1,len(pL))]
  val=abs(sum(c*2**i for i,c in enumerate(pL)))
  h=d-1-l
  out+=val*F(2)**(a-h)*sum(comb(h,m)*prod(range(3+l,3+l+h-m))*mx(a-d-2+m)*chi[m] for m in range(h+1))
 return out

def B_terms(z,N):
    q=[F(1)]
    for j in range(1,N+1):q.append(-q[-1]*z*z/(2*j*(23+2*j)))
    return q

def beta_mass(a,b):
    return 140*sum(F((-1)**j*comb(3,j),4+j)*(b**(4+j)-a**(4+j)) for j in range(4))

if __name__=='__main__':
    for z,lo in [(F(11),F(1,16)),(F(13),F(1,60)),(F(141,10),F(1,160))]:
        b=B_terms(z,61)
        assert sum(b)>lo and z*z<2*62*(23+2*62)
    coeff=B_terms(F(1),10)
    D=[sum(coeff[j]*coeff[k-j] for j in range(k+1)) for k in range(11)]
    P=[F(1)]
    for j in range(1,11):P.append(-sum(D[i]*P[j-i] for i in range(1,j+1)))
    pol=lambda z:sum(p*z**(2*j) for j,p in enumerate(P))
    m1=beta_mass(F(2,13),F(3,13)); m2=beta_mass(F(1,13),F(2,13))
    assert m1*20**4/16**2>20
    assert m2*32**4*pol(F(17,2))/60**2>20
    assert m2*40**4*pol(F(85,8))/160**2>20
    C=prod(range(1,24,2))**2*P[-1]; assert C>1600
    pdf=F(140)*F(64,13)*F(12,13)**3*F(1,13)**3
    assert pdf>F(6,25)
    assert 1-F(132)/(45*F(11,32))**2>F(7,16)>F(13,20)**2
    assert (F(7,16)-F(20,13)/(45*F(3,32)))/2>F(1,32)
    assert 1600*F(16,7)**4*F(6,25)*F(3,32)/32>20
    assert 140*F(64,13)*F(1,64)*2<32
    assert 140*F(64,13)**2/4<1024  # L1 of unnormalized beta second derivative <1/4
    assert 16*sum([F(1),F(0)])<32
    U=[[500,16,2,2,2,2],
       [11000,700,75,14,14,14],
       [400000,26000,4000,1000,300,300],
       [None,1500000,310000,100000,40000,18000],
       [None,None,46000000,21000000,9500000,4600000],
       [None,None,None,4100000000,2000000000,1000000000]]
    for d,row in enumerate(U):
        for a,upper in enumerate(row):
            if upper is not None: assert J(d,a)<upper, (d,a,float(J(d,a)))
    # h = A p(rho)-p derivative estimate; J(j,a) can be bounded by lower-a entries.
    def h(d,a):
        return F(U[d+1][a+1],4)+100*sum(U[j][a] for j in range(d+1))
    flat_rad=h(4,2)+4*h(3,1)
    cross=(4*(h(3,2)+2*h(2,1)+2*h(1,0))+2*(h(2,2)+2*h(1,1))
             +36*(h(2,2)+2*h(1,1))+4*h(1,2)+h(0,2))
    diag=[U[j][j]+1 for j in range(5)]
    anis=sum(a*b for a,b in zip([150,2000,1000,200,16],diag))
    assert 4*(flat_rad+cross+anis)<2**36
    assert 4*(h(0,2)+16*(U[0][2]+1))<2**20
    cap=F(9,400000)
    for d in range(5): assert cap**2*h(d,d)<1
    assert 2**20*cap**2<F(1,256)
    assert F(101,100)*2**35*cap**2<19
    assert (1-(cap/21)**2)**(-2)<F(101,100)
    assert 20/(1+5*cap/21)**4>19
    # General Cartesian curvature correction bounds (include volume 4*pi cancels p=pi*g).
    curve1=4*4*20000*25600*sum(diag)*cap**2
    curve2=4*4*20000*1600*sum(h(d,d) for d in range(5))*cap**2
    assert 4*(flat_rad+cross+anis)+curve1+curve2<2**35
    c5=F(2**15*8**5*21*5**4)*cap**3
    tailratio=8*cap*F(8,6)*F(6,5)**4
    assert 7+3240*cap+(4**4*24*2048)*cap**2+c5/(1-tailratio)<16
    print('K10 exact spectral, weighted derivative, and overlap budgets replayed.')
    g=exact_profile(10)
    assert g[1]<0
    assert sum(c*2**j for j,c in enumerate(g))==0
    assert sum(j*c*2**(j-1) for j,c in enumerate(g) if j)==0
    pole=4*sum(p*F(2**(j+1),j+1) for j,p in enumerate(g))
    assert pole==F(5372874714177540760015498375447,59863441286320708007812500000)
    assert pole<90 and (16+F(64,29)**2)/2<11
    print('Exact raw Newton pole cost and scale-moment bound replayed.')
