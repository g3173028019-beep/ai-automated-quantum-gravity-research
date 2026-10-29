"""Exact revision checks. Python >=3.10; standard library only.
Run: python verification/check_revision.py
No historical project, network, symbolic package or Monte Carlo required.
"""
from pathlib import Path
from fractions import Fraction
from math import comb,factorial
import json
rows=[]
for n in range(1,13):
    norms=[]
    for k in range(n+1):
        m=n-k
        conv=sum(comb(2*a,a)*comb(2*(k-a),k-a) for a in range(k+1))
        assert conv==4**k
        prior=Fraction(factorial(2*n)*comb(2*m,m)*conv,4**n)
        profile=Fraction(factorial(2*n)*comb(2*m,m),4**m)
        assert prior==profile
        norms.append(profile)
        rows.append({'n':n,'k':k,'convolution':conv,'sector_norm':str(profile)})
    assert sum(norms)==(2*n+1)*norms[0]
    if n==1:assert sum(norms[k] for k in range(1,n))==0
assert len(rows)==90

# Direct finite Gaussian expectation of the full Pfaffian polynomial, independent
# of the displayed principal-minor sum. Every Gaussian monomial is squarefree.
from collections import defaultdict
from itertools import combinations

def matchings(vertices):
    if not vertices:
        yield (),1
        return
    a=vertices[0]
    for j in range(1,len(vertices)):
        b=vertices[j]
        rest=vertices[1:j]+vertices[j+1:]
        for tail,sign in matchings(rest):
            yield ((a,b),)+tail,(-1)**(j+1)*sign

def fixed_pfaffian(M,vertices):
    ans=0
    for pairs,sgn in matchings(tuple(vertices)):
        term=sgn
        for a,b in pairs:term*=M[a][b]
        ans+=term
    return ans

def double_factorial_odd(r):
    ans=1
    for a in range(1,2*r,2):ans*=a
    return ans

convolution=[]
for n in [1,2,3]:
    N=2*n
    M=[[0]*N for _ in range(N)]
    for a,b in combinations(range(N),2):
        M[a][b]=((a+3)*(b+5)+a*b*b)%19-9
        M[b][a]=-M[a][b]
    coeff=defaultdict(int)
    for pairs,sgn in matchings(tuple(range(N))):
        for mask in range(1<<n):
            normal=[];value=sgn
            for j,(a,b) in enumerate(pairs):
                if mask>>j&1:normal.append((a,b))
                else:value*=M[a][b]
            coeff[tuple(sorted(normal))]+=value
    direct=[0]*(n+1)
    for monomial,value in coeff.items():direct[len(monomial)]+=value*value
    rhs=[]
    for r in range(n+1):
        rhs.append(double_factorial_odd(r)*sum(fixed_pfaffian(M,S)**2 for S in combinations(range(N),N-2*r)))
    assert direct==rhs
    convolution.append({'n':n,'coefficients_of_c_2r':direct,'passed':True})
out={'passed':True,'arithmetic':'exact integer/Fraction','prior_variance_reduction_cases':len(rows),'n_range':[1,12],'fourth_moment_sum_cases':12,'empty_intermediate_sector_n1':True,'cases':rows,'direct_pfaffian_convolution':convolution,'limits':['Finite arithmetic checks corroborate formulas; analytic proofs establish all-size statements.','Compact-density and operator-scope corrections are analytic, not numerically proved.','Unchanged historical Monte Carlo and Wolfram scripts were not rerun.']}
p=Path(__file__).resolve().with_name('revision-checks.json');p.write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
print(json.dumps({k:out[k] for k in ['passed','prior_variance_reduction_cases','fourth_moment_sum_cases','empty_intermediate_sector_n1']}))
