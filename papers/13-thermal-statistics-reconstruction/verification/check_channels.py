"""Portable exact convention checks and projection-channel diagnostics."""
from pathlib import Path
import json,math,sympy as s
ROOT=Path(__file__).resolve().parent
u,p=s.symbols('u p')
logmgf=s.log(1+p*(s.exp(u)-1))
cumulants=[s.simplify(s.diff(logmgf,u,k).subs(u,0)) for k in range(1,5)]
expected=[p,p-p*p,p-3*p*p+2*p**3,p-7*p*p+12*p**3-6*p**4]
assert all(s.expand(a-b)==0 for a,b in zip(cumulants,expected))
K=s.ones(2)/2; f=s.diag(10,1);g=s.diag(1,10)
lhs=s.trace(f*K*g*K);rhs=s.trace(f*g*K)
assert lhs==s.Rational(121,4) and rhs==10 and K*K==K
d=s.symbols('d',positive=True)
v=s.Matrix([s.sqrt(d),s.sqrt(1-d)])
P=v*v.T
assert all(s.simplify(z)==0 for z in P*P-P)
assert s.simplify(s.trace(P)-1)==0
diagnostics=[]
for beta in [10,100,1000,10000]:
    den=0.;num=0.
    for n in range(1,4000):
        dn=math.exp(-n*math.log(2))
        # All terms have been scaled by exp(-2 beta) to avoid overflow.
        den+=dn*math.exp(-2*beta/n)+(1-dn)*math.exp(-2*beta*(n+1))
        block=dn*math.exp(-beta/n)+(1-dn)*math.exp(-beta*(n+1))
        num+=block*block
    diagnostics.append({'beta':beta,'projection_two_channel_ratio':num/den,'Phi2':1-num/den,'terms':3999})
# Bernoulli unique law but two independent GNS vectors.
gram=s.diag(1,s.Rational(1,4)) # basis 1 and x-1/2
assert gram.det()>0
out={'status':'PASS','sympy':s.__version__,'bernoulli_cumulants':[str(s.expand(c)) for c in cumulants],'unordered_clock_counterexample':{'lhs':str(lhs),'rhs':str(rhs)},'rank_one_projection_block':'symbolic idempotency and trace1 PASS','pure_projection_condensation_diagnostics':diagnostics,'bernoulli_GNS_Gram_determinant':str(gram.det()),'scope':'Diagnostics illustrate the channel limit; the analytic proof establishes it. No full point-process reconstruction is inferred from a finite sample.'}
(ROOT/'channels.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
print('PASS: cumulants1-4; unordered-clock inequality; pure projection block; GNS distinction')
