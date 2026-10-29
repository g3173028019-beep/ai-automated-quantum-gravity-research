"""Exact finite checks and heat-trace diagnostics; run from any directory."""
from pathlib import Path
import json, math, sympy as s

ROOT=Path(__file__).resolve().parent
omega=s.symbols('Omega',positive=True)
rows=[]
for j in [s.Rational(1,2),s.Integer(1),s.Rational(3,2),s.Integer(2)]:
    d=int(2*j+1); eye=s.eye(d); jp=s.zeros(d)
    for a in range(d-1):
        m=j-a-1
        jp[a,a+1]=s.sqrt((j-m)*(j+m+1))
    js=[(jp+jp.T)/2,(jp-jp.T)/(2*s.I),s.diag(*[j-a for a in range(d)])]
    L=[s.kronecker_product(J,eye) for J in js]
    R=[s.kronecker_product(eye,J.T) for J in js]
    A=[L[a]-R[a] for a in range(3)]
    S=[L[a]+R[a] for a in range(3)]
    H=sum((A[a]**2+omega**2*S[a]**2 for a in range(2)),s.zeros(d*d))/j
    P=s.kronecker_product(eye,s.diag(*[s.exp(s.I*s.pi*(j-a)).expand(complex=True) for a in range(d)]))
    zero=lambda mat:all(s.simplify(v)==0 for v in mat)
    energy=s.diag(*[2*(2*j*(j+1)-(j-a)**2-(j-b)**2)/j for a in range(d) for b in range(d)])
    assert zero(H.subs(omega,1)-energy)
    assert zero(P*H*P.conjugate().T-omega**2*H.subs(omega,1/omega))
    assert zero(P**2-(-1)**(2*j)*s.eye(d*d))
    QY=sum((A[a]**2 for a in range(3)),s.zeros(d*d))+omega**2*A[2]**2+omega**2*sum((S[a]**2 for a in range(2)),s.zeros(d*d))
    assert zero(QY-j*H-(1+omega**2)*A[2]**2)
    for idx in [d-1,(d-1)*d]:
        e=s.eye(d*d)[:,idx]
        assert zero(H*e-2*(1+omega**2)*e)
    rows.append({'j':str(j),'matrix_dimension':d*d,'exact_spectrum':'PASS','symbolic_covariance':'PASS','projective_square':'PASS','cross_pole_states':'PASS','prior_operator_difference':'PASS'})

thermal=[]; t=1.;theta=1.
limit=4*math.exp(-4*t/theta)/(1-math.exp(-4*t/theta))**2
for j in [5,20,100,500]:
    slot=sum(math.exp(-2*t*(j*(j+1)-(j-n)**2)/(theta*j)) for n in range(2*j+1))
    thermal.append({'j':j,'heat_trace':slot**2,'analytic_limit':limit})

result={'status':'PASS','sympy':s.__version__,'normalization':'theta=1 for finite checks; general theta is an overall factor','finite_checks':rows,'heat_trace_diagnostics':thermal,'scope':'Finite checks corroborate conventions. Heat-trace convergence and Weyl-symbol transform are analytic proofs, not inferred from decimals.'}
(ROOT/'checks.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
print('PASS: 4 exact spin cases; symbolic Omega; axial comparison; heat diagnostics')
