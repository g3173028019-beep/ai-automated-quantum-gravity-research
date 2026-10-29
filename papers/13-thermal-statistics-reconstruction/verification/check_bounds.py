"""Bounded-statistic moment rescue and explicit assumption counterexamples."""
from pathlib import Path
import json, math
import sympy as sp
from scipy.integrate import quad
from scipy.special import ndtr

ROOT=Path(__file__).resolve().parent
n=sp.symbols("n",integer=True,nonnegative=True)
signed_moment=sp.exp(2*n*n-sp.pi**2/2)*sp.sin(2*sp.pi*n)
assert sp.simplify(signed_moment)==0

# Two deterministic positive, same-dimension spectra: equal one-clock data,
# unequal two-clock values. Exact rational Boltzmann weights suffice.
a=[sp.Rational(1,4),sp.Rational(3,4)]
b=[sp.Rational(1,3),sp.Rational(2,3)]
assert sum(a)==sum(b)==1
assert sum(x*x for x in a)==sp.Rational(5,8)
assert sum(x*x for x in b)==sp.Rational(5,9)

# Lambda is a variance-four lognormal Stieltjes class; conditional on Lambda
# use a Poisson process of intensity Lambda dE on [0,1], and append one fixed
# E=0 ground state. This is a finite, positive spectrum in every realization.
def laplace(s,eps,tol):
    c,cerr=quad(lambda E:1-math.exp(-s*math.exp(-E)),0,1,epsabs=tol,epsrel=tol)
    def integrand(x):
        density=math.exp(-x*x/8)/(2*math.sqrt(2*math.pi))
        return math.exp(-c*math.exp(x))*density*(1+eps*math.sin(math.pi*x/2))
    val,err=quad(integrand,-20,20,epsabs=tol,epsrel=tol,limit=200)
    # This error report is a quadrature diagnostic, not an interval certificate.
    return math.exp(-s)*val,math.exp(-s)*(err+cerr+3*ndtr(-10))

rows=[]
for s in [.5,1.,2.]:
    plus,e1=laplace(s,.5,2e-12);minus,e2=laplace(s,-.5,2e-12)
    plus2,e3=laplace(s,.5,2e-13);minus2,e4=laplace(s,-.5,2e-13)
    assert max(abs(plus-plus2),abs(minus-minus2))<1e-10
    rows.append({"s":s,"laplace_plus":plus2,"laplace_minus":minus2,
                 "difference":plus2-minus2,"reported_quadrature_error_sum":e3+e4})
assert any(abs(x["difference"])>100*x["reported_quadrature_error_sum"] for x in rows)

# A finite independent Bernoulli DPP checks the positive generating-function
# bound for two positive and one negative-but-uniformly-bounded energy.
energies=[-.3,.2,1.4];probabilities=[.2,.6,.8]
weights=[math.exp(-e) for e in energies];cap=max(weights)
mean=sum(p*w for p,w in zip(probabilities,weights))
mgf_checks=[]
for t in [.1,.5,1.,2.]:
    actual=math.prod(1-p+p*math.exp(t*w) for p,w in zip(probabilities,weights))
    bound=math.exp(mean/cap*math.expm1(t*cap))
    assert actual<=bound
    mgf_checks.append({"t":t,"actual":actual,"bound":bound})

report={
 "status":"PASS",
 "lognormal_signed_integer_moments":"zero for every nonnegative integer n",
 "same_clock_spectra":{"weights_A":[str(x) for x in a],"weights_B":[str(x) for x in b],
                       "Z_at_beta1":1,"Z_A_at_beta2":"5/8","Z_B_at_beta2":"5/9"},
 "positive_spectrum_Cox_counterexample":rows,
 "common_heat_mean_including_fixed_ground":1+math.exp(2)*(1-math.exp(-1)),
 "repulsion_failure_ratio_rho2_over_rho1_squared":math.exp(4),
 "bounded_DPP_MGF_examples":mgf_checks,
 "scope":"Exact integer-moment identity and finite examples; analytic general proofs in the accompanying manuscript. No black-hole dynamics or Airy density proof."
}
(ROOT/"checks.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8",newline="\n")
print(json.dumps(report,indent=2))
