"""Portable finite checks for resolved instruments and spectral references.

No historical imports or network. Exact identities use SymPy; the endpoint
covariance example uses NumPy and an explicit numerical tolerance.
"""
from pathlib import Path
from itertools import combinations, combinations_with_replacement, product
import json
import sys
import numpy as np
import sympy as S

OUT = Path(__file__).resolve().parents[1] / "data"

def operational_identities():
    chi=S.symbols('chi',positive=True)
    g=S.diag(S.exp(chi/2),S.exp(-chi/2)); k=S.exp(-chi/2)*g
    rejected=S.simplify(S.eye(2)-k.H*k)
    assert rejected==S.diag(0,1-S.exp(-2*chi))
    # Opposite boosts, locally maximal filters, coarse identity.
    invk=S.exp(-chi/2)*g.inv()
    assert S.simplify(invk*k)==S.exp(-chi)*S.eye(2)
    a,b=S.symbols('a b',real=True)
    assert S.expand(b*b*(1-a*a)+a*a*(1-b*b)+(1-a*a)*(1-b*b)-(1-a*a*b*b))==0
    one=S.eye(2); k0=one/S.sqrt(2); k1=one/S.sqrt(2)
    rho=S.diag(1,0)
    phase_plus=(k0+k1)/S.sqrt(2); phase_minus=(k0-k1)/S.sqrt(2)
    assert S.trace(phase_plus*rho*phase_plus.H)==1
    assert S.trace(phase_minus*rho*phase_minus.H)==0
    return dict(boost_complement=str(rejected),opposite_boost_product='exp(-chi) I',coarse_maximal_filter='I',jump_history_norm_identity=True,identical_diagonal_instrument_recombiner_probability=[1,0])

def centering():
    h1=S.diag(S.Rational(1,5),S.Rational(11,10)); h2=S.diag(S.Rational(7,10),S.Rational(3,2))
    center=lambda h:h-min(h.eigenvals())*S.eye(h.rows)
    union=S.diag(h1,h2); separate=S.diag(center(h1),center(h2)); together=center(union)
    assert together-separate==S.diag(0,0,S.Rational(1,2),S.Rational(1,2))
    return dict(separate_diagonal=[str(t) for t in separate.diagonal()],global_diagonal=[str(t) for t in together.diagonal()],spectator_shift='1/2',scalar_shift_blind=True,
                stronger_no_go='Unitary covariance + all finite direct sums + scalar shift invariance implies F_n(H)=a I_n; analytic proof in manuscript Section 7')

def endpoint_covariance():
    dagger=lambda m:m.conj().T
    def root(m,power):
        ev,u=np.linalg.eigh(m)
        assert min(ev)>0
        return (u*(ev**power))@dagger(u)
    def instrument(g,ux,uy):
        inv=np.linalg.inv(g)
        a=root(ux,-0.5)@inv@uy@dagger(inv)@root(ux,-0.5)
        c=float(np.linalg.eigvalsh(a)[-1])**(-0.5)
        k=c*dagger(inv)
        return c,dagger(k)@uy@k,ux-dagger(k)@uy@k
    ux=np.array([[7/5,1/7+1j/11],[1/7-1j/11,6/5]])
    uy=np.array([[9/7,-1/9+1j/13],[-1/9-1j/13,8/7]])
    g=np.array([[1,1j/3],[1/5,1+1j/15]])
    sx=np.array([[1,1/4],[0,1]],dtype=complex)
    sy=np.array([[5/4,0],[1j/5,4/5]],dtype=complex)
    c,e,n=instrument(g,ux,uy)
    cp,ep,np_=instrument(sy@g@np.linalg.inv(sx),sx@ux@dagger(sx),sy@uy@dagger(sy))
    state=np.array([[3/5,1/17],[1/17,2/5]],dtype=complex)
    state=state/np.trace(state@ux)
    invsx=np.linalg.inv(sx)
    statep=dagger(invsx)@state@invsx
    residuals=dict(scale=abs(c-cp),
                   selected=float(np.linalg.norm(ep-sx@e@dagger(sx))),
                   rejected=float(np.linalg.norm(np_-sx@n@dagger(sx))),
                   probability=float(abs(np.trace(state@e)-np.trace(statep@ep))))
    assert max(residuals.values())<1e-12
    assert np.linalg.eigvalsh(n)[0]>-1e-12
    return dict(residuals=residuals,rejected_minimum_eigenvalue=float(np.linalg.eigvalsh(n)[0]))


if __name__=='__main__':
    result=dict(python=sys.version.split()[0],sympy=S.__version__,numpy=np.__version__,
                operational=operational_identities(),centering=centering(),endpoint_covariance=endpoint_covariance())
    (OUT/'results.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result,indent=2))
