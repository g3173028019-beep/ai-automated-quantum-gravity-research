"""Portable checks for the energy-weighted JT recovery research note.
Run from packet root: python code/reproduce.py
No old-project input, network, or exact-JT simulation is used.
"""
from pathlib import Path
import datetime
import hashlib
import json
import math
import platform
import mpmath as mp
import numpy as np
import sympy as sp

BASE = Path(__file__).resolve().parents[1]
mp.mp.dps = 150
mp.iv.dps = 150
iv = mp.iv


def gaussian_closure():
    v,a,t=sp.symbols('v a t',real=True)
    mgf=sp.exp(a*a*v/2)
    cross=(mgf+2*sp.diff(mgf,a)+2*sp.diff(mgf,a,2)).subs(a,2)
    residual=sp.exp(8*v)-2*cross+1+8*v+12*v*v
    target=sp.exp(8*v)-2*sp.exp(2*v)*(1+6*v+8*v*v)+1+8*v+12*v*v
    assert sp.simplify(residual-target)==0
    F=mgf.subs(a,2)-2*mgf.subs(a,1)+1
    assert sp.limit(target/v**3,v,0)==sp.Rational(80,3)
    assert sp.limit(F/v,v,0)==1
    x,y=1-t*sp.tanh(t),sp.tanh(t)
    assert sp.simplify(sp.diff(x,t,2)+2*x/sp.cosh(t)**2)==0
    assert sp.simplify(sp.diff(y,t,2)+2*y/sp.cosh(t)**2)==0
    assert sp.simplify(x*sp.diff(y,t)-y*sp.diff(x,t)-1)==0
    return {'M_from_gaussian_mgf':'PASS','F_from_gaussian_mgf':'PASS',
            'M_leading_v_cubed':'80/3','F_leading_v':'1',
            'stability_ODE_and_Wronskian':'PASS'}


def interval(x):
    return {'lower': str(x.a), 'upper': str(x.b)}


def charge_polynomials(n):
    return ((n+1)*(n+2)**2*(n+3)/16,
            (n+1)**4/4,
            n*n*(n*n-1)/16)


def weighted(C, nu=1, N=None, T=1):
    if N is None:
        N = 0 if nu == 0 else int(mp.floor(4*mp.sqrt(2)*mp.sqrt(C)/3))
    c, v, energy = iv.mpf(C), 1/(8*iv.mpf(C)), iv.mpf(nu)
    M = iv.exp(8*v)-2*iv.exp(2*v)*(1+6*v+8*v*v)+1+8*v+12*v*v
    F = iv.exp(2*v)-2*iv.exp(v/2)+1
    grav = 2*c*iv.sqrt(M)
    k0 = iv.sqrt(1+(N+2)*energy)
    if N == 0:
        k0, charge = iv.mpf(1), iv.mpf(1)
    else:
        un, dn, ln = charge_polynomials(iv.mpf(N))
        charge = (iv.sqrt(iv.mpf(3)/4+(un-iv.mpf(3)/4)*energy/N)
                  +iv.sqrt(iv.mpf(1)/4+(dn-iv.mpf(1)/4)*energy/N)
                  +iv.sqrt(ln*energy/N))
    pieces = [2*T*grav, 2*T*k0*iv.sqrt(F), T*charge/(4*c),
              4*iv.sqrt(energy/(N+1))]
    total = sum(pieces)
    return {'C': C, 'nu': nu, 'N': N, 'T': T,
            'uncapped_bound': interval(total),
            'unhalved_trace_norm_bound_upper': str(min(mp.mpf(2), mp.mpf(total.b))),
            'pieces_gravity_redshift_charge_tail': [interval(x) for x in pieces],
            'C_quarter_scaled': interval(total*iv.sqrt(iv.sqrt(c)))}


def old(C, nu=1):
    N = int(mp.floor((8*mp.mpf(nu)*C)**(mp.mpf(1)/3)))
    c, v = mp.mpf(C), 1/(8*mp.mpf(C))
    M = mp.exp(8*v)-2*mp.exp(2*v)*(1+6*v+8*v*v)+1+8*v+12*v*v
    F = mp.exp(2*v)-2*mp.exp(v/2)+1
    return {'N': N, 'uncapped_bound': str(2*(2*c*mp.sqrt(M)+(N+1)*mp.sqrt(F)
                         +(N+1)*(N+2)/(8*c))+4*mp.sqrt(mp.mpf(nu)/(N+1)))}


def second_moment(C, mu=1):
    N = int(mp.floor(4*mp.sqrt(C)))
    c, m = iv.mpf(C), iv.mpf(mu)
    v, root = 1/(8*c), iv.sqrt(m)
    M = iv.exp(8*v)-2*iv.exp(2*v)*(1+6*v+8*v*v)+1+8*v+12*v*v
    F = iv.exp(2*v)-2*iv.exp(v/2)+1
    n = iv.mpf(N)
    charge = (iv.sqrt(iv.mpf(3)/4+7*root/4+(n*n/16+n/2+iv.mpf(23)/16)*m)
              +iv.sqrt(iv.mpf(1)/4+root+(n*n/4+n+iv.mpf(3)/2)*m)
              +n*root/4)
    total = 4*c*iv.sqrt(M)+2*(1+root)*iv.sqrt(F)+charge/(4*c)+4*root/(N+1)
    return {'C':C,'mu_second_excitation_moment':mu,'N':N,
            'bound':interval(total),'sqrt_C_scaled':interval(total*iv.sqrt(c))}


def algebra_checks():
    n = sp.symbols('n', nonnegative=True)
    un, dn, ln = charge_polynomials(n)
    # These quotients being increasing proves the three mean-weight inequalities.
    uq = sp.cancel((un-sp.Rational(3,4))/n)
    dq = sp.cancel((dn-sp.Rational(1,4))/n)
    lq = sp.cancel(ln/n)
    assert sp.expand(un-(n**4/16+n**3/2+23*n**2/16+7*n/4+sp.Rational(3,4)))==0
    assert sp.expand(dn-(n**4/4+n**3+3*n**2/2+n+sp.Rational(1,4)))==0
    assert all(a > 0 for a in sp.Poly(sp.diff(uq,n),n).all_coeffs())
    assert all(a > 0 for a in sp.Poly(sp.diff(dq,n),n).all_coeffs())
    assert sp.simplify(sp.diff(lq,n)-(3*n*n-1)/16) == 0
    # Last derivative is >=1/8 on n>=1, handling the lower shift.
    a = 4*sp.sqrt(2)/3
    coeff = sp.sqrt(a/2)+a**sp.Rational(3,2)/4+4/sp.sqrt(a)
    assert sp.simplify(1/sp.sqrt(2)+3*a/4-4/a)==0
    assert sp.simplify(coeff-7*sp.sqrt(2)*sp.sqrt(a)/3)==0
    return {'upper_shift_quotient': str(uq), 'diagonal_quotient': str(dq),
            'lower_shift_quotient': str(lq), 'optimized_a': str(a),
            'asymptotic_coefficient_exact': str(sp.simplify(coeff)),
            'asymptotic_coefficient_decimal': str(sp.N(coeff,40))}


def matrix_checks():
    rng = np.random.default_rng(20261028)
    worst_ratio, count = 0., 0
    for N in [1,2,3,8,31,100]:
        # N+4 basis is sufficient to retain the TWO upper-edge raising steps.
        dim = N+4
        K = np.zeros((dim,dim), complex)
        for n in range(dim-1):
            K[n+1,n] = np.sqrt((n+1)*(n+2))/(2j)
            K[n,n+1] = K[n+1,n].conjugate()
        squared = (K@K)[:, :N+1]
        for n in range(N+1):
            u,d,l = charge_polynomials(n)
            assert np.isclose(np.vdot(squared[:,n],squared[:,n]).real,u+d+l)
        for j in range(30):
            # Vector in matter tensor reference, with coherent tail budget.
            f = rng.normal(size=(N+1,3))+1j*rng.normal(size=(N+1,3))
            f[1:] *= 10**(-j/10)
            f /= np.linalg.norm(f)
            nu = float(np.sum(np.arange(N+1)[:,None]*np.abs(f)**2))
            mu = float(np.sum(np.arange(N+1)[:,None]**2*np.abs(f)**2))
            u,d,l = charge_polynomials(N)
            b = math.sqrt(.75+(u-.75)*nu/N)+math.sqrt(.25+(d-.25)*nu/N)+math.sqrt(l*nu/N)
            ratio = np.linalg.norm(squared@f)/b
            assert ratio <= 1+1e-12
            assert np.linalg.norm((np.arange(N+1)+1)[:,None]*f)**2 <= 1+(N+2)*nu+1e-10
            b2 = (math.sqrt(.75+7*math.sqrt(mu)/4+(N*N/16+N/2+23/16)*mu)
                  +math.sqrt(.25+math.sqrt(mu)+(N*N/4+N+1.5)*mu)+N*math.sqrt(mu)/4)
            assert np.linalg.norm(squared@f) <= b2*(1+1e-12)
            assert np.linalg.norm((np.arange(N+1)+1)[:,None]*f) <= 1+math.sqrt(mu)+1e-10
            worst_ratio = max(worst_ratio,ratio)
            count += 1
    return {'coherent_reference_vectors_checked': count, 'maximum_charge_ratio': worst_ratio,
            'operator_columns_checked_through_N':100, 'upper_edge_leakage_retained':True}


def hostile_controls():
    # A mean-energy constraint does NOT imply the moment required by raw Duhamel.
    # Finite coherent two-mode states already reject ||K0 psi||<=<K0>.
    rows=[]
    for N in [100,10000,1000000]:
        p=mp.mpf(1)/N
        rows.append({'N':N,'mean_n':1,'norm_K0_squared':str(1+(N+2)),
                     'mean_K0':2,'projection_error':str(mp.sqrt(4*p-3*p*p)),
                     'tail_probability':str(p)})
    return {'coherent_tail_controls':rows}


if __name__ == '__main__':
    algebra = algebra_checks()
    matrices = matrix_checks()
    examples=[]
    for C in [10**6,10**12,10**18,10**24]:
        row=weighted(C)
        row['inherited_C_minus_one_sixth_bound']=old(C)
        examples.append(row)
    report={'generated_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
            'python':platform.python_version(),
            'package_versions':{'mpmath':mp.__version__,'numpy':np.__version__,'sympy':sp.__version__},
            'gaussian_closure':gaussian_closure(),
            'algebra':algebra,'finite_matrix_checks':matrices,'examples':examples,
            'second_moment_new_input_class':[second_moment(C) for C in [10**6,10**12,10**18]],
            'zero_energy_example':weighted(10**12,0),'hostile_controls':hostile_controls(),
            'scope':'analytic bounds on exact full-space JT channel; finite matrices check charge algebra only',
            'status':'PASS', 'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (BASE/'results').mkdir(exist_ok=True)
    (BASE/'results/reproduction.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps({'status':report['status'],'algebra':algebra,'matrix_checks':matrices,
                      'examples':[{'C':r['C'],'N':r['N'],'new':r['unhalved_trace_norm_bound_upper'],
                                   'old':r['inherited_C_minus_one_sixth_bound']['uncapped_bound']} for r in examples]},indent=2))
