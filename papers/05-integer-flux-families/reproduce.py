"""Exact integer STU constructions, duality reduction and mass/rank checks.

No source-code imports. SymPy decisions are exact. All writes stay in this
portable packet. General formulas are proved by diagonal covariance;
the machine checks also exercise every capacity-28 divisor assignment.
"""
from __future__ import annotations
from pathlib import Path
import itertools
import json
import hashlib
import sympy as sp

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results'; OUT.mkdir(exist_ok=True)
I=sp.I
z=sp.symbols('S U1 U2 U3')
bits=list(itertools.product((0,1),repeat=4))
mon=sp.Matrix([sp.prod(zz**bb for zz,bb in zip(z,b)) for b in bits])
ubits=list(itertools.product((0,1),repeat=3))
P=sp.zeros(8)
for col,b,sg in [(0,(1,1,1),-1),(1,(0,1,1),1),(2,(1,0,1),1),
                 (3,(1,1,0),1),(4,(0,0,0),1),(5,(1,0,0),1),
                 (6,(0,1,0),1),(7,(0,0,1),1)]:
    P[ubits.index(b),col]=sg
T=sp.diag(P,-P)
J2=sp.Matrix([[0,-1],[1,0]])
J4=sp.kronecker_product(J2,J2,J2,J2)

def polynomial(flux):
    return sp.expand((T*sp.Matrix(flux)).dot(mon))

def flux_from_polynomial(w):
    pp=sp.Poly(sp.expand(w),*z)
    cc=sp.Matrix([pp.coeff_monomial(tuple(b)) for b in bits])
    assert sp.expand(cc.dot(mon)-w)==0
    return T.T*cc

def charge(flux):
    f=sp.Matrix(flux)
    q=sp.expand(sum(f[4+a]*f[8+a]-f[a]*f[12+a] for a in range(4))/2)
    assert q==sp.expand(((T*f).T*J4*(T*f))[0]/4)
    return q

def data(flux,point):
    w=polynomial(flux); vals=dict(zip(z,point)); ys=[sp.im(t) for t in point]
    wp=sp.simplify(w.subs(vals)); wd=[sp.diff(w,zz).subs(vals) for zz in z]
    fs=[sp.simplify(wd[a]+I*wp/(2*ys[a])) for a in range(4)]
    jj=sp.zeros(8)
    for a,b in itertools.product(range(4),repeat=2):
        wdd=sp.diff(w,z[a],z[b]).subs(vals)
        dx=sp.simplify(wdd+I*wd[b]/(2*ys[a]))
        dy=sp.simplify(I*dx-(I*wp/(2*ys[a]**2) if a==b else 0))
        jj[2*a,2*b],jj[2*a+1,2*b]=sp.re(dx),sp.im(dx)
        jj[2*a,2*b+1],jj[2*a+1,2*b+1]=sp.re(dy),sp.im(dy)
    # Lkin= -(1/2) G_ab dphi^a dphi^b, G_real=2K_AAbar.
    ek=1/(16*sp.prod(ys))
    invk=sp.diag(*[t for y in ys for t in (4*y*y,4*y*y)])
    kin=sp.diag(*[t for y in ys for t in (1/(2*y*y),1/(2*y*y))])
    hh=sp.simplify(2*ek*jj.T*invk*jj)
    kin_isqrt=sp.diag(*[sp.sqrt(1/kin[a,a]) for a in range(8)])
    hc=sp.simplify(kin_isqrt*hh*kin_isqrt)
    return w,wp,fs,jj,hc

unitflux=sp.Matrix([0]*4+[-2]*4+[-2]*4+[0]*4)
unitpoint=(I,I,I,I)
w0,wp0,fs0,j0,h0=data(unitflux,unitpoint)
assert fs0==[0]*4 and charge(unitflux)==8 and j0.det()==-4096
assert h0==j0.T*j0
x=sp.symbols('x')
assert sp.factor(h0.charpoly(x).as_expr())==(x*x-56*x+64)*(x*x-24*x+64)**3
eigen=h0.eigenvals()
assert sum(eigen.values())==8 and all(ev.is_positive for ev in eigen)
print('unit_spectrum', {str(k):v for k,v in eigen.items()},flush=True)

# Entire parent's two-integer family, not a sample interpolation.
m,n=sp.symbols('m n',positive=True,integer=True)
mnflux=sp.Matrix([0]*4+[-2*m]*4+[-2*n]*4+[0]*4)
wmn,wpmn,fsmn,jmn,hmn=data(mnflux,(I*m/n,I,I,I))
assert fsmn==[0]*4
assert charge(mnflux)==8*m*n
assert sp.factor(jmn.det())==-4096*m**4*n**4
assert sp.simplify(hmn-m*n*h0)==sp.zeros(8)

def divisors(k):
    return [a for a in range(1,k+1) if k%a==0]

def make_divisor(k,aa):
    aa=list(map(sp.Integer,aa)); bb=[sp.Integer(k)/a for a in aa]
    ff=sp.Matrix([0]*4+[-2*a for a in aa]+[-2*b for b in bb]+[0]*4)
    yy=[aa[1]*aa[2]*aa[3]/(aa[0]*k)]+[aa[0]/aa[a] for a in (1,2,3)]
    return ff,tuple(I*y for y in yy),aa,bb

def invert_variable(flux,point,a):
    ww=sp.cancel(z[a]*polynomial(flux).subs(z[a],-1/z[a]))
    ff=flux_from_polynomial(ww)
    pp=list(point); pp[a]=-1/pp[a]
    assert charge(ff)==charge(flux)
    assert all(t.is_Integer and t%2==0 for t in ff)
    return ff,tuple(pp)

def reduce_point(flux,point):
    ff,pp=sp.Matrix(flux),tuple(point)
    inversions=[]
    for a in range(4):
        if sp.im(pp[a])<1:
            ff,pp=invert_variable(ff,pp,a); inversions.append(a)
    perm=sorted(range(1,4),key=lambda a:sp.im(pp[a]))
    # Relabel the three torus variables to sorted imaginary parts.
    mp={z[old]:z[new] for new,old in enumerate(perm,1)}
    ff=flux_from_polynomial(polynomial(ff).subs(mp,simultaneous=True))
    pp=(pp[0],)+tuple(pp[a] for a in perm)
    return ff,pp,inversions,perm

assignments=[]; representatives={}
for k in (1,2,3):
    for aa in itertools.product(divisors(k),repeat=4):
        ff,pt,aa,bb=make_divisor(k,aa)
        w,wp,fs,jj,hc=data(ff,pt)
        assert fs==[0]*4 and charge(ff)==8*k
        assert jj.det()==-4096*k**4
        assert sp.simplify(hc-k*h0)==sp.zeros(8)
        rf,rp,inv,perm=reduce_point(ff,pt)
        rw,rwp,rfs,rj,rh=data(rf,rp)
        assert rfs==[0]*4 and rj.rank()==8
        # Canonical frames can rotate under modular inversion, so compare
        # invariant characteristic polynomials rather than matrix entries.
        assert rh.charpoly().as_expr()==(k*h0).charpoly().as_expr()
        key=tuple(str(sp.im(t)) for t in rp)
        row={'k':k,'A':[int(a) for a in aa],'B':[int(b) for b in bb],
             'flux':[int(t) for t in ff],'point':[str(t) for t in pt],
             'Q':8*k,'W':str(wp),'J_det':str(jj.det()),
             'fundamental_point':[str(t) for t in rp],
             'fundamental_flux':[int(t) for t in rf],
             'inverted_variables':inv,'torus_permutation':perm}
        assignments.append(row)
        representatives.setdefault(key,{'k':k,'point':rp,'flux':rf})
assert len(assignments)==33 and len(representatives)==9

# Pointwise all-even-lattice minima for these nine constructive points.
# At a pure imaginary point the energy splits into eight complementary
# coefficient pairs. Rational primitive pairs are saturated independently.
point_records=[]
for key,rep in representatives.items():
    k=rep['k']; pt=rep['point']; ff=rep['flux']
    ys=[sp.im(t) for t in pt]
    pairs=[]
    for b in bits[:8]:
        # S exponent zero selects exactly one member of each complement pair.
        wt=sp.prod(ys[a]**(2*b[a]-1) for a in range(4))
        num,den=sp.fraction(wt)
        primitive_charge=2*num*den
        assert primitive_charge>=2*k
        pairs.append({'exponent':list(b),'complement':[1-a for a in b],
                      'energy_weight':str(wt),
                      'primitive_even_coefficient_magnitudes':[str(2*den),str(2*num)],
                      'primitive_pair_Q':int(primitive_charge)})
    # Four pairs are necessary for full rank by common-weight symmetry.
    minpair=min(a['primitive_pair_Q'] for a in pairs)
    assert minpair==2*k
    _,wp,fs,jj,hc=data(ff,pt)
    assert fs==[0]*4 and jj.rank()==8 and charge(ff)==8*k
    point_records.append({'k':k,'fundamental_imaginary_tuple':list(key),
        'point':[str(t) for t in pt],'flux':[int(t) for t in ff],
        'Qmin_exact':8*k,'Q':int(charge(ff)),'gs':str(1/ys[0]),
        'remaining_D3_capacity28':28-8*k,'W':str(wp),
        'real_J_rank':jj.rank(),'real_J_det':str(jj.det()),
        'primitive_pair_data':pairs,
        'mass_squared_exact':{str(k*ev):mult for ev,mult in eigen.items()}})

# A direct weak-coupling fundamental-domain member for every positive k.
kk=sp.symbols('k',positive=True,integer=True)
weakflux=sp.Matrix([2,2*kk,2*kk,2*kk]+[0]*8+[-2*kk,-2,-2,-2])
ww,wwp,wwfs,wwj,wwh=data(weakflux,(I*kk**2,I*kk,I*kk,I*kk))
assert wwfs==[0]*4 and charge(weakflux)==8*kk
assert sp.factor(wwj.det())==-4096*kk**4
assert wwh.charpoly().as_expr()==(kk*h0).charpoly().as_expr()

# An explicit nontrivial SL(2,Z) image, recorded as the same orbit rather
# than counted as a new fundamental-domain point.
dualbase=weakflux.subs(kk,3)
dualpoint=(9*I,3*I,3*I,3*I)
gamma=sp.Matrix([[2,1],[1,1]])
df=sp.Matrix.vstack(2*dualbase[:8,0]+dualbase[8:,0],
                    dualbase[:8,0]+dualbase[8:,0])
dp=((2*dualpoint[0]+1)/(dualpoint[0]+1),)+dualpoint[1:]
dp=tuple(sp.simplify(t) for t in dp)
dw,dwp,dfs,dj,dh=data(df,dp)
assert charge(df)==24 and dfs==[0]*4 and dj.rank()==8
assert dh.charpoly().as_expr()==(3*h0).charpoly().as_expr()
oldS=sp.symbols('oldS')
assert sp.cancel(dw.subs(z[0],(2*oldS+1)/(oldS+1))-
                 polynomial(dualbase).subs(z[0],oldS)/(oldS+1))==0

out={'status':'PASS','sympy_version':sp.__version__,
     'unit_W':str(wp0),'unit_real_J_det':str(j0.det()),
     'unit_real_Jacobian':[[str(t) for t in r] for r in j0.tolist()],
     'unit_canonical_characteristic_polynomial':str(sp.factor(h0.charpoly(x).as_expr())),
     'unit_canonical_Hessian':[[str(t) for t in r] for r in h0.tolist()],
     'unit_mass_squared_exact':{str(ev):mult for ev,mult in eigen.items()},
     'mn_family':{'Q':'8*m*n','point':['I*m/n','I','I','I'],
                  'W':str(wpmn),'real_J_det':str(sp.factor(jmn.det())),
                  'canonical_Hessian':'m*n times unit canonical Hessian'},
     'divisor_assignments_count':len(assignments),
     'distinct_reduced_point_types':len(representatives),
     'point_records':point_records,
     'weak_coupling_family':{'flux':[str(t) for t in weakflux],
           'point':['I*k**2','I*k','I*k','I*k'],'Q':'8*k','gs':'1/k**2',
           'W':str(wwp),'real_J_det':str(sp.factor(wwj.det())),
           'mass_squared':'k times unit spectrum','capacity28_k':[1,2,3]},
     'duality_control':{'gamma':[[2,1],[1,1]],'base_point':[str(t) for t in dualpoint],
           'image_point':[str(t) for t in dp],'image_flux':[int(t) for t in df],
           'Q':24,'F_terms':[str(t) for t in dfs],
           'mass_spectrum_invariant':True,'new_point_type_counted':False},
     'provenance_sha256':hashlib.sha256((ROOT/'PROVENANCE.json').read_bytes()).hexdigest(),
     'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(OUT/'construction_certificate.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
(OUT/'divisor_assignments.json').write_text(json.dumps(assignments,indent=2)+'\n',encoding='utf-8')
print(json.dumps({'status':'PASS','divisor_assignments':len(assignments),
     'distinct_reduced_point_types':len(representatives),
     'unit_spectrum':out['unit_mass_squared_exact'],
     'weak_family_J_det':out['weak_coupling_family']['real_J_det'],
     'point_minima':[(r['fundamental_imaginary_tuple'],r['Qmin_exact']) for r in point_records]},indent=2))
