"""Constrained scalar/vector reduction and clock-foliated EFT repair.
Units tau=t/a0; action divided by a0**2/kappa; harmonics have unit L2 norm.
The old rounds are read only. Run directly with Python, SymPy, NumPy, SciPy.
"""
from pathlib import Path
import json
import sympy as s
import numpy as np
from scipy.linalg import eigvals, expm, sqrtm

OUT=Path(__file__).resolve().parents[1]/'results'
OUT.mkdir(exist_ok=True)
a,m,L=s.symbols('alpha m L', real=True)
b=(2-a)/4; om=s.sqrt(3); Lam=(6+a)/4
v,u,z,q,p=s.symbols('v u z q psi',real=True)
dv,du,dz,dq,dp=s.symbols('dv du dz dq dpsi',real=True)
A,B=s.symbols('A B',real=True)
x=s.Matrix([v,u,z,q,p]); dx=s.Matrix([dv,du,dz,dq,dp])
grav=-3*dp**2-2*L*dp*B-L*B**2+s.Rational(3,2)*(A+p)**2+2*L*p*(A+p)-L*p**2-Lam/2*(A+3*p)**2
kin=b/2*((dv-om*z)**2+(dz+om*v)**2+L*((du-om*q-B)**2+(dq+om*u)**2)-2*om*B*L*(z-q))
kin+=b*om*(3*p-A)*(dz+om*v)+(3*b+a)/4*(3*p-A)**2
spat=-b/2*((L+3)*(v*v+z*z)-4*L*(v*u+z*q)+L*(L-1)*(u*u+q*q))-b*(A+p)*(3*v-L*u)-3*b/4*(A+p)**2
lag=s.expand(grav+kin+spat-b*m*v*v/2)
cons=s.solve([s.diff(lag,A),s.diff(lag,B)],(A,B),simplify=False)
red=s.factor(lag.subs(cons))
K=s.simplify(s.hessian(red,dx))
D=s.simplify(s.Matrix([[s.diff(red,dx[i],x[j]) for j in range(5)] for i in range(5)]))
G=s.simplify(D-D.T); H=s.simplify(-s.hessian(red,x))
acc=s.simplify(-K.inv()*(G*dx+H*x))
def dt(expr):
    return s.simplify(sum(s.diff(expr,x[i])*dx[i]+s.diff(expr,dx[i])*acc[i] for i in range(5)))
clock=s.simplify(dt(cons[A])-3*dp-L*cons[B])
kg=s.simplify(s.Matrix([acc[0]-2*om*dz+(L+m)*v-2*L*u-6*(p-cons[A]),acc[1]-2*om*dq+(L-4)*u-2*v-dt(cons[B])-cons[A]-p,acc[2]+2*om*dv+L*z-2*L*q-om*clock,acc[3]+2*om*du+(L-4)*q-2*z-2*om*cons[B]]))
assert clock==0 and kg==s.zeros(4,1)

# Repair Lagrangian: -lambda_bar E^2 -zeta_bar Q^2 -V_Gram.
lam=s.Rational(1,10); zet=s.Integer(1); h=s.Integer(6)
strain=s.Rational(8,3)*b*b*L*(L-3)*u*u
compression=4*b*b*L*L*(u+p)**2
gram=b*L*h*q*q/2
delta=lam*strain+zet*compression+gram
Hr=s.simplify(H+s.hessian(delta,x))
params={a:s.Rational(1,10),m:-s.Rational(207,25)}
accr=s.simplify(-K.inv()*(G*dx+Hr*x))
# The repair has no linear clock dependence: K_0i=0 and the background
# clock slices have zero extrinsic curvature. Check the omitted clock EOM.
def dt_r(expr):
    return s.simplify(sum(s.diff(expr,x[i])*dx[i]+s.diff(expr,dx[i])*accr[i] for i in range(5)))
clock_r=s.simplify(dt_r(cons[A])-3*dp-L*cons[B])
kg_r=s.simplify(s.Matrix([accr[0]-2*om*dz+(L+m)*v-2*L*u-6*(p-cons[A]),accr[1]-2*om*dq+(L-4)*u-2*v-dt_r(cons[B])-cons[A]-p+s.diff(delta,u)/(b*L),accr[2]+2*om*dv+L*z-2*L*q-om*clock_r,accr[3]+2*om*du+(L-4)*q-2*z-2*om*cons[B]+s.diff(delta,q)/(b*L)]))
assert clock_r==0 and kg_r==s.zeros(4,1)
# Explicit exceptional-gauge identity: the dropped psi equation is implied
# by the other four at L=3, including compression and the symmetry-breaking
# Gram term. This is stronger than treating a singular K as invertible.
ac1=s.simplify(-K[:4,:4].subs(L,3).inv()*(G[:4,:4].subs(L,3)*dx[:4,0]+Hr[:4,:4].subs(L,3)*x[:4,0]))
all1=s.simplify((K*s.Matrix([*ac1,0])+G*dx+Hr*x).subs({L:3,p:0,dp:0}))
assert all1==s.zeros(5,1)
def exact_minor_tests(mat,offset):
    tests=[]
    for j in range(1,mat.rows+1):
        det=s.factor(mat[:j,:j].det())
        poly=s.Poly(s.expand(det.subs(L,s.Symbol('t')+offset)),s.Symbol('t'))
        assert all(c>=0 for c in poly.all_coeffs()) and poly.eval(0)>0
        tests.append({'minor_order':j,'polynomial':str(s.factor(det)),'shifted_coefficients':[str(c) for c in poly.all_coeffs()]})
    return tests
tests={}
for name,mat in [('kinetic',K),('real_stiffness',Hr.extract([0,1,4],[0,1,4])),('imaginary_stiffness',Hr.extract([2,3],[2,3]))]:
    tests[name]=exact_minor_tests(mat.subs(params),8)
# Exceptional ell=1: residual scalar spatial diffeo fixed by psi=0.
k1=K[:4,:4].subs(params).subs(L,3); hr1=Hr[:4,:4].subs(params).subs(L,3)
tests['ell1_K_minors']=[str(s.factor(k1[:j,:j].det())) for j in range(1,5)]
tests['ell1_H_minors']=[str(s.factor(hr1[:j,:j].det())) for j in range(1,5)]
assert all(s.Rational(t)>0 for t in tests['ell1_K_minors']+tests['ell1_H_minors'])
tests['original_scalar_ell_ge_4_real_H']=exact_minor_tests(H.extract([0,1,4],[0,1,4]).subs(params),24)
tests['original_scalar_ell_ge_4_imag_H']=exact_minor_tests(H.extract([2,3],[2,3]).subs(params),24)
root_cert={}
rr,ww=s.symbols('growth omega_squared')
for ell in [1,2,3]:
    sub={**params,L:ell*(ell+2)}
    n=4 if ell==1 else 5
    poly=s.Poly(s.expand((rr*rr*K[:n,:n]+rr*G[:n,:n]+H[:n,:n]).subs(sub).det(method='domain-ge')),rr)
    coeff={degree[0]//2:coef*(-1)**(degree[0]//2) for degree,coef in poly.terms()}
    assert all(degree[0]%2==0 for degree,coef in poly.terms())
    wp=s.Poly.from_dict({(degree,):coef for degree,coef in coeff.items()},ww).monic()
    nonzero=wp.exquo(s.Poly(ww,ww)) if wp.eval(0)==0 else wp
    root_cert[str(ell)]={'omega_squared_polynomial':str(wp.as_expr()),'positive_real_root_count':int(s.polys.polytools.count_roots(nonzero,0,s.oo)),'negative_real_root_count':int(s.polys.polytools.count_roots(nonzero,-s.oo,0)),'zero_root':wp.eval(0)==0}

def numeric_system(ell,repaired):
    sub={**params,L:ell*(ell+2)}
    kk=np.array(K.subs(sub),float);gg=np.array(G.subs(sub),float);hh=np.array((Hr if repaired else H).subs(sub),float)
    if ell==1:kk=kk[:4,:4];gg=gg[:4,:4];hh=hh[:4,:4]
    n=len(kk)
    return kk,gg,hh,np.block([[np.zeros((n,n)),np.eye(n)],[-np.linalg.solve(kk,hh),-np.linalg.solve(kk,gg)]])
spectra={}
for ell in [1,2,3,4,5,10]:
    spectra[str(ell)]={}
    for repaired in [False,True]:
        kk,gg,hh,aa=numeric_system(ell,repaired)
        ev=eigvals(aa)
        spectra[str(ell)]['repaired' if repaired else 'original']={'growth_rates':sorted(float(t.real) for t in ev),'positive_frequencies':sorted(float(t.imag) for t in ev if t.imag>1e-6),'min_K_eigenvalue':float(min(np.linalg.eigvalsh(kk))),'min_H_eigenvalue':float(min(np.linalg.eigvalsh(hh)))}
uv=[]
for ell in [30,100,300]:
    _,_,_,ax=numeric_system(ell,True)
    uv.append({'ell':ell,'omega_over_sqrtL':sorted(float(ev.imag/np.sqrt(ell*(ell+2))) for ev in eigvals(ax) if ev.imag>0)})

# Original ell=2 is a physical internal U(4) orientation zero mode, not a
# lapse/shift zero or floating point near-zero frequency.
H2=H.subs(params).subs(L,8)
zero=s.Matrix([0,0,1,s.Rational(1,2),0])
assert H2*zero==s.zeros(5,1)
assert H2.rank()==4

# Gaussian ground state and a bounded finite classical excitation at ell=1.
# Use canonical p=K xdot +D x after dropping a total derivative, D=G/2.
kk,gg,hh,aa=numeric_system(1,True); nn=len(kk); dd=gg/2
ki=np.linalg.inv(kk)
hc=np.block([[hh+dd.T@ki@dd,-dd.T@ki],[-ki@dd,ki]])
jj=np.block([[np.zeros((nn,nn)),np.eye(nn)],[-np.eye(nn),np.zeros((nn,nn))]])
hs=np.real_if_close(sqrtm(hc)); hi=np.linalg.inv(hs)
absq=np.real_if_close(sqrtm(-(hs@jj@hs)@(hs@jj@hs)))
cov=np.real_if_close(hi@absq@hi/2)
assert np.linalg.norm((jj@hc)@cov+cov@(jj@hc).T)<1e-9
assert min(np.linalg.eigvalsh(cov+0.5j*jj))>-1e-10
assert np.linalg.norm(cov@jj@cov-jj/4)<1e-9
init=np.array([1,0,0,0,0,0,0,0],float)*1e-4
times=np.linspace(0,100,1001)
states=np.array([expm(aa*t)@init for t in times])
energies=np.array([.5*st[nn:]@kk@st[nn:]+.5*st[:nn]@hh@st[:nn] for st in states])
E0=energies[0]; energy_metric=np.block([[hh,np.zeros_like(hh)],[np.zeros_like(kk),kk]])
bound=float(np.sqrt(2*E0/min(np.linalg.eigvalsh(energy_metric))))
reconstruct=s.lambdify([*x[:4,0],*dx[:4,0]],s.Matrix([cons[A],cons[B]]).subs({**params,L:3,p:0,dp:0}),'numpy')
metric=np.array([np.array(reconstruct(*st),float).ravel() for st in states])
np.savetxt(OUT/'controlled_ell1_evolution.csv',np.column_stack([times,states,metric,energies]),delimiter=',',header='tau,v,u,z,q,dv,du,dz,dq,A,B,E',comments='')

# Full vector family: d=n(n+2)-3, f=d/(d+2b). n=1 has V gauge.
vec=[];bv=19/40;ov=np.sqrt(3)
for n in [1,2,3,4,10]:
    d=n*(n+2)-3
    if n==1:
        vec.append({'n':n,'frequencies':[float(np.sqrt(12+6))],'gauge':'V is a Killing spatial diffeomorphism'})
    else:
        f=d/(d+2*bv); kv=np.diag([bv*f,bv]); gv=np.array([[0,-2*bv*ov*f],[2*bv*ov*f,0]])
        hv=np.diag([bv*d*(1+4*.1*bv),bv*(d+24*bv/(d+2*bv)+6)])
        av=np.block([[np.zeros((2,2)),np.eye(2)],[-np.linalg.solve(kv,hv),-np.linalg.solve(kv,gv)]])
        vec.append({'n':n,'frequencies':sorted(float(t.imag) for t in eigvals(av) if t.imag>0),'min_K':float(min(np.linalg.eigvalsh(kv))),'min_H':float(min(np.linalg.eigvalsh(hv)))})

output={'parameters':{'alpha':'1/10','m':'-207/25','g4':10,'lambda_bar':'1/10','zeta_bar':1,'h':6},'direct_original_covariant_checks':{'clock':str(clock),'four_KG':str(kg)},'all_tower_exact_positivity':tests,'spectra':spectra,'vector_repaired':vec,'ell2_physical_zero_vector':[str(t) for t in zero],'ell1_covariance':cov.tolist(),'ell1_state_checks':{'stationarity_residual':float(np.linalg.norm((jj@hc)@cov+cov@(jj@hc).T)),'purity_residual':float(np.linalg.norm(cov@jj@cov-jj/4)),'uncertainty_min_eigenvalue':float(min(np.linalg.eigvalsh(cov+.5j*jj))),'E0':E0,'max_energy_relative_drift':float(max(abs(energies-E0))/E0),'euclidean_phase_space_bound':bound,'max_sampled_phase_space_norm':float(max(np.linalg.norm(states,axis=1)))},'original_detH':str(s.factor(H.det()))}
output['direct_repaired_checks']={'clock':str(clock_r),'four_KG':str(kg_r),'exceptional_ell1_dropped_equation':str(all1)}
output['ell1_state_checks']['maximum_lapse_amplitude']=float(max(abs(metric[:,0])))
output['ell1_state_checks']['maximum_shift_amplitude']=float(max(abs(metric[:,1])))
output['original_low_mode_exact_root_counts']=root_cert
output['formal_large_ell_dispersion']=uv
# Actual finite-band scalar Gaussian: every harmonic of ell=1,2,3 can carry
# an independent copy of these matrices; no continuum microlocal claim.
band={}
for ell in [1,2,3]:
    bk,bg,bh,ba=numeric_system(ell,True);bn=len(bk);bd=bg/2;bki=np.linalg.inv(bk)
    bhc=np.block([[bh+bd.T@bki@bd,-bd.T@bki],[-bki@bd,bki]])
    bj=np.block([[np.zeros((bn,bn)),np.eye(bn)],[-np.eye(bn),np.zeros((bn,bn))]])
    bs=np.real_if_close(sqrtm(bhc));bsi=np.linalg.inv(bs)
    bc=np.real_if_close(bsi@sqrtm(-(bs@bj@bs)@(bs@bj@bs))@bsi/2)
    sr=float(np.linalg.norm((bj@bhc)@bc+bc@(bj@bhc).T))
    pr=float(np.linalg.norm(bc@bj@bc-bj/4))
    ue=float(min(np.linalg.eigvalsh(bc+.5j*bj)))
    assert sr<1e-8 and pr<1e-8 and ue>-1e-9
    band[str(ell)]={'spatial_degeneracy':(ell+1)**2,'canonical_covariance':bc.tolist(),'stationarity_residual':sr,'purity_residual':pr,'uncertainty_min_eigenvalue':ue}
output['finite_scalar_band_1_to_3']=band
# Tensor correction alone is new here; its original spectrum was already
# derived in round 02. EH kinetic coefficient is 1/8 for norm-one TT h.
output['tensor_repair']={'quadratic_E_squared':'b**2 h_TT**2','delta_omega_squared':'8*lambda_bar*b**2','delta_omega_squared_value':float(8*lam*s.Rational(19,40)**2),'repaired_omega_squared':'n*(n+2)+2261/2000'}
kin_inf=np.array([[2*1.9/6.1,2*1.9/6.1],[2*1.9/6.1,16/6.1]])
hinf=kin_inf+16/3*.1*(19/40)**2*np.diag([1,0])+8*(19/40)**2*np.ones((2,2))
output['principal_scalar_speeds']={'three_speeds':[1,1,1],'two_other_speeds':sorted(float(np.sqrt(t.real)) for t in eigvals(hinf,kin_inf)),'warning':'preferred-clock EFT cones differ from metric null cone'}
(OUT/'construction.json').write_text(json.dumps(output,indent=2),encoding='utf-8')
(OUT/'quadratic_action.txt').write_text('Unreduced L2 = '+str(lag)+'\nLapse = '+str(s.factor(cons[A]))+'\nShift = '+str(s.factor(cons[B]))+'\nK = '+str(K)+'\nG = '+str(G)+'\nH_original = '+str(H)+'\nH_repaired = '+str(Hr)+'\n',encoding='utf-8')
print(json.dumps({'exact_covariant_original_checks':'passed','all_scalar_harmonics_ell_ge_1_positive_K_H':'proved by exact principal minors','original_ell1':spectra['1']['original'],'original_ell3':spectra['3']['original'],'repaired_ell1':spectra['1']['repaired'],'state':output['ell1_state_checks']},indent=2))
