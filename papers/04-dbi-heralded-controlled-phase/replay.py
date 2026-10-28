"""Network-free frozen DBI gate replay; no original paths or refinement.

The saved full-half-line cell bounds are inherited evidence. This replay
checks their exact sums, initial mismatch, and new-time gate consequences;
it does not regenerate interval Taylor enclosures of the residual integrand.
"""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,math
from mpmath import iv
import sympy as s

ROOT=Path(__file__).resolve().parent
iv.dps=110
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def exact(x):
    f=F.from_float(float(x));return iv.mpf(f.numerator)/f.denominator
def hi(x):return math.nextafter(float(x.b),math.inf)
def lo(x):return math.nextafter(float(x.a),-math.inf)
def mag(z):return iv.sqrt(z.real**2+z.imag**2)
def box(x):return [lo(x),hi(x)]

# One exact symbolic check of the elementary tangent bound and sharp example.
e,t,a=s.symbols('eta t a',positive=True)
assert s.expand((t*t+a)-2*t*s.sqrt(a)-(t-s.sqrt(a))**2)==0
R=s.Matrix([[a,-e*s.sqrt(a)],[e*s.sqrt(a),a]])
assert ((R-s.eye(2)).T*(R-s.eye(2))-e**2*s.eye(2)).subs(a,1-e**2).applyfunc(s.simplify)==s.zeros(2)
assert (R.T*R-a*s.eye(2)).subs(a,1-e**2).applyfunc(s.simplify)==s.zeros(2)

prov=json.loads((ROOT/'PROVENANCE.json').read_text())
for name,h in prov['retained_file_sha256'].items():assert sha(ROOT/name)==h,(name,'hash mismatch')
trial=json.loads((ROOT/'inputs/continuum-trial.json').read_text())
tau=iv.mpf('107.61733846478977')
p0=iv.mpf('0.9177158589349711');p1=iv.mpf('0.69393706341036732')
targets=[iv.exp(iv.j*p0),iv.exp(iv.j*p1),-iv.exp(iv.j*(2*p1-p0))]
branch=[]
for row,d in zip(trial['branches'],targets):
    name=row['branch'];cert=json.loads((ROOT/f'inputs/residual-partition-{name}.json').read_text())
    assert cert['trial_sha256']==sha(ROOT/'inputs/continuum-trial.json')
    assert cert['integrator_sha256']==sha(ROOT/'reference/certify_trial_residual.py')
    assert cert['status']=='CERTIFIED_TARGET_MET'
    assert len(row['basis_vectors'])==16 and len(row['coefficients'])==8
    for j,r in enumerate(cert['residual_norm_upper_by_trial']):
        sq=F.from_float(cert['tail_norm_squared_upper_by_trial'][j])+sum((F.from_float(c[3][j]) for c in cert['cells']),F())
        assert F.from_float(r)**2>=sq
    sq=F()
    for n,v in enumerate(row['basis_vectors']):
        z=sum((F.from_float(x)*F.from_float(c) for x,c in zip(v,row['coefficients'])),F())-(n==0)
        sq+=z*z
    assert F.from_float(cert['initial_vector_error_upper'])**2>=sq
    amp=sum((exact(v)*exact(c)*iv.exp(iv.j*exact(lam)*tau) for v,c,lam in zip(row['basis_vectors'][0],row['coefficients'],row['energies'])),iv.mpc(0))
    slope=sum((abs(exact(c))*exact(r) for c,r in zip(row['coefficients'],cert['residual_norm_upper_by_trial'])),iv.mpf(0))
    err=exact(cert['initial_vector_error_upper'])+tau*slope
    dist=mag(amp-d)+err
    branch.append({'branch':name,'amplitude_real':box(amp.real),'amplitude_imag':box(amp.imag),
                   'residual_slope_upper':hi(slope),'Duhamel_error_upper':hi(err),'distance_upper':hi(dist)})

eta=iv.mpf('0.18972055982017294')
eta_saved=exact(max(r['distance_upper'] for r in branch))
assert eta.a>=eta_saved.b and eta.b<1
qmax=iv.mpf(4862)/729;M1=iv.sqrt(121+qmax*qmax*iv.mpf(4)/63)
extra=iv.mpf('.001')+iv.mpf('.001')+M1*iv.mpf('.001')+2*iv.mpf('.0001')
# The robust distance is certified directly from the branch bound plus extras;
# a second outward rounding of the nominal printed eta need not be added.
eta_rob=iv.mpf('0.20304819175570846')
assert (eta_saved+extra).b<=eta_rob.a and eta_rob.b<1
def guarantees(z):
    assert z.a>=0 and z.b<1
    return {'operator_distance_upper':hi(z),'all_input_success_lower':lo((1-z)**2),
            'all_reference_conditional_squared_fidelity_lower':lo(1-z*z),
            'conditional_Bell_negativity_lower':lo(iv.mpf('.5')-z*z),
            'successful_subnormalized_diamond_distance_upper':hi(2*z)}
out={'status':'PASS_FROZEN_HERALDED_GATE_SHARP_FIDELITY',
     'frozen_time':'107.61733846478977','phi0':'0.9177158589349711','phi1':'0.69393706341036732',
     'correction':'diag(1,exp(-i delta)) on each logical qubit; delta=phi1-phi0; standard Rz angle -delta',
     'sharp_lemma_exact_algebra':'PASS','input_hashes':'PASS',
     'inherited_residual_sums_and_initial_mismatch':'PASS','interval_digits':iv.dps,
     'branch_certificates':branch,'nominal':guarantees(eta),'robust':guarantees(eta_rob),
     'nominal_distance_bound_decimal':'0.18972055982017294',
     'robust_distance_bound_decimal':'0.20304819175570846',
     'M1_upper':hi(M1),'extra_operator_error_upper':hi(extra),
     'robust_box':{'source_normalized_pure_vector_norm':'.001','readout_normalized_pure_vector_norm':'.001',
                   'absolute_time_error':'.001','each_local_phase_error':'.0001'},
     'scope':'Common radial psi5 tensor arbitrary logical/reference input, exact invariant code and nominal Hamiltonians; complete P/Q instrument with failure flag.',
     'not_regenerated':'Old residual cellwise integrand and analytic-tail enclosures are retained certificates, not fresh continuum refinements.'}
(ROOT/'results/gate.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
print(json.dumps({k:out[k] for k in ['status','nominal','robust']},indent=2))
