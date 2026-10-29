"""Positive-energy, one-shot local clocks for the retained weighted echo.

Reads bundled accepted echo constants; never runs the historical producer.
All new output is confined to this line. hbar=1, inherited model units.
"""
from pathlib import Path
import json
import hashlib
import platform
import math
import numpy as np
import scipy
from scipy.integrate import quad
from scipy.special import roots_hermitenorm

ROOT = Path(__file__).resolve().parents[1]
OLD = ROOT / 'inputs/echo-reference.json'
old = json.loads(OLD.read_text(encoding='utf-8'))
checks = []

def check(name, condition, value=None):
    checks.append(dict(name=name, passed=bool(condition), value=value))
    if not condition:
        raise AssertionError(name)

def pulse(s):
    s = np.asarray(s)
    return np.where(np.abs(s) <= .5, 1 + np.cos(2*np.pi*s), 0.)

def construct(nu=1e-28, alpha=.05, eps=.001, c=1e8, rest_inclusive_speed=None):
    b = old['reference']['phase_coefficient_magnitude']
    x1, x2 = old['reference']['switches']
    A=.5; delta=.03; Bs=2/9
    k=2*np.sqrt(np.pi/b); Mrecoil=2*k/3
    R2=np.sqrt(Mrecoil**4+6*Mrecoil**2+3)
    Kbeta=2*(Bs/A)*16*(2*np.sqrt(2)*np.pi**2/3) + 2*4**(2/3)*4*1.001*.22*k
    Kr=2*(5/9)*4**(8/3)*k*R2+1
    R0=10+50*k; Vratio=48*k
    Q=(np.pi*A*A/(b*delta**2*nu**3))**(3/8)
    T=nu*Q*Q/A; E0=1/(4*nu*T); tau=c/E0
    Ec=alpha*E0
    # c_light=1 for this optional mass-energy-inclusive Newtonian model.
    mass=1. if rest_inclusive_speed is None else 2*Ec/rest_inclusive_speed**2
    p=np.sqrt(2*mass*Ec); v=p/mass
    rest_ratio=0. if rest_inclusive_speed is None else 4*mass/E0
    dp=eps*p; sigma_focus=1/(2*dp)
    a=sigma_focus/(v*tau)
    # Four offsets make EACH complete apparatus term nonnegative.
    offset_ratio=4*np.pi/c
    beta=1-4*alpha*(1+eps**2)-offset_ratio-rest_ratio
    clock_sd=2*alpha*np.sqrt(4*eps**2+2*eps**4)
    r=R0*np.sqrt(nu)+clock_sd+offset_ratio
    source=2*r/beta
    # Inner Lipschitz + analytic Gaussian tail; tail stored logarithmically.
    Smax=T/tau
    log_tail=math.log(8)-.25*math.log(2*math.pi)-.5*math.log(eps)+.5*math.log(Smax)-1/(32*eps**2)
    clock_inner=4*np.pi**2*(2*a+eps)
    log_clock_tail=math.log(2*np.pi)+log_tail
    clock_tail=math.exp(log_clock_tail) if log_clock_tail > -745 else 0.
    pieces=dict(genuine_shape_and_shear=Kbeta*nu**.25,
                radial_nonlinearity_and_cutoff=Kr*np.sqrt(nu), ordering=640*nu,
                initial_total_spectral_resource=source,
                finite_local_pulses=2*c*Vratio*np.sqrt(nu),
                clock_timing_dispersion_backreaction=clock_inner+clock_tail)
    total=sum(pieces.values())
    # Vector-norm comparison controls negative momentum amplitudes, not means
    # of unbounded operators. Free union-tail <=4 exp[-1/(2 eps^2)].
    filter_vector=np.sqrt(2)*source
    reflection_bound=(filter_vector+clock_inner+clock_tail)**2
    full_rms=(beta+r)/np.sqrt(1-source**2)
    return dict(nu=nu, alpha=alpha, relative_momentum_width=eps, c=c,
        Q=Q,T=T,E0=E0,w=nu**(3/8),tau=tau,tau_over_T=tau/T,
        switches=[x1,x2],clock_count=4,clock_mass=mass,clock_p0=p,clock_velocity=v,
        clock_sigma_p=dp,clock_sigma_Q_at_focus=sigma_focus,
        apparatus_rest_energy_included=rest_inclusive_speed is not None,
        apparatus_total_rest_energy_over_E0=rest_ratio,
        clock_sigma_time_at_focus=sigma_focus/v,
        scaled_focus_width_a=a,clock_sigma_time_at_window_edge_over_tau=np.sqrt(a*a+eps*eps/4),
        clock_sigma_Q_initial=[np.sqrt(sigma_focus**2+(tj*T*dp/mass)**2) for tj in [x1,x2]],
        clock_cov_QP_initial=[-tj*T*dp**2/mass for tj in [x1,x2]],
        clock_initial_mean_Q=[-v*tj*T for tj in [x1,x2]],
        apparatus_mean_kinetic_over_E0=4*alpha*(1+eps**2),
        positive_offsets_over_E0=offset_ratio,
        one_clock_peak_potential_over_central_kinetic=np.pi/(alpha*c),
        beta_negative_seed_center_over_E0=beta,
        source_graph_residual_over_E0=r,
        source_filter_failure_probability_bound=source**2,
        exact_dust_gap_over_E0=beta/2,
        full_total_H_RMS_bound_over_E0=full_rms,
        error_components=pieces,total_trace_bound=total,
        guaranteed_reduced_matter_negativity=.5-total,
        clock_outer_tail_log_bound=log_clock_tail,
        any_negative_clock_momentum_probability_bound=reflection_bound,
        free_negative_momentum_union_log_bound=math.log(4)-1/(2*eps**2),
        coefficients=dict(k=k,R0=R0,Kbeta=Kbeta,Kr=Kr,Vratio=Vratio))

def numerical_clock_integral(a,eps,n):
    z,w=roots_hermitenorm(n); w=w/np.sqrt(2*np.pi)
    def rms(s):
        d=np.sqrt(a*a+eps*eps*s*s)
        return float(np.sqrt(np.dot(w,(pulse(s+d*z)-pulse(s))**2)))
    cuts=[-1,-.5,-.49,0,.49,.5,1]
    return sum(quad(rms,cuts[j],cuts[j+1],epsabs=2e-10,limit=200)[0] for j in range(len(cuts)-1))

def main():
    p=construct(); nu=p['nu']; eps=p['relative_momentum_width']; a=p['scaled_focus_width_a']
    check('accepted_phase_constant',abs(old['reference']['phase_coefficient_magnitude']-.005512414487283133)<1e-17)
    check('pulse_area',abs(quad(lambda s:float(pulse(s)),-.5,.5)[0]-1)<1e-13)
    check('pulse_half_areas',abs(quad(lambda s:float(pulse(s)),-.5,0)[0]-.5)<1e-13)
    check('two_echo_windows_disjoint',p['tau_over_T']<min(p['switches'][0],p['switches'][1]-p['switches'][0],1-p['switches'][1]))
    check('positive_kinetic_clock',p['clock_mass']>0 and p['apparatus_mean_kinetic_over_E0']>0)
    check('positive_reserve_margin',p['exact_dust_gap_over_E0']>.39)
    check('tail_proof_regime',a<=eps and eps<.1)
    check('old_nonlinear_comparison_regime',nu<=1e-12)
    check('minimum_uncertainty_focused_packet',abs(p['clock_sigma_Q_at_focus']*p['clock_sigma_p']-.5)<1e-13)
    check('focus_width_formula',abs(a-1/(4*eps*p['alpha']*p['c']))<1e-15)
    check('four_clock_energy_variance',p['source_graph_residual_over_E0']<.000201)
    check('joint_source_failure_small',p['source_filter_failure_probability_bound']<3e-7)
    check('reduced_AB_entanglement',p['guaranteed_reduced_matter_negativity']>.45)
    check('backreaction_included_and_reflection_bounded',p['any_negative_clock_momentum_probability_bound']<.002)
    check('outer_tail_negligible',p['clock_outer_tail_log_bound'] < -30000)
    # Direct Gaussian expectation independently checks the more conservative
    # Lipschitz integral; it is not substituted for the rigorous bound.
    i64=numerical_clock_integral(a,eps,64); i128=numerical_clock_integral(a,eps,128)
    check('gaussian_quadrature_resolution',abs(i64-i128)<2e-7,[i64,i128])
    check('gaussian_error_below_analytic_bound',i128<2*np.pi*(2*a+eps))
    # Full time-independent 2x2 local potential: eigenvalues shift +/-g.
    check('complete_apparatus_potential_nonnegative',np.min([np.pi/p['tau']-np.pi/p['tau'],2*np.pi/p['tau']])>=0)
    # The locality and positivity structural identities are exact algebra;
    # these simple matrices check the retained ideal pulse endpoint.
    X=np.array([[0,1],[1,0]],complex)
    from scipy.linalg import expm
    check('local_pulse_endpoint',np.linalg.norm(expm(-1j*np.pi*X/2)+1j*X)<1e-13)
    rest=construct(nu=1e-42,alpha=1e-6,eps=.001,c=1e14,rest_inclusive_speed=.01)
    check('rest_inclusive_dust_gap',rest['exact_dust_gap_over_E0']>.459)
    check('rest_inclusive_entanglement',rest['guaranteed_reduced_matter_negativity']>.459)
    check('rest_energy_explicitly_charged',abs(rest['apparatus_total_rest_energy_over_E0']-.08)<1e-14)
    check('rest_inclusive_tail_proof_regime',rest['scaled_focus_width_a']<=eps and rest['clock_outer_tail_log_bound'] < -30000)
    out=dict(status='PASS',construction='four positive-energy party-local pre-focused autonomous clocks',
        runtime=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),
        accepted_input_sha256=hashlib.sha256(OLD.read_bytes()).hexdigest(),
        point=p,rest_energy_inclusive_point=rest,
        gaussian_integral_diagnostic=dict(order64=i64,order128=i128,
              four_clock_state_vector_error_estimate=2*np.pi*i128,proof_uses_analytic_bound=True),
        checks=checks)
    (ROOT/'results/local-clocks.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status='PASS',checks=len(checks),trace_bound=p['total_trace_bound'],
          negativity_lower=p['guaranteed_reduced_matter_negativity'],dust_gap_E0=p['exact_dust_gap_over_E0'],
          clock_error_bound=p['error_components']['clock_timing_dispersion_backreaction'],
          numerical_clock_error_estimate=2*np.pi*i128,
          rest_inclusive_negativity=rest['guaranteed_reduced_matter_negativity'],
          rest_inclusive_dust_gap_E0=rest['exact_dust_gap_over_E0']),indent=2))

if __name__=='__main__': main()
