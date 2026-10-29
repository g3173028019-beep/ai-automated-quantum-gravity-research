"""Exact consistency checks supporting the canonical finite-scalar manuscript.

Requires Python 3 and SymPy. No network, numerical integration, or PDE solver.
The abstract charged-adjoint calculation is dimension-symbolic. Tensor tests
use a pointwise orthonormal Lorentz frame; they do not prove a PDE theorem.
"""
import json
import sympy as s


def zero(x):
    return all(s.simplify(v) == 0 for v in x) if isinstance(x, s.MatrixBase) else s.simplify(x) == 0


# Schmid (2.25), (2.33): <h,k>_S2 = 2 h_ab k^ab, while
# B_phys = 2 <Ih,k>_S2 + 4 delta_AB phi^A psi^B.
g = s.diag(-1, 1, 1, 1)
def sym(prefix):
    v = {(i, j): s.Symbol(f'{prefix}{i}{j}', real=True) for i in range(4) for j in range(i, 4)}
    return s.Matrix(4, 4, lambda i, j: v[min(i, j), max(i, j)])
def tr(h):
    return s.trace(g*h)
def I(h):
    return h-g*tr(h)/2
def contract(h, k):
    return s.trace(g*h*g*k)
def pair_S2(h, k):
    return 2*contract(h, k)
def B_metric(h, k):
    return 2*pair_S2(I(h), k)

h, H = sym('h'), sym('H')
z = s.Symbol('z', real=True)
mixing_adjoint_residual = B_metric(h, 2*H*z)-4*pair_S2(H, I(h))*z
assert zero(mixing_adjoint_residual)
# Dropping the leading factor two in B_metric is detected.
wrong_normalization = pair_S2(I(h), 2*H*z)-4*pair_S2(H, I(h))*z
assert not zero(wrong_normalization)
# Criticality is tr Hess Phi = 0, while potential Hessian need not vanish.
H0 = H.subs(H[3,3], H[0,0]-H[1,1]-H[2,2])
assert zero(I(I(h))-h)
assert zero(I(H0)-H0)
assert not zero(I(g)-g)

# General charged-adjoint identity. Shapes stand only for well-typed maps;
# no entries, numerical spectra, or special matrix realization are supplied.
n, m, r = s.symbols('n m r', positive=True, integer=True)
Q1, Q2 = s.MatrixSymbol('Q1',m,m), s.MatrixSymbol('Q2',n,n)
K, R = s.MatrixSymbol('K',n,m), s.MatrixSymbol('R',r,n)
Kdag = s.Inverse(Q1)*K.adjoint()*Q2
Rdag = s.Inverse(Q2)*R.adjoint()
charged_residual = (Q1*Kdag*Rdag-(R*K).adjoint()).doit()
assert charged_residual == s.ZeroMatrix(m,r)

# Generic finite-rank defect: Q Hermitian; Pi K = Pi^dag K = 0.
# These assumptions give Q Pi K = Pi* Q K = Pi* Q Pi K = 0, term by term.
Pi = s.MatrixSymbol('Pi', n,n)
expanded_defect = (Q2-(s.Identity(n)-Pi.adjoint())*Q2*(s.Identity(n)-Pi)).expand()
expected_defect = Q2*Pi+Pi.adjoint()*Q2-Pi.adjoint()*Q2*Pi
assert (expanded_defect-expected_defect).doit() == s.ZeroMatrix(n,n)
# The spectral-part implication is a finite-dimensional theorem proved in
# ALGEBRA-NOTES.md, rather than asserted from a numerical example.

# Exact closed-FLRW hedgehog: a'=aH, f'=u, u'=-3Hu-3f/a^2.
a, f, u, hub, kap, lam = s.symbols('a f u H kappa Lambda', real=True, nonzero=True)
rho = (u*u+3*f*f/a**2)/2
p = (u*u-f*f/a**2)/2
udot = -3*hub*u-3*f/a**2
hdot = 1/a**2-kap*(rho+p)/2
flow = {a:a*hub, f:u, u:udot, hub:hdot}
def Dt(expr):
    return s.factor(sum(s.diff(expr,x)*v for x,v in flow.items()))
C = 3*(hub**2+a**(-2))-kap*rho-lam
assert zero(Dt(rho)+3*hub*(rho+p))
assert zero(Dt(C))
assert zero(-(2*hdot+3*hub**2+a**(-2))-kap*p+lam+C)
initial={a:1,f:1,u:3,hub:1,kap:1,lam:0}
assert zero(C.subs(initial))
assert [rho.subs(initial),p.subs(initial),hdot.subs(initial),udot.subs(initial)] == [6,4,-4,-12]
# At these data, Hess_ti X=(u-Hf) d_iY=2 d_iY, so mixing is nonzero.
assert (u-hub*f).subs(initial) == 2
# Radial/tangential dX decomposition gives |det dX|=|u f^3| times
# the oriented S3 volume density, nonzero wherever u*f != 0.
assert (u*f**3).subs(initial) == 3

# Optional fifth field psi=0 with V=mu^2 psi^2/2 leaves the background
# critical but gives a nonzero potential-Hessian entry.
psi, mu = s.symbols('psi mu', real=True)
V = mu**2*psi**2/2
assert s.diff(V,psi).subs(psi,0) == 0
assert s.diff(V,psi,2).subs(psi,0) == mu**2

print(json.dumps({
    'sympy_version':s.__version__,
    'normalization_mixing_adjoint':'PASS: raw metric/scalar weights 4/4',
    'normalization_negative_control':'PASS: missing metric factor detected',
    'critical_hessian_and_trace_reversal':'PASS',
    'charged_adjoint_identity':'PASS: dimension-symbolic',
    'physical_charge_defect_expansion':'PASS: dimension-symbolic',
    'flrw_continuity_constraint_propagation_acceleration':'PASS: exact symbolic',
    'flrw_initial_values':{'rho':6,'p':4,'Hdot':-4,'fddot':-12,'mixed_Hessian_coefficient':2,'oriented_rank_factor':3},
    'critical_nonzero_potential_hessian':'PASS: spectator mass squared retained',
    'scope':'algebraic consistency only; no projector/Fredholm/Hadamard PDE proof'
},indent=2))
