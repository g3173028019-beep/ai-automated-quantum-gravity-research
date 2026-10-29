"""Direct 4D verification of an exact non-pp Einstein--scalar Walker family.

Coordinates and metric conventions are those of the retained Walker
construction. The finite solution, real repair of a cubic obstruction and
coordinate-invariant non-pp witness are new in this continuation.
"""
from pathlib import Path
import json
import time
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
u,v,z,w=s.symbols('u v z w',real=True)
alpha,beta,kappa,eps=s.symbols('alpha beta kappa epsilon',real=True,nonzero=True)
X=[u,v,z,w]
checks=[]

def c(q):
    return s.cancel(s.expand(q))

def ck(expr,name):
    values=list(expr) if isinstance(expr,s.MatrixBase) else [expr]
    residual=[c(x) for x in values]
    if any(x!=0 for x in residual):
        raise AssertionError((name,[str(x) for x in residual if x!=0]))
    checks.append(name)

def main():
    started=time.time()
    H=alpha*w+beta
    theta=kappa*eps**2*(alpha*w**5/20+beta*w**4/12)+v*v/(2*H)
    psi=eps*u
    box=lambda f:-2*s.diff(f,u,v)+2*s.diff(f,z,w)
    det=s.diff(theta,v,2)*s.diff(theta,w,2)-s.diff(theta,v,w)**2
    E=c(box(theta)-det)
    ck(det-kappa*eps**2*w*w,'exact nonlinear Hessian determinant')
    ck(E+kappa*eps**2*w*w,'exact sourced heavenly equation')
    g=s.Matrix([[s.diff(theta,w,2),-1,s.diff(theta,v,w),0],
                [-1,0,0,0],
                [s.diff(theta,v,w),0,s.diff(theta,v,2),1],
                [0,0,1,0]]).applyfunc(c)
    inverse=g.inv().applyfunc(c)
    ck(g.det()-1,'constant nondegenerate determinant')
    ck(g*inverse-s.eye(4),'exact finite inverse')
    Gamma=[[[c(sum(inverse[r,l]*(s.diff(g[l,j],X[i])+s.diff(g[l,i],X[j])
                                 -s.diff(g[i,j],X[l])) for l in range(4))/2)
             for j in range(4)] for i in range(4)] for r in range(4)]
    Ric=s.Matrix(4,4,lambda i,j:c(sum(
        s.diff(Gamma[r][i][j],X[r])-s.diff(Gamma[r][i][r],X[j])
        +sum(Gamma[r][r][l]*Gamma[l][i][j]-Gamma[r][j][l]*Gamma[l][i][r]
             for l in range(4)) for r in range(4))))
    stress=s.zeros(4);stress[0,0]=kappa*eps**2
    ck(Ric-stress,'direct four-dimensional Einstein equation, all 16 components')
    dpsi=s.Matrix([s.diff(psi,x) for x in X])
    ck((dpsi.T*inverse*dpsi)[0],'scalar gradient exactly null')
    scalarwave=c(sum(s.diff(inverse[i,j]*dpsi[j],X[i])
                     for i in range(4) for j in range(4)))
    ck(scalarwave,'direct finite Klein--Gordon equation')
    ck(s.trace(inverse*Ric),'Ricci scalar vanishes')
    Hess=s.Matrix(4,4,lambda i,j:c(s.diff(psi,X[i],X[j])
                                -sum(Gamma[r][i][j]*dpsi[r] for r in range(4))))
    minor=c(Hess.extract([0,2],[0,2]).det())
    ck(minor+eps**2*alpha**2/(4*H**4),'non-pp rank-two scalar Hessian minor')
    assert minor!=0
    ck(Hess.extract([0,2],[0,2]).trace()+eps*alpha**2*v/H**3,
       'scalar Hessian uu normalization')
    ck((Hess*inverse*dpsi),'null gradient congruence is geodesic')
    # Vacuum subtraction leaves a finite matter-supported Kerr--Schild increment.
    vacuum=g.subs(eps,0)
    backreaction=s.zeros(4);backreaction[0,0]=kappa*eps**2*w**2*H
    ck(g-vacuum-backreaction,'exact quadratic scalar backreaction, no hidden higher orders')
    # Direct obstruction to the tempting homogeneous cubic square-source ansatz.
    q=s.symbols('q',real=True)
    cubic=q*w*w*v
    cubicdet=c(s.diff(cubic,v,2)*s.diff(cubic,w,2)-s.diff(cubic,v,w)**2)
    ck(cubicdet+4*q*q*w*w,'real cubic square-source sign obstruction')
    # Exact shifted metric determinant and second inverse coefficient are nilpotent.
    eta=s.Matrix([[0,-1,0,0],[-1,0,0,0],[0,0,0,1],[0,0,1,0]])
    h=g-eta
    ck((eta*h)**2,'Walker inclusion nilpotence is exact')
    report={'status':'PASS','sympy_version':s.__version__,'checks':checks,
       'elapsed_seconds':time.time()-started,
       'coordinates':['u','v','z','w'],'H':str(H),'theta':str(theta),'Psi':str(psi),
       'metric':[[str(x) for x in row] for row in g.tolist()],
       'Ricci':[[str(x) for x in row] for row in Ric.tolist()],
       'scalar_Hessian':[[str(x) for x in row] for row in Hess.tolist()],
       'rank_two_minor':str(minor),'exact_scalar_backreaction':str(backreaction[0,0]),
       'domain':'H=alpha*w+beta !=0; kappa>0, epsilon!=0, alpha!=0 for non-pp witness. Real coordinates give signature (2,2). Complex analytic continuation is allowed; a Lorentzian Real completion is not constructed.',
       'interpretation':'Exact finite Einstein--massless-scalar family beyond aligned-scalar pp waves. Algebraic coordinate patch, not asymptotic scattering, positive-energy Lorentzian spacetime, a global twistor geometry or interacting BV theory.'}
    (ROOT/'results/finite-backreaction.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'status':'PASS','checks':len(checks),'elapsed_seconds':report['elapsed_seconds'],
                      'Ricci_uu':str(Ric[0,0]),'Hessian_minor':str(minor)}))

if __name__=='__main__':
    main()
