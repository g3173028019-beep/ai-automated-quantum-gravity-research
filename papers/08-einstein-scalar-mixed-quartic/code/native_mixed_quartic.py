"""Native off-shell Einstein--scalar mixed quartic on the TT relative carrier.

The construction is a finite tensor formula at arbitrary off-shell momenta.
Direct Gamma--Gamma differentiation independently checks the quadratic
Einstein auxiliary source, including terms proportional to external waves.
No historical module is imported and no historical output is overwritten.
"""
from pathlib import Path
import itertools
import json
import time
import sympy as s

ROOT = Path(__file__).resolve().parents[1]
ETA = s.Matrix([[0,-1,0,0],[-1,0,0,0],[0,0,1,0],[0,0,0,1]])
CHECKS = []

def simp(x):
    return s.cancel(s.expand(x))

def check(expr, name):
    values = list(expr) if isinstance(expr, s.MatrixBase) else [expr]
    residuals = [simp(x) for x in values]
    if any(x != 0 for x in residuals):
        raise AssertionError((name, [str(x) for x in residuals if x != 0]))
    CHECKS.append(name)

def wave(p):
    return simp(-(p.T*ETA*p)[0])

def inclusion(p, f, g):
    """Metric factors f,g, before a canonical sqrt(2*kappa) insertion."""
    d,a,px,py = p
    C=s.Matrix([[f,g],[g,-f]])
    transverse=s.Matrix([px,py])
    B=C*transverse/a
    A=(transverse.T*C*transverse)[0]/a**2
    return s.Matrix([[A,0,B[0],B[1]],[0,0,0,0],
                     [B[0],0,f,g],[B[1],0,g,-f]])

def aux_basis():
    # trace coordinate is delta T = (t/2) I, not t I.
    t=s.zeros(4);t[2,2]=t[3,3]=s.Rational(1,2)
    A=s.zeros(4);A[0,0]=1
    Bx=s.zeros(4);Bx[0,2]=Bx[2,0]=1
    By=s.zeros(4);By[0,3]=By[3,0]=1
    return [t,A,Bx,By]

def lower_gamma(p,h):
    return [[[s.I*(p[i]*h[r,j]+p[j]*h[r,i]-p[r]*h[i,j])/2
              for j in range(4)] for i in range(4)] for r in range(4)]

def raise_gamma(lower,m):
    return [[[simp(sum(m[r,n]*lower[n][i][j] for n in range(4)))
              for j in range(4)] for i in range(4)] for r in range(4)]

def add_gamma(g,h):
    return [[[simp(g[r][i][j]+h[r][i][j]) for j in range(4)]
              for i in range(4)] for r in range(4)]

def mixed_ricci(p,q,h,j):
    lo1=lower_gamma(p,h);lo2=lower_gamma(q,j)
    ga1=raise_gamma(lo1,ETA);ga2=raise_gamma(lo2,ETA)
    ga12=add_gamma(raise_gamma(lo2,-ETA*h*ETA),
                   raise_gamma(lo1,-ETA*j*ETA))
    k=p+q
    return s.Matrix(4,4,lambda i,jj:simp(sum(
        s.I*k[r]*ga12[r][i][jj]-s.I*k[jj]*ga12[r][i][r]
        +sum(ga1[r][r][m]*ga2[m][i][jj]+ga2[r][r][m]*ga1[m][i][jj]
             -ga1[r][jj][m]*ga2[m][i][r]-ga2[r][jj][m]*ga1[m][i][r]
             for m in range(4)) for r in range(4))))

def einstein_density_current(p,q,h,j):
    """Coefficient [epsilon1 epsilon2] sqrt|g| G^{mu nu}.

    The two inclusions are trace free and transverse even off shell.
    The wave-dependent raising/density term MUST NOT be dropped off shell.
    """
    ric=mixed_ricci(p,q,h,j)
    tr=s.trace(ETA*ric)
    product=ETA*h*ETA*j*ETA+ETA*j*ETA*h*ETA
    trproduct=s.trace(ETA*h*ETA*j)
    correction=(wave(p)+wave(q))/2*(product-ETA*trproduct/2)
    density=ETA*ric*ETA-ETA*tr/2+correction
    return density.applyfunc(simp),correction.applyfunc(simp)

def gravitational_aux_source(p,q,h,j,kappa=s.Integer(1)):
    dens,corr=einstein_density_current(p,q,h,j)
    out=s.Matrix([-sum(dens[i,jj]*y[i,jj] for i in range(4)
                       for jj in range(4))/(2*kappa) for y in aux_basis()])
    return out.applyfunc(simp),dens,corr

def gamma_gamma_cubic(pp,hh,kappa=s.Integer(1)):
    """Independent coefficient of epsilon1 epsilon2 epsilon3 in S_EH.

    Integrated Gamma--Gamma representative, conserving all three momenta.
    This path does not use the Ricci/Euler source formula above.
    """
    low=[lower_gamma(p,h) for p,h in zip(pp,hh)]
    first=[raise_gamma(l,ETA) for l in low]
    second={}
    for i,j in itertools.combinations(range(3),2):
        second[i,j]=add_gamma(raise_gamma(low[j],-ETA*hh[i]*ETA),
                              raise_gamma(low[i],-ETA*hh[j]*ETA))
    def contract(m,g,h):
        return simp(sum(m[i,j]*(g[r][i][t]*h[t][j][r]
                                 -g[r][i][j]*h[t][r][t])
                        for i in range(4) for j in range(4)
                        for r in range(4) for t in range(4)))
    result=0
    for i in range(3):
        j,k=[n for n in range(3) if n!=i]
        second_jk=second[min(j,k),max(j,k)]
        dens1=s.trace(ETA*hh[i])*ETA/2-ETA*hh[i]*ETA
        result+=contract(ETA,first[i],second_jk)
        result+=contract(ETA,second_jk,first[i])
        result+=contract(dens1,first[j],first[k])
        result+=contract(dens1,first[k],first[j])
    return simp(result/(2*kappa))

def scalar_aux_source(p,q):
    d,a,px,py=p;e,b,qx,qy=q
    return s.Matrix([-(d*b+a*e)/2,-a*b,a*qx+b*px,a*qy+b*py])

def auxiliary_inverse(k,kappa=s.Integer(1)):
    d,a,px,py=k
    return kappa*s.Matrix([[0,-4/a**2,0,0],
        [-4/a**2,4*d/a**3+2*(px**2+py**2)/a**4,2*px/a**3,2*py/a**3],
        [0,2*px/a**3,2/a**2,0],[0,2*py/a**3,0,2/a**2]])

def gamma_gamma_quadratic(p,h,j,kappa=s.Integer(1)):
    g=raise_gamma(lower_gamma(p,h),ETA)
    q=raise_gamma(lower_gamma(-p,j),ETA)
    def contract(g,q):
        return sum(ETA[i,j]*(g[r][i][t]*q[t][j][r]-g[r][i][j]*q[t][r][t])
                   for i in range(4) for j in range(4)
                   for r in range(4) for t in range(4))
    return simp((contract(g,q)+contract(q,g))/(2*kappa))

def native_matter_contact(p,q,h,j):
    D2=ETA*h*ETA*j*ETA+ETA*j*ETA*h*ETA-ETA*s.trace(ETA*h*ETA*j)/2
    return simp((p.T*D2*q)[0])

def contact_from_local_density(p,q,h,j):
    # Polarized coefficient of the four explicit local terms in DERIVATION.md.
    C1=h[2:4,2:4];C2=j[2:4,2:4]
    B1=h[0,2:4].T;B2=j[0,2:4].T
    rho12=s.trace(C1*C2)  # [epsilon1 epsilon2] (f^2+g^2)
    dp,ap,px,py=p;dq,aq,qx,qy=q
    P=s.Matrix([px,py]);Q=s.Matrix([qx,qy])
    v1=rho12*(dp*aq+dq*ap)/2
    v2=-(ap*((B1.T*C2+B2.T*C1)*Q)[0]
          +aq*((B1.T*C2+B2.T*C1)*P)[0])
    v3=rho12*(P.T*Q)[0]/2
    v4=2*(B1.T*B2)[0]*ap*aq
    return simp(v1+v2+v3+v4)

def main():
    started=time.time()
    # Exact generic matrix identity for the native density coefficient.
    f,g,ep=s.symbols('f g ep')
    C=s.Matrix([[f,g],[g,-f]])
    check(C*C-(f*f+g*g)*s.eye(2),'native transverse inverse algebra')
    check(C.det()+f*f+g*g,'native transverse determinant')
    a,b,d,e,f1,g1,f2,g2=s.symbols('a b d e f1 g1 f2 g2',nonzero=True)
    sp=s.Matrix([d,a,0,0]);sq=s.Matrix([e,b,0,0])
    sh=inclusion(sp,f1,g1);sj=inclusion(sq,f2,g2)
    closedJA=(a*a+b*b+a*b)*(f1*f2+g1*g2)/2
    check(gamma_gamma_cubic([sp,sq,-sp-sq],[sh,sj,aux_basis()[1]])-closedJA,
          'generic longitudinal-frequency Einstein A source from direct action')
    jfull,_,_=gravitational_aux_source(sp,sq,sh,sj)
    check(jfull[1]-closedJA,'generic A source agrees with closed null constraint')
    examples=[]
    for seed in range(1,5):
        p=s.Matrix([seed+1,seed+2,2-seed,seed+3])
        q=s.Matrix([-seed-3,2*seed+1,seed+4,1-seed])
        r=s.Matrix([2*seed-1,-seed-2,seed+1,2*seed+3])
        t=-p-q-r;k=p+q
        h=inclusion(p,s.Rational(seed+1,seed+2),s.Rational(seed-2,seed+3))
        j=inclusion(q,s.Rational(2-seed,seed+4),s.Rational(seed+2,seed+1))
        check(s.trace(ETA*h),f'off-shell trace first {seed}')
        check(p.T*ETA*h,f'off-shell transversality first {seed}')
        check(q.T*ETA*j,f'off-shell transversality second {seed}')
        J,dens,corr=gravitational_aux_source(p,q,h,j)
        for index,y in enumerate(aux_basis()):
            gg=gamma_gamma_cubic([p,q,-k],[h,j,y])
            check(gg-J[index],f'actual off-shell EH source vs direct action {seed}:{index}')
        native=native_matter_contact(r,t,h,j)
        check(native-contact_from_local_density(r,t,h,j),
              f'native mixed quartic matrix vs four-term density {seed}')
        matter=scalar_aux_source(r,t)
        if seed==1:
            basis=aux_basis()
            Hactual=s.Matrix(4,4,lambda i,j:gamma_gamma_quadratic(k,basis[i],basis[j]))
            check(Hactual*auxiliary_inverse(k)-s.eye(4),
                  'retained inverse matches actual EH Hessian in declared trace/shift coordinates')
            for index,ybasis in enumerate(basis):
                check(gamma_gamma_quadratic(k,inclusion(k,2,3),ybasis),
                      f'physical inclusion orthogonal to auxiliary Hessian {index}')
        feedback=simp(-(J.T*auxiliary_inverse(k)*matter)[0])
        wrong=-(corr[2,2]+corr[3,3])/4
        assert wrong!=0 or any(corr)
        examples.append({'seed':seed,'external_waves':[str(wave(x)) for x in [p,q,r,t]],
            'momenta':[list(map(str,x)) for x in [p,q,r,t]],
            'metric_polarizations':[[str(h[2,2]),str(h[2,3])],[str(j[2,2]),str(j[2,3])]],
            'Einstein_auxiliary_source':list(map(str,J)),
            'omitted_offshell_density_term':[[str(x) for x in row] for row in corr.tolist()],
            'native_contact':str(native),'auxiliary_backreaction':str(feedback),
            'completed_mixed_contact':str(simp(native+feedback))})
    # Independent completion-of-square normalization, with arbitrary sources.
    H=s.Matrix([[2,1,0,0],[1,3,0,0],[0,0,4,1],[0,0,1,5]])
    jg=s.Matrix(s.symbols('g0:4'));jm=s.Matrix(s.symbols('m0:4'))
    x,z=s.symbols('x z');y=-H.inv()*(x*jg+z*jm)
    action=(y.T*H*y)[0]/2+(y.T*(x*jg+z*jm))[0]
    check(s.expand(action).coeff(x,1).coeff(z,1)+(jg.T*H.inv()*jm)[0],
          'mixed Schur complement has no extra factor one half')
    report={'status':'PASS','sympy_version':s.__version__,'checks':CHECKS,
            'elapsed_seconds':time.time()-started,'examples':examples,
            'closed_A_source':'(a1^2+a2^2+a1*a2)*trace(C1*C2)/(4*kappa)',
            'closed_trace_response':'t12_CC=(a1^2+a2^2+a1*a2)*trace(C1*C2)/(a1+a2)^2',
            'result':'Off-shell native C^2 phi^2 plus actual Einstein auxiliary backreaction',
            'scope':'Flat TT light-front inclusion, a1*a2*(a1+a2)!=0; arbitrary off-shell four momenta summing to zero. Metric factors explicit; kappa=1 in examples. Not full C^4, nonlinear chart descent, BV or a measure.'}
    (ROOT/'results/native-mixed-quartic.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'status':'PASS','checks':len(CHECKS),'examples':len(examples),
                      'elapsed_seconds':report['elapsed_seconds']}))

if __name__=='__main__':
    main()
