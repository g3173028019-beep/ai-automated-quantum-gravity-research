"""Symbolic replays of finite algebraic checks in the D5 manuscript.
Requires Python 3 and SymPy; do not run with -O.
"""

def check_1():
    import sympy as s
    z,h=s.symbols('z h', real=True)
    g=s.diag(-1,1,1)
    n=s.Matrix([1,0]); m=s.Matrix([z,h]); e=s.Matrix([0,1]); f=s.Matrix([-h,z])
    ks=[s.Matrix([1,*n]),s.Matrix([1,*(-n)]),s.Matrix([-1,*(-m)]),s.Matrix([-1,*m])]
    es=[s.Matrix([0,*x]) for x in [e,e,f,f]]
    red=lambda x:s.cancel(s.rem(s.fraction(s.cancel(x))[0],h*h+z*z-1,h)/s.rem(s.fraction(s.cancel(x))[1],h*h+z*z-1,h))
    dot=lambda a,b:red((a.T*g*b)[0])
    reps=[]
    for k,v in zip(ks,es):
        r=sum((dot(k,l)*l for l in ks),s.zeros(3,1)); D=dot(k,r)
        assert s.factor(D-2*(3+z*z))==0
        reps.append(v-k*dot(v,r)/D)
    gamma=z*h/(3+z*z)
    for i,j,a in [(0,1,1+2*gamma**2),(0,2,z+2*gamma*h+gamma**2*(1-z)),(0,3,z+2*gamma*h-gamma**2*(1+z))]:
        assert red(dot(reps[i],reps[j])-a)==0
    assert s.Rational(1,8)*4*sum([0,0,1,1,s.Rational(4,3)])==s.Rational(5,3)
    print('Exact four-leg projector pair contractions and wedge coefficient replayed.')
    from itertools import combinations
    V=[s.diag(1,-1,0)/s.sqrt(2),s.diag(1,1,-2)/s.sqrt(6)]
    for i,j in combinations(range(3),2):
        M=s.zeros(3,3); M[i,j]=M[j,i]=1/s.sqrt(2); V.append(M)
    avg=sum((e*e for e in V),s.zeros(3,3))/5
    assert avg==s.eye(3)/3
    approach=s.Matrix([1,0,0])
    gb=s.simplify(sum(1+2*s.trace(e**2*f**2)-2*(approach.T*(e**2+f**2)*approach)[0] for e in V for f in V)/25)
    assert gb==s.Rational(1,3)
    assert max(j for i in range(11) for j in range(11) if i+j==10)==10
    print('Exact D5 projected GB trace and ten-moment contact count replayed.')

def check_2():
    import sympy as s
    r,v,k,t,S,Theta=s.symbols('r v k t S Theta',positive=True)
    assert s.integrate(168*r*r*(1-r)**5,(r,0,1))==1
    assert s.expand(30*s.integrate((v-r)**2*(1-v)**2,(v,r,1))-(1-r)**5)==0
    H=s.integrate(r*(1-r)**2*s.sin(k*r)/k,(r,0,1))
    assert s.simplify(H-2*(k*(2+s.cos(k))-3*s.sin(k))/k**5)==0
    nu=s.symbols('nu')
    assert s.simplify(s.series(1/(S+t/2-nu)+1/(S+t/2+nu),nu,0,5).removeO().coeff(nu,4)-2/(S+t/2)**5)==0
    rho=s.sin(Theta*r/2)/s.sin(Theta/2)
    assert s.simplify(s.diff(rho,r).subs(r,0)-Theta/(2*s.sin(Theta/2)))==0
    print('Exact probability normalization, positive scale identity, radial transform, and ordinary k=4 kinematics replayed.')

def check_3():
    import sympy as s
    from math import comb
    x,h,t,mu,S,q=s.symbols('x h t mu S q')
    # Centered positive finite difference cancels lower monomials.
    for k in range(1,10):
        for d in range(2*k+1):
            C=sum((1-2*(j%2))*comb(2*k,k+j)*j**d for j in range(-k,k+1))
            assert s.sympify(C)==(0 if d<2*k else (-1)**k*s.factorial(2*k))
    # The exact D5 latitude moment identities.
    for J in range(1,12):
        U=s.chebyshevu(J,x)/(J+1)
        p1=s.diff(U,x).subs(x,1)
        p2=s.diff(U,x,2).subs(x,1)
        assert p1==s.Rational(J*(J+2),3)
        assert p2==s.Rational(J*(J-1)*(J+2)*(J+3),15)
        assert 0<=p2<=p1**2
        L=s.symbols('L')
        approx=1-L*(t**2/s.Integer(2)-t**4/s.Integer(24))
        F=s.series(U.subs(x,approx),t,0,5).removeO().expand()
        assert s.factor(-2*F.coeff(t,2)-L*p1)==0
        assert s.factor(24*F.coeff(t,4)-(L*p1+3*L*L*p2))==0
    # k=4 even-energy weight and the S^-5 to S^-3 bridge.
    nu=s.symbols('nu')
    res=s.series(1/(S+t/2-nu)+1/(S+t/2+nu),nu,0,5).removeO()
    assert s.simplify(res.coeff(nu,4)-2/(S+t/2)**5)==0
    print("Exact centered moment recovery, latitude second/fourth moments, and k=4 weight replayed.")

def check_4():
    import sympy as s
    a,h,t,x,S=s.symbols("a h t x S")
    N=4
    truncate=lambda p:s.series(p,h,0,N+1).removeO().expand()
    T=a
    for _ in range(N):
        T=truncate(a+h*T**3)
    assert truncate(T-a-h*T**3)==0
    f=(x**7+3*x**2+x**-1)*(1+2*t+t**4)
    coeff=lambda p,j:s.expand(p).coeff(x,j)
    series=truncate(f.subs(t,T))
    for n in range(1,N+1):
        lb=s.diff(a**(3*n)*s.diff(f.subs(t,a),a),a,n-1)/s.factorial(n)
        assert s.expand(series.coeff(h,n)-lb)==0
        assert coeff(s.expand(series.coeff(h,n)/x**n),1)==coeff(lb,n+1)
    
    # Fixed-transfer s/u pair resolvent and moving heavy x Jacobian.
    ss=s.symbols("ss")
    u=-ss-t
    inv=ss**2+ss*t+t*t
    XS=S*S+S*t+t*t
    assert s.simplify(1/(S-ss)+1/(S-u)-(2*S+t)/(XS-inv))==0
    XC=S**3/(S-a)
    assert s.factor(s.diff(XC,S)/XC**2-(2*S-3*a)/S**4)==0
    # Exact analytic inequalities in the chosen rescaled geometric proof
    epsilon=s.Rational(1,64)
    assert 4*(1+64*epsilon**2)**3<5
    assert 12*(1+64*epsilon**2)**2*epsilon**2<s.Rational(1,100)
    assert s.Rational(1,2)-(1+64*epsilon**2)*epsilon-(1+64*epsilon**2)**2*epsilon**2>s.Rational(2,5)
    print("Exact finite-order Lagrange series, crossing resolvent, moving Jacobian, and disc margins replayed.")

if __name__ == '__main__':
    for check in (check_1, check_2, check_3, check_4):
        check()
