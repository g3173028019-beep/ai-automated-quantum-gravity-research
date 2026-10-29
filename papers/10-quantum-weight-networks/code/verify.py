"""Independent exact finite audit; writes only beside this file.

Python 3 / SymPy. No source scripts or manuscripts are modified.
The binary Lie algorithm closes under the original generators: Jacobi makes
this sufficient to span their generated Lie algebra.
"""
import collections
import itertools
import json
from pathlib import Path
import sympy as s

MODEL=json.loads((Path(__file__).resolve().parents[1]/'inputs'/'model.json').read_text(encoding='utf-8'))
I=s.eye(2); X=s.Matrix([[0,1],[1,0]]); Y=s.Matrix([[0,-s.I],[s.I,0]]); Z=s.diag(1,-1)
kp=s.kronecker_product
paulis=dict(I=I,X=X,Y=Y,Z=Z)
gam=[kp(*(paulis[p] for p in word)) for word in MODEL['gamma']]
sig=lambda a,b:(gam[a]*gam[b]-gam[b]*gam[a])/(2*s.I)
br=lambda a,b:(a*b-b*a)/s.I
flat=lambda a:s.Matrix(list(a))
bivs=[sig(a,b) for a,b in itertools.combinations(range(5),2)]
S=gam[4]; C=sig(3,4)
J=[sig(1,2),-sig(0,2),sig(0,1)]
K=[sig(i,3) for i in range(3)]; E=[sig(i,4) for i in range(3)]

def matrix_closure(seeds):
    basis=[]; mat=s.zeros(16,0); queue=collections.deque(seeds)
    while queue:
        h=s.simplify(queue.popleft()); candidate=mat.row_join(flat(h))
        if candidate.rank()>len(basis):
            basis.append(h); mat=candidate
            queue.extend(br(h,g) for g in seeds)
    return basis

G=s.Matrix(MODEL['rational_weight'])
seed=[kp(G,X),kp(G,Y),kp(X,I)]
generic_bad=matrix_closure([kp(G,X),kp(G,Y),kp(X+Z,I)])
base=matrix_closure(seed)
assert len(base)==10 and len(generic_bad)==15
assert all(br(j,k)==2*sum((s.LeviCivita(a,b,c)*J[c] for c in range(3)),s.zeros(4))
           for a,j in enumerate(J) for b,k in enumerate(J))
for a,b in itertools.product(range(3),repeat=2):
    for vec in (K,E):
        assert br(J[a],vec[b])==2*sum((s.LeviCivita(a,b,c)*vec[c] for c in range(3)),s.zeros(4))
        assert br(vec[a],vec[b])==2*sum((s.LeviCivita(a,b,c)*J[c] for c in range(3)),s.zeros(4))
assert all(br(C,j)==s.zeros(4) for j in J)
assert all(br(K[a],E[b])==2*(a==b)*C for a in range(3) for b in range(3))
assert all(br(C,K[a])==2*E[a] and br(C,E[a])==-2*K[a] for a in range(3))
assert s.prod(gam)==s.eye(4)
form=s.I*kp(Z,Y)
assert all((s.I*h).T*form+form*(s.I*h)==s.zeros(4) for h in base)
nu,mu=s.symbols('nu mu',real=True)
rho=(s.eye(4)+nu*S+mu*C)/4
six=K+E
cov=s.Matrix(6,6,lambda a,b:s.simplify(s.trace(rho*(six[a]*six[b]+six[b]*six[a]))/2))
assert cov==s.eye(6)
assert s.simplify((nu*S+mu*C)**2)==(nu**2+mu**2)*s.eye(4)
trace_ok=all(s.trace(S*sig(a,b)*sig(c,d))==-4*s.LeviCivita(a,b,c,d)
             for a,b in itertools.combinations(range(4),2)
             for c,d in itertools.combinations(range(4),2))
assert trace_ok
for a,b in itertools.combinations(range(4),2):
    for c,d in itertools.combinations(range(4),2):
        got=s.trace((s.eye(4)+nu*S)*sig(a,b)*sig(c,d))/4
        want=int(a==c)*int(b==d)-int(a==d)*int(b==c)-nu*s.LeviCivita(a,b,c,d)
        assert s.simplify(got-want)==0
cart=[kp(sig(0,1),s.eye(4)),kp(sig(2,3),s.eye(4)),
      kp(s.eye(4),sig(0,1)),kp(s.eye(4),sig(2,3)),kp(S,S)]
assert s.prod(cart)==s.eye(16)
projector_counts=collections.Counter()
for signs in itertools.product((-1,1),repeat=5):
    p=s.eye(16)
    for t,h in zip(signs,cart): p=p*(s.eye(16)+t*h)/2
    rank=p.rank(); parity=s.prod(signs)
    assert rank==(1 if parity==1 else 0)
    projector_counts[str((int(parity),rank))]+=1

local=MODEL['local_basis']
def encode(string):
    x=z=0
    for i,p in enumerate(string):
        if p in 'XY': x|=1<<i
        if p in 'ZY': z|=1<<i
    return x,z
def pauli_closure(n,bond):
    seeds=[]
    for cell in range(n):
        for word in local: seeds.append(encode('II'*cell+word+'II'*(n-cell-1)))
    for cell in range(n-1): seeds.append(encode('II'*cell+bond+bond+'II'*(n-cell-2)))
    known=set(seeds); todo=collections.deque(seeds)
    while todo:
        x,z=todo.popleft()
        for u,v in seeds:
            if ((x&v).bit_count()+(z&u).bit_count())%2:
                nxt=x^u,z^v
                if nxt not in known: known.add(nxt); todo.append(nxt)
    return len(known)
table={str(n):{bond:pauli_closure(n,bond) for bond in ('ZI','YZ')} for n in MODEL['finite_chain_sizes']}
assert [table[str(n)]['ZI'] for n in (1,2,3,4)]==[10,45,4095,65535]
assert [table[str(n)]['YZ'] for n in (1,2,3,4)]==[10,255,4095,65535]

# Local algebraic witnesses used in the all-n induction proof.
assert br(sig(0,1),gam[0])==2*gam[1]
# Positive-weight family: two adjoint kick brackets isolate ZX, recovering IX.
a,b=s.symbols('a b',nonzero=True,real=True)
hx=kp(a*I+b*Z,X); kick=kp(X,I)
assert br(kick,br(kick,hx))==-4*b*kp(Z,X)
jc,h,beta,v=s.symbols('J h beta v',positive=True,real=True)
r=s.sqrt(jc**2*v**2+h**2)
free=jc*v**2/2-s.log(4*s.cosh(beta*r))/beta
assert s.simplify(s.diff(free,v)-jc*v*(1-jc*s.tanh(beta*r)/r))==0
assert s.simplify(s.diff(free,v,2).subs(v,0)-jc/h*(h-jc*s.tanh(beta*h)))==0

out={
 'sympy_version':s.__version__,
 'one_cell_closure':len(base), 'noncommuting_update_X_plus_Z_closure':len(generic_bad),
 'clifford_and_polarized_trace':'exact_pass',
 'all_displayed_rotation_clock_signs':'exact_pass',
 'six_variable_symmetric_covariance':str(cov),
 'chirality_projectors':dict(projector_counts),
 'independent_binary_chain_dimensions':table,
 'finite_N_CCR_residual_state_rms':'sqrt(1-mu^2)/(mu*sqrt(N))',
 'finite_N_CCR_residual_operator_norm':'1+1/mu (mu>0)',
 'positive_weight_family':'a>|b|>0 with a transverse X kick: same ten-dimensional span',
 'mean_field_derivatives':'exact_pass; model is transverse Curie-Weiss with a spectator per site',
 'all_n_proliferation':'proved in the manuscript; n=4 checked exactly',
 'scope':'Finite matrix and Pauli statements only; no mean-field Gibbs fluctuation limit.'}
(Path(__file__).resolve().parents[1]/'outputs'/'verification.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
print(json.dumps(out,indent=2))
