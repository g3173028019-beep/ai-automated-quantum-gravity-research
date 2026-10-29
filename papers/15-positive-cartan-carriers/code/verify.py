"""Bounded independent exact checks; writes only beside this script.

No import of the historical producer. Modular rank lower bounds plus an exact
nonzero null vector prove the rational invariant rank, without SVD tolerance.
"""
from pathlib import Path
from itertools import combinations, combinations_with_replacement, product
import json
import sys
import numpy as np
import sympy as S

OUT = Path(__file__).resolve().parents[1] / "data"

def eps(t):
    if len(set(t)) != len(t):
        return 0
    return (-1) ** sum(t[i] > t[j] for i in range(len(t)) for j in range(i+1, len(t)))

def rank_mod(a, p=1009):
    a = np.array(a, dtype=np.int64) % p
    row = 0
    for col in range(a.shape[1]):
        hits = np.flatnonzero(a[row:, col])
        if not len(hits):
            continue
        k = row + int(hits[0])
        a[[row, k]] = a[[k, row]]
        a[row] = a[row] * pow(int(a[row,col]), -1, p) % p
        idx = np.flatnonzero(a[row+1:,col]) + row+1
        a[idx] = (a[idx] - a[idx,col,None] * a[row]) % p
        row += 1
        if row == min(a.shape):
            break
    return row

def invariant_count(sign):
    pairs = list(combinations(range(5),2))
    syms = list(combinations_with_replacement(range(10),2))
    si = {p:i for i,p in enumerate(syms)}
    eta = [1,-1,-1,-1,sign]
    blocks = []
    for a,b in pairs:
        j = np.zeros((5,5),dtype=np.int64)
        j[a,b] = eta[b]
        j[b,a] = -eta[a]
        wedge = np.zeros((10,10),dtype=np.int64)
        for col,(c,d) in enumerate(pairs):
            for z in range(5):
                for term,coef in [((z,d),j[z,c]),((c,z),j[z,d])]:
                    if term[0] == term[1] or coef == 0:
                        continue
                    key = tuple(sorted(term))
                    wedge[pairs.index(key),col] += coef * (1 if term[0]<term[1] else -1)
        sym = np.zeros((55,55),dtype=np.int64)
        for col,(c,d) in enumerate(syms):
            for z in range(10):
                sym[si[tuple(sorted((z,d)))],col] += wedge[z,c]
                sym[si[tuple(sorted((c,z)))],col] += wedge[z,d]
        blocks.append(np.kron(j,np.eye(55,dtype=np.int64))+np.kron(np.eye(5,dtype=np.int64),sym))
    mat = np.vstack(blocks)
    v = np.array([eps((a,*pairs[c],*pairs[d])) for a in range(5) for c,d in syms],dtype=np.int64)
    null = bool(np.all(mat @ v == 0))
    assert null and np.any(v)
    lower = rank_mod(mat)
    assert lower == 274
    return dict(shape=list(mat.shape), modulus=1009, modular_rank=lower,
                exact_integer_null_vector_nonzero=True, exact_integer_null_residual=0,
                rational_rank=274, rational_nullity=1)

def variational_ranks():
    ps = list(combinations(range(4),2))
    ts = list(combinations(range(4),3))
    cols = list(product(range(4),range(6)))
    def tc(col,k,a,b):
        if a==b or col[0]!=k: return 0
        return (1 if a<b else -1)*int(col[1]==ps.index(tuple(sorted((a,b)))))
    tors = S.Matrix([[sum(eps((i,j,k,l))*sum(int(k==abc[z])*tc(col,l,abc[(z+1)%3],abc[(z+2)%3]) for z in range(3)) for k,l in product(range(4),repeat=2)) for col in cols] for (i,j),abc in product(ps,ts)])
    symps=list(combinations_with_replacement(range(6),2))
    bases=[]
    for i,j in symps:
        b=S.zeros(6); b[i,j]=1; b[j,i]=1; bases.append(b)
    bianchi=S.Matrix([[sum(eps((*ps[i],*ps[j]))*b[i,j] for i,j in product(range(6),repeat=2)) for b in bases]])
    ker=S.Matrix.hstack(*bianchi.nullspace())
    def rc(b,i,j,a,c):
        if i==j or a==c: return 0
        return (1 if i<j else -1)*(1 if a<c else -1)*b[ps.index(tuple(sorted((i,j)))),ps.index(tuple(sorted((a,c))))]
    emap=S.Matrix([[sum(eps((i,j,k,l))*sum(int(j==abc[z])*rc(b,k,l,abc[(z+1)%3],abc[(z+2)%3]) for z in range(3)) for j,k,l in product(range(4),repeat=3)) for b in bases] for i,abc in product(range(4),ts)])
    final=emap*ker
    assert tors.rank()==24 and tors.det()==16 and final.rank()==10
    return dict(torsion_rank=tors.rank(),torsion_determinant=int(tors.det()),
                bianchi_rank=bianchi.rank(),riemann_dimension=ker.cols,
                coframe_rank=final.rank(),weyl_kernel=ker.cols-final.rank(),arithmetic='SymPy rational exact')

def clifford_contact():
    one=S.eye(2); x=S.Matrix([[0,1],[1,0]]); y=S.Matrix([[0,-S.I],[S.I,0]]); z=S.diag(1,-1)
    eu=[S.kronecker_product(x,one),S.kronecker_product(y,one),S.kronecker_product(z,x),S.kronecker_product(z,y),S.kronecker_product(z,z)]
    out={}
    for sign in [-1,1]:
        gam=[eu[0],S.I*eu[1],S.I*eu[2],S.I*eu[3],(1 if sign==1 else S.I)*eu[4]]
        eta=[1,-1,-1,-1,sign]
        assert all(gam[a]*gam[b]+gam[b]*gam[a]==2*int(a==b)*eta[a]*S.eye(4) for a,b in product(range(5),repeat=2))
        cs=gam[0] if sign==-1 else S.I*gam[0]*gam[4]
        q=[(1 if sign==-1 else S.I)*g*cs for g in gam]
        spin={(a,b):(gam[a]*gam[b]-gam[b]*gam[a])/4 for a,b in combinations(range(5),2)}
        assert all(v==v.H for v in q)
        def gen(a,b):
            return S.zeros(4) if a==b else (spin[a,b] if a<b else -spin[b,a])
        for (a,b),(c,d) in product(spin,repeat=2):
            assert gen(a,b)*gen(c,d)-gen(c,d)*gen(a,b)==int(b==c)*eta[b]*gen(a,d)-int(a==c)*eta[a]*gen(b,d)-int(b==d)*eta[b]*gen(a,c)+int(a==d)*eta[a]*gen(b,c)
        for e in range(5):
            for (a,b),(c,d) in product(spin,repeat=2):
                assert S.trace(gam[e]*gen(a,b)*gen(c,d))==(-1 if sign==-1 else S.I)*eps((e,a,b,c,d))
        v=S.symbols('v0:5',real=True)
        contact=sum((v[a]*q[a] for a in range(5)),S.zeros(4))
        assert S.expand(contact.det()-(sum(eta[a]*v[a]**2 for a in range(5)))**2)==0
        for (b,c),v in spin.items():
            for a in range(5):
                assert v*q[a]+q[a]*v.H==int(c==a)*eta[c]*q[b]-int(b==a)*eta[b]*q[c]
        eig=q[4].eigenvals()
        assert eig=={-1:2,1:2}
        fixed=S.Matrix.vstack(*[S.kronecker_product(S.eye(4),v)+S.kronecker_product(v.conjugate(),S.eye(4)) for (a,b),v in spin.items() if b<4])
        # Lorentz-fixed Hermitian solutions X=C gamma0 +... may be signed;
        # rank is recorded without asserting all PSD boundary claims by numerics.
        out[str(sign)]=dict(contact_eigenvalues={'-1':2,'1':2},clifford_exact=True,all_spin_brackets_exact=True,all_epsilon_trace_components_exact=True,contact_determinant_identity_exact=True,all_contact_covariance_exact=True,lorentz_fixed_complex_matrix_dimension=16-fixed.rank())
    return out

if __name__=='__main__':
    result=dict(python=sys.version.split()[0],sympy=S.__version__,numpy=np.__version__,
                invariant={str(s):invariant_count(s) for s in [-1,1]},
                variational=variational_ranks(),clifford=clifford_contact())
    (OUT/'results.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result,indent=2))
