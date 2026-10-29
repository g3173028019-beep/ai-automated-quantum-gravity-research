"""Portable default: exact derivations, geometry, and frozen numerical evidence.

This does NOT recompute the complete Regge action Hessian. See run_full.py.
"""
from pathlib import Path
import hashlib,json,itertools,math
import numpy as np
import sympy as s
HERE=Path(__file__).resolve().parent
def decode(x):
    if isinstance(x,dict) and set(x)=={'real','imag'}: return complex(x['real'],x['imag'])
    if isinstance(x,list): return [decode(v) for v in x]
    return x
def geometry(m):
    q=np.array([[0,0,0,-1.5709397991101974],[.3762757793700152,0,0,-.42646145218382514],[.7702937141192545,0,0,1.1723475896580613],[-.7702937141192545,0,0,1.1723475896580613],[-.3762757793700152,0,0,-.42646145218382514]])
    qs=[q[0]]+[(1-j/m)*q[0]+j/m*q[1] for j in range(1,m)]+list(q[1:])+[(1-j/m)*q[4]+j/m*q[0] for j in range(1,m)]
    angles=np.arange(6)*math.pi/3
    rs=np.column_stack([np.zeros(6),np.cos(angles),np.sin(angles),np.zeros(6)])
    points=np.vstack([np.zeros((1,4)),qs,rs]); nq=len(qs)
    simplices=[(0,i+1,(i+1)%nq+1,nq+j+1,nq+(j+1)%6+1) for i in range(nq) for j in range(6)]
    ini=set(range(2,m+2)); fin={m+4}|set(range(m+5,2*m+4))
    return points,simplices,ini,fin,nq

def main():
    provenance=json.loads((HERE/'PROVENANCE.json').read_text())
    for row in provenance['frozen_files']:
        assert hashlib.sha256((HERE/row['path']).read_bytes()).hexdigest()==row['sha256'],row['path']
    O=s.Matrix([[0,1],[-1,0]])
    a,b,c=s.symbols('a b c',real=True,nonzero=True)
    M=s.Matrix([[-a/b,-1/b],[b-c*a/b,-c/b]])
    assert s.simplify(M.T*O*M-O)==s.zeros(2)
    lam=s.symbols('lambda')
    assert s.factor(M.charpoly(lam).as_expr()-(lam**2+(a+c)*lam/b+1))==0
    rot=s.Matrix([[-s.Rational(1,2),-s.sqrt(3)/2],[s.sqrt(3)/2,-s.Rational(1,2)]])
    refl=s.diag(1,-1); X=s.Matrix(2,2,s.symbols('x0:4'))
    cyclic=s.linsolve(list(rot.T*X*rot-X),list(X))
    dihedral=s.linsolve(list(rot.T*X*rot-X)+list(refl*X*refl-X),list(X))
    assert s.simplify(rot**6)==s.eye(2)
    shear=s.Matrix([[1,0],[-3,1]]); base=s.Matrix([[0,-1],[1,0]])
    changed=shear*base
    assert shear.T*O*shear==O and changed.trace()==3 and base.trace()==0
    assert shear.inv()*changed==base
    chi=s.symbols('chi',real=True)
    elliptic=s.Matrix([[chi,-1/b],[b*(1-chi**2),chi]])
    G=s.diag(b*b*(1-chi**2),1)
    assert s.simplify(elliptic.T*G*elliptic-G)==s.zeros(2)
    u=s.Matrix(s.symbols('u0:3')); v=s.Matrix(s.symbols('v0:3')); k=s.symbols('k',nonzero=True)
    assert (-u*v.T/k).rank()==1
    eta=np.diag([-1,1,1,1]); geoms=[]
    for level in range(1,6):
        pts,simplices,initial,final,nq=geometry(level)
        eig=[]; counts={'past_spacelike':0,'future_spacelike':0,'timelike_side':0}
        for simplex in simplices:
            assert not (set(simplex)&initial and set(simplex)&final)
            tang=pts[list(simplex[1:])]-pts[0]
            vals=np.linalg.eigvalsh(tang@eta@tang.T); eig.append(vals)
            assert vals[0]<0<vals[1]
            face=pts[list(simplex[1:])]; tangent=face[1:]-face[0]
            n=np.linalg.svd(tangent@eta)[2][-1]
            if n@eta@(pts[0]-face[0])>0:n=-n
            norm=n@eta@n
            key=('future_spacelike' if n[0]>0 else 'past_spacelike') if norm<0 else 'timelike_side'
            counts[key]+=1
        assert counts=={'past_spacelike':6*(level+1),'future_spacelike':6*(level+1),'timelike_side':6}
        assert abs(np.linalg.norm(pts[nq+1]-pts[nq+2])-1)<1e-14
        geoms.append({'level':level,'simplices':len(simplices),'counts':counts,'negative_gram_margin':float(min(abs(x[0]) for x in eig)),'positive_gram_margin':float(min(x[1] for x in eig)),'fixed_transverse_edge':1})
    saved=[]; repair=[]; target=.8150346718577313; target_b=.840003
    for level in [1,5]:
        raw=json.loads((HERE/f'data/level-{level}-reproduction.json').read_text())
        A=np.array(decode(raw['surviving_initial_self'])); B=np.array(decode(raw['surviving_transfer'])); C=np.array(decode(raw['surviving_final_self']))
        inv=np.linalg.inv(B); mat=np.block([[-inv@A,-inv],[B.T-C@inv@A,-C@inv]])
        omega=np.block([[np.zeros((2,2)),np.eye(2)],[-np.eye(2),np.zeros((2,2))]])
        eig=np.linalg.eigvals(mat); frozen=np.array(decode(raw['surviving_canonical_eigenvalues']))
        distance=max(min(abs(x-y) for y in frozen) for x in eig)
        assert distance<1e-10
        assert np.max(abs(mat.T@omega@mat-omega))<1e-12
        assert raw['spin_transfer_rank']==2
        if level==1: assert np.max(abs(abs(eig)-1))<1e-5
        else:
            assert sum(x.real>1.25 for x in eig)==2 and sum(.75<x.real<.8 for x in eig)==2
            assert max(abs(x.imag) for x in eig)<2e-6
        av=float(np.trace(A+C).real/4); bv=float(np.trace(B).real/2)
        delta=-bv*math.cos(target)-av; r=math.sqrt(target_b/bv)
        assert abs(r*r*bv-target_b)<1e-14
        assert abs(r*r*(av+delta)+target_b*math.cos(target))<1e-14
        repair.append({'level':level,'a':av,'b':bv,'chi':-av/bv,'delta_a':delta,'dilation':r})
        saved.append({'level':level,'eigenvalues':[[z.real,z.imag] for z in eig],'symplectic_residual':float(np.max(abs(mat.T@omega@mat-omega))),'frozen_spectrum_distance':float(distance),'original_full_producer_seconds':raw['elapsed_seconds']})
    result={'status':'pass','mode':'exact algebra + independent geometry + verification of frozen numerical evidence, not full Hessian replay','frozen_hashes_verified':len(provenance['frozen_files']),'commutant_C6':str(cyclic),'commutant_D6':str(dihedral),'endpoint_counterexample':{'old_trace':0,'new_trace':3},'structural_rank_bound':'at most two from unique internal E2 pair and zero direct cap cross support','geometry':geoms,'saved_numerics':saved,'real_scalar_calibrations':repair,'human_expert_review':False}
    (HERE/'results').mkdir(exist_ok=True)
    (HERE/'results/replay.json').write_text(json.dumps(result,indent=2))
    print(json.dumps({'status':'pass','frozen_hashes_verified':len(provenance['frozen_files']),'levels_checked':[1,5],'full_Hessian_recomputed':False}))
if __name__=='__main__':main()
