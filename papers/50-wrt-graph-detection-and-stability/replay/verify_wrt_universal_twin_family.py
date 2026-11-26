"""Exact count replay of the all-k WRT universal-twin family; no floating point."""
from collections import defaultdict
from functools import lru_cache
from math import factorial
from pathlib import Path
from subprocess import check_output,run
from tempfile import TemporaryDirectory

ROWS=((3,7,11,12),(3,5,14,14))

def add(*v):
    return tuple(map(sum,zip(*v))) if v else (0,)*5
def suffix(cols,A=15,B=15):
    """V_j = # old extensions with j vertices after the last surviving A."""
    full=A|(B<<4)
    def candidates(s):
        m=full^s
        while m:
            z=m&-m;m^=z;i=z.bit_length()-1
            if i<4 or not(cols[i-4]&A&~s):
                yield z
    @lru_cache(None)
    def f(s):
        if s==full:return 1
        return sum(f(s|z) for z in candidates(s))
    D=[0]*(1<<8);D[0]=1
    out=[0]*5
    for s,x in enumerate(D):
        if not x:continue
        for z in candidates(s):
            t=s|z;D[t]+=x
            if s&A!=A and t&A==A:
                j=(full^t).bit_count()
                out[j]+=x*f(t)
    return tuple(out)

def lower(rows):
    cols=tuple(sum(1<<i for i,r in enumerate(rows) if r>>j&1)for j in range(4))
    lead=suffix(cols)
    old,new,delet=[],[],[]
    for i,r in enumerate(rows):
        A=15^(1<<i)
        cc=tuple(c|A if r>>j&1 else c for j,c in enumerate(cols))
        new.append(suffix(cc,A))
        for j,c in enumerate(cols):
            if r>>j&1:
                cc=tuple(x|(c^(1<<i)) if r>>l&1 and l!=j else x
                         for l,x in enumerate(cols))
                old.append(suffix(cc,A,15^(1<<j)))
    for j,c in enumerate(cols):
        if c.bit_count()==2:
            delet.append(suffix(cols,15^c,15^(1<<j)))
    wo,wn=add(*old),add(*new)
    # k(k+j-1)!/j! = (k+j)!/j! - (k+j-1)!/(j-1)!; j=0 has wn_0=0.
    assert wn[0]==0
    both=tuple(wo[j]+wn[j]-(wn[j+1] if j<4 else 0) for j in range(5))
    return lead,both,add(*delet),wo,wn

def ints(s):
    return {(int(a),int(k),p):int(c) for a,k,p,c in map(str.split,s.splitlines())}
def lincomb(*cs):
    K={k for _,d in cs for k in d}
    return {k:v for k in K if (v:=sum(c*d.get(k,0) for c,d in cs))}

def main():
    H,G=map(lower,ROWS)
    assert H[:3]==G[:3]
    print('common lower data (V,W,D):',H[:3])
    print('left W ingredients',H[3:],'right W ingredients',G[3:])
    loc=Path(__file__).parent
    with TemporaryDirectory() as t:
        exe=str(Path(t)/'direct');tw=str(Path(t)/'twins')
        def build(src,dst):
            run(['g++','-O3','-std=c++17',str(loc/src),'-o',dst],check=True)
        build('verify_wrt_ordered_gap_signature.cpp',exe)
        build('verify_wrt_twin_signature.cpp',tw)
        def hist(a,b,rows):
            return ints(check_output([exe,str(a),str(b),*map(str,rows)],text=True))
        Q=lincomb((1,hist(2,2,(3,3))),(-2,hist(2,2,(1,3))),(1,hist(2,2,(1,2))))
        Q={p:c for p,c in Q.items()if p[0]==0}
        assert Q and all(k==3 for a,k,w in Q)
        P=defaultdict(dict)
        for d,j,a,k,w,c in map(str.split,check_output([tw],text=True).splitlines()):
            P[int(d),int(j)][int(a),int(k),w]=int(c)
        assert dict(P)=={(0,1):lincomb((-2,Q)),(0,2):lincomb((2,Q))}
        print('symbolic untouched-twin suffix j=1: −2Q; j=2: +2Q; all other prototype sectors vanish')
        for t in (1,2,3):
            m=((1<<t)-1)<<4
            rows=tuple(r|m for R in ROWS for r in R)
            assert hist(4,4+t,rows)==lincomb((factorial(t)*t*(t+1),Q))
        print('fourth signed-gap equality = k! k(k+1) times C4−2P4+matching, every k>=1')
if __name__=='__main__':main()
