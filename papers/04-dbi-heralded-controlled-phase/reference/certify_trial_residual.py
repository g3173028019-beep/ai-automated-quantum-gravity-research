"""Outward interval integration of actual continuum trial residuals.

The quadrature used to construct the trial plays NO role here. Each
trial's residual is integrated over the full half line with interval
Taylor cells and an analytic gamma-polynomial tail. Execution is bounded
by an internal 40-second slice; an incomplete slice saves a resumable
partition rather than pretending to be a completed certificate.
"""
from pathlib import Path
from datetime import datetime,timezone
from mpmath import iv
import argparse,hashlib,heapq,json,math,time
iv.dps=55
p=argparse.ArgumentParser();p.add_argument('--branch',choices=['00','01','11'],required=True)
p.add_argument('--seconds',type=float,default=40);p.add_argument('--target',type=float,default=.075)
args=p.parse_args();assert 1<=args.seconds<=50
HERE=Path(__file__).resolve().parent
trial=json.loads((HERE/'continuum-trial.json').read_text(encoding='utf-8'))
row=next(z for z in trial['branches'] if z['branch']==args.branch)
K,J=trial['K'],trial['J'];start=time.monotonic()
def hi(z):return math.nextafter(float(z.b),math.inf)
def mag(z):return max(abs(float(z.a)),abs(float(z.b)))
def absupper(z):return hi(abs(z))
def polyval(poly,z):
    ans=iv.mpf(0)
    for c in reversed(poly):ans=ans*z+c
    return ans
def derivative(poly):return [i*poly[i] for i in range(1,len(poly))] or [iv.mpf(0)]
def padded(poly,n):return poly+[iv.mpf(0)]*(n-len(poly))
def add(a,b):
    n=max(len(a),len(b));return [x+y for x,y in zip(padded(a,n),padded(b,n))]
def scale(a,c):return [v*c for v in a]
polys=[[iv.mpf(1)]]
for n in range(K-1):
    z=add([iv.mpf(0)]+polys[n],scale(polys[n],-(2*n+iv.mpf('5.5'))))
    if n:z=add(z,scale(polys[n-1],-iv.sqrt(n*(n+iv.mpf('4.5')))))
    polys.append(scale(z,1/iv.sqrt((n+1)*(n+iv.mpf('5.5')))))
trialp=[];trialh=[]
for j in range(J):
    pp=[iv.mpf(0)];hh=[iv.mpf(0)]
    for n in range(K):
        v=iv.mpf(row['basis_vectors'][n][j])
        pp=add(pp,scale(polys[n],v))
        hh=add(hh,scale(polys[n],v*(4*n+11-iv.mpf(row['energies'][j]))))
    trialp.append((pp,derivative(pp),derivative(derivative(pp))))
    trialh.append((hh,derivative(hh),derivative(derivative(hh))))
qbg=iv.mpf(121)/729;qlo=iv.mpf(121)/1458;qhi=iv.mpf(4741)/1458
qq={'00':[qbg,qlo,qlo],'01':[qbg,qlo,qhi],'11':[qbg,qhi,qhi]}[args.branch]
def fvalues(z):
    ff=iv.mpf(0);fp=iv.mpf(0);fpp=iv.mpf(0)
    for q in qq:
        root=iv.sqrt(z*z+2*q)
        ff+=2*q/(root+z)
        fp+=z/root-1
        fpp+=2*q/(root**3)
    return ff,fp,fpp
gamma=iv.mpf(945)/32*iv.sqrt(iv.pi)
def weight(z):return z**iv.mpf('4.5')*iv.exp(-z)/gamma
R=128.0
Fmax=sum(iv.sqrt(2*q) for q in qq)
tails=[]
for j in range(J):
    pp,_,_=trialp[j];hh,_,_=trialh[j]
    positive=[abs(h)+Fmax*abs(a) for a,h in zip(pp,hh)]
    bound=weight(iv.mpf(R))*polyval(positive,iv.mpf(R))**2/(1-(iv.mpf('4.5')+2*(K-1))/R)
    tails.append(hi(bound))
def cell(a,b):
    mid=(a+b)/2;half=(b-a)/2;point=iv.mpf(mid);box=iv.mpf([a,b])
    f0,f1,_=fvalues(point);fi,fpi,fppi=fvalues(box)
    wmax=max(hi(weight(iv.mpf(a))),hi(weight(iv.mpf(b))))
    if a<=4.5<=b:wmax=max(wmax,hi(weight(iv.mpf('4.5'))))
    uppers=[]
    for j in range(J):
        pp,p1,p2=trialp[j];hh,h1,h2=trialh[j]
        pm=polyval(pp,point)
        r0=polyval(hh,point)-f0*pm
        r1=polyval(h1,point)-f1*pm-f0*polyval(p1,point)
        r2=polyval(h2,box)-fppi*polyval(pp,box)-2*fpi*polyval(p1,box)-fi*polyval(p2,box)
        rb=abs(r0)+abs(r1)*half+abs(r2)*(half*half/2)
        uppers.append(hi(iv.mpf(wmax)*(b-a)*rb**2))
    # This score is only an adaptive scheduling heuristic, not evidence.
    score=sum(abs(row['coefficients'][j])*math.sqrt(uppers[j]) for j in range(J))
    return [-score,a,b,uppers]
statefile=HERE/f'residual-partition-{args.branch}.json'
identity=hashlib.sha256((HERE/'continuum-trial.json').read_bytes()).hexdigest()
if statefile.exists():
    old=json.loads(statefile.read_text(encoding='utf-8'))
    assert old['trial_sha256']==identity
    assert old['integrator_sha256']==hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    heap=old['cells']
else:
    heap=[cell(float(a),float(a+2)) for a in range(0,int(R),2)]
heapq.heapify(heap)
initialsq=iv.mpf(0)
for n in range(K):
    cn=sum(iv.mpf(row['basis_vectors'][n][j])*iv.mpf(row['coefficients'][j]) for j in range(J))-(1 if n==0 else 0)
    initialsq+=cn*cn
initial=iv.sqrt(initialsq)
def estimate():
    sq=[iv.mpf(z) for z in tails]
    for entry in heap:
        for j,v in enumerate(entry[3]):sq[j]+=iv.mpf(v)
    residuals=[iv.sqrt(z) for z in sq]
    e=initial+(27*iv.pi/2)*sum(abs(iv.mpf(row['coefficients'][j]))*residuals[j] for j in range(J))
    return hi(e),[hi(v) for v in residuals]
bound,resnorm=estimate();iterations=0
while bound>args.target and time.monotonic()-start<args.seconds:
    item=heapq.heappop(heap);a,b=item[1],item[2];mid=(a+b)/2
    heapq.heappush(heap,cell(a,mid));heapq.heappush(heap,cell(mid,b));iterations+=1
    if iterations%25==0:bound,resnorm=estimate()
bound,resnorm=estimate()
out={'schema':'e20-continuum-residual-interval-v1','utc':datetime.now(timezone.utc).isoformat(),
    'branch':args.branch,'trial_sha256':identity,
    'integrator_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'status':'CERTIFIED_TARGET_MET' if bound<=args.target else 'PARTIAL_RIGOROUS_BOUND_RESUMABLE',
    'vector_error_upper':bound,'initial_vector_error_upper':hi(initial),
    'residual_norm_upper_by_trial':resnorm,'tail_norm_squared_upper_by_trial':tails,
    'new_refinements':iterations,'cell_count':len(heap),'elapsed_seconds':time.monotonic()-start,
    'proof_rule':'initial mismatch + T sum_j |c_j| ||(A_D-lambda_j)f_j||; interval Taylor cell integrals plus gamma-polynomial tail.',
    'cells':heap}
statefile.write_text(json.dumps(out,separators=(',',':'))+'\n',encoding='utf-8')
print(json.dumps({k:v for k,v in out.items() if k!='cells'},indent=2))
