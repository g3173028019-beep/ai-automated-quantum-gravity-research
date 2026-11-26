"""Check the branch normalizations in equations (7)–(8), not a new topological example."""
import numpy as np
from scipy.integrate import quad

def bump(t):
    t=np.asarray(t)
    a=np.zeros_like(t,dtype=float)
    keep=np.abs(t)<1
    a[keep]=np.exp(-t[keep]**2/(1-t[keep]**2))
    return a

def evaluate(rank):
    center=.44
    width=.13
    q=np.pi*np.array([[1,.2],[.2,-.6]]) if rank==2 else np.pi*np.ones((2,2))
    u=np.array([1,-1])/np.sqrt(2)
    nu=0 if rank==2 else 1.35
    def f(a,b): return bump((a-center)/width)*bump((b-center)/width)
    eig=np.linalg.eigvalsh(q)
    eig=eig[np.abs(eig)>1e-8]
    const=(2*np.pi)**(rank/2)*np.exp(1j*np.pi*np.sign(eig).sum()/4)/np.sqrt(abs(eig).prod())
    null=1 if rank==2 else 2*quad(
        lambda t: float(f(center+t*u[0],center+t*u[1]))*np.cos(nu*t),
        0,width*np.sqrt(2))[0]
    target=const*null
    for r in [1500,6000,18000]:
        a=np.arange(int((center-width)*r),int((center+width)*r)+1)
        x,y=np.meshgrid(a/r-center,a/r-center,sparse=True)
        # Constant exp(-i r x_0 Q x_0/2) removed explicitly.
        phase=r*(q[0,0]*x*x/2+q[0,1]*x*y+q[1,1]*y*y/2)+nu*(u[0]*x+u[1]*y)
        out=(f(x+center,y+center)*np.exp(1j*phase)).sum()/r**(2-rank/2)
        print(rank,r,out,target,abs(out-target))
    assert abs(out-target)<0.015

if __name__=='__main__':
    evaluate(2)
    evaluate(1)
