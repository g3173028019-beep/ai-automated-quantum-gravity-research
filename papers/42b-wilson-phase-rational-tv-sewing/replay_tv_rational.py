"""Checks exact rational kernel against cyclic modular words and a genuine Euler-zero example."""
from itertools import product
import numpy as np

def candidate(M,N):
    a,b,c,d=M.flatten()
    n,m=np.arange(N)[None,:],np.arange(N)[:,None]
    out=np.zeros((N,N),complex)
    for l in range(abs(c)):
        nn=n+N*l
        num=(a*m*m+d*nn*nn-2*m*nn)%(2*N*abs(c))
        out+=np.exp(1j*np.pi*num/(N*c))
    return out/np.sqrt(N*abs(c))

def replay_kernel():
    S=np.array([[0,-1],[1,0]])
    T=np.array([[1,1],[0,1]])
    for N in (6,8,10,12,18,24,30):
        n=np.arange(N)
        F=np.exp(-2j*np.pi*((n[:,None]*n[None,:])%N)/N)/np.sqrt(N)
        D=np.diag(np.exp(1j*np.pi*((n*n)%(2*N))/N))
        gen=(S,T,np.array([[1,0],[6,1]]))
        ops=(F,D,F@np.linalg.matrix_power(D,-6)@F.conj().T)
        for word in product(range(3), repeat=5):
            M=np.eye(2,dtype=int); W=np.eye(N,dtype=complex)
            for j in word:M=M@gen[j]; W=W@ops[j]
            if M[1,0]==0:continue
            K=candidate(M,N)
            ph=np.vdot(W,K)/N
            assert abs(abs(ph)-1)<2e-12
            assert np.max(np.abs(K-ph*W))<1e-10, (M,N)
    print('Kernel agrees with cyclic S/T for every length-5 test word at N=6,8,10,12,18,24,30.')

def rational_benchmark():
    for r in (41,81,201,601):
        a=np.arange(1,r)
        S=np.sqrt(2/r)*np.sin(np.pi*((a[:,None]*a[None,:])%(2*r))/r)
        T=np.exp(1j*np.pi*((a*a-1)%(4*r))/(2*r))
        id2=S@S
        e=lambda b: id2[:,b-1]
        g1=S@(T*T*(S@(T**-2*(S@e(1))))) # S_chiral is projectively equal to S_chiral^-1
        g2=S@(T**10*(S@e(2)))
        g3=S@(T**-2*(S@e(2)))
        prod=(np.sqrt(r)**3)*g1*g2*g3
        x=a/r
        target=(a%2==1)*np.sin(np.pi*((((4*r+1)//5)*a)%(2*r))/r)**2*np.sin(np.pi*x)
        C=np.vdot(target,prod)/np.vdot(target,target)
        Z=np.dot(np.sin(np.pi/r)/np.sin(np.pi*x),prod)/C
        assert np.max(np.abs(prod-C*target))<3e-10
        assert abs(abs(C)-4*np.sqrt(2))<3e-10
        assert abs(Z-(r/4)*np.sin(np.pi/r))<3e-10
        print(r, 'closure units',Z,'common magnitude',abs(C))

if __name__=='__main__':
    replay_kernel()
    rational_benchmark()
