"""Numerical checks for the analytic growing-color and mesoscopic theorems."""
from math import pi, sqrt
import numpy as np
from scipy.integrate import quad

def stationary_witness():
    xstar=131/1386
    kappa=3*pi/5
    expected=pi/(2*sqrt(6)*np.sin(pi*xstar))
    for m in (20,100,500):
        r=2541+4620*m
        a=np.arange(max(1,int((xstar-.01)*r)),min(r,int((xstar+.01)*r)+1))
        x=a/r; t=(x-xstar)/.01
        sels=np.abs(t)<1
        a,x,t=a[sels],x[sels],t[sels]
        chi=np.exp(-t*t/(1-t*t))
        cap=np.ones(len(a),complex)
        def G(y):
            return sum(np.exp(2j*pi*((2*(r%5)*l*l+y*l)%5)/5) for l in range(5))
        for q in (3,7,11):
            b=r//q
            cap*=(G(2*b-a)*np.exp(-1j*pi*x*b/5)
                    -G(-2*b-a)*np.exp(1j*pi*x*b/5))/sqrt(10)
        Z=np.sum(np.sin(pi/r)/np.sin(pi*x)*chi*cap*np.exp(-1j*kappa*a.astype(float)**2/r))
        print('stationary',r,sqrt(r)*abs(Z),'predicted',expected)
    assert abs(sqrt(r)*abs(Z)-expected)<0.025

def fresnel_mesoscopic():
    kappa=pi/2
    theta=pi*.37
    H=lambda th: .5*np.exp(-1j*pi/4)*sqrt(pi/kappa)*(
        quad(lambda u:np.cos(u*u/(4*kappa)),0,th)[0]
        +1j*quad(lambda u:np.sin(u*u/(4*kappa)),0,th)[0])
    for r in (20001,100001,560001):
        a=np.arange(1,r,dtype=float)
        b=int(round(theta*sqrt(r)/pi))
        if b%2==0:b+=1
        v=np.sin(pi/r)/np.sin(pi*a/r)*np.exp(-1j*kappa*a*a/r)
        Z=sqrt(2)*np.sum(v*np.sin(pi*b*a/r))
        print('nonzero Euler mesoscopic',r,b,Z,'predicted',sqrt(2)*H(theta))
    assert abs(Z-sqrt(2)*H(theta))<0.02
    for r in (30001,150001):
        a=np.arange(1,r,dtype=float)
        b=int(r**.42)|1
        Z=sqrt(2)*np.sum(np.sin(pi/r)/np.sin(pi*a/r)*np.sin(pi*b*a/r))
        assert abs(Z-sqrt(2)*np.sin(pi/r)*(r-b))<1e-10
        print('zero Euler mesoscopic',r,b,Z,'limiting',pi*sqrt(2))

if __name__=='__main__':
    stationary_witness()
    fresnel_mesoscopic()
