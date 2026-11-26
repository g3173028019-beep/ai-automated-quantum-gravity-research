"""Finite-level check; the trigonometric identity is proved in the article."""
from cmath import exp
from math import cos, pi, sin

for r in range(3, 501):
    a = range(1, r)
    v = [exp(1j*pi*(j*j-1)/(2*r))/sin(pi*j/r) for j in a]
    d = sum(abs(z)**2 for z in v)
    ribbon = sum(x.conjugate()*y for x,y in zip(v, reversed(v)))/d
    # B X acts on the reflected input, with the B eigenvalue at the output.
    quadrature = sum(1j*cos(pi*j/r)*x.conjugate()*y
                     for j,x,y in zip(a,v,reversed(v)))/d
    exact = -(-1)**(r//2)*(r*r+2)/(2*(r*r-1)) if r%2 == 0 else 0
    assert abs(ribbon-exact) < 2e-11, (r,ribbon,exact)
    if r%2:
        assert abs(quadrature+1j*exp(1j*pi*r/2)/2) < 2e-11
print("498 levels: fixed-cap ribbon formula and odd-level quadrature passed.")
