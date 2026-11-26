#!/usr/bin/env python3
"""Independent floating-point regressions of exact finite sums, not proof inputs."""
from math import pi, sqrt
import numpy as np
from scipy.fft import dst
from scipy.special import zeta


def S(f):
    return dst(f, type=1, norm="ortho")


def data(r):
    a = np.arange(1, r, dtype=float)
    return a, sqrt(2 / r) * np.sin(pi * a / r)


def cap(r, p):
    a, v = data(r)
    return S(v * np.exp(1j * pi * p * (a*a - 1) / (2*r)))


def cap_recip(r, p):
    a, _ = data(r)
    b = a[:, None] + 2*r*np.arange(abs(p))[None, :]
    return (1j * np.exp(1j * pi * np.sign(p)/4 - 1j*pi*(p+1/p)/(2*r))
            * sqrt(2/(r*abs(p)))
            * np.sum(np.sin(pi*b/(r*p))*np.exp(-1j*pi*b*b/(2*r*p)), axis=1))


def kappa(r):
    return np.exp(1j*pi*(3/4 - 3/(2*r)))


def H(r):
    a, v = data(r)
    f = abs(cap(r, 2))**2 / v
    finite = np.dot(f, S(f))
    d = a.astype(int)[a.astype(int) % 2 == 1]
    cotangent = (-1)**r * sqrt(2/r) * np.sum(1/np.tan(pi*d/(2*r)))
    return finite, cotangent


def closed_star(r):
    _, v = data(r)
    outer = S(cap(r, 2)**4 / v**3)
    middle = cap(r, 4)**2 * cap(r, -2) / v**3
    return kappa(r)**(-8) * np.dot(middle, outer**2)


def triangle(r, paired=True):
    a, v = data(r)
    if paired:
        f = abs(cap(r, 2))**2 / v**2
        k = kappa(r)**(-1)
    else:
        # Smooth positive-corner Wilson cutoff, no reflection or parity selection.
        cutoff = np.zeros_like(a)
        x = a/(.028*r)
        cutoff[x <= .55] = 1
        mid = (.55 < x) & (x < 1)
        u = (x[mid] - .55)/.45
        e1, e2 = np.exp(-1/u), np.exp(-1/(1-u))
        cutoff[mid] = e2/(e1 + e2)
        f = cutoff * cap(r, 3) * cap(r, -5) / v**2
        k = kappa(r)**(-1)
    keep = f != 0
    a, f = a[keep], f[keep]
    M = f[:, None] * sqrt(2/r) * np.sin(pi*np.outer(a,a)/r)
    return k * np.sum((M @ M) * M.T)


print("cap reciprocity, maximum direct discrepancy")
error = max(np.max(abs(cap(r,p)-cap_recip(r,p)))
            for r in (13, 20, 31, 50, 81) for p in (-5,-2,2,3,4))
print(f"{error:.3g}")
assert error < 1e-11
print("H exact finite cotangent, relative discrepancies")
for r in (41, 82, 161, 320):
    actual, formula = H(r)
    discrepancy = abs(actual-formula)/abs(formula)
    print(r, f"{discrepancy:.3g}")
    assert discrepancy < 1e-10
print("star Z/r^5 (even reference, odd leading coefficient 0)")
star = (np.exp(5j*pi/4)/(2*sqrt(2)*pi**2)
        + 3*zeta(3)*np.exp(-1j*pi/4)/(sqrt(2)*pi**4))
print("even reference", star)
for r in (480, 481, 1200, 1201, 2400, 2401):
    print(r, closed_star(r)/r**5)
print("opposite-cap triangle, dephase by (-1)^r exp(-3*pi*i/4)")
for r in (240, 480, 960):
    print(r, triangle(r)/(r**3*(-1)**r * np.exp(-3j*pi/4)))
print("reference", 4/pi**3)
print("marked triangle, positive corner; full complex coefficient")
t = 671/3375
print("reference", -15**(-4.5) * (1/(2*sqrt(2*(1-t))) + 1j/(2*sqrt(2*(1+t)))))
for r in (3840, 7680, 15360, 30720):
    print(r, triangle(r,paired=False))
