#!/usr/bin/env python3
"""Exact checks supplementing the causal area-memory paper (Python 3.8+, SymPy)."""
from sympy import (Symbol, symbols, Rational, Poly, legendre, expand, prod, cancel,
                   simplify, diff)


def main():
    t, s, r, k, gamma = symbols("t s r k gamma", positive=True)
    xi = (t**2 + s**2 - r**2)/(2*t*s)
    assert simplify(diff(xi, t)**2 - diff(xi, r)**2 - (xi**2-1)/t**2) == 0
    assert simplify(diff(xi, t, 2)-diff(xi, r, 2)-2*xi/t**2) == 0
    print("Legendre retarded-layer chain-rule identities passed.")

    X, p = symbols("X p")
    for n in range(11):
        laurent = expand(legendre(n, (X+1/X)/2))
        coeffs = Poly(expand(X**n*laurent), X).terms()
        assert all(c > 0 for _, c in coeffs)
        lt = sum(c/(p-(power-n)) for (power,), c in coeffs)
        num = prod(p-(n-1-2*j) for j in range(n))
        den = prod(p-(n-2*j) for j in range(n+1))
        assert cancel(lt - num/den) == 0
        lt = cancel(lt)
        assert cancel((p*p-p-n*(n+1))*lt.subs(p,p-1)
                      -(p*p+p-n*(n+1))*lt.subs(p,p+1)) == 0
        residues = [cancel(den/diff(num,p)).subs(p,n-1-2*j) for j in range(n)]
        assert all(x.is_negative for x in residues)
        assert sum(residues) == -Rational(n*(n+1),2)
    print("Inverse-square ladder n=0..10: exact transform, shift, inverse poles, residues passed.")

    q, sig, mpl2 = symbols("q sig mpl2", positive=True)
    k2 = 2*mpl2*q-sig
    potT2 = mpl2*(3*q*q-q)-sig
    assert expand(k2/2+3*sig/2+potT2-3*mpl2*q*q) == 0
    assert expand(k2+sig-2*mpl2*q) == 0
    assert expand(k2*(3*q-1)+3*sig*(q-1)-2*potT2) == 0
    p0 = q/(1-q)
    prescribed_sigma = mpl2*(q*(2*q-1)-gamma*(1-q)**2)/2
    assert simplify(p0*(p0-1)-2*prescribed_sigma/(mpl2*(1-q)**2)-gamma) == 0
    print("Power-law Einstein, weighted clock and target shear equations passed.")


if __name__ == "__main__":
    main()
