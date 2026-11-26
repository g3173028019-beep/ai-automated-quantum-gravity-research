#!/usr/bin/env python3
"""Check marked Dehn-filling supports over Z[X]/(Phi_{4r}(X)).

This optional check uses integer arithmetic and the Python standard library.
The paper proves the counts without this finite verification.
"""
from functools import lru_cache
from math import gcd


def trim(poly):
    while poly and poly[-1] == 0:
        poly.pop()
    return poly


def divide_monic(poly, divisor):
    rem = trim(list(poly))
    quotient = [0] * max(0, len(rem) - len(divisor) + 1)
    terms = [(i, x) for i, x in enumerate(divisor) if x]
    while len(rem) >= len(divisor):
        const = rem[-1]
        off = len(rem) - len(divisor)
        quotient[off] = const
        for i, x in terms:
            rem[off + i] -= const * x
        trim(rem)
    return trim(quotient), rem


@lru_cache(None)
def cyclotomic(n):
    poly = [-1] + [0] * (n - 1) + [1]
    for d in range(1, n):
        if n % d == 0:
            poly, rem = divide_monic(poly, cyclotomic(d))
            assert not rem
    return tuple(poly)


def marked_support(r, n, b=1):
    modulus = 4 * r
    cyclo = cyclotomic(modulus)
    result = set()
    for a in range(1, r):
        coefficients = [0] * modulus
        for j in range(1, r):
            twist = n * (j * j - 1)
            for u, v in ((1, 1), (1, -1), (-1, 1), (-1, -1)):
                coefficients[(twist + 2 * j * (u * a + v * b)) % modulus] += u * v
        _, rem = divide_monic(coefficients, cyclo)
        if rem:
            result.add(a)
    return result


def odd_support_formula(r, m, b=1):
    d = gcd(m, r)
    assert r % 2 and (b == 1 or (d > 1 and b % d))
    return {
        a for a in range(1, r)
        if (a + b + m) % 2 == 0 and
        (d == 1 or (a - b) % d == 0 or (a + b) % d == 0)
    }


if __name__ == "__main__":
    for r in range(4, 22):
        assert marked_support(r, 2) == {a for a in range(1, r) if (a + r) % 2}
        assert marked_support(r, 2 * r) == {r - 1}
        if r % 2:
            for m in range(r + 1):
                assert marked_support(r, 2 * m) == odd_support_formula(r, m)

    a, b = (marked_support(15, 10, t) for t in (1, 3))
    assert a == {4, 6, 14} and b == {2, 8, 12}
    assert not a & b

    for b, expected_joint in ((1, 3), (3, 6)):
        first = marked_support(135, 18, b)
        second = marked_support(135, 30, b)
        assert (len(first), len(second), len(first & second)) == (15, 9, expected_joint)
    print("Exact cyclotomic checks passed.")
