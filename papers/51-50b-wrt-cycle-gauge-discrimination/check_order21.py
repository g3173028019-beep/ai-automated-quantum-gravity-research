#!/usr/bin/env python3
"""Exact rank certificates for the odd-order DW separator; Python 3 only."""
from collections import Counter
from itertools import combinations, product

P = 7
POW = (1, 2, 4)
FREE = tuple(product(range(1, 4), repeat=2))
N = 9

def signs(case):
    return [[1 if not i or not j else
             (-1 if (i == j) == (case == 0) else 1)
             for j in range(4)] for i in range(4)]

def determinant3(a):
    return (a[0][0]*(a[1][1]*a[2][2]-a[1][2]*a[2][1])
            -a[0][1]*(a[1][0]*a[2][2]-a[1][2]*a[2][0])
            +a[0][2]*(a[1][0]*a[2][1]-a[1][1]*a[2][0]))

def rank_eliminate(a, p=P):
    a = [[x % p for x in row] for row in a]
    rank = 0
    for j in range(len(a[0])):
        piv = next((k for k in range(rank, len(a)) if a[k][j]), None)
        if piv is None:
            continue
        a[rank], a[piv] = a[piv], a[rank]
        ur = pow(a[rank][j], -1, p)
        for k in range(rank+1, len(a)):
            c = a[k][j] * ur
            if c:
                a[k] = [(x-c*y) % p for x, y in zip(a[k], a[rank])]
        rank += 1
    return rank

def rank_small_minors(a):
    if determinant3(a) % P:
        return 3
    for i, k in combinations(range(3), 2):
        for j, l in combinations(range(3), 2):
            if (a[i][j]*a[k][l]-a[i][l]*a[k][j]) % P:
                return 2
    return int(any(x % P for row in a for x in row))

def ranks(sign, exponents):
    plus = [row.copy() for row in sign]
    minus_t = [row.copy() for row in sign]  # transpose of reverse block
    for (i,j), e in zip(FREE, exponents):
        plus[i][j] *= POW[e]
        minus_t[i][j] *= POW[-e % 3]
    general = (rank_eliminate(plus), rank_eliminate(minus_t))
    schur = tuple(1+rank_small_minors([[m[i][j]-1
                                      for j in range(1,4)] for i in range(1,4)])
                  for m in (plus, minus_t))
    assert general == schur
    return general

def analyze(case):
    sign = signs(case)
    assert rank_eliminate(sign, 3) == 4
    hist = Counter()
    for exponents in product(range(3), repeat=N):
        rp, rm = ranks(sign, exponents)
        hist[8-rp-rm] += 1
    assert sum(hist.values()) == 3**N
    weighted = sum(P**nullity * count for nullity, count in hist.items())
    return dict(sorted(hist.items())), weighted, P**N * weighted

def main():
    actual = [analyze(case) for case in range(2)]
    expected = [
        ({0: 14415, 1: 5070, 2: 198}, 59607, 2405357452449),
        ({0: 13687, 1: 5712, 2: 284}, 67587, 2727379236309),
    ]
    assert actual == expected
    print("case | nullity histogram | weighted count | vacuum-divided DW")
    for case, data in enumerate(actual):
        print(f"D{case} | {data[0]} | {data[1]} | {data[2]}")
    print("All 39,366 rank pairs verified by two exact methods.")

if __name__ == "__main__":
    main()
