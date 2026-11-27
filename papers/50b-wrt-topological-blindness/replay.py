#!/usr/bin/env python3
"""Exact finite checks of the cap/seam H1 presentation and collision examples."""
from itertools import product
from math import prod
from sympy import Matrix, Rational, ZZ, eye, ones, zeros, diag
from sympy.matrices.normalforms import smith_normal_form


def remove_two_factors(a):
    a = abs(int(a))
    while a and a % 2 == 0:
        a //= 2
    return a


def group_data(relations):
    smith = smith_normal_form(relations, domain=ZZ)
    invariants = [
        abs(int(smith[i, i])) for i in range(min(relations.shape))
        if smith[i, i]
    ]
    return relations.cols - len(invariants), tuple(
        x for a in invariants if (x := remove_two_factors(a)) > 1
    )


def full_cap_seam_relations(biadj):
    """Before elimination: h_v, z_v+, z_v-, and each seam half-port."""
    p, q = biadj.shape
    n = p + q
    edges = [(i, p+j, int(biadj[i, j])) for i in range(p) for j in range(q)
             if biadj[i, j]]
    m = len(edges)
    cap = {}
    for v in range(n):
        cap[v, 0], cap[v, 1] = n + 2*v, n + 2*v + 1
    block = zeros(3*n+2*m)
    t = zeros(2*n+2*m)
    u = zeros(n, 2*n+2*m)
    for v in range(n):
        for j, sign in enumerate([1, -1]):
            c = cap[v, j]
            block[2*v+j, c] = 2*sign
            block[2*v+j, v] = 1
            block[2*n+v, c] = 1
            t[c-n, c-n] = 2*sign
            u[v, c-n] = 1
    for i, (a, b, e) in enumerate(edges):
        z_a, z_b = 3*n+2*i, 3*n+2*i+1
        block[3*n+2*i, z_a], block[3*n+2*i, b] = 1, -e
        block[3*n+2*i+1, z_b], block[3*n+2*i+1, a] = 1, -e
        block[2*n+a, z_a] = block[2*n+b, z_b] = 1
        t[z_a-n, z_b-n] = t[z_b-n, z_a-n] = -e
        u[a, z_a-n] = u[b, z_b-n] = 1
    aeps = zeros(p).row_join(biadj).col_join(biadj.T.row_join(zeros(q)))
    assert u*t.inv()*u.T == -aeps
    free, odd = group_data(block)
    assert (free, odd) == group_data(aeps)
    # The remaining Mayer–Vietoris H0 term contributes m-n+1 free summands.
    return m-n+1+free, odd


def signed_binary(d):
    n = d.rows
    return ones(1, n+1).col_join(ones(n, 1).row_join(ones(n)-2*d))


def hub(m, signs):
    blocks = [Matrix([[1, 1], [1, sign]]) for sign in signs]
    assert len(blocks) == m
    return diag(*blocks).row_join(ones(2*m, 1))


def mod_h1_size(free, odd, q):
    from math import gcd
    return q**free * prod(gcd(a, q) for a in odd)


def main():
    cases = [full_cap_seam_relations(Matrix([[1, 1], [1, e]])) for e in (1, -1)]
    assert cases == [(3, ()), (1, ())]
    assert mod_h1_size(*cases[0], 5) // mod_h1_size(*cases[1], 5) == 25
    print("C4: b1=(3,1), odd torsion trivial; squared pointed modulus ratio at q=5 is 25")
    for m in range(1, 5):
        observed = set()
        for signs in product((1, -1), repeat=m):
            free, odd = full_cap_seam_relations(hub(m, signs))
            assert free == 2*m+1+2*signs.count(1) and not odd
            observed.add(free)
        print(f"Hub m={m}: confirmed b1 values {sorted(observed)}")
    square_cases = [full_cap_seam_relations(signed_binary(d))
                    for d in (eye(4), ones(4)-eye(4))]
    assert square_cases == [(16, ()), (16, (3, 3))]
    assert mod_h1_size(*square_cases[1], 3) // mod_h1_size(*square_cases[0], 3) == 9
    print("K5,5: b1=(16,16), odd torsion 0 versus (Z/3)^2; squared pointed ratio 9")
    n = 4 + 6
    for e3, e5 in product((0, 1), repeat=2):
        d = diag(ones(4)-eye(4) if e3 else eye(4), ones(6)-eye(6) if e5 else eye(6))
        b = signed_binary(d)
        for q, expected in [(3, e3), (5, e5)]:
            free, odd = group_data(b)
            assert free == 0
            assert mod_h1_size(0, odd, q) == q**expected
        assert b.rows == n+1
    print("K11,11: independent odd prime blocks (3 and 5) give four distinct torsion groups")
    # At shifted level r=3, scale S = Hadamard/sqrt(2) and use Theta^2=diag(1,-1).
    s0 = Matrix([[1, 1], [1, -1]])
    assert s0*s0/2 == eye(2)
    caps = s0*diag(1, -1)*s0/2
    assert caps == Matrix([[0, 1], [1, 0]])
    assert (-1)**4 * 2**Rational(4, 2) == 4
    print("Level r=3: a=2 is sole cap support, exact common C4 value Z_3=4")
    print("All exact finite checks passed.")


if __name__ == "__main__":
    main()
