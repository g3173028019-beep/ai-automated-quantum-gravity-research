"""Bounded exact checks of omitted long-manuscript supporting claims."""
import itertools
import json
import math
from fractions import Fraction
from pathlib import Path


def splits(vertices, q):
    vertices = tuple(vertices)
    first = vertices[0]
    for rest in itertools.combinations(vertices[1:], q - 1):
        a = (first,) + rest
        aset = set(a)
        b = tuple(x for x in vertices if x not in aset)
        yield (a, b)


def sign_of_partition(blocks):
    seq = tuple(x for block in sorted(blocks) for x in block)
    parity = sum(seq[i] > seq[j] for i in range(len(seq))
                 for j in range(i + 1, len(seq))) % 2
    return 1 - 2 * parity


def four_block_sector(nvertices, q):
    assert nvertices == 4 * q
    universe = tuple(range(nvertices))
    mr = math.comb(2 * q, q) // 2
    h2h2 = 4 * mr * mr * (math.comb(nvertices, q)
                          * math.comb(nvertices-q, q) // 2)
    four_j = 0
    completion_sums = {}
    for support in itertools.combinations(universe, 2 * q):
        complement = tuple(x for x in universe if x not in set(support))
        xx, rr = list(splits(support, q)), list(splits(complement, q))
        signs = [[sign_of_partition(x+r) for r in rr] for x in xx]
        for i in range(len(xx)):
            for j in range(i + 1, len(xx)):
                s = sum(a*b for a,b in zip(signs[i], signs[j]))
                four_j += 4*s*s
                completion_sums[str(s)] = completion_sums.get(str(s),0)+1
    return {"N":nvertices,"q":q,"H2H2":h2h2,"four_J":four_j,
            "total":h2h2+four_j,"completion_sum_counts":completion_sums}


def quadrature():
    # Work in Q[s]/(s^2-13) so exact radicals need no external dependency.
    def add(x,y): return (x[0]+y[0], x[1]+y[1])
    def mul(x,y): return (x[0]*y[0]+13*x[1]*y[1], x[0]*y[1]+x[1]*y[0])
    def power(x,k):
        out=(Fraction(1),Fraction(0))
        for _ in range(k): out=mul(out,x)
        return out
    lo=(Fraction(20),Fraction(-5))
    hi=(Fraction(20),Fraction(5))
    wa=(Fraction(1,4),Fraction(17,260))
    wb=(Fraction(1,4),Fraction(-17,260))
    got=[]
    expected=[1,3,45,1575]
    for k in range(4):
        val=mul((Fraction(2),Fraction(0)),add(mul(wa,power(lo,k)),mul(wb,power(hi,k))))
        assert val==(Fraction(expected[k]),Fraction(0))
        got.append(str(val[0]))
    return {"p4":"x^4-40x^2+75","squared_nodes":"20 +/- 5 sqrt(13)",
            "individual_weights":"1/4 +/- 17/(20 sqrt(13))",
            "even_moments_0_2_4_6":got,"odd_moments_through_7":0,
            "all_exact":True}


results={"odd_q":four_block_sector(12,3),
         "even_q_reference":four_block_sector(8,2),
         "quadrature":quadrature()}
assert results["odd_q"]["total"]==17580480
assert results["even_q_reference"]["total"]==15120
out=Path(__file__).with_name("chaos-extras-results.json")
out.write_text(json.dumps(results,indent=2)+"\n",encoding="utf-8")
print(json.dumps(results,indent=2))
