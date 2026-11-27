"""Exact, small independent checks for the K_{4,4} collision in 50B.

Uses only SymPy and the standard library.  The 24-generator presentation
retains fiber and both cap-section generators; it is compared with the
separately derived 8-generator matrix 4A.  No floating point is used.
Run: python check_signed_homology.py
"""

from __future__ import annotations

import json
from math import gcd, prod
from pathlib import Path

import sympy as sp
from sympy.matrices.normalforms import smith_normal_form


def smith_data(matrix: sp.Matrix) -> dict:
    diagonal = smith_normal_form(matrix, domain=sp.ZZ)
    factors = [abs(int(diagonal[i, i])) for i in range(min(matrix.shape))]
    return {
        "all_diagonal_entries": factors,
        "nontrivial_torsion_factors": [d for d in factors if d > 1],
        "free_rank": factors.count(0),
    }


def hom_count(free_rank: int, factors: list[int], q: int) -> int:
    return q**free_rank * prod(gcd(d, q) for d in factors)


def check_case(name: str, d_matrix: sp.Matrix) -> dict:
    part_size = 4
    vertex_count = 8
    edge_count = 16
    graph_cycle_rank = edge_count - vertex_count + 1
    b_matrix = sp.ones(part_size)
    b_matrix[1:, 1:] = sp.ones(3) - 2 * d_matrix
    adjacency = sp.zeros(vertex_count)
    adjacency[:part_size, part_size:] = b_matrix
    adjacency[part_size:, :part_size] = b_matrix.T

    identity = sp.eye(vertex_count)
    zero = sp.zeros(vertex_count)
    # Generator order: h_v, a_v, b_v.
    # Relations: A h + a + b = 0; h + 2a = 0; h - 2b = 0.
    full_relations = sp.BlockMatrix([
        [adjacency, identity, identity],
        [identity, 2 * identity, zero],
        [identity, zero, -2 * identity],
    ]).as_explicit()
    full = smith_data(full_relations)
    reduced = smith_data(4 * adjacency)
    assert full["nontrivial_torsion_factors"] == reduced["nontrivial_torsion_factors"]
    assert full["free_rank"] == reduced["free_rank"]
    betti = graph_cycle_rank + full["free_rank"]
    assert betti == 9
    factors = full["nontrivial_torsion_factors"]
    counts = {
        str(q): hom_count(betti, factors, q)
        for q in (2, 4, 8, 16, 3, 5, 9)
    }
    for q in (3, 5, 9):
        assert counts[str(q)] == q**9
    return {
        "name": name,
        "B": [[int(x) for x in row] for row in b_matrix.tolist()],
        "det_D": int(d_matrix.det()),
        "det_B": int(b_matrix.det()),
        "D_smith": smith_data(d_matrix),
        "full_24_by_24_smith": full,
        "reduced_8_by_8_smith": reduced,
        "graph_cycle_rank": graph_cycle_rank,
        "first_betti_number": betti,
        "cyclic_Hom_counts": counts,
    }


def main() -> None:
    cases = [
        check_case("D0 = I3", sp.eye(3)),
        check_case("D1 = J3 - I3", sp.ones(3) - sp.eye(3)),
    ]
    expected_factors = [
        [4, 4, 8, 8, 8, 8, 8, 8],
        [4, 4, 8, 8, 8, 8, 16, 16],
    ]
    for case, expected in zip(cases, expected_factors):
        assert case["full_24_by_24_smith"]["nontrivial_torsion_factors"] == expected
    for q in (2, 4, 8, 3, 5, 9):
        assert cases[0]["cyclic_Hom_counts"][str(q)] == cases[1]["cyclic_Hom_counts"][str(q)]
    assert cases[0]["cyclic_Hom_counts"]["16"] == 2**58
    assert cases[1]["cyclic_Hom_counts"]["16"] == 2**60
    output = {
        "scope": "Exact integral relation matrices and finite cyclic Hom counts; no asymptotics.",
        "sympy_version": sp.__version__,
        "cases": cases,
        "all_assertions_passed": True,
        "Z16_ratio": 4,
    }
    target = Path(__file__).with_name("signed_homology_check.json")
    target.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(target), "all_assertions_passed": True, "Z16_ratio": 4}))


if __name__ == "__main__":
    main()
