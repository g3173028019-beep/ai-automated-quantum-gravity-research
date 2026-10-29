"""Exact SymPy rank audit for fixed-boundary simplex-metric matching."""

from itertools import combinations

import sympy as sp
from sympy.polys.matrices import DomainMatrix


POINTS = [
    (1, -3, 2, 7),
    (8, 11, -4, 2),
    (-5, 4, 13, -2),
    (12, -7, 8, 15),
    (-9, 14, -11, 6),
    (3, -12, 17, -8),
]
COMPONENTS = [(i, i) for i in range(4)] + list(combinations(range(4), 2))


def subtract(first, second):
    return tuple(a - b for a, b in zip(first, second))


def contraction_row(first, second):
    return [
        first[a] * second[b]
        if a == b
        else first[a] * second[b] + first[b] * second[a]
        for a, b in COMPONENTS
    ]


def exact_matrices(side):
    omitted = list(range(2, 6)) if side == "initial" else [0, 1]
    simplices = {k: [v for v in range(6) if v != k] for k in omitted}
    matching_rows = []
    for first_position, first in enumerate(omitted):
        for second_position in range(first_position + 1, len(omitted)):
            second = omitted[second_position]
            shared = sorted(set(simplices[first]) & set(simplices[second]))
            base = shared[0]
            tangents = [subtract(POINTS[v], POINTS[base]) for v in shared[1:]]
            for a in range(3):
                for b in range(a, 3):
                    row = [0] * (10 * len(omitted))
                    block = contraction_row(tangents[a], tangents[b])
                    row[10 * first_position : 10 * (first_position + 1)] = block
                    row[10 * second_position : 10 * (second_position + 1)] = [
                        -value for value in block
                    ]
                    matching_rows.append(row)

    boundary_rows = []
    for first, second in combinations(range(6), 2):
        if (first, second) == (0, 1):
            continue
        containing = [k for k in omitted if k not in (first, second)]
        position = omitted.index(containing[0])
        tangent = subtract(POINTS[second], POINTS[first])
        row = [0] * (10 * len(omitted))
        row[10 * position : 10 * (position + 1)] = contraction_row(tangent, tangent)
        boundary_rows.append(row)
    matching = sp.Matrix(matching_rows)
    boundary = sp.Matrix(boundary_rows)
    return matching, boundary


if __name__ == "__main__":
    for cluster_side in ("initial", "final"):
        matching_matrix, boundary_matrix = exact_matrices(cluster_side)
        combined_matrix = matching_matrix.col_join(boundary_matrix)
        matching_rank = DomainMatrix.from_Matrix(matching_matrix).rank()
        combined_rank = DomainMatrix.from_Matrix(combined_matrix).rank()
        print(
            f"{cluster_side}: matching_rank={matching_rank} "
            f"combined_rank={combined_rank} "
            f"combined_nullity={combined_matrix.cols - combined_rank}"
        )
