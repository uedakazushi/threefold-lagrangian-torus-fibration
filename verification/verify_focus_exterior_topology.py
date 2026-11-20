#!/usr/bin/env python3
"""Finite exact-integer checks for the marked focus-exterior boundary maps.

Standard library only.  This checks the displayed homology presentations and
their Smith normal forms for finitely many degrees and finite cyclic covers.
It does not prove the geometric identifications, the general formulas, a
homotopy equivalence, or the existence of a global Lagrangian fibration.
The proofs remain in the current threefold manuscript.

Run from any directory:
    python3 verification/verify_focus_exterior_topology.py
    python3 verification/verify_focus_exterior_topology.py --json

The optional compact report is written to build/focus-exterior-topology.json.
No report is written by default.  Updated 2026-11-21.
"""

from argparse import ArgumentParser
from itertools import combinations, permutations
from math import gcd, prod
from pathlib import Path
import json


REPO = Path(__file__).resolve().parents[1]


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def zero_matrix(rows, columns):
    return [[0] * columns for _ in range(rows)]


def smith_diagonal(matrix):
    """Return positive nonzero Smith entries using integer elementary moves.

    Euclidean row/column operations clear each pivot.  If the pivot does not
    divide the remaining block, adding an offending row exposes a smaller
    gcd.  All operations are unimodular; no floating-point arithmetic is used.
    The input is not mutated.  Zero diagonal entries are omitted.
    """
    a = [list(row) for row in matrix]
    if not a:
        return []
    rows, columns = len(a), len(a[0])
    require(all(len(row) == columns for row in a), "Ragged matrix")
    require(all(isinstance(value, int) for row in a for value in row),
            "Smith input must consist of integers")
    pivot = 0
    while pivot < min(rows, columns):
        location = next(((i, j) for i in range(pivot, rows)
                         for j in range(pivot, columns) if a[i][j]), None)
        if location is None:
            break
        i, j = location
        a[pivot], a[i] = a[i], a[pivot]
        for row in a:
            row[pivot], row[j] = row[j], row[pivot]

        while True:
            below = next((i for i in range(pivot + 1, rows)
                          if a[i][pivot]), None)
            if below is not None:
                quotient = a[below][pivot] // a[pivot][pivot]
                a[below] = [x - quotient * y
                            for x, y in zip(a[below], a[pivot])]
                if a[below][pivot]:
                    a[pivot], a[below] = a[below], a[pivot]
                continue

            right = next((j for j in range(pivot + 1, columns)
                          if a[pivot][j]), None)
            if right is not None:
                quotient = a[pivot][right] // a[pivot][pivot]
                for row in a:
                    row[right] -= quotient * row[pivot]
                if a[pivot][right]:
                    for row in a:
                        row[pivot], row[right] = row[right], row[pivot]
                continue

            # A unit divides the whole block without any additional search.
            offending = None if abs(a[pivot][pivot]) == 1 else next(
                ((i, j) for i in range(pivot + 1, rows)
                 for j in range(pivot + 1, columns)
                 if a[i][j] % a[pivot][pivot]), None)
            if offending is None:
                break
            i, _ = offending
            a[pivot] = [x + y for x, y in zip(a[pivot], a[i])]

        if a[pivot][pivot] < 0:
            a[pivot] = [-x for x in a[pivot]]
        pivot += 1

    require(all(value == 0 for i, row in enumerate(a)
                for j, value in enumerate(row) if i != j),
            "Smith reduction left an off-diagonal entry")
    diagonal = [a[i][i] for i in range(pivot)]
    require(all(diagonal[i + 1] % diagonal[i] == 0
                for i in range(len(diagonal) - 1)),
            "Smith diagonal does not form a divisibility chain")
    return diagonal


def small_determinant(matrix):
    """Independent Leibniz expansion, used only for minors of order <= 3."""
    size = len(matrix)
    total = 0
    for permutation in permutations(range(size)):
        inversions = sum(permutation[i] > permutation[j]
                         for i in range(size) for j in range(i + 1, size))
        total += (-1) ** inversions * prod(
            matrix[i][permutation[i]] for i in range(size))
    return total


def check_smith_implementation():
    """Compare the algorithm with exact gcds of all small minors."""
    examples = [
        [[0, 0], [0, 0]], [[2, 0], [0, 3]], [[0, 6], [10, 0]],
        [[4, 6], [3, 9]], [[-3, 9, 0], [6, -18, 0]],
    ]
    examples.extend(
        [[((i + 1) * (j + 2) * (seed + 1) + i * i + j * j * seed) % 13 - 6
          for j in range(4)] for i in range(3)]
        for seed in range(12)
    )
    for matrix in examples:
        diagonal = smith_diagonal(matrix)
        for order in range(1, min(len(matrix), len(matrix[0])) + 1):
            divisor = 0
            for rr in combinations(range(len(matrix)), order):
                for cc in combinations(range(len(matrix[0])), order):
                    minor = [[matrix[i][j] for j in cc] for i in rr]
                    divisor = gcd(divisor, small_determinant(minor))
            expected = prod(diagonal[:order]) if order <= len(diagonal) else 0
            require(divisor == expected,
                    f"Determinantal-divisor check failed: {matrix}, order {order}")
    return len(examples)


def check_boundary_degree(d):
    e = d - 2
    genus = (d - 3) * (d - 4) // 2
    rank_curve = 2 * genus
    r = rank_curve

    # H2(boundary): u, x_i*b, x_i*t, F, T.
    # H2(V): h, x_i*mu, nu.  The image of u is e*h.
    v2 = zero_matrix(r + 2, 2 * r + 3)
    v2[0][0] = e
    for i in range(r):
        v2[1 + i][1 + r + i] = 1
    v2[r + 1][2 * r + 2] = 1
    require(smith_diagonal(v2) == [1] * (r + 1) + [e],
            f"H2 boundary cokernel failed for d={d}")

    # H2(U): u, x_i*beta, sigma.  The full MV map has no torsion cokernel.
    u2 = zero_matrix(r + 2, 2 * r + 3)
    u2[0][0] = 1
    for i in range(r):
        u2[1 + i][1 + i] = 1
    u2[r + 1][2 * r + 1] = 1
    mv2 = u2 + [[-entry for entry in row] for row in v2]
    require(smith_diagonal(mv2) == [1] * (2 * r + 3),
            f"Full H2 MV map is not primitive for d={d}")

    # H1(boundary): x_i, b, t; only t maps to H1(V)=Z.
    v1 = [[0] * (r + 1) + [1]]
    require(smith_diagonal(v1) == [1], f"H1 surjectivity failed for d={d}")

    # H4(boundary): u*F, u*T, x_i*[Y].  Images are 0, u*nu, ell_i.
    v4 = zero_matrix(r + 1, r + 2)
    v4[0][1] = 1
    for i in range(r):
        v4[i + 1][i + 2] = 1
    require(smith_diagonal(v4) == [1] * (r + 1),
            f"H4 surjectivity failed for d={d}")
    return {"degree": d, "curve_genus": genus,
            "h2_boundary_index": e, "boundary_matrices_checked": 4}


def check_cover_degree(d, sheets):
    e = d - 2
    genus = (d - 3) * (d - 4) // 2
    m = sheets

    # Generators u,h_0,...,h_(m-1); column j is e*h_j-u.
    relation = zero_matrix(m + 1, m)
    for j in range(m):
        relation[0][j] = -1
        relation[j + 1][j] = e
    require(smith_diagonal(relation) == [1] + [e] * (m - 1),
            f"Cyclic-cover relation matrix failed for d={d}, m={m}")

    # Independent presentation of Z[t^+-1]/(e(t-1),t^m-1): e(P-I).
    cyclic = zero_matrix(m, m)
    for j in range(m):
        cyclic[j][j] -= e
        cyclic[(j + 1) % m][j] += e
    require(smith_diagonal(cyclic) == [e] * (m - 1),
            f"Laurent-module coinvariant matrix failed for d={d}, m={m}")

    return {"degree": d, "sheets": m,
            "h2_free_rank": m * (2 * genus + 1) + 1,
            "h2_torsion_order_each": e, "h2_torsion_summands": m - 1}


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--degree-min", type=int, default=6)
    parser.add_argument("--degree-max", type=int, default=20)
    parser.add_argument("--cover-max", type=int, default=12)
    parser.add_argument("--json", action="store_true",
                        help="Write a compact report under build/")
    args = parser.parse_args()
    require(6 <= args.degree_min <= args.degree_max,
            "Degree range must satisfy 6 <= min <= max")
    require(args.cover_max >= 1, "Cover bound must be positive")

    self_checks = check_smith_implementation()
    boundary = [check_boundary_degree(d)
                for d in range(args.degree_min, args.degree_max + 1)]
    covers = [check_cover_degree(d, m)
              for d in range(args.degree_min, args.degree_max + 1)
              for m in range(1, args.cover_max + 1)]
    report = {
        "status": "passed", "arithmetic": "exact integer; standard library only",
        "scope": "Finite presentation checks; not a substitute for the general proofs.",
        "not_checked": ["geometric identifications", "homotopy equivalences",
                        "general theorem proofs", "global Lagrangian gluing"],
        "degree_range": [args.degree_min, args.degree_max],
        "cover_range": [1, args.cover_max],
        "smith_self_checks_by_determinantal_divisors": self_checks,
        "boundary_matrices_checked": 4 * len(boundary),
        "cover_relation_matrices_checked": len(covers),
        "coinvariant_matrices_checked": len(covers),
        "sextic_examples": [item for item in covers
                            if item["degree"] == 6 and item["sheets"] in (1, 2, 12)],
    }
    print(f"PASS: {self_checks} Smith-algorithm checks by determinantal divisors.")
    print(f"PASS: {4 * len(boundary)} boundary matrices; "
          f"{len(covers)} cyclic-cover presentations and "
          f"{len(covers)} Laurent-module coinvariant matrices.")
    print(f"Range: d={args.degree_min}..{args.degree_max}, "
          f"m=1..{args.cover_max}. Exact integer arithmetic only.")
    print("Finite checks are not a substitute for the general proofs or geometric constructions.")
    if args.json:
        destination = REPO / "build" / "focus-exterior-topology.json"
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n",
                               encoding="utf-8")
        print(f"Report: {destination.relative_to(REPO)}")


if __name__ == "__main__":
    main()
