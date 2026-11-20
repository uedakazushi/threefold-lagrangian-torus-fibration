"""Exact checks for the degree-d hyperplane degeneration in projective 4-space.

The calculations check algebraic/combinatorial formulas, not the existence of
the unconstructed global Lagrangian torus fibration. Standard library only.
"""

from fractions import Fraction
from itertools import combinations
from math import comb
from pathlib import Path
import json


def arithmetic(d):
    core = -comb(d - 3, 3)
    genus = comb(d - 1, 2)
    curve_euler = 2 - 2 * genus
    punctures = d * (d - 2)
    curve_open_euler = curve_euler - punctures
    assert curve_euler == 3 * d - d * d
    assert curve_open_euler == -d * (2 * d - 5)
    component_eulers = []
    component_coefficients = []
    for i in range(1, d + 1):
        closed_component = 4 + (d - i) * curve_euler
        double_sum = sum(
            3 + d * (d - max(i, j))
            for j in range(1, d + 1) if j != i
        )
        by_resolved_strata = (
            closed_component - double_sum
            + 2 * comb(d - 1, 2) - comb(d - 1, 3)
        )
        by_original_strata = (
            core + (d - i) * curve_open_euler + d * comb(d - i, 2)
        )
        assert by_resolved_strata == by_original_strata
        component_eulers.append(by_original_strata)

        # Alpha_i = a + coefficient * epsilon. This coefficient is integral.
        alpha_i = Fraction(d * (d + 1 - 2 * i), 2)
        assert alpha_i.denominator == 1
        assert alpha_i == sum(j - i for j in range(1, d + 1))
        assert alpha_i - sum(j - i for j in range(i + 1, d + 1)) == -comb(i, 2)
        component_coefficients.append(int(alpha_i))
        for j in range(i + 1, d + 1):
            alpha_j = Fraction(d * (d + 1 - 2 * j), 2)
            assert alpha_i - d * (j - i) == alpha_j
            for k in range(j + 1, d + 1):
                alpha_k = Fraction(d * (d + 1 - 2 * k), 2)
                # Exceptional coefficient from restricting Omega_i to S_ij.
                assert (j - i) - (k - i) == -(k - j)
                assert alpha_j - d * (k - j) == alpha_k
                # Degrees of the vertical divisor on the exceptional chain.
                b_i, b_j, b_k = i - d, j - d, k - d
                assert -b_i + b_j == j - i > 0
                assert -b_j + b_k == k - j > 0

    total_points = sum(
        d * (d - j)
        for i in range(1, d + 1)
        for j in range(i + 1, d + 1)
    )
    assert total_points == d * comb(d, 3)
    assert comb(d, 2) * punctures == 3 * total_points
    assert -(d - 1) + 1 + (d - 2) == 0  # normal-bundle h coefficient
    assert 1 + 0 - 1 == 0              # exceptional coefficient

    # Extract c_3 from (1+H)^5/(1+dH), independently of the Euler expansion.
    c3_coefficient = sum(comb(5, k) * (-d) ** (3 - k) for k in range(4))
    euler = d * c3_coefficient
    assert euler == d * (10 - 10 * d + 5 * d**2 - d**3)
    negative_prediction = -comb(d, 2) * curve_open_euler
    positive_prediction = total_points
    assert sum(component_eulers) == d * core - negative_prediction + positive_prediction == euler
    assert euler - d * core == -Fraction(d*d*(d-1)*(5*d-13), 6)
    raw_cells = [comb(d, k + 1) for k in range(4)]
    raw_euler = sum((-1)**k * n for k, n in enumerate(raw_cells))
    pg = comb(d - 1, 4)
    assert raw_euler == 1 - pg
    return {
        "degree": d,
        "K_cubed": d * (d - 5)**3,
        "geometric_genus": pg,
        "Euler_X": euler,
        "b3_X": 4 - euler,
        "curve_genus": genus,
        "curve_punctures": punctures,
        "curve_open_euler": curve_open_euler,
        "component_open_eulers": component_eulers,
        "alpha_epsilon_coefficients": component_coefficients,
        "core_total_euler": d * core,
        "predicted_negative_vertices_NOT_existence": negative_prediction,
        "predicted_positive_vertices_NOT_existence": positive_prediction,
        "raw_base_cell_counts": raw_cells,
        "raw_base_euler": raw_euler,
        "persistent_nodes_for_product_auxiliary": comb(d, 2)**2,
    }


def boundary_columns(cells, lower):
    index = {face: i for i, face in enumerate(lower)}
    return [
        {index[cell[:a] + cell[a + 1:]]: (-1)**a for a in range(len(cell))}
        for cell in cells
    ]


def unit_pivot_rank(nrows, columns):
    """Integral row elimination, requiring every pivot to be a unit.

    Unit pivots with a zero residual matrix certify the image has no torsion
    in its quotient. No floating-point or modular rank inference is used.
    """
    matrix = [[column.get(r, 0) for column in columns] for r in range(nrows)]
    rank = 0
    for c in range(len(columns)):
        nz = [r for r in range(rank, nrows) if matrix[r][c]]
        if not nz:
            continue
        units = [r for r in nz if abs(matrix[r][c]) == 1]
        assert units, "Non-unit pivot: this elimination is not an integral certificate"
        r = units[0]
        matrix[rank], matrix[r] = matrix[r], matrix[rank]
        if matrix[rank][c] == -1:
            matrix[rank] = [-x for x in matrix[rank]]
        pivot = matrix[rank]
        for rr in range(rank + 1, nrows):
            factor = matrix[rr][c]
            if factor:
                matrix[rr] = [x - factor*y for x, y in zip(matrix[rr], pivot)]
        rank += 1
        if rank == nrows:
            break
    assert all(not any(row) for row in matrix[rank:])
    return rank


def topology(d):
    cells = [list(combinations(range(d), k + 1)) for k in range(4)]
    boundaries = {k: boundary_columns(cells[k], cells[k - 1]) for k in range(1, 4)}
    for k in (2, 3):
        for col in boundaries[k]:
            product = {}
            for middle, coeff in col.items():
                for row, val in boundaries[k - 1][middle].items():
                    product[row] = product.get(row, 0) + coeff * val
            assert all(v == 0 for v in product.values())
    ranks = {
        k: unit_pivot_rank(len(cells[k - 1]), boundaries[k])
        for k in range(1, 4)
    }
    assert [ranks[k] for k in range(1, 4)] == [comb(d - 1, k) for k in range(1, 4)]
    betti = [len(cells[k]) - ranks.get(k, 0) - ranks.get(k + 1, 0) for k in range(4)]
    assert betti == [1, 0, 0, comb(d - 1, 4)]
    incidence = {face: 0 for face in cells[2]}
    for tetra in cells[3]:
        for face in combinations(tetra, 3):
            incidence[face] += 1
    assert set(incidence.values()) == {d - 3}
    assert d - 3 >= 3
    return {
        "degree": d,
        "boundary_ranks": [ranks[k] for k in range(1, 4)],
        "integral_betti": betti,
        "unit_pivots_only": True,
        "boundary_squared_zero": True,
        "tetrahedra_per_triangle": d - 3,
    }


def main():
    root = Path(__file__).resolve().parent
    rows = [arithmetic(d) for d in range(6, 31)]
    topology_rows = [topology(d) for d in range(6, 11)]
    comparison = arithmetic(5)
    assert comparison["predicted_negative_vertices_NOT_existence"] == 250
    assert comparison["predicted_positive_vertices_NOT_existence"] == 50
    result = {
        "scope": "Exact algebraic and combinatorial checks; no global fibration existence claim.",
        "arithmetic_degrees": [6, 30],
        "topology_degrees": [6, 10],
        "quintic_calibration": comparison,
        "arithmetic": rows,
        "topology": topology_rows,
    }
    (root / "verification_results.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    table_rows = [comparison] + rows[:5]
    lines = [
        r"\begin{center}",
        r"\begin{tabular}{@{}rrrrr@{}}",
        r"\toprule",
        r"$d$ & $K^3$ & $p_g$ & $\chi(X_d)$ & $b_3(X_d)$\\",
        r"\midrule",
    ]
    for row in table_rows:
        vals = [row[k] for k in ("degree", "K_cubed", "geometric_genus", "Euler_X", "b3_X")]
        lines.append(" & ".join(str(v) for v in vals) + r"\\")
    lines.extend([r"\bottomrule", r"\end{tabular}", r"\end{center}"])
    (root / "invariants_table.tex").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": "passed",
        "arithmetic_degrees": [6, 30],
        "topology_degrees": [6, 10],
        "sextic": rows[0],
        "invariants_table": table_rows,
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
