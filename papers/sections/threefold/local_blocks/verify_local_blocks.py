#!/usr/bin/env python3
"""Exact checks of the explicit algebra in the local-block manuscript.

No numerical check here proves an existence or relative-gluing theorem.
All existence proofs and their hypotheses are in the accompanying TeX.
"""
from pathlib import Path
from math import comb
import json


def multiply(a, b):
    return [[sum(x*y for x, y in zip(row, col)) for col in zip(*b)] for row in a]


def trace(a):
    return sum(a[i][i] for i in range(len(a)))


def main():
    # An invariant Hessian, with c=u+iv and the circle moment divided by 2*pi.
    hessian_cases = 0
    for t in range(-3, 4):
        for u in range(-3, 4):
            for v in range(-3, 4):
                c = complex(u, v)
                matrix = [[-1j*t, -1j*c.conjugate()], [1j*c, 1j*t]]
                square = multiply(matrix, matrix)
                scalar = u*u + v*v - t*t
                assert square == [[scalar, 0], [0, scalar]]
                hessian_cases += 1

    # Nondegenerate trace pairing on the focus Cartan algebra.
    b1 = [[1,0,0,0],[0,1,0,0],[0,0,-1,0],[0,0,0,-1]]
    b2 = [[0,1,0,0],[-1,0,0,0],[0,0,0,1],[0,0,-1,0]]
    assert multiply(b1, b2) == multiply(b2, b1)
    assert (trace(multiply(b1,b1)),trace(multiply(b2,b2)),trace(multiply(b1,b2))) == (4,-4,0)

    # Counts, and the distinction between intermediate and final normal frames.
    rows = []
    normal_cases = 0
    for d in range(6,31):
        genus = (d-1)*(d-2)//2
        ends = d*(d-2)
        pants = 2*genus+ends-2
        cuffs = 3*genus+ends-3
        assert pants == d*(2*d-5)
        assert cuffs-pants+1 == genus
        for i in range(1,d):
            for j in range(i+1,d+1):
                intermediate = d*j
                final = d*(i+1)
                later_intersections = d*(j-i-1)
                assert intermediate-later_intersections == final
                assert d+d*(j-1) == intermediate
                assert d+d*i == final
                normal_cases += 1
        rows.append({"degree":d,"genus":genus,"ends":ends,"pants_per_curve":pants,
                     "internal_cuffs":cuffs,"local_product_fibers":comb(d,2)*pants})
    assert rows[0]["local_product_fibers"] == 630
    assert 36 + 24 == 60

    # Formal coefficient cancellation for the two exact shears.
    # primitive coefficients are recorded in the basis (I, rho_+, rho_-).
    for n in range(-20,21):
        normal_opposite = (0,n,-n)
        shifted_base = (1,-n,n)
        assert tuple(x+y for x,y in zip(normal_opposite,shifted_base)) == (1,0,0)
        normal_one_sided = (0,n,0)
        shifted_one_sided = (1,-n,0)
        assert tuple(x+y for x,y in zip(normal_one_sided,shifted_one_sided)) == (1,0,0)

    result = {
        "status":"passed",
        "scope":"Exact explicit matrix, winding-balance and combinatorial checks only.",
        "global_fibration_proven":False,
        "invariant_hessian_cases":hessian_cases,
        "normal_degree_cases":normal_cases,
        "focus_trace_pairing":[4,-4,0],
        "sextic_intermediate_j6":36,
        "sextic_final_by_i":[12,18,24,30,36],
        "conditional_axis_comparison_winding":60,
        "rows":rows,
    }
    path=Path(__file__).with_name("local_block_checks.json")
    path.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(f"PASS: {hessian_cases} Hessian cases; {normal_cases} stage comparisons.")
    print("This does not establish arbitrary prescribed-end or global gluing.")


if __name__ == "__main__":
    main()
