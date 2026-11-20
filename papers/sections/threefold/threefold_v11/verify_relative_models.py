#!/usr/bin/env python3
"""Exact arithmetic audit of the relative curve blow-up class formulas.

Run with Python 3 (standard library only):
    python3 verify_relative_models.py

The default audit covers every d = 6,...,30 and every 1 <= i < j <= d.
All affine polynomials in a, epsilon, mu use integer coefficient vectors;
the intersection pairing is the integral diagonal form (1,-1,...,-1).

This verifies arithmetic consistency of specified formulas.  It does not
prove symplectic positivity, a relative Moser comparison, or existence of
any Lagrangian fibration.  The curve-rank lemmas are not tested here.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path
from typing import Dict, Mapping, Tuple


# Coefficient order: constant, a, epsilon, mu.  No floating point is used.
Poly = Tuple[int, int, int, int]
Class = Dict[str, Poly]
IntegralClass = Dict[str, int]
ZERO: Poly = (0, 0, 0, 0)
MU: Poly = (0, 0, 0, 1)


def add(left: Poly, right: Poly) -> Poly:
    return tuple(x + y for x, y in zip(left, right))  # type: ignore[return-value]


def sub(left: Poly, right: Poly) -> Poly:
    return tuple(x - y for x, y in zip(left, right))  # type: ignore[return-value]


def scale(coefficient: int, value: Poly) -> Poly:
    return tuple(coefficient * x for x in value)  # type: ignore[return-value]


def at_mu(value: Poly, epsilon_multiple: int) -> Poly:
    """Substitute mu = epsilon_multiple * epsilon exactly."""
    constant, a, epsilon, mu = value
    return (constant, a, epsilon + epsilon_multiple * mu, 0)


def alpha(d: int, index: int) -> Poly:
    numerator = d * (d + 1 - 2 * index)
    if numerator % 2:
        raise AssertionError("The coefficient of epsilon in alpha is not integral.")
    return (0, 1, numerator // 2, 0)


def polynomial_text(value: Poly) -> str:
    names = ("", "a", "epsilon", "mu")
    pieces = []
    for coefficient, name in zip(value, names):
        if coefficient == 0:
            continue
        magnitude = abs(coefficient)
        term = str(magnitude) if not name else (
            name if magnitude == 1 else f"{magnitude}*{name}"
        )
        if not pieces:
            pieces.append(("-" if coefficient < 0 else "") + term)
        else:
            pieces.append((" - " if coefficient < 0 else " + ") + term)
    return "".join(pieces) or "0"


def poly_record(value: Poly) -> dict:
    return {"expression": polynomial_text(value), "coefficients": list(value)}


def intersection(left: Mapping[str, int], right: Mapping[str, int]) -> int:
    """Integral diagonal pairing; omitted coordinates are zero."""
    if len(left) > len(right):
        left, right = right, left
    return sum(
        (1 if key == "h" else -1) * coefficient * right.get(key, 0)
        for key, coefficient in left.items()
        if coefficient
    )


class Audit:
    def __init__(self) -> None:
        self.counts: Counter[str] = Counter()

    def require(self, condition: bool, category: str, context: str) -> None:
        if not condition:
            raise AssertionError(f"{category}: {context}")
        self.counts[category] += 1

    def class_equal(
        self, actual: Mapping[str, Poly], expected: Mapping[str, Poly],
        category: str, context: str,
    ) -> None:
        if actual.keys() != expected.keys():
            raise AssertionError(f"{category}: mismatched bases at {context}")
        for key in actual:
            self.require(
                actual[key] == expected[key],
                category,
                f"{context}, basis {key}: {actual[key]} != {expected[key]}",
            )


def verify_pair(d: int, i: int, j: int, audit: Audit) -> dict:
    context = f"d={d}, i={i}, j={j}"
    exceptional = [
        (k, p, f"e_{k}_{p}")
        for k in range(j + 1, d + 1)
        for p in range(1, d + 1)
    ]
    basis = ["h"] + [key for _, _, key in exceptional]

    # Independent input formulas for the stage before and after the blow-up.
    beta: Class = {"h": alpha(d, i)}
    expected_post: Class = {"h": alpha(d, j)}
    curve: IntegralClass = {"h": d}
    normal_before: IntegralClass = {"h": 1}
    for k, _, key in exceptional:
        beta[key] = (0, 0, -(k - i), 0)
        expected_post[key] = (0, 0, -(k - j), 0)
        curve[key] = -1
        normal_before[key] = 0

    delta_multiple = j - i
    delta: Poly = (0, 0, delta_multiple, 0)
    post = {
        key: sub(beta[key], scale(curve[key], delta))
        for key in basis
    }
    audit.class_equal(
        post, expected_post, "beta_minus_delta_C_coefficients", context
    )

    # Expand the two branches of beta - mu*L - max(0, delta-mu)*C.
    below = {
        key: sub(
            sub(beta[key], scale(normal_before[key], MU)),
            scale(curve[key], sub(delta, MU)),
        )
        for key in basis
    }
    above = {
        key: sub(beta[key], scale(normal_before[key], MU))
        for key in basis
    }

    # The lower branch independently follows from the post-blow-up normal L-C.
    normal_after = {
        key: normal_before[key] - curve[key] for key in basis
    }
    expected_below = {
        key: sub(expected_post[key], scale(normal_after[key], MU))
        for key in basis
    }
    audit.class_equal(
        below, expected_below, "DH_lower_branch_coefficients", context
    )
    audit.class_equal(
        {key: at_mu(below[key], 0) for key in basis},
        expected_post,
        "DH_zero_level_coefficients",
        context,
    )

    below_at_wall = {
        key: at_mu(below[key], delta_multiple) for key in basis
    }
    above_at_wall = {
        key: at_mu(above[key], delta_multiple) for key in basis
    }
    expected_wall = {
        key: sub(beta[key], scale(normal_before[key], delta))
        for key in basis
    }
    audit.class_equal(
        below_at_wall, above_at_wall, "DH_continuity_coefficients", context
    )
    audit.class_equal(
        below_at_wall, expected_wall, "DH_wall_value_coefficients", context
    )

    # Polynomial differentiation is the coefficient of mu.
    for key in basis:
        audit.require(
            below[key][3] == curve[key] - normal_before[key],
            "DH_lower_slope_coefficients",
            f"{context}, basis {key}",
        )
        audit.require(
            above[key][3] == -normal_before[key],
            "DH_upper_slope_coefficients",
            f"{context}, basis {key}",
        )
        audit.require(
            above[key][3] - below[key][3] == -curve[key],
            "DH_slope_jump_coefficients",
            f"{context}, basis {key}",
        )

    curve_square = intersection(curve, curve)
    positive_degree = intersection(normal_before, curve)
    negative_normal = {
        key: curve[key] - normal_before[key] for key in basis
    }
    negative_degree = intersection(negative_normal, curve)
    audit.require(curve_square == d * j, "curve_self_intersection", context)
    audit.require(positive_degree == d, "positive_weight_normal_degree", context)
    audit.require(
        negative_degree == d * (j - 1),
        "negative_weight_normal_degree",
        context,
    )
    audit.require(
        positive_degree + negative_degree == curve_square,
        "normal_degree_sum",
        context,
    )

    # The old strict triple line is h - sum_{p=1}^d e_{k,p}.
    for k in range(j + 1, d + 1):
        strict_line = {"h": 1}
        for p in range(1, d + 1):
            strict_line[f"e_{k}_{p}"] = -1
        audit.require(
            intersection(curve, strict_line) == 0,
            "old_strict_line_disjointness",
            f"{context}, k={k}",
        )
    for k, p, key in exceptional:
        exceptional_curve = {key: 1}
        audit.require(
            intersection(curve, exceptional_curve) == 1,
            "old_exceptional_curve_intersection",
            f"{context}, k={k}, p={p}",
        )

    remaining_indices = [k for k in range(1, j) if k != i]
    audit.require(
        len(remaining_indices) == j - 2,
        "remaining_strict_line_count",
        context,
    )
    for k in remaining_indices:
        audit.require(
            intersection(curve, normal_before) == d,
            "remaining_strict_line_intersection",
            f"{context}, k={k}",
        )
    strict_marks = d * len(remaining_indices)
    exceptional_marks = len(exceptional)
    audit.require(
        exceptional_marks == d * (d - j),
        "old_exceptional_mark_count",
        context,
    )
    audit.require(
        strict_marks + exceptional_marks == d * (d - 2),
        "marked_end_split",
        context,
    )

    return {
        "i": i,
        "j": j,
        "basis_size": len(basis),
        "delta": poly_record(delta),
        "beta_h": poly_record(beta["h"]),
        "post_h": poly_record(post["h"]),
        "DH_below_h": poly_record(below["h"]),
        "DH_above_h": poly_record(above["h"]),
        "DH_wall_h": poly_record(below_at_wall["h"]),
        "all_requested_identities_passed": True,
    }


def sextic_summary(j: int, pairs: list[dict]) -> dict:
    d = 6
    return {
        "j": j,
        "number_of_pairs_i_less_than_j": j - 1,
        "reduced_surface": f"Bl_{d * (d - j)} P^2",
        "number_of_old_exceptional_curves": d * (d - j),
        "curve_class": {
            "h_coefficient": d,
            "each_old_exceptional_coefficient": -1,
        },
        "curve_self_intersection": d * j,
        "fixed_curve_normal_degrees": {
            "positive_weight": d,
            "negative_weight": d * (j - 1),
        },
        "intersection_with_each_old_strict_triple_line": (
            0 if j < d else None
        ),
        "intersection_with_each_old_exceptional_curve": (
            1 if j < d else None
        ),
        "marked_ends": {
            "remaining_strict_triple_lines": d * (j - 2),
            "old_exceptional_curves": d * (d - j),
            "total": d * (d - 2),
        },
        "post_reduced_class": {
            "h": poly_record(alpha(d, j)),
            "coefficient_of_each_e_kp_by_k": {
                str(k): poly_record((0, 0, -(k - j), 0))
                for k in range(j + 1, d + 1)
            },
        },
        "DH_slopes": {
            "below_wall": {"h": d - 1, "each_old_exceptional": -1},
            "above_wall": {"h": -1, "each_old_exceptional": 0},
            "jump_above_minus_below": {
                "h": -d, "each_old_exceptional": 1
            },
        },
        "pair_details": pairs,
    }


def run(min_degree: int, max_degree: int) -> dict:
    if min_degree < 6 or max_degree < min_degree:
        raise ValueError("Use 6 <= min_degree <= max_degree.")
    audit = Audit()
    degree_summaries = []
    sextic_pairs: Dict[int, list[dict]] = {}
    total_pairs = 0
    for d in range(min_degree, max_degree + 1):
        pair_count = 0
        coefficient_basis_count = 0
        before = sum(audit.counts.values())
        for j in range(2, d + 1):
            for i in range(1, j):
                result = verify_pair(d, i, j, audit)
                pair_count += 1
                coefficient_basis_count += result["basis_size"]
                if d == 6:
                    sextic_pairs.setdefault(j, []).append(result)
        audit.require(
            pair_count == d * (d - 1) // 2,
            "all_ordered_pairs_covered",
            f"d={d}",
        )
        total_pairs += pair_count
        degree_summaries.append({
            "d": d,
            "pairs_checked": pair_count,
            "coefficient_basis_instances": coefficient_basis_count,
            "exact_equalities_checked": sum(audit.counts.values()) - before,
            "status": "passed",
        })
    return {
        "schema_version": 1,
        "status": "passed",
        "scope": (
            "Exact arithmetic consistency of the specified curve blow-up "
            "cohomology, Duistermaat-Heckman, intersection, and marking formulas."
        ),
        "limitations": [
            "This is not a proof or test of symplectic positivity.",
            "This does not establish a relative Moser comparison or a singular "
            "reduced family.",
            "This does not establish existence or gluing of Lagrangian fibrations.",
            "Curve-rank lemmas and relative framing obstruction theorems are "
            "not tested by arithmetic.",
        ],
        "arithmetic": {
            "polynomial_variables": ["constant", "a", "epsilon", "mu"],
            "coefficient_type": "exact integers",
            "intersection_matrix": "diag(1,-1,...,-1) in basis h,e_{k,p}",
            "exceptional_index_range": "j < k <= d, 1 <= p <= d",
            "floating_point_operations": 0,
        },
        "formulas": {
            "alpha_i": "a + d*(d+1-2*i)/2 * epsilon",
            "beta": "alpha_i*h - sum_{k>j,p}(k-i)*epsilon*e_{k,p}",
            "delta": "(j-i)*epsilon",
            "C": "d*h - sum_{k>j,p}e_{k,p}",
            "L": "h",
            "post_class": "alpha_j*h - sum_{k>j,p}(k-j)*epsilon*e_{k,p}",
            "DH_below": "beta - mu*L - (delta-mu)*C",
            "DH_above": "beta - mu*L",
            "old_strict_line": "h - sum_{p=1}^d e_{k,p}, k>j",
        },
        "degree_range": {"minimum": min_degree, "maximum": max_degree},
        "degrees_checked": max_degree - min_degree + 1,
        "pairs_checked": total_pairs,
        "exact_equalities_checked": sum(audit.counts.values()),
        "checks_by_category": dict(sorted(audit.counts.items())),
        "degree_summaries": degree_summaries,
        "sextic_by_j": {
            str(j): sextic_summary(j, pairs)
            for j, pairs in sorted(sextic_pairs.items())
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--min-degree", type=int, default=6)
    parser.add_argument("--max-degree", type=int, default=30)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).resolve().with_name("relative_model_results.json"),
    )
    args = parser.parse_args()
    result = run(args.min_degree, args.max_degree)
    args.output.write_text(
        json.dumps(result, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(
        "PASS: "
        f"d={args.min_degree}..{args.max_degree}; "
        f"{result['pairs_checked']} pairs; "
        f"{result['exact_equalities_checked']} exact equalities."
    )
    print(f"Results: {args.output}")
    print("Arithmetic consistency only; no symplectic positivity/existence claim.")


if __name__ == "__main__":
    main()
