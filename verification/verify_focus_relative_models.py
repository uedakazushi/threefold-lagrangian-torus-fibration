#!/usr/bin/env python3
"""Exact finite checks for the relative focus models (2026-11-21).

This checks algebraic consequences of the geometric statements in the paper.
It does not establish the round-handle decomposition, relative Hurewicz input,
or any symplectic/Lagrangian construction.  Standard library only; no files
are written.  Run with --json for a machine-readable report on stdout.
"""

from __future__ import annotations

import argparse
import itertools
import json


Laurent = dict[int, int]
Betti = dict[int, int]


def clean(coefficients: Laurent) -> Laurent:
    return {power: value for power, value in coefficients.items() if value}


def add(left: Laurent, right: Laurent) -> Laurent:
    result = dict(left)
    for power, value in right.items():
        result[power] = result.get(power, 0) + value
    return clean(result)


def scale(poly: Laurent, scalar: int) -> Laurent:
    return clean({power: scalar * value for power, value in poly.items()})


def shift(poly: Laurent, amount: int) -> Laurent:
    return {power + amount: value for power, value in poly.items()}


def multiply(left: Laurent, right: Laurent) -> Laurent:
    result: Laurent = {}
    for p, a in left.items():
        for q, b in right.items():
            result[p + q] = result.get(p + q, 0) + a * b
    return clean(result)


def boundary(poly: Laurent) -> tuple[int, int]:
    """Coefficients of (a,b): delta(f B) = f'(1) a + f(1) b."""
    return sum(p * c for p, c in poly.items()), sum(poly.values())


def divide_square_difference(poly: Laurent) -> Laurent:
    """Exact Laurent division by (tau-1)^2 via monic long division.

    Raises AssertionError if the input is not divisible over Z[tau,tau^-1].
    Negative powers are first cleared; no rational or floating-point values
    enter the calculation.
    """
    if not poly:
        return {}
    offset = max(0, -min(poly))
    remainder = shift(poly, offset)
    quotient: Laurent = {}
    divisor = {0: 1, 1: -2, 2: 1}
    while remainder and max(remainder) >= 2:
        highest = max(remainder)
        coefficient = remainder[highest]
        term = {highest - 2: coefficient}
        quotient = add(quotient, term)
        remainder = add(remainder, scale(multiply(term, divisor), -1))
    assert not remainder, (poly, remainder)
    result = shift(quotient, -offset)
    assert multiply(result, divisor) == clean(poly)
    return result


def convolve(left: Betti, right: Betti) -> Betti:
    """Ranks from the torsion-free Kunneth direct sum."""
    answer: Betti = {}
    for i, a in left.items():
        for j, b in right.items():
            answer[i + j] = answer.get(i + j, 0) + a * b
    return clean(answer)


def check_laurent_model() -> dict:
    exponents = tuple(range(-2, 3))
    samples = kernel_samples = action_checks = 0
    square_difference = {0: 1, 1: -2, 2: 1}

    # These two classes give the independently visible a and b generators.
    assert boundary({0: 1}) == (0, 1)
    assert boundary({0: -1, 1: 1}) == (1, 0)
    assert boundary(square_difference) == (0, 0)

    for coefficients in itertools.product((-1, 0, 1), repeat=len(exponents)):
        poly = clean(dict(zip(exponents, coefficients)))
        a, b = boundary(poly)
        # Independent Euclidean division checks the evaluation/derivative
        # description of the quotient Lambda/(tau-1)^2.
        canonical_remainder = clean({0: b - a, 1: a})
        assert boundary(canonical_remainder) == (a, b)
        residual = add(poly, scale(canonical_remainder, -1))
        quotient = divide_square_difference(residual)
        assert add(multiply(square_difference, quotient), canonical_remainder) == poly

        for n in range(-3, 4):
            assert boundary(shift(poly, n)) == (a + n * b, b)
            action_checks += 1
        if (a, b) == (0, 0):
            kernel_samples += 1
            # For a kernel sample supported in [-2,2], the quotient must be
            # supported in [-2,0]; these are the three obvious kernel shifts.
            assert all(-2 <= p <= 0 for p in quotient)
        samples += 1

    assert samples == 3**5
    return {
        "coefficient_values": [-1, 0, 1],
        "exponent_range": [-2, 2],
        "samples_checked": samples,
        "kernel_samples": kernel_samples,
        "deck_action_checks": action_checks,
        "deck_power_range": [-3, 3],
        "exact_divisions_after_canonical_remainder": samples,
        "boundary_convention": "(a,b)=(f'(1),f(1)); tau acts by (a+b,b)",
    }


def check_product_pairs() -> dict:
    closed_cases = punctured_cases = 0
    sextic: dict = {}
    # These are topological inputs proved in the manuscript, not computed
    # geometric facts: H_*(E,S^3) and H_*(E,Y), respectively.
    normal_out = {1: 1, 2: 1}
    normal_in = {2: 1, 3: 1}

    for d in range(6, 21):
        g = (d - 3) * (d - 4) // 2
        curve = {0: 1, 1: 2 * g, 2: 1}
        from_outer = convolve(curve, normal_out)
        from_inner = convolve(curve, normal_in)
        expected_outer = {1: 1, 2: 2 * g + 1, 3: 2 * g + 1, 4: 1}
        expected_inner = {2: 1, 3: 2 * g + 1, 4: 2 * g + 1, 5: 1}
        assert from_outer == expected_outer
        assert from_inner == expected_inner
        # Separate duality and product computations must agree in dimension6.
        assert from_inner == {6 - k: rank for k, rank in from_outer.items()}
        assert from_outer[4] == 1 and from_inner[5] == 1
        assert sum((-1) ** k * rank for k, rank in from_outer.items()) == 0
        closed_cases += 1

        for punctures in range(1, 5):
            # A connected orientable genus-g surface with r boundary circles
            # has the homology of a graph with 2g+r-1 independent cycles.
            graph = {0: 1, 1: 2 * g + punctures - 1}
            punctured = convolve(graph, normal_out)
            expected = {1: 1, 2: 2 * g + punctures, 3: 2 * g + punctures - 1}
            assert punctured == expected
            assert all(k <= 3 for k in punctured)
            punctured_cases += 1
            if d == 6 and punctures == 1:
                sextic["one_boundary_circle_outer_relative_betti"] = punctured

        if d == 6:
            sextic["closed_outer_relative_betti"] = from_outer
            sextic["closed_inner_relative_betti"] = from_inner

    return {
        "degree_range": [6, 20],
        "closed_curve_cases": closed_cases,
        "punctured_curve_cases": punctured_cases,
        "puncture_range": [1, 4],
        "sextic_examples": sextic,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", action="store_true", help="emit only JSON on stdout")
    args = parser.parse_args()
    report = {
        "status": "passed",
        "arithmetic": "exact integer; Python standard library only",
        "scope": "Finite algebra checks of stated models; not proofs of geometric identifications.",
        "laurent_boundary": check_laurent_model(),
        "relative_product_betti": check_product_pairs(),
        "not_checked": [
            "round-handle decomposition or attaching framings",
            "relative Hurewicz hypotheses and crossed-module topology",
            "symplectic or Lagrangian disk/handle realization",
            "existence of a Weinstein structure after puncturing",
            "global Lagrangian fibration or boundary-free completion",
        ],
    }
    if args.json:
        print(json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True))
    else:
        laurent = report["laurent_boundary"]
        pairs = report["relative_product_betti"]
        print("PASS: focus relative models (exact integer arithmetic)")
        print(f"Laurent samples and exact divisions: {laurent['samples_checked']}")
        print(f"Deck-action checks: {laurent['deck_action_checks']}")
        print(f"Closed-curve product/duality cases: {pairs['closed_curve_cases']}")
        print(f"Punctured-curve product cases: {pairs['punctured_curve_cases']}")
        print(report["scope"])
        print("No geometric identification, framing, or symplectic existence is verified here.")


if __name__ == "__main__":
    main()
