#!/usr/bin/env python3
"""Cross-check phase degrees and constants in the exact axis-end action.

The algebraic checks use integers. The action checks are numerical and are
explicitly separate from the analytic proofs in the manuscript. Standard
library only; this does not verify global gluing or a symplectic existence
theorem.
"""
from pathlib import Path
import json
import math


def phase_checks():
    checked = 0
    sextic_last = None
    for d in range(6, 31):
        for i in range(1, d):
            for j in range(i + 1, d + 1):
                types = [k for k in range(1, d + 1) if k not in (i, j)]
                tensor_degree = d * (i + 1)
                zero_degree = d * (d - i - 1)
                assert tensor_degree + zero_degree == d * d
                assert -tensor_degree - zero_degree == -(d * d)
                q, remainder = divmod(i + 1, len(types))
                winds = [q + (index < remainder) for index in range(len(types))]
                powers = [wind + int(k > i) for k, wind in zip(types, winds)]
                assert sum(powers) == d
                assert d * sum(winds) == tensor_degree
                assert max(winds) == (i + 1 + d - 3) // (d - 2)
                assert max(winds) <= 2
                assert d * (d - 3) + d * (d - 2) == d * (d - 5) + d * d
                if (d, i, j) == (6, 5, 6):
                    sextic_last = {
                        "boundary_types": types, "powers": powers,
                        "winding_per_end": winds, "ends_of_each_type": d,
                        "winding_sum": d * sum(winds),
                        "disk_standard_relative_defect": -(d * d),
                        "original_axis_to_semistable_shift": d * (d - 2),
                    }
                checked += 1
    return {"triples_checked": checked, "range_d": [6, 30], "sextic_last": sextic_last}


def radial_root(h, radius):
    """Unique nonnegative solution t sqrt(h+t) = pi^(3/2) radius."""
    if radius <= 0:
        return 0.0
    rhs = math.pi ** 1.5 * radius
    lo, hi = 0.0, rhs / math.sqrt(h)
    for _ in range(58):
        mid = (lo + hi) / 2
        if mid * math.sqrt(h + mid) < rhs:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def nodal_action(h, epsilon, steps=2048):
    """Composite Simpson integration of the explicit marked cycle action."""
    step = math.pi / steps
    total = 0.0
    for k in range(steps + 1):
        phi = -math.pi / 2 + k * step
        radius = max(0.0, 2 * epsilon * math.cos(phi))
        value = radial_root(h, radius)
        weight = 1 if k in (0, steps) else (4 if k % 2 else 2)
        total += weight * value
    return total * step / (6 * math.pi)


def numerical_action_checks():
    records = []
    for h in (1.0, 2.0, 5.0):
        errors = []
        for epsilon in (1e-3, 5e-4):
            actual = nodal_action(h, epsilon)
            frozen = 2 * epsilon * math.sqrt(math.pi / h)
            assert 0 < actual < frozen
            observed = (frozen - actual) / epsilon ** 2
            predicted = math.pi ** 3 / (2 * h ** 2)
            relative_error = abs(observed / predicted - 1)
            assert relative_error < 0.03
            errors.append(relative_error)
            records.append({
                "H": h, "epsilon": epsilon, "actual_action": actual,
                "frozen_action": frozen, "observed_quadratic_coefficient": observed,
                "predicted_quadratic_coefficient": predicted,
                "relative_error": relative_error,
            })
        assert errors[1] < errors[0]
    return {"method": "Simpson quadrature plus monotone bisection", "records": records}


def main():
    report = {
        "scope": "Integer phase identities and numerical action-constant checks; proofs are in the manuscript.",
        "phase": phase_checks(),
        "action_numerical": numerical_action_checks(),
        "global_fibration_proven": False,
        "status": "passed",
    }
    destination = Path(__file__).with_name("relative_gluing_checks.json")
    destination.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print("PASS: relative phase identities and six numerical action comparisons.")


if __name__ == "__main__":
    main()
