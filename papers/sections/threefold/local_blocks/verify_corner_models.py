#!/usr/bin/env python3
"""Small reproducible checks of corner-model signs and normalizations.

Standard library only. Exact tests use rational samples, not symbolic proofs.
Pi is kept as a separate affine coefficient or absorbed into area radii;
the exact tests never replace pi by a rational approximation. Quadrature and
finite weight enumeration are supplementary checks, not existence or
nonexistence proofs. No global or relative gluing theorem is verified here.
"""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import json
import math


def graph_chart_checks():
    values = (F(0), F(1, 2), F(-2, 3))
    counts = {"b_nonzero": 0, "a_d_nonzero": 0, "a_c_nonzero": 0}
    for parameters in product(values, repeat=4):
        t, X, y, z = parameters
        chart_one = (t * X, y, z, t, X * y * z, X, F(1), X * y, F(1))
        x, T, Y, z = parameters
        chart_two = (x, T * Y, z, x * T, Y * z, F(1), T, Y, F(1))
        x, y, U, w = parameters
        chart_three = (x, y, U * w, x * y * U, w, F(1), y * U, F(1), U)
        for name, chart in zip(counts, (chart_one, chart_two, chart_three)):
            x, y, z, t, w, a, b, c, d = chart
            assert x * y * z == t * w
            assert x * b == t * a
            assert b * c == a * y * d
            assert c * z == d * w
            assert (a, b) != (0, 0) and (c, d) != (0, 0)
            counts[name] += 1
    return {
        "arithmetic": "Fraction; real rational coordinates, including zero strata",
        "samples_per_chart": counts,
        "equations_per_sample": 4,
    }


def moments_from_graph(x, y, z, a, b, c, d, e1, e2):
    """Use the two homogeneous graph coordinates; return pi-affine pairs."""
    A1 = e1 * a * a / (a * a + b * b)
    A2 = e2 * c * c / (c * c + d * d)
    return ((x * x - y * y, A1), (y * y - z * z, A2))


def central_moment_checks():
    """Compare coefficient pairs (coefficient of pi, constant term)."""
    values = (F(0), F(1, 2), F(-2, 3))
    area_pairs = ((F(2, 3), F(5, 4)), (F(3, 2), F(2, 5)))
    checked = {"Y_chart": 0, "U_chart": 0}
    for (e1, e2), (x, coordinate, transverse) in product(
        area_pairs, product(values, repeat=3)
    ):
        # D_5 in the middle graph chart: T = y = t = 0.
        Y, z = coordinate, transverse
        M1, M2 = moments_from_graph(x, F(0), z, F(1), F(0), Y, F(1), e1, e2)
        J = (z * z, e2 / (1 + Y * Y))
        assert M1 == (x * x, e1)
        assert M2 == (-J[0], e2 - J[1])
        checked["Y_chart"] += 1

        # D_5 in the last graph chart: y = 0; z = U w.
        # This includes U = 0, hence the actual strict triple divisor.
        U, w = coordinate, transverse
        z = U * w
        M1, M2 = moments_from_graph(x, F(0), z, F(1), F(0), F(1), U, e1, e2)
        J = (z * z, e2 * U * U / (1 + U * U))
        assert M1 == (x * x, e1)
        assert M2 == (-J[0], e2 - J[1])
        checked["U_chart"] += 1
    return {
        "arithmetic": "Fraction pairs for affine expressions in the exact symbol pi",
        "samples_per_chart": checked,
        "identities": ["M1 = e1 + rho_x", "M2 = e2 - J"],
    }


def surface_identity_checks():
    areas = (F(1, 2), F(2), F(7, 3))
    norms = (F(0), F(1, 4), F(3))
    radii = (F(0), F(2, 5), F(4, 3))
    fractions = (F(0), F(1, 4), F(3, 5))
    radial_count = slant_count = admissible = excluded = 0
    for e, u, rho_z in product(areas, norms, radii):
        J = rho_z + e / (1 + u)
        K = u * rho_z + e * u / (1 + u)
        a, b = rho_z * (1 + u), u * J
        j = J - e
        pi_R_squared = u * rho_z
        assert K == b and J == e + a - b
        assert (1 + u) * J == e + a
        # Inverse radial coordinates, and the squared reciprocal factor.
        assert b / J == u
        assert a * J / (e + a) == rho_z
        assert a * (a - j) == pi_R_squared * (e + a)
        assert a * b / (e + a) == pi_R_squared
        # Primitive coefficients in (theta_z, theta_Y).
        assert (rho_z + u * rho_z, K) == (a, b)
        assert (J - e + K, K) == (a, b)
        radial_count += 1

        for normal_fraction in fractions:
            scale = 1 - normal_fraction
            a_star, b_star = a / scale, b / scale
            J_star, K_star = e + a_star - b_star, b_star
            assert J_star == (J - e * normal_fraction) / scale
            assert K_star == K / scale
            assert scale * (J_star - e) == J - e
            assert scale * K_star == K
            assert J == scale * J_star + e * normal_fraction
            assert (J_star > 0) == (J > e * normal_fraction)
            if J_star > 0:
                admissible += 1
            else:
                excluded += 1
            slant_count += 1

    residue_count = 0
    for e, s, rho_w in product(areas, norms, radii):
        # U-chart: s = |U|^2, rho_z = s rho_w.
        J = s * (rho_w + e / (1 + s))
        K = rho_w + e / (1 + s)
        assert J + K - e == (1 + s) * rho_w
        assert (J - e) + e == s * (rho_w + e / (1 + s))
        if s == 0:
            assert J - e == -e
        # The formal target face pulls back to J = e mu / w_x.
        for normal_fraction in fractions:
            face_J = e * normal_fraction
            assert e + (face_J - e) / (1 - normal_fraction) == 0
        residue_count += 1

    # Weights of the same J-circle in the two blow-up charts.
    weight_Y, weight_z = -1, 1
    weight_U, weight_w = -weight_Y, weight_Y + weight_z
    assert (weight_U, weight_w) == (1, 0)
    return {
        "arithmetic": "Fraction; u and s are squared norms, rho variables are areas",
        "radial_and_primitive_samples": radial_count,
        "strict_slant_samples": slant_count,
        "slant_domain_admissible": admissible,
        "slant_domain_excluded_or_boundary": excluded,
        "residue_samples": residue_count,
        "positive_meridian_residue": "-e",
        "weights_at_p0_Y_z": [weight_Y, weight_z],
        "weights_at_p_infinity_U_w": [weight_U, weight_w],
    }


def radial_root_by_bisection(e, radius):
    """Solve a^2 = pi R^2 (e+a), without using its quadratic formula."""
    if radius == 0.0:
        return 0.0
    coefficient = math.pi * radius * radius
    lo, hi = 0.0, coefficient + math.sqrt(math.pi * e) * radius
    for _ in range(64):
        mid = (lo + hi) / 2
        if mid * mid < coefficient * (e + mid):
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def nodal_action_simpson(e, epsilon, steps=2048):
    """Use W = epsilon (1 + exp(i theta)), with theta in [-pi, pi].

    Away from its endpoint node, d(arg W)/dtheta = 1/2. The endpoint
    integrand has limit zero. Bisection supplies the radial action.
    """
    assert steps % 2 == 0
    step = 2 * math.pi / steps
    total = 0.0
    for index in range(steps + 1):
        theta = -math.pi + index * step
        if index in (0, steps):
            integrand = 0.0
        else:
            radius = epsilon * math.hypot(1 + math.cos(theta), math.sin(theta))
            integrand = radial_root_by_bisection(e, radius) / 2
        weight = 1 if index in (0, steps) else (4 if index % 2 else 2)
        total += weight * integrand
    return total * step / (6 * math.pi)


def nodal_action_checks():
    records = []
    for e, epsilon in ((1.0, 0.1), (0.3, 0.2), (2.5, 0.7), (0.8, 1.1)):
        closed = (
            math.pi * epsilon ** 2 / 2
            + epsilon * math.sqrt(e / math.pi)
            + (e / math.pi + epsilon ** 2)
            * math.atan(epsilon * math.sqrt(math.pi / e))
        )
        quadrature = nodal_action_simpson(e, epsilon)
        error = abs(closed - quadrature)
        assert error <= 2e-10 * max(1.0, abs(closed))
        frozen = 2 * epsilon * math.sqrt(e / math.pi)
        lower_bound = frozen + math.pi * epsilon ** 2 / 2
        assert closed > lower_bound
        records.append({
            "e": e,
            "epsilon": epsilon,
            "closed_action": closed,
            "simpson_action": quadrature,
            "absolute_error": error,
            "positive_remainder_after_quadratic_term": closed - lower_bound,
        })
    return {
        "method": "64-step monotone bisection plus 2048-panel composite Simpson quadrature",
        "tolerance": "2e-10 * max(1, absolute closed action)",
        "records": records,
        "is_analytic_proof": False,
    }


def determinant(u, v):
    return u[0] * v[1] - u[1] * v[0]


def finite_weight_checks():
    bound = 2
    vectors = [
        (a, b) for a, b in product(range(-bound, bound + 1), repeat=2)
        if (a, b) != (0, 0)
    ]
    satisfying = effective = 0
    indices = {}
    for x, y, z in product(vectors, repeat=3):
        dxy, dyz, dzx = determinant(x, y), determinant(y, z), determinant(z, x)
        if dxy == 0 or not (dxy == dyz == dzx):
            continue
        satisfying += 1
        assert (x[0] + y[0] + z[0], x[1] + y[1] + z[1]) == (0, 0)
        lattice_index = math.gcd(math.gcd(abs(dxy), abs(dyz)), abs(dzx))
        indices[lattice_index] = indices.get(lattice_index, 0) + 1
        if lattice_index == 1:
            effective += 1
            assert abs(dxy) == abs(dyz) == abs(dzx) == 1
    assert satisfying > effective > 0
    return {
        "weight_coordinate_range": [-bound, bound],
        "ordered_triples_examined": len(vectors) ** 3,
        "opposite_nonzero_stabilizer_weight_conditions_satisfied": satisfying,
        "effective_triples": effective,
        "counts_by_character_lattice_index": {
            str(index): indices[index] for index in sorted(indices)
        },
        "sum_zero_failures_in_this_finite_box": 0,
        "finite_enumeration_is_proof": False,
    }


def main():
    report = {
        "scope": (
            "Exact rational sample checks and small independent numerical checks of local "
            "corner formulas. These are not substitutes for the manuscript proofs, and "
            "do not verify global existence, nonexistence, or relative gluing."
        ),
        "graph_chart_equations": graph_chart_checks(),
        "central_FS_moments": central_moment_checks(),
        "surface_radial_slant_residue": surface_identity_checks(),
        "nodal_action": nodal_action_checks(),
        "torus_weight_enumeration": finite_weight_checks(),
        "global_fibration_proven": False,
        "status": "passed",
    }
    destination = Path(__file__).with_name("corner_model_checks.json")
    destination.write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print("PASS: corner graph/moments, surface identities, four nodal actions, finite weight checks.")


if __name__ == "__main__":
    main()
