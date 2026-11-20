#!/usr/bin/env python3
"""Finite consistency checks for upper caps, cycle covectors, and SU(2) gauges.

Only the Python standard library is used. These checks do not establish an
ambient embedding, a relative cap, or a global Lagrangian fibration.
"""

from __future__ import annotations

import cmath
from fractions import Fraction
import json
import math
from pathlib import Path
import sys


PI = math.pi
NODES = {}


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def close(actual, expected, *, absolute=1.0e-8, relative=2.0e-6, label=""):
    error = abs(actual - expected)
    limit = absolute + relative * max(abs(actual), abs(expected))
    require(error <= limit, f"{label}: {actual!r} != {expected!r}; error={error}")
    return error


def radial_maximum(h, j):
    u1 = (h + j + math.sqrt(h * h - h * j + j * j)) / 3.0
    f1 = math.sqrt((h - u1) * u1 * (u1 - j) / PI**3)
    return u1, f1


def trig_nodes(number):
    if number not in NODES:
        NODES[number] = [
            (math.cos(2.0 * PI * t / number), math.sin(2.0 * PI * t / number) ** 2)
            for t in range(number)
        ]
    return NODES[number]


def upper_action(h, j, c, q, number=256):
    """Integrate U+, U_h, U_j, U_q with h,j,c,q treated as independent."""
    u1, f1 = radial_maximum(h, j)
    require(c * c < q < (f1 - c) ** 2, "Parameters leave the upper regular band")
    totals = [0.0, 0.0, 0.0, 0.0]
    for cosine, sine_squared in trig_nodes(number):
        root = math.sqrt(q - c * c * sine_squared)
        radius = c * cosine + root
        target = PI**3 * radius * radius
        low, high = u1, h
        # G is strictly decreasing on this interval.
        for _ in range(57):
            u = (low + high) / 2.0
            g = (h - u) * u * (u - j)
            if g > target:
                low = u
            else:
                high = u
        u = (low + high) / 2.0
        gu = -3.0 * u * u + 2.0 * (h + j) * u - h * j
        gh = u * (u - j)
        gj = -(h - u) * u
        uh = -gh / gu
        uj = -gj / gu
        uq = PI**3 * radius / (root * gu)
        require(gu < 0.0 and uq < 0.0, "Wrong upper-branch sign")
        close(uh + uj - 1.0, -(h - u) * (u - j) / gu,
              absolute=2.0e-13, relative=2.0e-12, label="pointwise derivative identity")
        for index, value in enumerate((u, uh, uj, uq)):
            totals[index] += value
    return tuple(value / number for value in totals)


def solve_upper_q(h, j, c, target_action, initial_q, number=256):
    """Safeguarded Newton solve of I(h,j,q)=target_action."""
    _, f1 = radial_maximum(h, j)
    q_low = c * c
    q_high = (f1 - c) ** 2
    margin = (q_high - q_low) * 1.0e-10
    low, high = q_low + margin, q_high - margin
    q = min(max(initial_q, low), high)
    for _ in range(14):
        value, _, _, iq = upper_action(h, j, c, q, number)
        error = value - target_action
        if abs(error) < 2.0e-13:
            return q
        if error > 0.0:  # I_q<0, so the root is to the right.
            low = q
        else:
            high = q
        candidate = q - error / iq
        q = candidate if low < candidate < high else (low + high) / 2.0
    raise AssertionError("The upper inverse did not converge")


def check_upper_geometry():
    h = 3.0
    examples = []
    maxima = {name: 0.0 for name in ("I_h", "I_j", "I_q", "cut_jacobian", "fixed_A_q_M")}
    sign_margins = {"minus_I_q": math.inf, "cut_jacobian": math.inf, "fixed_A_q_M": math.inf}
    for j in (-0.1, 0.0, 0.1):
        _, f1 = radial_maximum(h, j)
        for ratio in (0.08, 0.2, 0.4):
            c = ratio * f1
            for fraction in (0.2, 0.5, 0.8):
                q = c * c + fraction * ((f1 - c) ** 2 - c * c)
                action, ih, ij, iq = upper_action(h, j, c, q)
                dh, dj, dq = 2.0e-5 * h, 2.0e-5, 1.0e-6 * f1 * f1
                ih_fd = (upper_action(h + dh, j, c, q)[0] - upper_action(h - dh, j, c, q)[0]) / (2.0 * dh)
                ij_fd = (upper_action(h, j + dj, c, q)[0] - upper_action(h, j - dj, c, q)[0]) / (2.0 * dj)
                iq_fd = (upper_action(h, j, c, q + dq)[0] - upper_action(h, j, c, q - dq)[0]) / (2.0 * dq)
                determinant = ih + ij - 1.0
                determinant_fd = (1.0 - ih_fd) * (ij_fd - 1.0) + ij_fd * ih_fd
                q_m = (1.0 - ih - ij) / iq
                ax, az = h - action, action - j
                dm = 2.0e-5
                plus = solve_upper_q(ax + action + dm, action + dm - az, c, action + dm, q)
                minus = solve_upper_q(ax + action - dm, action - dm - az, c, action - dm, q)
                q_m_fd = (plus - minus) / (2.0 * dm)
                require(iq < 0.0 and determinant > 0.0 and q_m > 0.0, "Projection/orientation sign failed")
                require(ax > 0.0 and az > 0.0, "A noncollapsing toric action is nonpositive")
                values = (("I_h", ih, ih_fd), ("I_j", ij, ij_fd), ("I_q", iq, iq_fd),
                          ("cut_jacobian", determinant, determinant_fd), ("fixed_A_q_M", q_m, q_m_fd))
                errors = {}
                for name, analytic, numerical in values:
                    error = close(analytic, numerical, absolute=2.0e-8, relative=3.0e-6, label=name)
                    errors[name] = error
                    maxima[name] = max(maxima[name], error)
                sign_margins["minus_I_q"] = min(sign_margins["minus_I_q"], -iq)
                sign_margins["cut_jacobian"] = min(sign_margins["cut_jacobian"], determinant)
                sign_margins["fixed_A_q_M"] = min(sign_margins["fixed_A_q_M"], q_m)
                examples.append({"h": h, "j": j, "c_over_F": ratio, "q_fraction": fraction,
                                 "F": f1, "c": c, "q": q, "I": action,
                                 "I_h": ih, "I_j": ij, "I_q": iq,
                                 "A_x": ax, "A_z": az, "cut_jacobian": determinant,
                                 "fixed_A_q_M": q_m, "finite_difference_errors": errors})
    # Higher quadrature resolution checks the three c/F extremes at the most
    # difficult selected q-level, independently of the finite differences.
    refinement = []
    _, f1 = radial_maximum(h, 0.0)
    for ratio in (0.08, 0.2, 0.4):
        c = ratio * f1
        q = c * c + 0.8 * ((f1 - c) ** 2 - c * c)
        coarse = upper_action(h, 0.0, c, q, 256)
        fine = upper_action(h, 0.0, c, q, 512)
        errors = [close(x, y, absolute=1.0e-10, relative=1.0e-10, label="quadrature refinement") for x, y in zip(coarse, fine)]
        refinement.append({"c_over_F": ratio, "component_errors": errors})
    return {"examples_checked": len(examples), "quadrature_nodes": 256,
            "derivatives_hold_c_and_q_fixed": True,
            "maximum_absolute_finite_difference_errors": maxima,
            "minimum_positive_sign_margins": sign_margins,
            "quadrature_refinement": refinement, "examples": examples}


def check_cycle_lattice():
    sphere_checks = 0
    threshold_checks = 0
    raw_checks = 0
    a = Fraction(1)
    for ab in (Fraction(1, 10), Fraction(1, 3), Fraction(3, 4)):
        for portion in (Fraction(1, 4), Fraction(1, 2), Fraction(4, 5)):
            a0 = ab * portion
            delta = ab - a0
            for k in (1, 2):
                for p in range(-3, 4):
                    for ell in range(-3, 4):
                        for r in range(-2, 3):
                            local_area = p * a0 + ell * ab
                            local_mu = 2 * ((k - 1) * p + k * ell)
                            invariant = local_area + a * Fraction(local_mu, 2)
                            for n in range(-2, 3):
                                area = local_area + a * n
                                mu = local_mu - 2 * n
                                require(area + a * Fraction(mu, 2) == invariant, "Sphere cancellation failed")
                                sphere_checks += 1
                            forced_n = (k - 1) * p + k * ell - 1
                            require(local_mu - 2 * forced_n == 2, "Wrong forced sphere degree")
                            expression = p * (a0 + a * (k - 1)) + ell * (ab + a * k) - a
                            require(local_area + a * forced_n == expression == invariant - a,
                                    "Maslov-two area expressions disagree")
                for p in (-1, 1):
                    lower = 0 if (k == 2 and p == 1) else 1
                    for ell in range(-4, 5):
                        for r in range(-2, 3):
                            area = p * (a0 + a * (k - 1)) + ell * (ab + a * k) - a
                            require((area > 0) == (ell >= lower), "Integer threshold failed")
                            require(math.gcd(math.gcd(abs(p), abs(ell)), abs(r)) == 1, "Complement not primitive")
                            threshold_checks += 1
                # raw upper meridian: p=-1, ell=1, r=0
                require(2 * (-(k - 1) + k) == 2, "Raw cycle Maslov index failed")
                require(-a0 + ab == delta > 0, "Raw cycle area failed")
                require(-(k - 1) + k - 1 == 0, "Raw cycle sphere degree failed")
                # Same-orientation, zero-base-degree tube candidate.
                require(2 * ((k - 1) - k) == -2, "Tube candidate Maslov index failed")
                require(a0 - ab == -delta < 0, "Tube candidate area failed")
                raw_checks += 1
            require(a0 - a < 0, "k=1 natural/only-J-twisted cycle must be excluded")
            require(a0 + ab > 0, "k=1 minimal positive-orientation candidate failed")
    return {"arithmetic": "fractions.Fraction; exact comparisons",
            "area_pairs": 9, "sphere_cancellation_checks": sphere_checks,
            "unimodular_threshold_checks": threshold_checks, "raw_cycle_and_tube_checks": raw_checks,
            "thresholds": [
                {"k": 1, "p": 1, "minimum_ell": 1},
                {"k": 1, "p": -1, "minimum_ell": 1},
                {"k": 2, "p": 1, "minimum_ell": 0},
                {"k": 2, "p": -1, "minimum_ell": 1}],
            "j": 0, "r_does_not_affect_area_or_index_on_this_slice": True}


def matrix_product(a, b):
    return tuple(tuple(sum(a[i][k] * b[k][j] for k in range(2)) for j in range(2)) for i in range(2))


def adjoint(a):
    return tuple(tuple(a[j][i].conjugate() for j in range(2)) for i in range(2))


def apply_matrix(a, z):
    return tuple(sum(a[i][j] * z[j] for j in range(2)) for i in range(2))


def inner(z, w):
    return sum(x.conjugate() * y for x, y in zip(z, w))


def matrix_distance(a, b):
    return max(abs(a[i][j] - b[i][j]) for i in range(2) for j in range(2))


def homotopy(m, time, theta):
    """A continuous piecewise smooth periodic contraction; flatten time at joins."""
    zero = ((0j, 0j), (0j, 0j))
    if time <= 0.5:
        angle = PI * time
        cosine, sine = math.cos(angle), math.sin(angle)
        return ((complex(cosine), complex(-sine)), (complex(sine), complex(cosine))), zero, 0.0
    xi = PI * (time - 0.5)
    sine, cosine = math.sin(xi), math.cos(xi)
    phase = cmath.exp(2j * PI * m * theta)
    matrix = ((sine * phase, complex(-cosine)), (complex(cosine), sine * phase.conjugate()))
    derivative = ((2j * PI * m * sine * phase, 0j), (0j, -2j * PI * m * sine * phase.conjugate()))
    return matrix, derivative, sine


def check_su2_gauge():
    identity = ((1 + 0j, 0j), (0j, 1 + 0j))
    maxima = {name: 0.0 for name in ("unitarity", "determinant", "periodicity", "theta_derivative", "antihermitian", "generator_square", "norm_preservation", "primitive_identity", "endpoint_shift")}
    matrices = 0
    vectors = 0
    endpoints = 0
    max_bound_ratio = 0.0
    test_vectors = [(1 + 0j, 0j), (0j, 1 + 0j),
                    (complex(1 / math.sqrt(2)), 1j / math.sqrt(2)),
                    (0.3 + 0.4j, -0.2 + 0.7j), (0j, 0j)]
    tangent = (0.2 - 0.3j, 0.4 + 0.1j)
    theta_speed = 0.37
    moment = 0.23
    for m in (-2, -1, 0, 1, 2):
        for time in (0.0, 0.125, 0.25, 0.5, 0.625, 0.75, 0.875, 1.0):
            for theta in (0.0, 0.137, 0.25, 0.5, 0.917, 1.0):
                g, dg, sine = homotopy(m, time, theta)
                star_g = adjoint(g)
                maxima["unitarity"] = max(maxima["unitarity"], matrix_distance(matrix_product(star_g, g), identity))
                determinant = g[0][0] * g[1][1] - g[0][1] * g[1][0]
                maxima["determinant"] = max(maxima["determinant"], abs(determinant - 1))
                maxima["periodicity"] = max(maxima["periodicity"], matrix_distance(g, homotopy(m, time, theta + 1)[0]))
                step = 1.0e-6
                plus = homotopy(m, time, theta + step)[0]
                minus = homotopy(m, time, theta - step)[0]
                fd = tuple(tuple((plus[i][j] - minus[i][j]) / (2 * step) for j in range(2)) for i in range(2))
                maxima["theta_derivative"] = max(maxima["theta_derivative"], matrix_distance(dg, fd))
                logarithmic = matrix_product(star_g, dg)
                skew_error = max(abs(logarithmic[i][j] + logarithmic[j][i].conjugate()) for i in range(2) for j in range(2))
                maxima["antihermitian"] = max(maxima["antihermitian"], skew_error)
                generator = tuple(tuple((-1j / (2 * PI)) * logarithmic[i][j] for j in range(2)) for i in range(2))
                square = matrix_product(generator, generator)
                target = tuple(tuple((m * m * sine * sine if i == j else 0.0) + 0j for j in range(2)) for i in range(2))
                maxima["generator_square"] = max(maxima["generator_square"], matrix_distance(square, target))
                matrices += 1
                if time == 0.0:
                    require(matrix_distance(g, identity) < 1.0e-14, "Homotopy does not start at identity")
                if time == 1.0:
                    phase = cmath.exp(2j * PI * m * theta)
                    require(matrix_distance(g, ((phase, 0j), (0j, phase.conjugate()))) < 1.0e-14,
                            "Wrong endpoint gauge")
                for z in test_vectors:
                    norm = PI * inner(z, z).real
                    gz = apply_matrix(g, z)
                    maxima["norm_preservation"] = max(maxima["norm_preservation"], abs(PI * inner(gz, gz).real - norm))
                    correction = -0.5 * inner(z, apply_matrix(logarithmic, z)).imag
                    bound = abs(m) * norm
                    require(abs(correction) <= bound + 2.0e-13, "Gauge action shift exceeds |m|N")
                    if bound > 0:
                        max_bound_ratio = max(max_bound_ratio, abs(correction) / bound)
                    gw = apply_matrix(g, tangent)
                    dgz = apply_matrix(dg, z)
                    image_tangent = tuple(gw[i] + dgz[i] * theta_speed for i in range(2))
                    source_primitive = moment * theta_speed + 0.5 * inner(z, tangent).imag
                    image_primitive = (moment + correction) * theta_speed + 0.5 * inner(gz, image_tangent).imag
                    maxima["primitive_identity"] = max(maxima["primitive_identity"], abs(image_primitive - source_primitive))
                    if time == 1.0:
                        expected = -m * PI * (abs(z[0]) ** 2 - abs(z[1]) ** 2)
                        maxima["endpoint_shift"] = max(maxima["endpoint_shift"], abs(correction - expected))
                        endpoints += 1
                    vectors += 1
    # Check the join from both formulas, not just the piecewise dispatch.
    rotation = ((0j, -1 + 0j), (1 + 0j, 0j))
    for m in (-2, -1, 0, 1, 2):
        for theta in (0.0, 0.17, 0.5, 1.0):
            phase = cmath.exp(2j * PI * m * theta)
            second_start = ((0 * phase, -1 + 0j), (1 + 0j, 0 * phase.conjugate()))
            require(matrix_distance(rotation, second_start) == 0.0, "Homotopy join is not continuous")
    for name, error in maxima.items():
        tolerance = 5.0e-8 if name == "theta_derivative" else 2.0e-12
        require(error < tolerance, f"SU(2) check {name} failed: {error}")
    return {"windings": [-2, -1, 0, 1, 2], "matrix_samples": matrices,
            "vector_and_lift_samples": vectors, "endpoint_shift_samples": endpoints,
            "maximum_errors": maxima, "maximum_shift_to_bound_ratio": max_bound_ratio,
            "normal_primitive": "(1/2) Im <z,dz>; rho_j=pi |z_j|^2",
            "base_action_shift": "Delta I=-(1/2) Im <z,G* partial_theta G z>",
            "endpoint_shift": "-m(rho_1-rho_2)",
            "bound": "|Delta I| <= |m| (rho_1+rho_2)",
            "time_smoothing": "The displayed path is continuous and piecewise smooth; smooth time flattening is an analytic step, not numerically verified."}


def main():
    report = {"status": "passed", "scope": "Finite consistency checks only; not a proof of cap embeddings, relative gluings, or a global fibration.",
              "upper_projection_and_orientation": check_upper_geometry(),
              "maslov_area_lattice": check_cycle_lattice(),
              "su2_periodic_gauge": check_su2_gauge()}
    path = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).with_name("cap_extension_checks.json")
    path.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print("PASS: 27 upper-action cases, exact rational lattice identities, and periodic SU(2) gauges.")
    print("Upper finite-difference maximum errors:", report["upper_projection_and_orientation"]["maximum_absolute_finite_difference_errors"])
    print("SU(2) maximum errors:", report["su2_periodic_gauge"]["maximum_errors"])
    print("Report:", path.resolve())


if __name__ == "__main__":
    main()
