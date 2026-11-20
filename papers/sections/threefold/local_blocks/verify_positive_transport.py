#!/usr/bin/env python3
"""Independent algebraic samples for toric transport and annular insertion.

These checks supplement the proofs. They do not certify smooth existence,
global gluing, or the positivity of an arbitrary prescribed global form.
Only the Python standard library is used.
"""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import json
import math


def determinant3(m):
    return (m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]))


def lower_root(h, radius):
    left, right = 0.0, 2.0 * h / 3.0
    for _ in range(100):
        mid = (left + right) / 2
        value = mid * math.sqrt(h - mid) / math.pi**1.5
        if value < radius:
            left = mid
        else:
            right = mid
    return (left + right) / 2


def main():
    # pi is an arbitrary positive scalar in these rational identities.
    # Using a rational surrogate checks exact identities, not pi numerically.
    pi = F(22, 7)
    values = (F(1, 5), F(2, 3), F(7, 4))
    e1, e2 = F(2, 5), F(3, 7)
    toric_count = 0
    for rx, ry, rz, a in product(values, repeat=4):
        pw = pi * rx * ry * rz / a
        A1 = e1 * rx / (rx + a)
        A2 = e2 * rx * ry / (rx * ry + a)
        P = (pi * rx + pw + A1 + A2,
             pi * ry + pw + A2, pi * rz + pw)
        q = (1 + ry * rz / a + e1 / (pi * (rx + a))
             + e2 * ry / (pi * (rx * ry + a)),
             1 + rx * rz / a + e2 * rx / (pi * (rx * ry + a)),
             1 + rx * ry / a)
        assert all(qi > 0 for qi in q)
        assert tuple(pi * r * qi for r, qi in zip((rx, ry, rz), q)) == P
        hessian = [[pw for _ in range(3)] for _ in range(3)]
        for i, r in enumerate((rx, ry, rz)):
            hessian[i][i] += pi * r
        hessian[0][0] += e1 * a * rx / (rx + a)**2
        cross = e2 * a * rx * ry / (rx * ry + a)**2
        for i in (0, 1):
            for j in (0, 1):
                hessian[i][j] += cross
        assert hessian[0][0] > 0
        assert hessian[0][0] * hessian[1][1] - hessian[0][1]**2 > 0
        assert determinant3(hessian) > 0
        toric_count += 1

    chart_count = 0
    for r, v, z2, a in product(values, repeat=4):
        ry = a * v / r
        pw = pi * r * ry * z2 / a
        A1 = e1 * r / (r + a)
        A2 = e2 * r * ry / (r * ry + a)
        original = (pi * r + pw + A1 + A2,
                    pi * ry + pw + A2, pi * z2 + pw)
        semistable = (pi * r + pi * v * z2 + e1 * r / (r + a)
                      + e2 * v / (1 + v),
                      pi * a * v / r + pi * v * z2 + e2 * v / (1 + v),
                      pi * z2 * (1 + v))
        assert original == semistable
        xi = pi * z2 * (1 + v)
        zeta = pi * v * z2 + e2 * v / (1 + v)
        central = (pi * r + e1 + zeta, zeta, xi)
        assert central[0] - central[1] == pi * r + e1
        assert central[1] - central[2] < e2
        # (theta_1,theta_2,theta_3)=(theta_x,theta_zeta-theta_x,theta_xi).
        assert (central[0] - central[1], central[1], central[2]) == (
            pi * r + e1, zeta, xi)
        # B=e^{-2i arg x} zeta, H=rho_x+e1+2rho_zeta.
        H = pi * r + e1 + 2 * zeta
        assert H - 2 * zeta == pi * r + e1
        assert H - zeta == central[0]
        chart_count += 1

    insertion_count = 0
    for k, scale, rx, ry, rz in product(
            range(-3, 5), (F(1, 100), F(1, 17)), values, values, values):
        n, b = k - 1, F(3, 2)
        I, u, v = b + scale * (rx - n * ry), scale * ry, scale * rz
        assert I > 0
        assert I + n * u == b + scale * rx
        assert I - b + k * u == scale * (rx + ry)
        assert u - v == scale * (ry - rz)
        assert I - b + n * u == scale * rx > 0
        assert -k + n == -1
        insertion_count += 1

    chambers = []
    for h, E1, c in ((3.0, 1.0, 0.02), (1.1, 0.2, 0.003),
                      (4.0, 1.5, 0.03), (0.8, 0.1, 0.001)):
        b = (h - E1) / 2
        fb = b * math.sqrt(h - b) / math.pi**1.5
        maximum = 2 * h**1.5 / (3 * math.sqrt(3) * math.pi**1.5)
        assert 2 * c < min(fb, maximum - fb)
        focus_u = lower_root(h, 2 * c)
        saddle_u = lower_root(h, maximum - 2 * c)
        assert focus_u < b < saddle_u
        wall_min, wall_max = (fb - c)**2, (fb + c)**2
        qs = (maximum - c)**2
        assert c*c < wall_min < wall_max < qs
        chambers.append(dict(h=h, e1=E1, c=c, wall_u=b,
                             focus_u_max=focus_u, saddle_u_min=saddle_u,
                             focus_q=c*c, wall_q=[wall_min, wall_max], saddle_q=qs))

    saddle_count = 0
    for h, j in product((0.5, 1.0, 3.0), (-3.0, -0.1, 0.0, 0.1, 0.4)):
        if h <= max(0.0, j):
            continue
        u = (h + j + math.sqrt(h*h - h*j + j*j)) / 3
        assert h - u < u and h - u < u - j
        saddle_count += 1

    distributed_count = 0
    for alpha, delta, scale, rx, ry, rz in product(
            range(-2, 4), range(-2, 4), (F(1, 100), F(1, 17)),
            values, values, values):
        b = F(-2, 3)
        k = alpha + delta
        u, v = scale * ry, scale * rz
        I = b + scale * (rx - (alpha - 1) * ry - delta * rz)
        assert I + (alpha - 1) * u + delta * v == b + scale * rx
        assert I - b + alpha * u + delta * v == scale * (rx + ry)
        assert I - b + (alpha - 1) * u + delta * v == scale * rx > 0
        assert -k + (alpha - 1) + delta == -1
        distributed_count += 1
    twists = [(1, 1), (0, 2), (0, 1), (0, 1)]
    degree_plus = 6 * sum(a for a, d in twists)
    degree_minus = 6 * sum(d for a, d in twists)
    assert (degree_plus, degree_minus) == (6, 30)
    assert [a+d for a, d in twists] == [2, 2, 1, 1]
    assert degree_plus + degree_minus == 36

    result = {
        "scope": "Rational identities and selected numerical inequalities; proofs are in the manuscript.",
        "toric_radial_and_convex_hessian_samples": toric_count,
        "semistable_and_central_chart_samples": chart_count,
        "annular_primitive_action_and_phase_samples": insertion_count,
        "annular_windings_checked": list(range(-3, 5)),
        "distributed_normal_twist_samples": distributed_count,
        "C56_twists_by_end_type": twists,
        "C56_individual_degrees": [degree_plus, degree_minus],
        "C56_total_winding": degree_plus + degree_minus,
        "whole_fiber_chamber_examples": chambers,
        "saddle_chamber_samples": saddle_count,
        "status": "passed",
    }
    path = Path(__file__).with_name("positive_transport_checks.json")
    path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: v for k, v in result.items()
                      if k not in ("whole_fiber_chamber_examples",)}, indent=2))


if __name__ == "__main__":
    main()
