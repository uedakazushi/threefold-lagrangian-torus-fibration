#!/usr/bin/env python3
"""Finite checks for the closed-curve shell, delta budget, and membrane families.

Only the Python standard library is used.  The exact arithmetic and selected
finite-difference/quadrature checks below do not prove geometric existence,
global embeddings, collar compatibility, or a boundary-free fibration.  Those
claims and their hypotheses must be read in the accompanying manuscript.
"""
from collections import deque
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import cmath
import json
import math


def determinant3(a):
    return (a[0][0] * (a[1][1]*a[2][2]-a[1][2]*a[2][1])
            - a[0][1] * (a[1][0]*a[2][2]-a[1][2]*a[2][0])
            + a[0][2] * (a[1][0]*a[2][1]-a[1][1]*a[2][0]))


def matrix_product(a, b):
    return [[sum(a[i][k]*b[k][j] for k in range(len(b)))
             for j in range(len(b[0]))] for i in range(len(a))]


def simpson(fn, start, end, n):
    assert n > 0 and n % 2 == 0
    step = (end-start)/n
    return step/3*(fn(start)+fn(end)+sum(
        (4 if i % 2 else 2)*fn(start+i*step) for i in range(1, n)))


def cyclic_split_merge_graph(genus):
    """Retain parallel edges, including the three-edge genus-two case."""
    assert genus >= 2
    blocks = genus-1
    vertices = [(kind, j) for j in range(blocks) for kind in ('split', 'merge')]
    edges = []
    for j in range(blocks):
        edges.extend([(('split', j), ('merge', j)) for _ in range(2)])
        edges.append((('merge', j), ('split', (j+1) % blocks)))
    adjacent = {v: [] for v in vertices}
    incoming = {v: 0 for v in vertices}
    outgoing = {v: 0 for v in vertices}
    for a, b in edges:
        adjacent[a].append(b)
        adjacent[b].append(a)
        outgoing[a] += 1
        incoming[b] += 1
    seen = {vertices[0]}
    queue = deque(seen)
    while queue:
        for neighbor in adjacent[queue.popleft()]:
            if neighbor not in seen:
                seen.add(neighbor)
                queue.append(neighbor)
    assert len(seen) == len(vertices)
    assert all(len(adjacent[v]) == 3 for v in vertices)
    for v in vertices:
        assert (incoming[v], outgoing[v]) == ((1, 2) if v[0] == 'split' else (2, 1))
    critical_values = [F(3*j+offset, 3*blocks)
                       for j in range(blocks) for offset in (1, 2)]
    assert len(set(critical_values)) == 2*blocks
    assert all(0 < x < 1 for x in critical_values)
    actual_v, actual_e = len(vertices), len(edges)
    first_betti = actual_e-actual_v+1
    assert actual_v == 2*genus-2
    assert actual_e == 3*genus-3
    assert first_betti == genus
    assert -actual_v == 2-2*genus  # Euler sum of the glued pants.
    return {'genus': genus, 'vertices': actual_v, 'edges': actual_e,
            'first_betti_number': first_betti, 'connected': True,
            'all_vertices_trivalent': True}


def graph_checks():
    records = []
    for degree in range(6, 21):
        genus = (degree-1)*(degree-2)//2
        record = cyclic_split_merge_graph(genus)
        assert record['vertices'] == degree*(degree-3)
        assert record['edges'] == 3*degree*(degree-3)//2
        assert 2-2*genus == degree*(3-degree)
        records.append({'degree': degree, **record})
    low_genus = [cyclic_split_merge_graph(g) for g in range(2, 9)]
    return {'degrees_checked': [6, 20], 'degree_examples': records,
            'additional_low_genus_examples': low_genus,
            'scope': 'Explicit finite multigraph construction; not a smooth-map existence test.'}


def closed_budget_checks():
    count = 0
    examples = []
    fractions = [F(1, 40), F(1, 9), F(2, 11), F(3, 8)]
    for degree, area_scale, f1, f2 in product(
            range(6, 21), [F(1), F(7, 5), F(11, 8)], fractions, fractions):
        p1, p2 = area_scale*f1, area_scale*f2
        assert 0 < p1+p2 < area_scale
        nu = area_scale/F(degree-5)
        reduced_area = degree*(area_scale-p1-p2)
        chi = degree*(3-degree)
        theta_chern = degree*(p1+nu)+degree*(p2+nu)
        assert reduced_area > 0
        assert reduced_area+nu*chi+theta_chern == 0
        assert reduced_area+nu*chi-theta_chern != 0
        assert reduced_area-(-nu*chi) == -theta_chern
        assert reduced_area/(degree*area_scale) == 1-(p1+p2)/area_scale
        count += 1
        if degree in (6, 10, 20) and area_scale == 1 and f1 == fractions[0] and f2 == fractions[1]:
            examples.append({'degree': degree, 'nu': str(nu), 'p1': str(p1), 'p2': str(p2),
                             'area': str(reduced_area), 'chi': chi,
                             'chern_vector': [degree, degree],
                             'theta_chern_correction': str(theta_chern), 'budget_residual': '0'})
    return {'exact_rational_cases': count, 'examples': examples,
            'tested_identity': 'Area + nu*chi + Theta_Gamma(c) = 0, c=(d,d)',
            'scope': 'Checks the Euler correction and its sign; no bundle embedding is inferred.'}


def delta_budget_checks():
    tube_count = 0
    wide_count = 0
    positive_examples = []
    for ambient, ratio, divisor, sfraction, k in product(
            [F(1), F(7, 5)], [F(1, 4), F(1, 2), F(3, 4)],
            [4, 5, 7], [F(1, 4), F(1, 2), F(3, 4)], [1, 2]):
        base_area = ambient*ratio
        delta = base_area/divisor
        a0 = base_area-delta
        total_inner = delta*sfraction
        assert 0 < total_inner < delta < base_area < ambient
        for genus, boundary_count in product(range(6), range(1, 6)):
            defect_x = 0
            defect_z = -2*(genus+boundary_count-1)
            stokes_area = total_inner+defect_z*delta
            signed_budget = (total_inner+2*ambient*(genus+boundary_count-1)
                             + defect_x*(a0+ambient*(k-1))
                             + defect_z*(delta+ambient))
            assert stokes_area == signed_budget
            if genus+boundary_count >= 2:
                assert stokes_area < 0
            else:
                assert stokes_area == total_inner > 0
            tube_count += 1
        for defect_x in range(-6, 7):
            defect_z = -2-(k-1)*defect_x
            stokes_area = total_inner+defect_x*a0+defect_z*delta
            signed_budget = (total_inner+2*ambient+defect_x*(a0+ambient*(k-1))
                             + defect_z*(delta+ambient))
            formula = total_inner-2*delta+defect_x*(a0-(k-1)*delta)
            assert stokes_area == signed_budget == formula
            coefficient = a0-(k-1)*delta
            assert coefficient > 0
            assert (formula > 0) == (defect_x > (2*delta-total_inner)/coefficient)
            wide_count += 1
    for k, base_area, sfraction in product([1, 2], [F(1, 5), F(2, 3)], [F(1, 3), F(2, 3)]):
        delta = base_area/(k+3)
        total_inner = delta*sfraction
        a0 = base_area-delta
        defect_x, defect_z = 1, -(k+1)
        candidate_area = total_inner+base_area-(k+2)*delta
        assert delta < base_area/(k+2)
        assert candidate_area == total_inner+defect_x*a0+defect_z*delta > 0
        assert candidate_area < 2  # In these a=1 examples, no universal extra 2a follows.
        positive_examples.append({'k': k, 'base_area': str(base_area), 'delta': str(delta),
                                  'sum_inner_areas': str(total_inner),
                                  'relative_defect': [defect_x, defect_z],
                                  'candidate_area': str(candidate_area)})
    return {'tube_genus_boundary_cases': tube_count,
            'wide_exact_defect_cases': wide_count,
            'positive_wide_numerical_examples': positive_examples,
            'scope': 'These positive budgets do not construct an annulus or coisotropic connection.'}


def chart_image(coords, e1=0.2):
    r, theta, u, v, s, t = coords
    zeta = complex(u, v)
    phase = cmath.exp(2j*math.pi*theta)
    square = r+e1+math.pi*abs(zeta)**2
    assert square > 0
    x = phase*math.sqrt(square/math.pi)
    y = phase.conjugate()*zeta
    return [x.real, x.imag, y.real, y.imag, s, t]


def six_dimensional_chart_checks():
    points = [(r, theta, zeta.real, zeta.imag, xi.real, xi.imag)
              for r, theta, zeta, xi in product(
                  [-0.15, -0.02, 0.0, 0.02, 0.3], [0.0, 0.17, 0.63],
                  [0j, complex(0.1, 0.07), complex(-0.13, 0.22)],
                  [0j, complex(0.04, 0.1), complex(-0.17, 0.03)])]
    # Also sample the non-product part A<0, where p> -A is essential.
    points += [(-0.3, theta, 0.3, -0.2, 0.08, 0.03) for theta in [0.0, 0.27, 0.71]]
    maximum_form_error = 0.0
    maximum_primitive_error = 0.0
    maximum_inverse_error = 0.0
    for point in points:
        target = chart_image(point)
        columns = []
        for j in range(6):
            step = 1e-6
            hi, lo = list(point), list(point)
            hi[j] += step
            lo[j] -= step
            columns.append([(a-b)/(2*step) for a, b in zip(chart_image(hi), chart_image(lo))])
        for i in range(6):
            for j in range(6):
                pullback = sum(columns[i][k]*columns[j][k+1]
                               - columns[i][k+1]*columns[j][k] for k in [0, 2, 4])
                expected = (1 if (i, j) in [(0, 1), (2, 3), (4, 5)] else
                            -1 if (j, i) in [(0, 1), (2, 3), (4, 5)] else 0)
                maximum_form_error = max(maximum_form_error, abs(pullback-expected))
        r, theta, u, v, s, t = point
        expected_primitive = [0, r+0.2, -v/2, u/2, -t/2, s/2]
        for j in range(6):
            pullback = sum((target[k]*columns[j][k+1]
                            - target[k+1]*columns[j][k])/2 for k in [0, 2, 4])
            maximum_primitive_error = max(maximum_primitive_error, abs(pullback-expected_primitive[j]))
        x, y = complex(*target[0:2]), complex(*target[2:4])
        recovered_r = math.pi*(abs(x)**2-abs(y)**2)-0.2
        recovered_zeta = x/abs(x)*y
        maximum_inverse_error = max(maximum_inverse_error,
                                    abs(recovered_r-r), abs(recovered_zeta-complex(u, v)))
    assert maximum_form_error < 2e-8
    assert maximum_primitive_error < 2e-9
    assert maximum_inverse_error < 2e-14
    return {'six_dimensional_samples': len(points),
            'samples_with_r_negative': sum(p[0] < 0 for p in points),
            'samples_with_r_zero': sum(p[0] == 0 for p in points),
            'maximum_symplectic_pullback_error': maximum_form_error,
            'maximum_primitive_pullback_error': maximum_primitive_error,
            'maximum_inverse_error': maximum_inverse_error,
            'scope': 'Finite differences in the explicit target chart; Phi_tau exactness is proved separately.'}


def focus_radii(a, w):
    total = math.hypot(a, 2*math.pi*abs(w))
    return (total+a)/2, (total-a)/2


def focus_radial_graph(a, c, q, psi):
    angle = 2*math.pi*psi
    root = math.sqrt(q-c*c*math.sin(angle)**2)
    radius = c*math.cos(angle)+root
    radical = math.hypot(a, 2*math.pi*radius)
    return (radical-a)/2, math.pi**2*radius/(root*radical)


def winding_number(c, radius, n=1024):
    samples = [c+radius*cmath.exp(2j*math.pi*j/n) for j in range(n)]
    return sum(cmath.phase(samples[(j+1) % n]/samples[j]) for j in range(n))/(2*math.pi)


def focus_checks():
    maximum_radius_error = 0.0
    radius_count = 0
    for a, w in product([-1.0, -0.2, 0.0, 0.3, 1.2],
                        [0j, complex(0.04, 0.03), complex(-0.3, 0.1), complex(1.2, -0.8)]):
        px, py = focus_radii(a, w)
        assert px >= 0 and py >= 0
        if px > 0:
            x = math.sqrt(px/math.pi)
            y = w/x
        else:
            y = math.sqrt(py/math.pi)
            x = 0j
        maximum_radius_error = max(maximum_radius_error,
                                   abs(px-py-a), abs(px*py-math.pi**2*abs(w)**2),
                                   abs(x*y-w), abs(math.pi*abs(x)**2-px),
                                   abs(math.pi*abs(y)**2-py))
        assert (px == 0 and py == 0) == (a == 0 and w == 0)
        radius_count += 1
    assert maximum_radius_error < 2e-13
    proper_box_samples = 0
    amax, c, qmax = 0.7, 0.2, 0.3
    uniform_sum_bound = math.hypot(amax, 2*math.pi*(c+math.sqrt(qmax)))
    for a, q, angle in product([-amax, 0, amax], [0.01, c*c, qmax], range(64)):
        w = c+math.sqrt(q)*cmath.exp(2j*math.pi*angle/64)
        px, py = focus_radii(a, w)
        assert px+py <= uniform_sum_bound+1e-14
        proper_box_samples += 1
    refinement_error = 0.0
    derivative_relative_error = 0.0
    minimum_derivative = float('inf')
    action_examples = []
    for a, c, ratio in product([-0.2, 0.0, 0.05, 0.3, 1.0],
                                [0.03, 0.1, 0.25], [1.1, 1.6, 3.0]):
        q = ratio*c*c
        action = lambda qq, n: simpson(lambda psi: focus_radial_graph(a, c, qq, psi)[0], 0, 1, n)
        coarse, fine = action(q, 512), action(q, 1024)
        derivative = simpson(lambda psi: focus_radial_graph(a, c, q, psi)[1], 0, 1, 1024)
        step = q*1e-5
        fd = (action(q+step, 1024)-action(q-step, 1024))/(2*step)
        assert derivative > 0 and fd > 0
        refinement_error = max(refinement_error, abs(fine-coarse))
        derivative_relative_error = max(derivative_relative_error, abs(fd-derivative)/derivative)
        minimum_derivative = min(minimum_derivative, derivative)
        assert fine+a > 0  # The other section action is J+A.
        action_examples.append({'A': a, 'c': c, 'q_over_c_squared': ratio,
                                'normal_action': fine, 'derivative': derivative})
    assert refinement_error < 1e-10
    assert derivative_relative_error < 2e-6
    clockwise = [[1, 0, 1], [0, 1, 0], [0, 0, 1]]
    counterclockwise = [[1, 0, -1], [0, 1, 0], [0, 0, 1]]
    assert determinant3(clockwise) == determinant3(counterclockwise) == 1
    assert matrix_product(clockwise, counterclockwise) == [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    winding_records = []
    for c in [0.03, 0.2, 1.0]:
        below = winding_number(c, c*math.sqrt(0.4))
        above = winding_number(c, c*math.sqrt(2.0))
        assert abs(below) < 1e-12 and abs(above-1) < 1e-12
        winding_records.append({'c': c, 'section_degree_below': below, 'section_degree_above': above})
    return {'radius_inverse_samples': radius_count, 'maximum_radius_inverse_error': maximum_radius_error,
            'compact_box_bound_samples': proper_box_samples,
            'action_quadrature_cases': len(action_examples),
            'maximum_512_to_1024_refinement_difference': refinement_error,
            'maximum_derivative_relative_error': derivative_relative_error,
            'minimum_action_derivative_sampled': minimum_derivative,
            'action_examples': action_examples,
            'clockwise_monodromy_in_basis_a_b_gamma': clockwise,
            'section_transition_winding_checks': winding_records,
            'scope': 'Numerical radii/action/section-degree checks; properness and monodromy are proved in the text.'}


def lattice_checks():
    # Rows give (a,b,m) in the old basis (gamma_U,Gamma_H,Gamma_J).
    transition = [[2, -1, 0], [-1, 1, -1], [-1, 1, 0]]
    assert determinant3(transition) == 1
    count = 0
    for u, v, w, k, a0, base_area in product(
            range(-3, 4), range(-3, 4), range(-3, 4), [1, 2],
            [F(1, 5), F(1, 3)], [F(1, 2), F(4, 5)]):
        old_p, old_ell, old_r = 2*u-v-w, -u+v+w, -v
        toric = (old_p+old_ell, old_ell+old_r, -old_r)
        assert toric == (u, -u+w, v)
        old_area = old_p*a0+old_ell*base_area
        new_area = u*(2*a0-base_area)+(v+w)*(base_area-a0)
        old_maslov = 2*((k-1)*old_p+k*old_ell)
        new_maslov = 2*((k-2)*u+v+w)
        assert old_area == new_area
        assert old_maslov == new_maslov
        source_maslov = 2*sum(toric)
        assert source_maslov == 2*(v+w)
        assert old_maslov-source_maslov == 2*(k-2)*u
        count += 1
    return {'transition_matrix': transition, 'determinant': 1,
            'exact_area_maslov_cases': count,
            'scope': 'Conditional marking compatibility; it does not identify two actual charts by itself.'}


def main():
    report = {
        'scope': ('Finite exact-arithmetic, combinatorial, and selected numerical checks only; '
                  'no existence proof for geometric embeddings, port connections, or a global fibration.'),
        'closed_curve_graphs': graph_checks(),
        'closed_reduction_euler_budget': closed_budget_checks(),
        'delta_port_budget': delta_budget_checks(),
        'six_dimensional_membrane_chart': six_dimensional_chart_checks(),
        'focus_block': focus_checks(),
        'integral_lattice_transition': lattice_checks(),
        'status': 'passed',
    }
    target = Path(__file__).with_name('closed_curve_family_checks.json')
    temporary = target.with_suffix('.json.tmp')
    temporary.write_text(json.dumps(report, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    temporary.replace(target)
    print('Closed-curve and membrane-family checks passed:', target.name)


if __name__ == '__main__':
    main()
