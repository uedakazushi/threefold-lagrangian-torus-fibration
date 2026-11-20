#!/usr/bin/env python3
"""Finite exact checks for normal focus bundles and product regions.

Only the Python standard library is used.  These checks verify displayed
linear algebra, multigraphs, integer degrees, Euler identities and polynomial
jets.  They do not prove geometric existence, smoothness of a hypersurface,
symplectic embeddings, nonlinear gluing, or a boundary-free global fibration.
Those assertions and their hypotheses are established in the manuscript.
"""
from collections import deque
from fractions import Fraction as F
from itertools import product
from math import comb
from pathlib import Path
import json
import os


def eye(n):
    return [[int(i == j) for j in range(n)] for i in range(n)]


def zero(n, m=None):
    return [[0 for _ in range(m if m is not None else n)] for _ in range(n)]


def transpose(a):
    return [list(row) for row in zip(*a)]


def multiply(a, b):
    return [[sum(a[i][k] * b[k][j] for k in range(len(b)))
             for j in range(len(b[0]))] for i in range(len(a))]


def scale(c, a):
    return [[c * entry for entry in row] for row in a]


def add(*matrices):
    return [[sum(a[i][j] for a in matrices)
             for j in range(len(matrices[0][0]))]
            for i in range(len(matrices[0]))]


def block_diag(*matrices):
    sizes = [len(a) for a in matrices]
    out = zero(sum(sizes))
    offset = 0
    for a in matrices:
        for i, row in enumerate(a):
            for j, entry in enumerate(row):
                out[offset + i][offset + j] = entry
        offset += len(a)
    return out


def canonical_form(n):
    out = zero(2*n)
    for i in range(n):
        out[i][i+n] = 1
        out[i+n][i] = -1
    return out


def determinant(matrix):
    a = [[F(x) for x in row] for row in matrix]
    value = F(1)
    for i in range(len(a)):
        pivot = next((k for k in range(i, len(a)) if a[k][i]), None)
        if pivot is None:
            return F(0)
        if pivot != i:
            a[i], a[pivot] = a[pivot], a[i]
            value = -value
        entry = a[i][i]
        value *= entry
        for k in range(i+1, len(a)):
            ratio = a[k][i]/entry
            for j in range(i+1, len(a)):
                a[k][j] -= ratio*a[i][j]
            a[k][i] = 0
    return value


def rank(matrix):
    a = [[F(x) for x in row] for row in matrix]
    row = 0
    for column in range(len(a[0])):
        pivot = next((k for k in range(row, len(a)) if a[k][column]), None)
        if pivot is None:
            continue
        a[row], a[pivot] = a[pivot], a[row]
        divisor = a[row][column]
        a[row] = [x/divisor for x in a[row]]
        for k in range(len(a)):
            if k != row:
                factor = a[k][column]
                a[k] = [x-factor*y for x, y in zip(a[k], a[row])]
        row += 1
        if row == len(a):
            break
    return row


def flatten(a):
    return [entry for row in a for entry in row]


def commutator(a, b):
    return add(multiply(a, b), scale(-1, multiply(b, a)))


def graph(genus, punctures):
    """Explicit cyclic split-merge graph; subdivide edges and add leaf edges."""
    assert genus >= 2 and punctures >= 0
    blocks = genus-1
    vertices = [('split', j) for j in range(blocks)]
    vertices += [('merge', j) for j in range(blocks)]
    edges = []
    for j in range(blocks):
        edges += [(('split', j), ('merge', j))]*2
        edges.append((('merge', j), ('split', (j+1) % blocks)))
    for j in range(punctures):
        old_start, old_end = edges.pop(0)
        vertex, leaf = ('puncture_saddle', j), ('leaf', j)
        vertices += [vertex, leaf]
        edges += [(old_start, vertex), (vertex, old_end), (leaf, vertex)]
    adjacency = {v: [] for v in vertices}
    for a, b in edges:
        adjacency[a].append(b)
        adjacency[b].append(a)
    seen = {vertices[0]}
    todo = deque(seen)
    while todo:
        for other in adjacency[todo.popleft()]:
            if other not in seen:
                seen.add(other)
                todo.append(other)
    assert len(seen) == len(vertices)
    internal = [v for v in vertices if v[0] != 'leaf']
    leaves = [v for v in vertices if v[0] == 'leaf']
    assert all(len(adjacency[v]) == 3 for v in internal)
    assert all(len(adjacency[v]) == 1 for v in leaves)
    assert len(internal) == 2*genus+punctures-2
    assert len(leaves) == punctures
    assert len(edges) == 3*genus+2*punctures-3
    first_betti = len(edges)-len(vertices)+1
    assert first_betti == genus
    assert -len(internal) == 2-2*genus-punctures
    return {'genus': genus, 'punctures': punctures,
            'trivalent_vertices': len(internal), 'leaves': len(leaves),
            'edges': len(edges), 'first_betti_number': first_betti,
            'product_region_euler': -len(internal)}


def graph_and_degree_checks():
    punctured_cases = 0
    closed_cases = 0
    budget_cases = 0
    examples = []
    closed_examples = []
    for d in range(6, 21):
        genus = (d-1)*(d-2)//2
        assert 2-2*genus == d*(3-d)
        generic_normal_degree = (5-d)*d-(2-2*genus)
        assert generic_normal_degree == 2*d
        for r in range(1, 2*d+1):
            record = graph(genus, r)
            assert record['trivalent_vertices'] == d*(d-3)+r
            assert record['product_region_euler'] == -d*(d-3)-r
            punctured_cases += 1
            if d == 6 and r in (1, 6, 12):
                examples.append({'degree': d, **record})

        curve_degree = d-2
        closed_genus = (curve_degree-1)*(curve_degree-2)//2
        closed = graph(closed_genus, 0)
        normal_degree = (5-d)*curve_degree-(2-2*closed_genus)
        determinant_formula = curve_degree*(curve_degree+2-d)
        assert normal_degree == determinant_formula == 0
        node_count = (d-2)*(d-5)
        assert closed['trivalent_vertices'] == node_count
        assert closed['product_region_euler'] == -node_count
        chi_x = d*(10-10*d+5*d*d-d*d*d)
        complement_euler = chi_x-closed['product_region_euler']
        assert complement_euler == chi_x+node_count
        if d == 6:
            assert closed_genus == 3 and node_count == 4
            assert chi_x == -516 and complement_euler == -512
        for area_scale in (F(1), F(3, 2), F(11, 7)):
            nu = area_scale/F(d-5)
            curve_area = curve_degree*area_scale
            assert curve_area+nu*(2-2*closed_genus) == 0
            budget_cases += 1
        if d in (6, 7, 10, 20):
            closed_examples.append({
                'ambient_degree': d, 'curve_degree': curve_degree,
                **closed, 'normal_degree': normal_degree,
                'ambient_euler': chi_x, 'complement_euler': complement_euler,
            })
        closed_cases += 1
    return {'degrees_checked': [6, 20],
            'explicit_punctured_multigraphs': punctured_cases,
            'explicit_closed_multigraphs': closed_cases,
            'exact_closed_curve_budget_cases': budget_cases,
            'sextic_punctured_examples': examples,
            'closed_normal_trivial_examples': closed_examples,
            'scope': 'Combinatorics, adjunction numbers and Euler identities only; no geometric existence test.'}


def cartan_checks():
    ident2, ident4 = eye(2), eye(4)
    omega4 = canonical_form(2)
    rotation = [[0, -1], [1, 0]]
    reflection = [[1, 0], [0, -1]]
    p = block_diag(ident2, scale(-1, ident2))
    r = block_diag(rotation, rotation)
    s = scale(-1, omega4)
    k = block_diag(reflection, reflection)
    assert commutator(p, r) == zero(4)
    for a in (p, r):
        assert add(multiply(transpose(a), omega4), multiply(omega4, a)) == zero(4)
    assert multiply(p, p) == ident4
    assert multiply(r, r) == scale(-1, ident4)
    assert rank([flatten(p), flatten(r)]) == 2
    assert multiply(multiply(transpose(p), omega4), p) == scale(-1, omega4)
    assert multiply(s, s) == scale(-1, ident4)
    assert multiply(k, k) == ident4
    for a in (s, k):
        assert multiply(multiply(transpose(a), omega4), a) == omega4
        assert determinant(a) == 1
    s_inverse = scale(-1, s)
    assert multiply(multiply(s, p), s_inverse) == scale(-1, p)
    assert multiply(multiply(s, r), s_inverse) == r
    assert multiply(multiply(k, p), k) == p
    assert multiply(multiply(k, r), k) == scale(-1, r)
    assert commutator(s, k) == zero(4)

    # A compatible J exchanging the standard Lagrangian pair.
    j = scale(-1, omega4)
    assert multiply(j, j) == scale(-1, ident4)
    assert multiply(omega4, j) == ident4
    assert add(multiply(p, j), multiply(j, p)) == zero(4)

    rotation_rational = [[F(3, 5), F(-4, 5)], [F(4, 5), F(3, 5)]]
    centralizer_cases = 0
    for radius in (F(1, 3), F(1), F(5, 2)):
        a = block_diag(scale(radius, rotation_rational),
                       scale(1/radius, rotation_rational))
        assert multiply(multiply(transpose(a), omega4), a) == omega4
        assert commutator(a, p) == commutator(a, r) == zero(4)
        centralizer_cases += 1

    characteristic_cases = 0
    for a, b in product(range(-3, 4), repeat=2):
        if a == b == 0:
            continue
        matrix = add(scale(a, p), scale(b, r))
        square = multiply(matrix, matrix)
        fourth = multiply(square, square)
        polynomial = add(fourth, scale(2*(b*b-a*a), square),
                         scale((a*a+b*b)**2, ident4))
        assert polynomial == zero(4)
        assert determinant(matrix) == (a*a+b*b)**2
        characteristic_cases += 1

    omega6 = block_diag(omega4, canonical_form(1))
    p6, r6 = block_diag(p, zero(2)), block_diag(r, zero(2))
    hyperbolic = block_diag(zero(4), [[1, 0], [0, -1]])
    triple = [p6, r6, hyperbolic]
    assert rank([flatten(a) for a in triple]) == 3
    for a, b in product(triple, repeat=2):
        assert commutator(a, b) == zero(6)
    for a in triple:
        assert add(multiply(transpose(a), omega6), multiply(omega6, a)) == zero(6)
    direct_sum_cases = 0
    for a, b, h in product((1, 2), (1, 3), (-2, 1, 4)):
        matrix = add(scale(a, p6), scale(b, r6), scale(h, hyperbolic))
        assert determinant(matrix) == -h*h*(a*a+b*b)**2
        square = multiply(matrix, matrix)
        fourth = multiply(square, square)
        focus_polynomial = add(fourth, scale(2*(b*b-a*a), square),
                               scale((a*a+b*b)**2, eye(6)))
        hyperbolic_polynomial = add(square, scale(-h*h, eye(6)))
        assert multiply(focus_polynomial, hyperbolic_polynomial) == zero(6)
        direct_sum_cases += 1
    return {'focus_characteristic_polynomial_cases': characteristic_cases,
            'centralizer_samples_exact': centralizer_cases,
            'normalizer_generators': {'S_P_Sinv': '-P', 'S_R_Sinv': 'R',
                                     'K_P_Kinv': 'P', 'K_R_Kinv': '-R',
                                     'S_squared': '-I', 'K_squared': 'I'},
            'focus_hyperbolic_direct_sum_cases': direct_sum_cases,
            'sp_and_commutator_checks': 'passed',
            'scope': 'Exact matrices for the claimed normal forms; not nonlinear normal-form verification.'}


def complex_multiply(a, b):
    return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])


def complex_subtract(a, b):
    return (a[0]-b[0], a[1]-b[1])


def determinant_phase_checks():
    phase_cases = 0
    for b in ((F(1), F(0)), (F(0), F(1)), (F(-1), F(0)),
              (F(3, 5), F(4, 5)), (F(5, 13), F(-12, 13))):
        assert b[0]*b[0]+b[1]*b[1] == 1
        conjugate = (b[0], -b[1])
        i_conjugate = complex_multiply((0, 1), conjugate)
        # Columns (1, i*conj(b)) and (i, conj(b)) span an unstable plane.
        det = complex_subtract(conjugate, complex_multiply((0, 1), i_conjugate))
        assert det == (2*conjugate[0], 2*conjugate[1])
        assert complex_multiply(b, conjugate) == (1, 0)
        squared_phase = complex_multiply(conjugate, conjugate)
        assert complex_multiply(complex_multiply(b, b), squared_phase) == (1, 0)
        phase_cases += 1
    return {'exact_unit_complex_samples': phase_cases,
            'determinant_phase': 'b^(-1)', 'squared_phase': 'b^(-2)',
            'scope': 'Checks the phase and factor of two in the fixed unitary model.'}


def clutching_action_checks():
    cases = 0
    omega6 = canonical_form(3)
    max_order = 0
    for a, b in product(range(21), repeat=2):
        m = a+b
        if not m:
            continue
        angle_matrix = [[1, 0, 0], [a, 1, 0], [b, 0, 1]]
        angle_inverse = [[1, 0, 0], [-a, 1, 0], [-b, 0, 1]]
        action_matrix = transpose(angle_inverse)
        transform = block_diag(action_matrix, angle_matrix)
        assert multiply(angle_matrix, angle_inverse) == eye(3)
        assert multiply(transpose(action_matrix), angle_matrix) == eye(3)
        assert multiply(multiply(transpose(transform), omega6), transform) == omega6
        assert determinant(transform) == 1
        assert determinant([[a, 1], [b, -1]]) == -m
        stabilizer = {(F(j, m), F((-a*j) % m, m)) for j in range(m)}
        assert len(stabilizer) == m
        for u, v in stabilizer:
            assert (a*u+v).denominator == 1
            assert (b*u-v).denominator == 1
        assert (F(1 % m, m), F((-a) % m, m)) in stabilizer
        max_order = max(max_order, m)
        cases += 1
    return {'nonnegative_integer_weight_pairs': cases,
            'weights_in_each_coordinate': [0, 20],
            'maximum_isotropy_order': max_order,
            'exact_primitive_and_symplectic_matrix_checks': 'passed',
            'scope': 'Action-angle algebra and stabilizer arithmetic; no collar embedding or cap existence test.'}


def leading_polynomial(m, real_k, imag_k, c):
    """The exact polynomial -2c Re(K (u-i v)^m)."""
    out = {}
    powers_minus_i = ((1, 0), (0, -1), (-1, 0), (0, 1))
    for j in range(m+1):
        real, imag = powers_minus_i[j % 4]
        coefficient = -2*c*comb(m, j)*(real_k*real-imag_k*imag)
        if coefficient:
            out[(m-j, j)] = coefficient
    return out


def derivative(poly, variable):
    out = {}
    for powers, coefficient in poly.items():
        exponent = powers[variable]
        if exponent:
            new = list(powers)
            new[variable] -= 1
            out[tuple(new)] = coefficient*exponent
    return out


def polynomial_jet_checks():
    records = []
    real_k, imag_k, c = F(3, 5), F(4, 5), F(7, 3)
    for m in range(1, 41):
        poly = leading_polynomial(m, real_k, imag_k, c)
        assert poly and all(sum(exponents) == m for exponents in poly)
        gradient = [derivative(poly, i).get((0, 0), 0) for i in range(2)]
        hessian = [[derivative(derivative(poly, i), j).get((0, 0), 0)
                    for j in range(2)] for i in range(2)]
        invariant_linear_dimension = sum((2*i-1) % m == 0 for i in range(2))
        invariant_quadratic_dimension = sum((2*i-2) % m == 0 for i in range(3))
        if m == 1:
            assert gradient != [0, 0]
            assert invariant_linear_dimension == 2
            assert invariant_quadratic_dimension == 3
            kind = 'regular'
        elif m == 2:
            assert gradient == [0, 0]
            assert determinant(hessian) == -16*c*c*(real_k*real_k+imag_k*imag_k) < 0
            assert invariant_linear_dimension == 0
            assert invariant_quadratic_dimension == 3
            kind = 'nondegenerate_hyperbolic'
        else:
            assert gradient == [0, 0] and hessian == zero(2)
            assert invariant_linear_dimension == 0
            assert invariant_quadratic_dimension == 1
            # The sole invariant quadratic monomial is u^2+v^2 (z*conj(z)).
            assert [i for i in range(3) if (2*i-2) % m == 0] == [1]
            kind = 'order_m_degenerate'
        # For K=1 the leading zero rays have normalized angles (1+2j)/(4m).
        rays = {F(1+2*j, 4*m) for j in range(2*m)}
        assert len(rays) == 2*m and all(0 < angle < 1 for angle in rays)
        if m in (1, 2, 3, 12, 40):
            records.append({'isotropy_order': m, 'generic_center_type': kind,
                            'leading_degree': m, 'leading_half_rays': 2*m,
                            'real_invariant_linear_dimension': invariant_linear_dimension,
                            'real_invariant_quadratic_dimension': invariant_quadratic_dimension})
    return {'orders_checked': [1, 40], 'examples': records,
            'invariant_quadratic_rule': 'For m>2 only the radial quadratic form survives.',
            'scope': 'Exact polynomial jets and cyclic characters; not a nonlinear repair construction.'}


def degree_distribution_checks():
    count = 0
    examples = []
    for d in range(6, 21):
        for r in range(1, 2*d+1):
            # A basis of the integer boundary-evaluation image, sum coordinates=0.
            basis = []
            for i in range(r-1):
                vector = [0]*r
                vector[i], vector[-1] = 1, -1
                basis.append(vector)
            if basis:
                assert rank(basis) == r-1
                assert all(sum(v) == 0 for v in basis)
            target = [((3*i+2) % 7)-3 for i in range(r-1)]
            target.append(-sum(target))
            rebuilt = [sum(target[i]*basis[i][j] for i in range(r-1)) for j in range(r)]
            assert rebuilt == target
            initial = [0]*r
            initial[0] = d
            a = [initial[i]+target[i] for i in range(r)]
            b = [initial[i]-target[i] for i in range(r)]
            k = [x+y for x, y in zip(a, b)]
            assert sum(a) == sum(b) == d and sum(k) == 2*d
            mixed_windings = [-x for x in k]
            squared_phase_windings = [-2*x for x in mixed_windings]
            assert sum(mixed_windings) == -2*d
            assert sum(squared_phase_windings) == 4*d
            count += 1

        for r, assignments, bound in (
            (d, [(1, 1)]*d, 2),
            (2*d, [(1, 0)]*d+[(0, 1)]*d, 1),
        ):
            assert len(assignments) == r
            assert sum(x for x, _ in assignments) == d
            assert sum(y for _, y in assignments) == d
            charges = [x+y for x, y in assignments]
            assert all(0 < abs(k) <= bound for k in charges)
            assert r == (2*d+bound-1)//bound
            nodes = d*(d-3)+r
            assert nodes == (d*(d-2) if bound == 2 else d*(d-1))
            if d == 6:
                examples.append({'degree': d, 'punctures': r,
                                 'charge_per_puncture': bound,
                                 'focus_hyperbolic_whole_fibers_in_region': nodes})
    return {'integer_boundary_distribution_cases': count,
            'boundary_evaluation_image': 'kernel of the coordinate sum map Z^r -> Z',
            'mixed_coefficient_total_winding': '-2d',
            'squared_determinant_phase_total_winding': '4d',
            'sextic_small_charge_allocations': examples,
            'scope': 'Integer degree arithmetic and conditional bounds; not completed caps or global fiber counts.'}


def main():
    result = {
        'scope': (
            'Finite exact-arithmetic, multigraph, matrix, degree and polynomial-jet checks only. '
            'No computational claim of geometric existence, hypersurface smoothness, actual embedding, '
            'cap gluing, or a boundary-free global fibration.'
        ),
        'graphs_and_curve_degrees': graph_and_degree_checks(),
        'focus_cartan_linear_algebra': cartan_checks(),
        'determinant_phase': determinant_phase_checks(),
        'weighted_clutching_actions': clutching_action_checks(),
        'weighted_center_jets': polynomial_jet_checks(),
        'boundary_degree_distribution': degree_distribution_checks(),
        'status': 'passed',
    }
    destination = Path(__file__).with_name('focus_normal_region_checks.json')
    temporary = destination.with_suffix('.json.tmp')
    temporary.write_text(json.dumps(result, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    os.replace(temporary, destination)
    print('PASS: normal focus linear algebra, explicit region graphs, degree budgets and weighted jets.')
    print('These finite checks do not prove geometric existence or global gluing.')
    print(destination.name)


if __name__ == '__main__':
    main()
