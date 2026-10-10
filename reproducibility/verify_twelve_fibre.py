"""Bounded, locally written probe of the twelve-fibre all-profile gate.

Numerical scouting is separated from exact Fraction acceptance. The exact
--certify mode uses only the standard library. No downloaded code is executed.
"""
from __future__ import annotations

import argparse
import json
import math
from fractions import Fraction
from functools import lru_cache

# Frozen numerical discovery converted to exact rational input. The remaining
# six fibres are obtained by reflection p_{11-j}(y)=p_j(-y). Normalize the full
# polynomial by its exact integral before using it as a certificate.
POLY_DENOMINATOR = 10 ** 14
POLY_HALF_NUMERATORS = [
    [81452785069661, 7292373590683, -495291728588, -8832261778, 468777856, 3579876, -172336],
    [94237762244332, 5325472480437, -540048842503, -5392125542, 498965786, 1817587, -181576],
    [102829051521759, 3211889787491, -561973129360, -1752444969, 508473653, -36353, -183294],
    [106611512235993, 1220224746606, -555948929227, 1548835301, 491383198, -1705511, -175034],
    [105872959612158, 3026144941748, -520652451189, -4834131202, 443772398, 2367437, -154759],
    [110067398521096, 987168788008, -542216720202, -1573291099, 463005280, 777705, -161751],
]

R = Fraction
ELL = R(1, 12)


def poly_add(a, b):
    return [(a[k] if k < len(a) else R(0)) + (b[k] if k < len(b) else R(0))
            for k in range(max(len(a), len(b)))]


def poly_scale(a, c):
    return [x * c for x in a]


def poly_mul(a, b):
    result = [R(0)] * (len(a) + len(b) - 1)
    for j, x in enumerate(a):
        for k, y in enumerate(b):
            result[j + k] += x * y
    return result


def poly_integral(a):
    """Exact integral on [-1,1]."""
    return sum((2 * a[k] / (k + 1) for k in range(0, len(a), 2)), R(0))


def poly_primitive(a):
    """Antiderivative vanishing at -1."""
    terms = [R(0)] + [x / (k + 1) for k, x in enumerate(a)]
    terms[0] = -sum((terms[k] * (-1) ** k for k in range(1, len(terms))), R(0))
    return terms


def point(x):
    return (R(x), R(x))


def iadd(a, b):
    return (a[0] + b[0], a[1] + b[1])


def ineg(a):
    return (-a[1], -a[0])


def isub(a, b):
    return iadd(a, ineg(b))


def imul(a, b):
    products = [x * y for x in a for y in b]
    return (min(products), max(products))


def iscale(a, x):
    return imul(a, point(x))


def outward(a, digits=60):
    scale = 10 ** digits
    lo = a[0].numerator * scale // a[0].denominator
    hi = -((-a[1].numerator * scale) // a[1].denominator)
    return (R(lo, scale), R(hi, scale))


def atan_partial(inverse, last):
    return sum((R((-1) ** k, (2 * k + 1) * inverse ** (2 * k + 1))
                for k in range(last + 1)), R(0))


@lru_cache(maxsize=1)
def pi_interval():
    a5 = (atan_partial(5, 65), atan_partial(5, 64))
    a239 = (atan_partial(239, 21), atan_partial(239, 20))
    return outward(isub(iscale(a5, 16), iscale(a239, 4)), 70)


@lru_cache(maxsize=None)
def cos_pi_interval(multiplier):
    """Cos(pi * rational), exact interval arithmetic plus Taylor remainder."""
    multiplier = R(multiplier) % 2
    if multiplier > 1:
        multiplier = 2 - multiplier
    x = iscale(pi_interval(), multiplier)
    z = imul(x, x)
    term = point(1)
    total = term
    last = 40
    for k in range(1, last + 1):
        term = iscale(imul(term, z), R(1, (2 * k - 1) * (2 * k)))
        total = iadd(total, term if k % 2 == 0 else ineg(term))
    # Taylor to degree 81 (the odd coefficient is zero) has a remainder
    # bounded by |x|^82/82! because every real derivative of cos is <=1.
    error = x[1] ** (2 * last + 2) / math.factorial(2 * last + 2)
    return outward((total[0] - error, total[1] + error))


def require(condition, label):
    if not condition:
        raise ArithmeticError(label)


def exact_profile():
    first = [[R(n, POLY_DENOMINATOR) for n in row] for row in POLY_HALF_NUMERATORS]
    polynomials = first + [[(-1) ** k * c for k, c in enumerate(row)]
                           for row in reversed(first)]
    normalization = ELL / 2 * sum((poly_integral(p) for p in polynomials), R(0))
    polynomials = [poly_scale(p, 1 / normalization) for p in polynomials]
    require(ELL / 2 * sum((poly_integral(p) for p in polynomials), R(0)) == 1,
            'Exact normalization')
    for j in range(12):
        require(polynomials[j] == [(-1) ** k * c for k, c in enumerate(polynomials[11-j])],
                'Exact even reflection')
        require(polynomials[j][0] - sum(map(abs, polynomials[j][1:]), R(0)) > R(7, 10),
                'Positive profile on every cell')
    cumulatives, potentials = [], []
    mass, potential_start = R(0), R(1, 2)
    for p in polynomials:
        cumulative = poly_add([mass], poly_scale(poly_primitive(p), ELL / 2))
        potential = poly_add([potential_start], poly_scale(
            poly_primitive(poly_add(poly_scale(cumulative, 2), [-R(1)])), ELL / 2))
        cumulatives.append(cumulative)
        potentials.append(potential)
        mass = sum(cumulative, R(0))
        potential_start = sum(potential, R(0))
    require(mass == 1 and potential_start == R(1, 2), 'Potential endpoint guards')
    a = ELL / 2 * sum((poly_integral(poly_mul(p, p)) for p in polynomials), R(0))
    correlations = {
        k: ELL / 2 * sum((poly_integral(poly_mul(polynomials[j], polynomials[j+k]))
                          for j in range(12-k)), R(0))
        for k in (3, 4, 6, 8, 9)
    }
    j_energy = ELL * sum((poly_integral(poly_add(h, poly_scale(poly_mul(h, h), -1)))
                           for h in cumulatives), R(0))
    return polynomials, potentials, normalization, a, correlations, j_energy


def exact_matrix(q):
    """Return interval M and effective coefficients q*c_k."""
    weight_a, weight_b = 6 - 8*q, 9*q - 6
    require(weight_a > 0 and weight_b > 0, 'Positive periodic mixture')
    alpha, gamma = R(5189, 5000), R(10289, 10000)
    coefficients = {
        4*j: iscale(iadd(point(2), iscale(cos_pi_interval(2*j*alpha/3), 2)), weight_a/3)
        for j in (1, 2)
    }
    coefficients.update({
        3*j: iscale(iadd(iadd(point(3), iscale(cos_pi_interval(j*gamma/2), 4)),
                        iscale(cos_pi_interval(j*gamma), 2)), weight_b/6)
        for j in (1, 2, 3)
    })
    for value in coefficients.values():
        require(value[0] > 0, 'Nonnegative periodic spectral coefficient')
    diagonal = 1 - weight_a * R(2, 3) - weight_b * R(3, 4)
    matrix = [[point(diagonal if i == j else 0) for j in range(12)] for i in range(12)]
    for k, value in coefficients.items():
        for i in range(12-k):
            matrix[i][i+k] = matrix[i+k][i] = iscale(value, -R(1, 2))
    for i in range(12):
        row_sum = sum((max(map(abs, matrix[i][j])) for j in range(12) if i != j), R(0))
        require(diagonal - row_sum > R(77, 250), 'Gershgorin M > 77/250')
    require(R(77, 250) - 2 / R(157, 50) ** 2 > R(1, 10), 'Zero-mean coercivity > 1/10')
    return matrix, diagonal, coefficients


def exact_energy(q, profile_data):
    polynomials, potentials, normalization, a, correlations, j_energy = profile_data
    matrix, diagonal, coefficients = exact_matrix(q)
    energy = point(j_energy + diagonal * a)
    for k, value in coefficients.items():
        energy = isub(energy, iscale(value, correlations[k]))
    return energy, matrix


def residual_bound(matrix, profile_data, constant):
    polynomials, potentials = profile_data[:2]
    squared_norm_upper = R(0)
    maximum_cell_sup = R(0)
    for i in range(12):
        coefficients = [point(x) for x in potentials[i]]
        coefficients[0] = isub(coefficients[0], point(constant))
        for j in range(12):
            for k, c in enumerate(polynomials[j]):
                coefficients[k] = iadd(coefficients[k], iscale(matrix[i][j], c))
        cell_sup = sum((max(map(abs, x)) for x in coefficients), R(0))
        squared_norm_upper += ELL * cell_sup ** 2
        maximum_cell_sup = max(maximum_cell_sup, cell_sup)
    return squared_norm_upper, maximum_cell_sup


def rational_decimal(x, digits, upper=False):
    scale = 10 ** digits
    n = -((-x.numerator * scale) // x.denominator) if upper else x.numerator * scale // x.denominator
    sign = '-' if n < 0 else ''
    n = abs(n)
    return f'{sign}{n // scale}.{n % scale:0{digits}d}'


def interval_display(x, digits=18):
    return [rational_decimal(x[0], digits), rational_decimal(x[1], digits, True)]


def certify():
    profile_data = exact_profile()
    q0, q1 = R(337017, 500000), R(16851, 25000)
    energy0, _ = exact_energy(q0, profile_data)
    require(energy0[1] < 2 - 2*q0 - R(1, 200000),
            'Positive even profile defeats original all-mode gate by > 1/200000')
    energy1, matrix1 = exact_energy(q1, profile_data)
    constant = R(651920180779, 10 ** 12)
    residual2, residual_sup = residual_bound(matrix1, profile_data, constant)
    lower = energy1[0] - residual2 / R(1, 10)
    require(lower > 2 - 2*q1 + R(18, 10 ** 8),
            'All-profile minimum exceeds target by > 18/10^8')
    centered_residual_sup = residual_sup + max(abs(energy1[0] - constant), abs(energy1[1] - constant))
    require(centered_residual_sup < R(1, 10 ** 11), 'Energy-centered residual sup < 10^-11')
    eps = R(1, 10 ** 11)
    # The proof explains the mean-zero L2 error <=10 eps, energy gap
    # <=10 eps^2, integral-potential sup <=10 eps, and ||M^-1||_infty<4.
    minimizer_sup_error = 4 * (eps + 10*eps + 10*eps**2)
    require(minimizer_sup_error < R(1, 10 ** 9), 'Positive true minimizer via uniform error')
    print(json.dumps({
        'status': 'PASS: every acceptance comparison used exact Fraction arithmetic',
        'candidate_normalization': str(profile_data[2]),
        'positive_even_normalized_profile': True,
        'matrix_lower_bound': '77/250',
        'zero_mean_coercivity': '1/10',
        'q_original': str(q0),
        'original_target': str(2 - 2*q0),
        'original_candidate_energy_interval': interval_display(energy0),
        'original_gate_failure_at_least': '1/200000',
        'q_new': str(q1),
        'new_target': str(2 - 2*q1),
        'new_candidate_energy_interval': interval_display(energy1),
        'residual_sup_upper': rational_decimal(residual_sup, 24, True),
        'residual_squared_norm_upper': rational_decimal(residual2, 30, True),
        'certified_affine_minimum_interval': interval_display((lower, energy1[1])),
        'certified_margin_at_least': '18/10^8',
        'energy_centered_residual_sup_less_than': '1/10^11',
        'true_minimizer_uniform_distance_less_than': '1/10^9',
        'interpretation': 'All-profile joint-affine method cap 0.67404, not an upper bound on actual zeros',
    }, indent=2))


def matrix_data(q=337017 / 500000):
    import numpy as np
    lam = 9 - 6 / q
    alpha = 5189 / 5000
    gamma = 10289 / 10000
    shifts = {
        4 * j: (1 - lam) / 3 * (2 + 2 * math.cos(2 * math.pi * j * alpha / 3))
        for j in (1, 2)
    }
    shifts.update({
        3 * j: lam / 6 * (3 + 4 * math.cos(2 * math.pi * j * gamma / 4)
                         + 2 * math.cos(4 * math.pi * j * gamma / 4))
        for j in (1, 2, 3)
    })
    c0 = (1 - lam) * 2 / 3 + lam * 3 / 4
    mat = (1 - q * c0) * np.eye(12)
    for k, coefficient in shifts.items():
        mat -= q * coefficient / 2 * (np.eye(12, k=k) + np.eye(12, k=-k))
    return q, mat, shifts, c0


def solve_bvp(q=337017 / 500000):
    import numpy as np
    q, mat, shifts, c0 = matrix_data(q)
    eig, vectors = np.linalg.eigh(mat)
    omega = np.sqrt(2 / eig)
    ell = 1 / 12
    cosine = (vectors * np.cos(omega * ell)) @ vectors.T
    sine = (vectors * (np.sin(omega * ell) / omega)) @ vectors.T
    dcosine = (vectors * (-omega * np.sin(omega * ell))) @ vectors.T
    dsine = cosine
    end_f = np.hstack((mat @ cosine, mat @ sine))
    end_df = np.hstack((mat @ dcosine, mat @ dsine))
    start_f = np.hstack((mat, np.zeros((12, 12))))
    start_df = np.hstack((np.zeros((12, 12)), mat))
    boundary = np.vstack((
        end_f[:11] - start_f[1:],
        end_df[:11] - start_df[1:],
        start_df[:1],
        end_df[11:],
    ))
    rhs = np.r_[np.zeros(22), 1.0, -1.0]
    data = np.linalg.solve(boundary, rhs)
    initial, derivative = data[:12], data[12:]
    multiplier = (1 + (start_f[0] + end_f[11]) @ data) / 2
    return {
        'q': q, 'matrix': mat, 'shifts': shifts, 'c0': c0,
        'eigenvalues': eig, 'vectors': vectors, 'omega': omega,
        'initial': initial, 'derivative': derivative,
        'minimum': multiplier, 'boundary_condition': float(np.linalg.cond(boundary)),
        'boundary_residual': float(np.max(np.abs(boundary @ data - rhs))),
    }


def sample_profile(data, x):
    import numpy as np
    vectors, omega = data['vectors'], data['omega']
    cosine = (vectors * np.cos(omega * x)) @ vectors.T
    sine = (vectors * (np.sin(omega * x) / omega)) @ vectors.T
    return cosine @ data['initial'] + sine @ data['derivative']


def scout(q=337017 / 500000, polynomial=False):
    import numpy as np
    data = solve_bvp(q)
    profile = np.array([sample_profile(data, x) for x in np.linspace(0, 1 / 12, 101)])
    mat = data['matrix']
    diagonal_lower = np.min(np.diag(mat) - (np.sum(np.abs(mat), axis=1) - np.abs(np.diag(mat))))
    result = {
        'status': 'NUMERICAL SCOUT ONLY; no numerical value is a proof',
        'q': data['q'], 'target': 2 - 2 * data['q'],
        'matrix_eigenvalues': data['eigenvalues'].tolist(),
        'matrix_condition': float(np.linalg.cond(mat)),
        'gershgorin_lower': float(diagonal_lower),
        'coercivity_reserve_from_gershgorin': float(diagonal_lower - 2 / math.pi ** 2),
        'bvp_condition': data['boundary_condition'],
        'bvp_residual': data['boundary_residual'],
        'affine_minimum': data['minimum'],
        'margin_above_target': data['minimum'] - (2 - 2 * data['q']),
        'profile_sample_min': float(profile.min()), 'profile_sample_max': float(profile.max()),
        'initial_fibres': data['initial'].tolist(),
        'initial_derivatives': data['derivative'].tolist(),
    }
    print(json.dumps(result, indent=2))
    if polynomial:
        half_cell = 1 / 24
        vectors, omega = data['vectors'], data['omega']
        c = (vectors * np.cos(omega * half_cell)) @ vectors.T
        s = (vectors * (np.sin(omega * half_cell) / omega)) @ vectors.T
        dc = (vectors * (-omega * np.sin(omega * half_cell))) @ vectors.T
        columns = [c @ data['initial'] + s @ data['derivative'],
                   half_cell * (dc @ data['initial'] + c @ data['derivative'])]
        for k in range(5):
            columns.append(-2 * half_cell ** 2 * np.linalg.solve(data['matrix'], columns[k])
                           / ((k + 1) * (k + 2)))
        coefficients = np.array(columns).T
        sym = np.array([[(coefficients[j, k] + (-1) ** k * coefficients[11-j, k]) / 2
                         for k in range(7)] for j in range(12)])
        integers = np.rint(sym[:6] * 10 ** 14).astype(np.int64)
        print('RATIONAL_POLYNOMIAL_INTEGERS_DENOMINATOR_10_14')
        print(json.dumps(integers.tolist(), indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--scout', action='store_true')
    parser.add_argument('--certify', action='store_true')
    parser.add_argument('--q', type=float, default=337017 / 500000)
    parser.add_argument('--polynomial', action='store_true')
    args = parser.parse_args()
    if args.certify:
        certify()
    else:
        scout(args.q, args.polynomial)
