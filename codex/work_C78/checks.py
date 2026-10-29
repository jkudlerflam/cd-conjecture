"""Small algebra and quadrature controls for C78. Not an XXZ RG calculation."""

from fractions import Fraction as F
import itertools
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import quad


def matchings(items):
    if not items:
        yield ()
        return
    first = items[0]
    for j in range(1, len(items)):
        rest = items[1:j] + items[j+1:]
        for pairs in matchings(rest):
            yield ((first, items[j]),) + pairs


def ring_count_checks():
    rows = []
    for n in (2, 4, 6, 8):
        pairs = list(matchings(tuple(range(n))))
        connected = 0
        for left in pairs:
            for right in pairs:
                edges = left + right
                seen = {0}
                while True:
                    new = seen | {v for a, b in edges for u, v in ((a, b), (b, a))
                                  if u in seen}
                    if new == seen:
                        break
                    seen = new
                connected += len(seen) == n
        assert connected == math.factorial(n-1)
        coefficient = F(connected, math.factorial(n))
        assert coefficient == F(1, n)
        rows.append(dict(vertices=n, connected_pairings=connected,
                         exponential_series_coefficient=str(coefficient)))
    return rows


def perm_sign(p):
    return (-1)**sum(p[i] > p[j] for i in range(len(p)) for j in range(i+1, len(p)))


def determinant(A, indices):
    return sum((F(perm_sign(p))*math.prod(A[indices[i]][indices[p[i]]]
                                        for i in range(len(indices)))
                for p in itertools.permutations(range(len(indices)))), F(0))


def poly_product(left, right):
    out = {}
    for a, x in left.items():
        for b, y in right.items():
            if not a & b:
                out[a | b] = out.get(a | b, F(0)) + x*y
    return out


def fourth_cumulant_check():
    # Rational antisymmetric high covariance samples. The determinant is
    # the exact density moment of the finite Grassmann Gaussian algebra.
    data = ((F(1, 2), F(2, 5), F(1, 7), F(3, 8), F(2, 9), F(1, 3)),
            (F(2, 3), F(-1, 4), F(3, 7), F(1, 5), F(-2, 7), F(4, 9)))
    mats = []
    for values in data:
        A = [[F(0) for _ in range(4)] for _ in range(4)]
        for (i, j), value in zip(itertools.combinations(range(4), 2), values):
            A[i][j], A[j][i] = value, -value
        mats.append(A)
    moments = [{mask: determinant(A, [i for i in range(4) if mask >> i & 1])
                for mask in range(1, 16)} for A in mats]
    Zminus1 = {mask: moments[0][mask]*moments[1][mask] for mask in range(1, 16)}
    power, log_coefficient = {0: F(1)}, F(0)
    for order in range(1, 5):
        power = poly_product(power, Zminus1)
        log_coefficient += F((-1)**(order+1), order)*power.get(15, F(0))

    partitions = list(matchings(tuple(range(4))))
    pair_products, cumulants = [], []
    for A in mats:
        def bubble(pair):
            i, j = pair
            return -A[i][j]*A[j][i]
        pairs = [math.prod(bubble(pair) for pair in partition) for partition in partitions]
        c4 = -sum(math.prod(A[cycle[j]][cycle[(j+1) % 4]] for j in range(4))
                  for tail in itertools.permutations((1, 2, 3))
                  for cycle in [(0,)+tail])
        assert determinant(A, list(range(4))) == c4 + sum(pairs)
        pair_products.append(pairs)
        cumulants.append(c4)
    ring = sum(pair_products[0][i]*pair_products[1][j]
               for i in range(3) for j in range(3) if i != j)
    nonring = (cumulants[0]*cumulants[1]
               +cumulants[0]*sum(pair_products[1])
               +sum(pair_products[0])*cumulants[1])
    assert log_coefficient == ring + nonring
    assert nonring != 0
    return dict(connected_four_vertex_cumulant=str(log_coefficient),
                ring_sector=str(ring), nonring_sector=str(nonring),
                high_density_four_cumulants=list(map(str, cumulants)),
                scope='Exact finite Gaussian algebra, not an integrated continuum coefficient.')


def cutoff(p):
    s = abs(p)
    if s <= 1:
        return 1.0
    if s >= 2:
        return 0.0
    z = s-1
    a, b = math.exp(-1/z), math.exp(-1/(1-z))
    return b/(a+b)


def thermal_coefficient_checks():
    A, g = 1.3, .2
    predicted = -g*g/(48*math.pi*A**3)
    rows = []
    for T in (.4, .2, .1, .05):
        def integrand(x):
            bose = x*math.exp(-x)/(-math.expm1(-x)) if x else 1.0
            return bose*cutoff(T*x/A)**2
        bound = 2*A/T
        integral = quad(integrand, 0, bound, points=[A/T],
                        epsabs=2e-13, epsrel=2e-13)[0]
        coefficient = -g*g*integral/(8*math.pi**3*A**3)
        rows.append(dict(T=T, thermal_coefficient_over_T_squared=coefficient,
                         ratio_to_limit=coefficient/predicted))
    assert abs(rows[-1]['ratio_to_limit']-1) < 2e-10
    assert all(a['ratio_to_limit'] < b['ratio_to_limit']
               for a, b in zip(rows, rows[1:]))
    frequencies = []
    B = g/(2*math.pi)
    for p, chi in ((.3, 1.), (.8, .7), (1.4, .2)):
        u = math.sqrt(A*A-B*B*chi*chi)
        direct = quad(lambda omega: math.log1p(-B*B*chi*chi*p*p
                                               /(omega*omega+A*A*p*p)),
                      0, np.inf, epsabs=2e-13, epsrel=2e-12)[0]/(2*math.pi)
        target = .5*p*(u-A)
        residual = abs(direct-target)
        assert residual < 5e-16
        frequencies.append(dict(p=p, residual=residual))
    return dict(predicted_coefficient=predicted, thermal_integrals=rows,
                frequency_logarithm_checks=frequencies)


def shell_checks():
    A, B = 1.3, .12
    def annulus(r):
        if r <= 1 or r >= 2:
            return 0.0
        z = r-1
        return math.exp(4-1/(z*(1-z)))
    def symbol(omega, p, q):
        radius2 = omega*omega+p*p
        eta = annulus(math.sqrt(radius2)/q)
        if eta == 0:
            return 0.0
        denominator = omega*omega+(A*A-B*B*math.exp(-2*radius2))*p*p
        ratio = B*B*p*p*math.expm1(-2*radius2)/denominator
        return .5*eta*math.log1p(ratio)

    theta, weights = np.polynomial.legendre.leggauss(96)
    theta, weights = (theta+1)*math.pi, weights*math.pi
    rows = []
    for q in (.125, .0625, .03125):
        def radial(r):
            vals = np.array([symbol(q*r*math.cos(t), q*r*math.sin(t), q)
                             for t in theta])
            return q*q*r*float(vals@weights)/(2*math.pi)**2
        zero = quad(radial, 1, 2, epsabs=2e-15, epsrel=2e-12)[0]
        thermal_rows = []
        for ratio in (.2, .1, .05, .025):
            T = ratio*q
            total = 0.
            for n in range(-math.floor(q/(math.pi*T)), math.floor(q/(math.pi*T))+1):
                omega = 2*math.pi*n*T
                outer = math.sqrt(max(0., 4*q*q-omega*omega))
                inner = math.sqrt(max(0., q*q-omega*omega))
                if outer:
                    total += 2*quad(lambda p: symbol(omega, p, q), inner, outer,
                                    epsabs=2e-15, epsrel=2e-12)[0]
            thermal = T*total/(2*math.pi)
            difference = abs(thermal-zero)
            normalized = difference/(q**4*ratio**3)
            assert math.isfinite(normalized)
            thermal_rows.append(dict(T_over_q=ratio, difference=difference,
                                     normalized_power_three=normalized))
        assert thermal_rows[-1]['difference'] < thermal_rows[0]['difference']/100
        rows.append(dict(q=q, zero_temperature_shell=zero,
                         shell_over_q_four=zero/q**4, thermal=thermal_rows))
    scaled = [row['shell_over_q_four'] for row in rows]
    assert abs(scaled[-1]/scaled[-2]-1) < .02
    return rows


if __name__ == '__main__':
    result = dict(scope='Mixed ring and finite coefficient controls. Complete TP remains unresolved.',
                  pairings=ring_count_checks(), fourth_cumulant=fourth_cumulant_check(),
                  thermal=thermal_coefficient_checks(), gaussian_shells=shell_checks())
    Path(__file__).with_name('checks.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(dict(status='passed', ring_orders=[r['vertices'] for r in result['pairings']],
                         exact_fourth_cumulant_decomposition=True,
                         thermal_ratio=result['thermal']['thermal_integrals'][-1]['ratio_to_limit'],
                         largest_log_integral_residual=max(r['residual'] for r in
                                                          result['thermal']['frequency_logarithm_checks'])),
                     indent=2))
