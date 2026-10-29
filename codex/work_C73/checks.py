"""Small coefficient controls for C73. No interacting XXZ flow computation."""

import itertools
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import quad


def smooth_r(s):
    s = abs(s)
    if s <= 1:
        return 1.0
    if s >= 2:
        return 0.0
    u = s - 1
    left, right = math.exp(-1/u), math.exp(-1/(1-u))
    return right/(left+right)


def cycle_sum(points, covariance):
    terms = []
    for perm in itertools.permutations(range(1, len(points))):
        order = (0,) + perm
        terms.append(np.prod([covariance(points[order[j]]-points[order[(j+1) % len(order)]])
                              for j in range(len(order))]))
    return sum(terms)


def cycle_checks():
    c = lambda z: 1/(2*math.pi*complex(z[0], -z[1]))
    a = lambda z: math.exp(-float(np.dot(z, z))/4)*c(z)
    points = np.array([[0., 0.], [1., 0.], [2., 0.], [3., 0.]])
    full = cycle_sum(points, c)
    high = cycle_sum(points, a)
    formula = 2/(2*math.pi)**4 * (-math.exp(-3)/3
                                 +math.exp(-2.5)/4+math.exp(-4.5)/12)
    assert abs(full) < 1e-18
    assert high.real > 0 and abs(high-formula) < 1e-18
    odd_rows = []
    for size in (3, 5):
        pts = np.array([[0.31*j, 0.13*j*j] for j in range(size)])
        value = cycle_sum(pts, a)
        assert abs(value) < 2e-17
        odd_rows.append(dict(length=size, residual=float(abs(value))))
    return dict(full_four_cycle_residual=float(abs(full)),
                high_four_cycle=float(high.real), formula=formula,
                formula_residual=float(abs(high-formula)), odd_cycles=odd_rows)


def ward_checks():
    A = 1.3
    residuals, omitted = [], []
    for k, q in ((np.array([.3, .7]), np.array([.8, -.2])),
                 (np.array([-.5, .4]), np.array([.2, .9]))):
        D = lambda z: A*z[1]-1j*z[0]
        f = lambda z: math.exp(-float(z@z))
        a = lambda z: (1-f(z))/D(z)
        lhs = D(q)*a(k)*a(k+q)
        remainder = f(k)*a(k+q)-f(k+q)*a(k)
        rhs = a(k)-a(k+q)+remainder
        residuals.append(float(abs(lhs-rhs)))
        omitted.append(float(abs(remainder)))
    assert max(residuals) < 3e-16
    assert min(omitted) > .01
    return dict(maximum_residual=max(residuals), retained_remainder_magnitudes=omitted)


def loop_checks():
    A = 1.2
    rows = []
    maximum_quadrature_error = 0.
    maximum_plateau_error = 0.
    for size in (2, 3, 4, 5):
        shifts = np.array([.37*j + .04*j*j for j in range(size)])
        nus = np.array([.23*j + .07*j*j for j in range(size)])
        differences = A*shifts-1j*nus
        coeffs = np.array([1/np.prod([differences[l]-differences[j]
                                     for l in range(size) if l != j])
                           for j in range(size)])

        def residue_value(k):
            return np.dot(coeffs, np.sign(k+shifts))/2

        sample_k = [-shifts[-1]-.6, .7]
        sample_k += [-(shifts[j]+shifts[j+1])/2 for j in range(size-1)]
        for k in sample_k:
            def integrand(omega):
                return 1/np.prod(A*(k+shifts)-1j*(omega+nus))/(2*math.pi)
            direct = (quad(lambda w: integrand(w).real, -np.inf, np.inf,
                           epsabs=2e-11, epsrel=2e-11)[0]
                      +1j*quad(lambda w: integrand(w).imag, -np.inf, np.inf,
                               epsabs=2e-11, epsrel=2e-11)[0])
            error = float(abs(direct-residue_value(k)))
            maximum_quadrature_error = max(maximum_quadrature_error, error)
            assert error < 2e-10
        assert abs(residue_value(sample_k[0])) < 5e-15
        assert abs(residue_value(sample_k[1])) < 5e-15

        breaks = sorted(-shifts)
        expected = sum((right-left)*residue_value((left+right)/2)
                       for left, right in zip(breaks, breaks[1:]))/(2*math.pi)
        cutoff_rows = []
        for cutoff in (shifts[-1]+.2, 2*shifts[-1]+1., 4*shifts[-1]+2.):
            actual = 0j
            for left, right in zip(breaks, breaks[1:]):
                weight = quad(lambda k: math.prod(smooth_r((k+s)/cutoff)**2
                                                  for s in shifts),
                              left, right, epsabs=1e-13)[0]
                actual += weight*residue_value((left+right)/2)/(2*math.pi)
            error = float(abs(actual-expected))
            maximum_plateau_error = max(maximum_plateau_error, error)
            assert error < 2e-14
            cutoff_rows.append(dict(cutoff=float(cutoff), plateau_residual=error))
        rows.append(dict(cycle_length=size, momentum_span=float(shifts[-1]),
                         plateau_checks=cutoff_rows))
    return dict(maximum_frequency_quadrature_error=maximum_quadrature_error,
                maximum_plateau_error=maximum_plateau_error, cases=rows)


def winding_checks():
    def a(t, x):
        return math.exp(-(t*t+x*x)/4)/(2*math.pi*complex(t, -x))
    rows = []
    for size in (3, 4):
        points = np.array([[.11*j, .17*j] for j in range(size)])
        for beta in (2., 3.):
            thermal = 1+0j
            for j in range(size):
                dt, dx = points[j]-points[(j+1) % size]
                thermal *= sum((-1.)**n*a(dt+n*beta, dx) for n in range(-16, 17))
            periodized, nonzero_windings = 0j, 0j
            for lifts in itertools.product(range(-8, 9), repeat=size-1):
                k = (0,)+lifts
                tree = 1+0j
                for j in range(size-1):
                    dt, dx = points[j]-points[j+1]
                    tree *= a(dt+(k[j]-k[j+1])*beta, dx)
                dt, dx = points[-1]-points[0]
                for ell in range(-16, 17):
                    term = (-1.)**ell*tree*a(dt+(k[-1]+ell)*beta, dx)
                    if ell == 0:
                        periodized += term
                    else:
                        nonzero_windings += term
            error = float(abs(thermal-periodized-nonzero_windings))
            assert error < 2e-12, (size, beta, error, thermal, periodized, nonzero_windings)
            rows.append(dict(cycle_length=size, beta=beta,
                             unfolding_residual=error,
                             thermal_defect_magnitude=float(abs(thermal-periodized))))
    return rows


if __name__ == '__main__':
    result = dict(scope='Chiral cycle algebra and fixed-order controls. Full TP remains unproved.',
                  cycles=cycle_checks(), ward=ward_checks(),
                  loops=loop_checks(), winding=winding_checks())
    Path(__file__).with_name('checks.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(dict(status='passed',
                         high_four_cycle=result['cycles']['high_four_cycle'],
                         contour_residual=result['loops']['maximum_frequency_quadrature_error'],
                         plateau_residual=result['loops']['maximum_plateau_error'],
                         winding_residual=max(row['unfolding_residual']
                                              for row in result['winding'])), indent=2))
