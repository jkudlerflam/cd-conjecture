"""Finite controls of C71's ultraviolet lemmas, not an XXZ RG calculation."""

import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.special import expit


def r(s):
    s = abs(s)
    if s <= 1:
        return 1.0
    if s >= 2:
        return 0.0
    u = s - 1
    left, right = math.exp(-1 / u), math.exp(-1 / (1 - u))
    return right / (left + right)


def logarithm_check():
    A, delta = 1.3, 0.7
    rows = []
    for cutoff in (4., 8., 16., 32., 64., 128.):
        # Split the plateau and the smooth transition, so the largest scale
        # does not hide the integrable k=0 boundary layer from quadrature.
        plateau = quad(lambda k: -math.expm1(-2 * A * delta * k) / k,
                       0, cutoff, epsabs=1e-12, epsrel=1e-12)[0]
        edge = quad(lambda s: r(s)**4 * (-math.expm1(-2*A*delta*cutoff*s)) / s,
                    1, 2, epsabs=1e-12, epsrel=1e-12)[0]
        value = (plateau + edge) / (4 * math.pi * A)
        rows.append(dict(cutoff=cutoff, absolute_loop=value,
                         subtracted_constant=value-math.log(cutoff)/(4*math.pi*A)))
    expected = math.log(2) / (4 * math.pi * A)
    increments = [b['absolute_loop']-a['absolute_loop']
                  for a, b in zip(rows, rows[1:])]
    assert abs(increments[-1]-expected) < 2e-13
    assert abs(rows[-1]['subtracted_constant']-rows[-2]['subtracted_constant']) < 2e-13
    return dict(rows=rows, doubling_increments=increments,
                predicted_doubling_increment=expected)


def dynamic_checks():
    A = 1.3
    plateau_rows = []
    for cutoff in (1.5, 3., 6.):
        for p in (-0.8, -0.2, 0.3, 1.0):
            numerator = math.copysign(1., p) * quad(
                lambda k: r(k/cutoff)**2*r((k+p)/cutoff)**2,
                min(-p, 0), max(-p, 0), epsabs=1e-13)[0]/(2*math.pi)
            for omega in (0., 0.7, 3.):
                value = numerator/(A*p-1j*omega)
                target = p/(2*math.pi*(A*p-1j*omega))
                assert abs(value-target) < 1e-14
                plateau_rows.append(abs(value-target))

    # Independent Matsubara sum at fixed internal momentum. Symmetric
    # sums have a leading 1/N tail; Richardson removes it.
    rows = []
    for T, k, p, mode in ((.3, -.17, .4, 1), (.2, .5, -.3, 2),
                          (.4, -.4, .2, 0)):
        omega = 2*math.pi*mode*T
        values = []
        for size in (2048, 4096):
            frequencies = (2*np.arange(-size, size)+1)*math.pi*T
            value = -T*np.sum(1/((A*k-1j*frequencies)
                                 *(A*(k+p)-1j*(frequencies+omega))))
            values.append(value)
        extrapolated = 2*values[1]-values[0]
        expected = (expit(-A*k/T)-expit(-A*(k+p)/T))/(A*p-1j*omega)
        residual = float(abs(extrapolated-expected))
        assert residual < 3e-11
        rows.append(dict(T=T, k=k, p=p, bosonic_mode=mode,
                         extrapolated_residual=residual))
    return dict(plateau_max_residual=max(plateau_rows),
                plateau_cases=len(plateau_rows), matsubara=rows)


def chiral_derivative_check():
    # Evaluate the distributional identity away from its singular point
    # by a numerical derivative, independently of a Fourier bubble formula.
    A = 1.3
    rows = []
    for t, x in ((.23, .41), (-.7, .11), (.13, -1.2)):
        def c(y):
            return 1/(2*math.pi*(A*t-1j*y))
        h = 1e-4
        derivative = (c(x-2*h)-8*c(x-h)+8*c(x+h)-c(x+2*h))/(12*h)
        residual = abs((-1j/(2*math.pi))*derivative-c(x)**2)
        assert residual < 2e-12
        rows.append(dict(t=t, x=x, derivative_residual=residual))
    return rows


def image_checks():
    # Explicit continuum high covariance for Gaussian radial subtraction,
    # A=q0=1. This is a control of the singular-bound/image argument.
    # The compact plateau identity above uses the actual smooth regulator.
    def a(t, x):
        return math.exp(-(t*t+x*x)/4)/(2*math.pi*(t-1j*x))
    identity_residual = 0.
    for beta in (2., 3., 5.):
        for t, x in ((.13, .23), (-.2, .8), (.3, -.4)):
            terms = np.array([(-1.)**n*a(t+n*beta, x) for n in range(-8, 9)])
            lhs = terms.sum()**2-np.sum(terms**2)
            rhs = sum(terms[i]*terms[j] for i in range(len(terms))
                      for j in range(len(terms)) if i != j)
            identity_residual = max(identity_residual, float(abs(lhs-rhs)))
    assert identity_residual < 1e-14

    # Absolute single-winding integral, without the polynomial output
    # weight. Center coordinates exploit reflection symmetry. The only
    # two singularities are integrable; split at their projected position.
    rows = []
    for beta in (2., 3., 4., 6., 8.):
        half = beta/2
        def outer(u):
            def integrand(x):
                denominator = math.hypot(u-half, x)*math.hypot(u+half, x)
                return math.exp(-(u*u+x*x)/2)/denominator
            return quad(integrand, 0, 12, epsabs=2e-10, epsrel=2e-9)[0]
        value = sum(quad(outer, left, right, epsabs=3e-9, epsrel=3e-8)[0]
                    for left, right in ((0, half), (half, half+12)))
        value *= math.exp(-beta*beta/8)/(math.pi**2)
        rows.append(dict(beta=beta, absolute_single_winding=value,
                         normalized_power_three=beta**3*value))
    assert all(b['absolute_single_winding'] < arow['absolute_single_winding']
               for arow, b in zip(rows, rows[1:]))
    assert rows[-1]['normalized_power_three'] < rows[0]['normalized_power_three']/100
    return dict(finite_sum_identity_residual=identity_residual,
                gaussian_high_covariance=rows)


if __name__ == '__main__':
    result = dict(scope='Ultraviolet coefficient and image controls. TP remains unproved.',
                  logarithm=logarithm_check(), dynamic=dynamic_checks(),
                  chiral_derivative=chiral_derivative_check(), images=image_checks())
    Path(__file__).with_name('checks.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(dict(status='passed',
                         logarithm_increment=result['logarithm']['doubling_increments'][-1],
                         predicted_increment=result['logarithm']['predicted_doubling_increment'],
                         matsubara_max_residual=max(x['extrapolated_residual']
                                                   for x in result['dynamic']['matsubara']),
                         image_identity_residual=result['images']['finite_sum_identity_residual']),
                     indent=2))
