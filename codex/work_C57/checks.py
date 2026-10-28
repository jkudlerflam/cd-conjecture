"""Small independent controls of C57's partial results, not an RG simulation."""

import json
import math
from decimal import Decimal, localcontext
from pathlib import Path

from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import expit


def moments(t, mu, temperature):
    p = math.acos(-mu / (2 * t))
    def occupation(k):
        return expit((2 * t * math.cos(k) + mu) / temperature)
    n = quad(lambda k: occupation(k) / math.pi, 0, math.pi,
             points=[p], epsabs=2e-13, epsrel=2e-13)[0]
    m = quad(lambda k: math.cos(k) * occupation(k) / math.pi,
             0, math.pi, points=[p], epsabs=2e-13, epsrel=2e-13)[0]
    return n, m


def first_order_checks():
    out = []
    for t, p in [(0.5, math.acos(-1 / 3)), (0.8, 0.9),
                 (0.5, math.pi / 2)]:
        mu = -2 * t * math.cos(p)
        rho = p / math.pi
        m = math.sin(p) / math.pi
        a = rho - 0.5
        v = 2 * t * math.sin(p)
        ground = a*a - m*m
        predicted_gc = 2 * math.pi * t / (3 * v**3) * (m-a*math.cos(p))
        predicted_c = 1 / (12*t*t)
        u1_gc = 4*t/v*(m-a*math.cos(p))
        assert abs(predicted_gc - math.pi*u1_gc/(6*v*v)) < 1e-14
        rows = []
        errors_gc, errors_c = [], []
        for temperature in [0.02*t, 0.01*t, 0.005*t, 0.0025*t]:
            n_t, m_t = moments(t, mu, temperature)
            gc = ((n_t-0.5)**2-m_t*m_t-ground)/temperature**2
            mu_t = brentq(lambda z: moments(t, z, temperature)[0]-rho,
                          mu-0.05*t, mu+0.05*t, xtol=5e-15)
            n_c, m_c = moments(t, mu_t, temperature)
            canonical = ((n_c-0.5)**2-m_c*m_c-ground)/temperature**2
            errors_gc.append(abs(gc-predicted_gc))
            errors_c.append(abs(canonical-predicted_c))
            rows.append(dict(T=temperature, grand_canonical=gc,
                             canonical=canonical, mu_thermal=mu_t))
        # The leading error is T^2 after division by T^2 in both ensembles.
        assert errors_gc[-1] < errors_gc[0] / 8
        assert errors_c[-1] < errors_c[0] / 8
        assert errors_gc[-1] < 5e-5, (t, p, errors_gc, rows)
        assert errors_c[-1] < 5e-5, (t, p, errors_c, rows)
        out.append(dict(t=t, p=p, predicted_grand_canonical=predicted_gc,
                        predicted_canonical=predicted_c,
                        error_grand_canonical=errors_gc[-1],
                        error_canonical=errors_c[-1], rows=rows))
    return out


def four_history_check():
    # A nonlinear analytic flow with nonzero reference thermal forcing.
    # It checks the mixed-difference identity and suppression, not XXZ.
    theta, p = 0.5, 3
    D = Decimal
    initial = D('0.01')
    r, b, c, d = map(D, ['0.6', '0.04', '0.003', '0.02'])
    worst_identity, worst_ratio, terminal_ratios = 0.0, 0.0, []
    for H in [6, 12, 24, 40]:
        temperature = D(2)**(-H)
        x0, y0, xt, yt = initial, D(0), initial, D(0)
        for h in range(H):
            q = D(2)**(-h)
            tau = (temperature/q)**p
            eta = float(initial)*2.0**(-theta*h)
            z, a = x0-y0, yt-y0
            mixed = (xt-yt)-z
            def R(w, thermal):
                return r*w+b*w*w+thermal*(c+d*w)
            nx0, ny0 = R(x0, 0), R(y0, 0)
            nxt, nyt = R(xt, tau), R(yt, tau)
            # Use the polynomial form to avoid subtracting nearly equal R's.
            rhs = ((r+2*b*(y0+a+z)+d*tau)*mixed+b*mixed*mixed
                   +2*b*a*z+d*tau*z)
            lhs = (nxt-nyt)-(nx0-ny0)
            worst_identity = max(worst_identity, float(abs(lhs-rhs)))
            assert float(abs(z)) <= eta*(1+1e-12)
            # F_T/q^2 = w+w^2+tau*(0.07+0.13*w).
            # Its exact rectangular difference is used at very small T.
            scalar = ((1+2*(y0+a+z)+D('0.13')*tau)*mixed+mixed*mixed
                      +2*a*z+D('0.13')*tau*z)
            worst_ratio = max(worst_ratio, float(abs(scalar))/(eta*float(tau)))
            x0, y0, xt, yt = nx0, ny0, nxt, nyt
        terminal_eta = float(initial)*2.0**(-theta*H)
        terminal_ratios.append(float(abs(xt-yt))/terminal_eta)
    assert worst_identity < 5e-18
    assert worst_ratio < 2
    assert max(terminal_ratios) < 2
    return dict(max_mixed_identity_error=worst_identity,
                max_normalized_scalar_ratio=worst_ratio,
                terminal_defect_ratios=terminal_ratios)


if __name__ == '__main__':
    with localcontext() as ctx:
        ctx.prec = 90
        results = dict(first_order=first_order_checks(),
                       four_history=four_history_check(),
                       scope='First Taylor coefficient and abstract algebra only. No all-order RG claim.')
    Path(__file__).with_name('checks.json').write_text(json.dumps(results, indent=2)+'\n')
    print(json.dumps(results, indent=2))
