"""Lightweight independent checks of the analytic formulas in PROOF.md."""
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad

from smooth_checks import gradient, profile


xx = []
for q in [.4, .2, .1, .05]:
    s = q/2
    integral = quad(lambda u: np.sin(u)**2/np.cos(u), -s, s,
                    epsabs=1e-15, epsrel=1e-12)[0]/(8*np.pi*np.sin(s))
    closed = (np.arctanh(np.sin(s))-np.sin(s))/(4*np.pi*np.sin(s))
    assert abs(integral-closed) < 1e-15
    xx.append(dict(q=q, second_order_deficit_per_site=integral,
                   kappa_eff=integral/(q*q/2),
                   kappa_eff_over_kappa=integral/(q*q/(48*np.pi)),
                   closed_form_error=float(abs(integral-closed))))

ferro = []
for length in [128, 512, 2048, 8192]:
    eps, t, b = .1, 1., .5
    q = 2*np.pi/length
    x = np.arange(length)
    w = profile(length, 'cosine', eps, 1)
    nz = b*np.sin(q*x)
    nx = np.sqrt(1-nz*nz)
    cost = np.sum(w*((nx-np.roll(nx, -1))**2+(nz-np.roll(nz, -1))**2))/8
    field = np.sum(t*(w-np.roll(w, 1))*nz)/2
    lower = float(-cost-field)
    conservative = t*b*eps*length*np.sin(q/2)/2-(1+eps)*length*b*b*q*q/(8*(1-b*b))
    assert lower >= conservative-1e-14
    ferro.append(dict(L=length, epsilon=eps, t=t, b=b, M=.5,
                      trial_deficit_lower_bound=lower,
                      analytic_deficit_lower_bound=float(conservative),
                      ratio_lower_bound=lower/(.5*gradient(w)),
                      G=gradient(w)))

# Independently verify that the numerical optimized representatives leave
# uniform spectra unchanged in the stored full-space random-bond runs.
data = json.loads(Path('small_results.json').read_text())
changes = []
for d in [2, 3, 4]:
    for draw in range(3):
        rows = [r for r in data['random'] if r['d'] == d and r['draw'] == draw and r['shape'] == 'cosine']
        change = abs(rows[0]['uniform_energy']-rows[1]['uniform_energy'])
        assert change < 1e-10
        changes.append(dict(d=d, draw=draw, uniform_energy_change=change))

result = dict(status='numerical', xx_susceptibility=xx,
              ferromagnetic_trial_bounds=ferro, gauge_spectrum_check=changes)
Path('scalar_results.json').write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps(result, indent=2))
