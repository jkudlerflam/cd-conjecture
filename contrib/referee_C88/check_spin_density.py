"""(5b) Pin cost versus vacancy density on one-cell and two-cell rings (K = infinity),
using check_spin.run. Single cell: squeezed ring of one cell of p data sites plus
marker, with v vacancies up to and beyond half density. Compares the pin cost with
J * cbar * 2 cos(pi/(p+2)), where cbar is the mean contact probability measured on
the pinned cell, and with the free-fermion half-density reference alpha_contact."""
import sys, os, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from check_spin import run

out = []
alpha = 0.5 - 2 / np.pi ** 2
for p, vmax in ((2, 6), (3, 6), (4, 4)):
    for v in range(0, vmax + 1):
        for t in (0.3, 1., 3.):
            if v == 0 and t != 1.:
                continue
            r = run(p, v, 1., t, 0., ncell=1)
            cb = float(np.mean(r['pinned_cell_contacts']))
            r['mean_contact'] = cb
            r['cost_over_J_cbar_2cos'] = r['pin_cost_1'] / (cb * 2 * np.cos(np.pi / (p + 2))) if cb > 0 else None
            r['material_density'] = (p + 1) / (p + 1 + v)
            out.append(r)
json.dump(dict(alpha_contact=alpha, runs=out), open(os.path.join(HERE, 'check_spin_density.json'), 'w'), indent=1)
