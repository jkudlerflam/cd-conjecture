"""Referee ED cross-check (not done by the C139/C140 authors at these Delta).
(1) Bethe (C139/C140 solver) vs ED sector minima at minimizing flux L*pi, L=2..12, all 1<=n<=L/2,
    Delta in {0.25,0.5,0.9,0.98}, both parities; also a flux scan confirming that flux pi per bond minimizes.
(2) L|F_L(n)-L e_bulk| against K_rep(gamma) of C139.
(3) End-to-end L(E_2L-2E_L)/M with field and DM angle, L=2..6 (2L<=12), against 345600.
Floating point. Output referee_ed.json."""
import json
import numpy as np
from bethe_scan import solve
from check_c129 import energy as ed, bond_M
from continuum_rep import Model, e_bulk
import sys
sys.path.insert(0, '../work_C139')
from constants import K as Krep

out = dict(rows=[], flux=[], defect=[])
worst = 0.
for D in [0.25, 0.5, 0.9, 0.98]:
    M = Model(D); Kc = Krep(M.g)
    for L in range(2, 13):
        for n in range(1, L//2+1):
            E, dg = solve(L, n, D)
            Eed = ed(L, n, D)
            err = abs(float(E)-Eed); worst = max(worst, err)
            eb = e_bulk(M, n/L)
            out['rows'].append(dict(Delta=D, L=L, n=n, bethe=float(E), ed=Eed, err=err,
                                    LdF=float(L*(Eed-L*eb)), Krep=Kc, ratio=float(abs(L*(Eed-L*eb))/Kc)))
            if L <= 11:
                fl = [ed(L, n, D, a=np.pi+phi/L) for phi in np.linspace(0, 2*np.pi, 13)[1:-1]]
                out['flux'].append(dict(Delta=D, L=L, n=n, min_over_flux_minus_pi=float(min(fl)-Eed)))
    print(D, 'max Bethe-ED err', worst,
          'max ratio', max(r['ratio'] for r in out['rows'] if r['Delta'] == D),
          'min flux gap', min(r['min_over_flux_minus_pi'] for r in out['flux'] if r['Delta'] == D), flush=True)

for D in [0.25, 0.5, 0.9, 0.98]:
    for dm in [0.0, 0.3, 1.0]:
        J = float(np.hypot(1, dm)); th = float(np.arctan2(dm, 1)); jz = D*J
        for L in range(2, 7):
            small = np.array([ed(L, n, D, th, J) for n in range(L//2+1)])
            large = np.array([ed(2*L, n, D, th, J) for n in range(L+1)])
            for B in np.linspace(0, 2.2*J, 23):
                EL = np.min(small-B*(L/2-np.arange(len(small))))
                E2 = np.min(large-B*(L-np.arange(len(large))))
                Mb = bond_M(J, jz, B)
                out['defect'].append(dict(Delta=D, dm=dm, L=L, B=float(B), Q=float(L*(E2-2*EL)/Mb)))
    print(D, 'max L(E2L-2EL)/M', max(r['Q'] for r in out['defect'] if r['Delta'] == D), flush=True)

out['max_bethe_ed_err'] = worst
out['max_ratio_to_Krep'] = max(r['ratio'] for r in out['rows'])
out['min_flux_gap'] = min(r['min_over_flux_minus_pi'] for r in out['flux'])
out['max_defect_Q'] = max(r['Q'] for r in out['defect'])
print({k: out[k] for k in ['max_bethe_ed_err', 'max_ratio_to_Krep', 'min_flux_gap', 'max_defect_Q']})
json.dump(out, open('referee_ed.json', 'w'), indent=1)
