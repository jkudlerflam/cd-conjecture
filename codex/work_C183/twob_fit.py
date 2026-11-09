#!/usr/bin/env python3
"""C177: fit effective two-body interactions to the full-chain two-magnon ring levels (twob_full.json), then compute
the infinite-ring (L = 400) two-body ground state at K = pi relative to 2 eps_min with the fitted interactions.
Models: (A) V2, V3 free; (B) V2, V3, V4 free; (C) V2 free, V3 bare; (D) V3 free, V2 bare.  Data: the two lowest bosonic
levels at K = pi for L = 14, 16, 18, 20.  Also evaluates C175's resonance function
D = 2 V2 V3 + 4 V2 + V3 u^2 + 4 V3 + 8 (units |t2| = 1, u = |w|/|t2| = 1; D < 0 bound side) for the fitted V2, V3."""
import json, numpy as np
from scipy.optimize import least_squares
import twob as T
R = json.load(open('twob_full.json'))
out = []
for r in sorted(set(x['r'] for x in R)):
    recs = [x for x in R if x['r'] == r]
    t2 = recs[0]['t2']; f = recs[0]['fit']; tm = T.tm_from_fit(f); c0 = f['c0']; em = T.band_min(tm, c0)
    data = [(x['L'], np.array(x['E2_full'][:2])) for x in recs]
    def model(p, keys):
        V = {2: -2.0 * t2, 3: -2 / 3 * t2}
        for k, v in zip(keys, p):
            V[k] = v * t2
        return V
    def resid(p, keys):
        V = model(p, keys)
        res = []
        for L, ef in data:
            e = T.eff2(L, np.pi, tm, c0, V, nev=2)
            res += list((e - ef) / t2)
        return res
    line = f"r={r:5.1f}"
    for name, keys, p0 in (('A', (2, 3), [-2, -2 / 3]), ('B', (2, 3, 4), [-2, -2 / 3, 0]), ('C', (2,), [-2]), ('D', (3,), [-2 / 3])):
        sol = least_squares(resid, p0, args=(keys,))
        V = model(sol.x, keys)
        Einf = T.eff2(400, np.pi, tm, c0, V, nev=1)[0]
        v2, v3 = V[2] / t2, V[3] / t2
        Dres = 2 * v2 * v3 + 4 * v2 + v3 + 4 * v3 + 8
        rms = float(np.sqrt(np.mean(np.square(sol.fun))))
        out.append(dict(r=r, model=name, V={k: V[k] / t2 for k in V}, rms=rms, E_K_pi_inf=(Einf - 2 * em) / t2, D=Dres))
        line += f" | {name}: V/|t2|={ {k: round(V[k] / t2, 4) for k in V} } rms={rms:.1e} (E-2emin)/|t2|(L=400)={(Einf - 2 * em) / t2:+.5f} D={Dres:+.3f}"
    print(line, flush=True)
json.dump(out, open('twob_fit.json', 'w'), indent=1)
