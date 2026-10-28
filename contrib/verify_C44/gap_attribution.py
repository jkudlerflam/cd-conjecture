"""Attribute C44's n=2 residual: rerun the scale-separated bond with the p=2
tunnelling gap set to C44's saved value (its plain eigensolve differs from the
reflection-split value by ~1e-7 relative) and compare again."""
import json
import numpy as np
import precise
from wrap import Charge

ref = np.load('../../codex/work_C44/source_bonds.npz')
cases = [c for c in json.load(open('../../codex/work_C44/source_results.json'))['cases'] if c['n'] == 2]
orig = Charge.__init__
out = []
for c in cases:
    gaps_c44 = c['charge_gaps']

    def patched(self, k, p, U, wvac, _g=gaps_c44):
        orig(self, k, p, U, wvac)
        E, V = self.eig
        E = E.copy()
        E[1] = E[0] + _g[p - 1]          # impose C44's tunnelling gap
        self.eig = (E, V)
        self.gap = _g[p - 1]
    res = {}
    for tag, init in (('own_gap', orig), ('C44_gap', patched)):
        Charge.__init__ = init
        B, _ = precise.precise_bonds(2, c['wvac'], s=c['s'])
        Charge.__init__ = orig
        res[tag] = max(float(np.abs(B[k] - ref['n2_w%s_t%s_%s' % (c['wvac'], c['ratio'], k)]).max()) for k in B)
    rec = dict(wvac=c['wvac'], ratio=c['ratio'], **res)
    out.append(rec)
    print(json.dumps(rec), flush=True)
json.dump(out, open('gap_attribution.json', 'w'), indent=1)
