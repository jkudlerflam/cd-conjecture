"""Roundoff floor of the scale-separated rebuild: repeat it in randomly
rotated core bases (an exact symmetry of (11)) and report the spread."""
import json
import numpy as np
import precise
from core46 import Core

orig_sector = Core.sector
orig_source = Core.source


def run_rotated(n, wvac, s, seed):
    rng = np.random.default_rng(seed)
    cache = {}

    def rot(self, bus):
        key = (id(self), bus)
        if key not in cache:
            idx = orig_sector(self, bus)
            M = rng.standard_normal((len(idx), len(idx))) + 1j * rng.standard_normal((len(idx), len(idx)))
            Q, _ = np.linalg.qr(M)
            cache[key] = (idx, Q)
        return cache[key]
    # rotate bonds, ends and sources inside each bus sector by the same unitary
    orig_init = precise.PreciseCell.__init__

    def init(self, core, ch, s_, bus, groups=(1, 1)):
        idx, Q = rot(core, bus)
        class RotCore:
            pass
        rc = RotCore()
        full = core.dim
        U = np.eye(full, dtype=complex)
        U[np.ix_(idx, idx)] = Q
        import scipy.sparse as sp
        Us = sp.csr_matrix(U)
        rc.bonds = [Us.conj().T @ b @ Us for b in core.bonds]
        rc.aL, rc.bR = Us.conj().T @ core.aL @ Us, Us.conj().T @ core.bR @ Us
        rc.sector = lambda b: orig_sector(core, b)
        rc.source = lambda end, b, sg: Us.conj().T @ orig_source(core, end, b, sg)
        orig_init(self, rc, ch, s_, bus, groups)
        self.core = rc
    precise.PreciseCell.__init__ = init
    try:
        B, ex = precise.precise_bonds(n, wvac, s=s)
    finally:
        precise.PreciseCell.__init__ = orig_init
    return B


if __name__ == '__main__':
    cases = json.load(open('../../codex/work_C44/source_results.json'))['cases']
    ref = np.load('../../codex/work_C44/source_bonds.npz')
    out = []
    for c in cases:
        B0, _ = precise.precise_bonds(c['n'], c['wvac'], s=c['s'])
        spreads = []
        for seed in (1, 2):
            B = run_rotated(c['n'], c['wvac'], c['s'], seed)
            spreads.append(max(float(np.abs(B[k] - B0[k]).max()) for k in B))
        dC = max(float(np.abs(B0[k] - ref['n%d_w%s_t%s_%s' % (c['n'], c['wvac'], c['ratio'], k)]).max()) for k in B0)
        rec = dict(n=c['n'], wvac=c['wvac'], ratio=c['ratio'], rotation_spread=max(spreads), C44_minus_rebuild=dC)
        out.append(rec)
        print(json.dumps(rec), flush=True)
    json.dump(out, open('roundoff_check.json', 'w'), indent=1)
