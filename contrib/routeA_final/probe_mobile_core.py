"""Floating-point probe of the C32 mobile (28-letter positive) core beyond n = 2.

Imports ../../codex/work_C32/mobile.py and band_checks.py read-only (run with
python3 -B so no bytecode is written there). For each (n, D, loop) it builds
the exact one-cell source eigensystem, the complete second-order two-cell bond
at the given conversion g, fits its quadratic part (hopping T, onsite mu) and
reports the eigenvalues of the Hermitian part of the phase-rotated hopping in
units of max|T|.

Precision diagnostic. At loop = 0 the onsite matrix mu vanishes identically
(C33 PROOF section 8, C32 PROOFS section 5), so the reported |mu|/|T| at
loop = 0 measures the double-precision floor of the full-in-g extraction,
which subtracts an O(1e-2) scalar background from an O(|T|) signal.
A row is called positive (indefinite) only when its smallest eigenvalue
is above (below minus) ten times this floor, and unresolved otherwise. The quadratic fit is a
Hilbert-Schmidt projection, so a large non-Gaussian remainder at g = 0.05
does not by itself invalidate it, but where the floor is large the
remainder is also noise.

numpy and scipy only. Runtime about one minute (n = 4 rows take ~20 s each).
Output: probe_mobile_core.json in this folder.
"""
import json
import os
import sys
import time

sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'codex', 'work_C32'))
import numpy as np  # noqa: E402
import mobile as m  # noqa: E402
import band_checks as bc  # noqa: E402

G, MASS = 0.05, 0.2
RUNS = [(2, 1.0, 0.0), (3, 1.0, 0.0), (3, 3.0, 0.0), (3, 10.0, 0.0),
        (3, 3.0, 1e-6), (4, 2.0, 0.0), (4, 3.0, 0.0), (4, 5.0, 0.0)]


def probe(n, D, loop):
    t0 = time.time()
    c = m.core(n, g=G, mass=MASS, express=D, loop=loop)
    me = m.measure(c)
    H = m.bond(me, me)
    T, mu, ng = bc.fit(H, n)
    R = T * np.exp(-1j * np.angle(T[0, 0]))
    ev = np.linalg.eigvalsh((R + R.conj().T) / 2)
    sc = float(np.max(abs(T)))
    return dict(n=n, D=D, loop=loop, g=G, mass=MASS, core_dim=int(c['d']),
                T_scale=sc, Rs_eigenvalues_over_T=(ev / sc).tolist(),
                onsite_over_T=float(np.max(abs(mu)) / sc),
                nonGaussian_over_T=float(ng / sc),
                scalar_per_state=float(abs(np.trace(H)) / H.shape[0]),
                seconds=round(time.time() - t0, 2))


if __name__ == '__main__':
    out = []
    for n, D, loop in RUNS:
        r = probe(n, D, loop)
        if loop == 0.0:
            floor = r['onsite_over_T']
            lo = min(r['Rs_eigenvalues_over_T'])
            r['floor_estimate_over_T'] = floor
            # trust a sign only when it exceeds ten times the floor estimate
            r['verdict'] = ('positive' if lo > 10 * floor else
                            'indefinite' if lo < -10 * floor else 'unresolved')
        print(json.dumps(r))
        out.append(r)
    with open(os.path.join(HERE, 'probe_mobile_core.json'), 'w') as f:
        json.dump({'status': 'numerical, double precision', 'runs': out}, f, indent=1)
        f.write('\n')
