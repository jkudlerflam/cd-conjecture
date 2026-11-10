"""C195 local check of the offload bookkeeping (Claude Opus 5.5 agent, 2026-11-10).  numpy only.

Usage: python3 validate_off.py [L]

Runs p153c at small L with the unmodified C190g solver (ref/z2dmrg_xp.py, numpy) and with the C195 copy in
OFFLOAD_TEST mode (offloaded sweeps with numpy as the device and boxed host copies, so that any arithmetic on a
host copy fails).  GPU_MIN_CHI = 16 in the C195 run, so the first sweeps run in memory and the later ones
offloaded, which tests the hand-over.  A second test restarts both codes from the same mid-run state with an
offloaded first sweep (the site-by-site canonicalization and environment build of a resume).  Prints the largest
differences per sweep in E, in the entropy at any cut, and in <Z_j>.
"""
import os
import sys
import importlib.util
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import z2dmrg_xp as ZN            # C195 copy
import models as Mo


def load_ref():
    spec = importlib.util.spec_from_file_location('z2dmrg_xp_ref', os.path.join(HERE, 'ref', 'z2dmrg_xp.py'))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


ZR = load_ref()
ZR.USE_GPU = False
ZN.USE_GPU = False
ZN.OFFLOAD = True
ZN.OFFLOAD_TEST = True
ZN.GPU_MIN_CHI = 16

L = int(sys.argv[1]) if len(sys.argv) > 1 else 20
mpo, gph, shift, info = Mo.candidate('p153', L, -19.97589125)


def sweeps(spec):
    out = []
    for chi, alpha, rt, ns in spec:
        out += [dict(chi=chi, alpha=alpha, res_tol=rt, maxit=300, mmax=32)] * ns
    return out


def compare(h1, h2, label):
    worst = dict(E=0.0, S=0.0, Z=0.0)
    for a, b in zip(h1, h2):
        dE = abs(a['E'] - b['E'])
        dZ = float(np.max(np.abs(np.array(a['Zprof']) - np.array(b['Zprof']))))
        dS = float(np.max(np.abs(np.array(a['Sprof']) - np.array(b['Sprof'])))) if 'Sprof' in a else 0.0
        print(f"  {label} sweep {a['sweep']} chi {a['chi']} alpha {a['alpha']}: E {a['E']:.12f} dE {dE:.1e} "
              f"dS {dS:.1e} dZ {dZ:.1e} mv {a['mv_mean']:.2f}/{b['mv_mean']:.2f}", flush=True)
        worst['E'] = max(worst['E'], dE); worst['S'] = max(worst['S'], dS); worst['Z'] = max(worst['Z'], dZ)
    assert len(h1) == len(h2)
    return worst


sp1 = sweeps([[8, 1e-3, 1e-7, 2], [16, 1e-4, 1e-8, 2], [24, 1e-5, 1e-8, 1], [32, 0.0, 1e-9, 2]])
p0 = ZR.random_mps(L, 0, D0=8, seed=1, dtype=complex)
pr, hr = ZR.dmrg(ZR.MPS([m.copy() for m in p0.M], list(p0.n0), 0), mpo, sp1)
pn, hn = ZN.dmrg(ZN.MPS([m.copy() for m in p0.M], list(p0.n0), 0), mpo, sp1)
w1 = compare(hr, hn, 'run')
dM = max(float(np.max(np.abs(a - b))) for a, b in zip(pr.M, pn.M))
print('test 1 (in memory then offloaded):', w1, 'final tensors max diff', dM, 'n0 equal', pr.n0 == pn.n0)

# test 2: restart both from the final state of test 1, first sweep offloaded (resume path)
sp2 = sweeps([[32, 0.0, 1e-9, 1], [40, 0.0, 1e-9, 1]])
pr2, hr2 = ZR.dmrg(ZR.MPS([m.copy() for m in pr.M], list(pr.n0), 0), mpo, sp2)
pn2, hn2 = ZN.dmrg(ZN.MPS([m.copy() for m in pr.M], list(pr.n0), 0), mpo, sp2)
w2 = compare(hr2, hn2, 'resume')
dM2 = max(float(np.max(np.abs(a - b))) for a, b in zip(pr2.M, pn2.M))
print('test 2 (resume, offloaded from the start):', w2, 'final tensors max diff', dM2)

# odd sector, product start
p1 = ZR.product_mps(L, [3, 8, 13], dtype=complex)
sp3 = sweeps([[8, 1e-3, 1e-7, 1], [16, 1e-4, 1e-8, 1], [32, 0.0, 1e-9, 2]])
pr3, hr3 = ZR.dmrg(ZR.MPS([m.copy() for m in p1.M], list(p1.n0), 1), mpo, sp3)
pn3, hn3 = ZN.dmrg(ZN.MPS([m.copy() for m in p1.M], list(p1.n0), 1), mpo, sp3)
w3 = compare(hr3, hn3, 'odd')
print('test 3 (odd sector, product start):', w3)
ok = all(w['E'] <= 1e-10 and w['S'] <= 1e-8 for w in (w1, w2, w3))
print('PASS' if ok else 'FAIL')
