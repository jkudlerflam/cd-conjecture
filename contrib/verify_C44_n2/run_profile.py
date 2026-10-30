"""Open-chain DMRG entanglement profile of the projected n = 2 chain (block-sparse, bsdmrg.py).
python3 -B run_profile.py MODEL L CHI [SWEEPS]  -> results/MODEL_L{L}_chi{CHI}.json
Fixed charges: payload half filling N_p = L; length number N_l = round(nu L) (0 for C0, C1).
For M0 and C0 the exact fixed-N free-fermion profile is stored alongside."""
import sys
import os
import json
import time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from model2 import build
from bsdmrg import local_charges, mpo_charged, finite_dmrg, QSHIFT


def exact_fixedN(name, L, info, T, Nl):
    def prof(H, N, nm):
        e, v = np.linalg.eigh(H)
        occ = v[:, :N]
        C = occ @ occ.conj().T
        out = []
        for l in range(1, L):
            ev = np.linalg.eigvalsh(C[:nm * l, :nm * l])
            ev = ev[(ev > 1e-15) & (ev < 1 - 1e-15)]
            out.append(float(-(ev * np.log(ev) + (1 - ev) * np.log(1 - ev)).sum()))
        return np.array(out), float(e[:N].sum())
    H = np.zeros((2 * L, 2 * L), complex)
    for i in range(L - 1):
        H[2 * i:2 * i + 2, 2 * i + 2:2 * i + 4] = T
        H[2 * i + 2:2 * i + 4, 2 * i:2 * i + 2] = T.conj().T
    S, E = prof(H, L, 2)
    if name == 'M0':
        Hl = np.diag([-info['kappa']] * L) - info['jl'] * (np.eye(L, k=1) + np.eye(L, k=-1))
        Sl, El = prof(Hl, Nl, 1)
        S, E = S + Sl, E + El
    return S, E


if __name__ == '__main__':
    name, L, chi = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    nsw = int(sys.argv[4]) if len(sys.argv) > 4 else 6
    h2, h1, info, ops, T = build(name)
    qs = local_charges()
    W, delta, err, viol = mpo_charged(h2, h1, 8, qs)
    Nl = 0 if name in ('C0', 'C1') else int(round(info['nu'] * L))
    sched = [c for c in (64, 128, 256, 512) if c < chi] + [chi] * nsw
    t0 = time.time()
    o = finite_dmrg(W, delta, L, L * QSHIFT + Nl, sched, qs, lanczos_k=100,
                    diag_ops={'nl': np.diag(ops['nl']).real, 'nz': np.diag(ops['nz']).real}, verbose=True)
    o['S'] = list(o['S'])
    if name in ('M0', 'C0'):
        Sx, Ex = exact_fixedN(name, L, info, T, Nl)
        o['S_exact'] = list(Sx)
        o['E_exact'] = Ex
        o['max_dev_from_exact'] = float(np.abs(np.array(o['S']) - Sx).max())
        print('max dev from exact', o['max_dev_from_exact'], 'dE', o['E'] - Ex)
    o.update(info=info, Nl=Nl, mpo_D=int(W.shape[0]), mpo_err=err, charge_violation=viol, seconds=time.time() - t0)
    os.makedirs('results', exist_ok=True)
    fn = 'results/%s_L%d_chi%d.json' % (name, L, chi)
    json.dump(o, open(fn, 'w'), indent=1, default=float)
    print(fn, o['E'], o['seconds'])
