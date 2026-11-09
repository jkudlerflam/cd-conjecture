#!/usr/bin/env python3
"""C183 (Claude Opus 5.5 agent, 2026-11-09): free-fermion controls for the entropy estimators, same DMRG (dmrg2.py)
and the same estimators (S-bar over the three middle cuts, dmrg_c.profile_fits) as fc.py.

  xx  : the XX chain h = XX + YY + (hz/2)(Z1 + 1Z), hz = -2 cos(0.17 pi) * 2 (density 0.17), c = 1.
  leg2: two decoupled XX legs as a d = 4 chain (site = (a, b)), each leg the xx chain above, c = 2.
        A chain of L/2 rungs is the zigzag image of an L-site qubit chain with two decoupled NNN legs, so the tag
        length L means L/2 rungs.
usage: ctrl183.py xx|leg2 L --chis ... --outdir DIR
"""
import os, sys, json, time, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import dmrg2 as D2
import dmrg_c as DC

SX = np.array([[0, 1], [1, 0]], float); SY = np.array([[0, -1j], [1j, 0]]); SZ = np.diag([1.0, -1.0]); I2 = np.eye(2)
HZ = 2 * np.cos(0.83 * np.pi)  # band 4 cos k - 2 hz, filled fraction 1 - arccos(hz/2)/pi = 0.17 (or 0.83, same entropy)


def xx_bond(hz=HZ):
    k = np.kron
    return (k(SX, SX) + k(SY, SY)).real + 0.5 * hz * (k(SZ, I2) + k(I2, SZ))


def leg2_bond(hz=HZ):
    """sites ordered (a_j, b_j, a_{j+1}, b_{j+1}); legs a and b each carry xx_bond between j and j+1."""
    hx = xx_bond(hz).reshape(2, 2, 2, 2)           # (s_a s_a' ; t_a t_a')
    I = np.eye(2)
    # h on a-leg: acts on (a_j, a_{j+1}), identity on (b_j, b_{j+1})
    ha = np.einsum('prPR,qQ,sS->pqrsPQRS', hx, I, I)
    hb = np.einsum('qsQS,pP,rR->pqrsPQRS', hx, I, I)
    return (ha + hb).reshape(16, 16)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('kind'); ap.add_argument('L', type=int)
    ap.add_argument('--chis', type=int, nargs='+', default=[32, 64, 96, 128, 128])
    ap.add_argument('--outdir', default='ctrl')
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    t0 = time.time()
    if a.kind == 'xx':
        h = xx_bond(); Ls = a.L
        edge = 0.5 * HZ * SZ
    else:
        h = leg2_bond(); Ls = a.L // 2
        edge = 0.5 * HZ * (np.kron(SZ, I2) + np.kron(I2, SZ))
    r = D2.dmrg(h, Ls, chis=tuple(a.chis), cutoff=1e-12, seed=0, verbose=True, res_tol=1e-8, edge_field=edge)
    S = np.asarray(r['S'])
    cpf, csm, mid = DC.profile_fits(S, Ls)
    c = Ls // 2
    out = dict(kind=a.kind, L=a.L, sites=Ls, E=r['E'], S=S.tolist(), Sbar=float(np.mean(S[c - 2:c + 1])), c_prof=cpf,
               c_sm=csm, trunc=r['hist'][-1][2], hist=r['hist'], time_s=time.time() - t0)
    fn = os.path.join(a.outdir, f"{a.kind}_{a.L}.json")
    json.dump(out, open(fn, 'w'), default=float)
    print(f"done {fn} Sbar={out['Sbar']:.4f} c_prof={cpf:.3f} trunc={out['trunc']:.1e} t={out['time_s']:.0f}s", flush=True)


if __name__ == '__main__':
    main()
