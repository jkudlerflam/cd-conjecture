#!/usr/bin/env python3
"""Open-chain DMRG at filling 1/3 for the micro chain and its logical chains.

For each run: ground state in the charge sector Q0 = M/6 (down-spin density n = 1/3), its entropy profile and
left-block charge moments, the ground energies at Q0 + 1 and Q0 - 1, and the first excited state at Q0.
Micro runs start from the encoded logical ground state (chain.encode).

usage: run_obc.py NAME --model {micro,logical,logical_x} --M 48 --lam 0.1 --rho 2 --chis 32 64 128
                  [--chis_log 32 64] [--out data/obc/NAME.json]
Claude (Opus 5.5), 2026-10-28.
"""
import os, sys, json, time, argparse
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np
import chain as CH
import dmrg as DM


def build(model, M, lam, rho):
    Jp, Jz = CH.couplings(rho)
    Jxy, Jzz, Delta = CH.second_order(Jp, Jz)
    if model == 'micro':
        ch = CH.micro(M, lam, Jp, Jz)
    elif model == 'logical':
        ch = CH.logical(M, Jxy, Jzz, lam=lam)
    elif model == 'logical_x':
        Jx, Jzx, _ = CH.exact_two_cell(lam, Jp, Jz)
        ch = CH.logical(M, Jx / lam ** 2, Jzx / lam ** 2, lam=lam)
    J = 2 * abs(Jxy) * lam ** 2
    return ch, J, dict(Jp=Jp, Jz=Jz, Jxy=Jxy, Jzz=Jzz, Delta=Delta, J=J)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('name')
    ap.add_argument('--model', default='micro')
    ap.add_argument('--M', type=int, default=24)
    ap.add_argument('--lam', type=float, default=0.1)
    ap.add_argument('--rho', type=float, default=2.0)
    ap.add_argument('--chis', type=int, nargs='+', default=[32, 64, 128])
    ap.add_argument('--chis_log', type=int, nargs='+', default=[16, 32, 64])
    ap.add_argument('--etol', type=float, default=1e-10)
    ap.add_argument('--kmax', type=int, default=40)
    ap.add_argument('--sweeps_max', type=int, default=30)
    ap.add_argument('--sweeps_min', type=int, default=3)
    ap.add_argument('--inter_max', type=int, default=None, help='cap on sweeps at intermediate chi (micro: 6)')
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    assert a.M % 6 == 0
    ch, J, info = build(a.model, a.M, a.lam, a.rho)
    Q0 = a.M / 6
    out = a.out or os.path.join(HERE, 'data', 'obc', a.name + '.json')
    os.makedirs(os.path.dirname(out), exist_ok=True)
    data = dict(name=a.name, args=vars(a), info=info, Q0=Q0, L=ch.L, cell_len=ch.cell_len, sectors={})
    t0 = time.time()
    log = lambda s: print(f"[{time.time() - t0:7.0f}s] {s}", flush=True)
    w_pen = 3 * J
    # logical seeds (for micro) built on the plain logical chain with the same Jxy, Jzz
    if a.model == 'micro':
        # seed: the logical chain with the exact two-cell couplings (closest logical chain to the micro one)
        Jx, Jzx, _ = CH.exact_two_cell(a.lam, info['Jp'], info['Jz'])
        lch = CH.logical(a.M, Jx / a.lam ** 2, Jzx / a.lam ** 2, lam=a.lam)

    def seed(Q, ortho=None):
        if a.model != 'micro':
            return None
        Wl = lch.mpo(Q0=Q, w=w_pen)
        o = [] if ortho is None else [(2 * J, ortho)]
        dl, _, _ = DM.run_ground(Wl, lch.dims, a.chis_log, log=None, etol=1e-14, sweeps_min=4, ortho=o)
        return dl.M, CH.encode(dl.M)

    chis = a.chis
    gs = {}
    for dq in (0, 1, -1):
        Q = Q0 + dq
        Ws = ch.mpo(Q0=Q, w=w_pen)
        s = seed(Q)
        log(f"sector Q={Q}")
        dm, E, hist = DM.run_ground(Ws, ch.dims, chis, mps=None if s is None else s[1], log=log, etol=a.etol, sweeps_min=a.sweeps_min, inter_max=a.inter_max if a.inter_max else (6 if a.model == 'micro' else None),
                                    kmax=a.kmax, sweeps_max=a.sweeps_max)
        rec = dict(E=E, hist=hist)
        if dq == 0:
            gs['dm'] = dm; gs['seed'] = s
            rec['S'] = dm.S.tolist()
            m1, m2 = dm.block_moments(ch.charge)
            rec['QA'] = m1.tolist(); rec['QA2'] = m2.tolist()
            rec['q'] = dm.local(ch.charge).tolist()
        data['sectors'][str(dq)] = rec
        json.dump(data, open(out + '.tmp', 'w')); os.replace(out + '.tmp', out)
    # first excited state at Q0
    Ws = ch.mpo(Q0=Q0, w=w_pen)
    log("excited state at Q0")
    init = None
    if a.model == 'micro':
        s1 = seed(Q0, ortho=gs['seed'][0])
        init = s1[1]
    dm1, E1, hist1 = DM.run_ground(Ws, ch.dims, chis, mps=init, log=log, etol=a.etol, sweeps_min=a.sweeps_min, inter_max=a.inter_max if a.inter_max else (6 if a.model == 'micro' else None), kmax=a.kmax,
                                   sweeps_max=a.sweeps_max, ortho=[(2 * J, gs['dm'].M)], seed=7)
    data['sectors']['ex'] = dict(E=E1, hist=hist1, overlap=dm1.overlap(gs['dm'].M))
    E0 = data['sectors']['0']['E']
    Ep, Em = data['sectors']['1']['E'], data['sectors']['-1']['E']
    dc = Ep + Em - 2 * E0
    data['summary'] = dict(E0=E0, charge_gap=dc, E1=E1 - E0, K_gap=(E1 - E0) / dc, J=J, time_s=time.time() - t0)
    log(f"summary {data['summary']}")
    json.dump(data, open(out + '.tmp', 'w')); os.replace(out + '.tmp', out)


if __name__ == '__main__':
    main()
