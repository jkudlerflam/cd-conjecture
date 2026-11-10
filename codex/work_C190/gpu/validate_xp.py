"""C190g gates G1 and G2 (Claude Opus 5.5 agent, 2026-11-09).

Usage: python3 validate_xp.py <numpy|cupy> [L] [--svd-dev] [--sweep32] [--only32] [--min-chi N]

For every case of validate.py (p153, rat, sym40 with and without the staggered gauge, xx1, xx2) at L = 10, both
parity sectors, random and product starts, runs the original z2dmrg.dmrg and z2dmrg_xp.dmrg on the chosen
backend with the same schedule, and prints (i) the errors of the xp run against sparse ED (as validate.py) and
(ii) the differences between the xp run and the original code.  With --sweep32 it also runs p153 at the
collaborator's field at L = 32 with both codes and compares every sweep.
"""
import os
import sys
import time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
if not os.path.exists(os.path.join(HERE, 'z2dmrg.py')):
    sys.path.insert(0, os.path.dirname(HERE))
import z2dmrg as Z0
import z2dmrg_xp as ZX
import models as Mo
import validate as V

backend = sys.argv[1]
L = int(sys.argv[2]) if len(sys.argv) > 2 and sys.argv[2].isdigit() else 10
ZX.USE_GPU = backend == 'cupy'
ZX.GPU_MIN_CHI = int(sys.argv[sys.argv.index('--min-chi') + 1]) if '--min-chi' in sys.argv else 0
ZX.SVD_HOST = '--svd-dev' not in sys.argv


def run(mod, mpo, gph, L, P, init, chi=64):
    dt = complex if mpo.dtype is complex else float
    if init == 'rand':
        psi = mod.random_mps(L, P, D0=4, seed=7 + P, dtype=dt)
    else:
        occ = [1, 5] if P == 0 else [3]
        psi = mod.product_mps(L, occ, dtype=dt)
    sw = [dict(chi=16, alpha=1e-3, res_tol=1e-8), dict(chi=32, alpha=1e-5, res_tol=1e-9),
          dict(chi=chi, alpha=0.0, res_tol=1e-11), dict(chi=chi, alpha=0.0, res_tol=1e-12),
          dict(chi=chi, alpha=0.0, res_tol=1e-12)]
    psi, hist = mod.dmrg(psi, mpo, sw)
    if mod is ZX:
        psi = psi.host()
    assert psi.check(), 'parity structure broken'
    E = Z0.energy(psi, mpo) + mpo.const
    var, _ = Z0.variance(psi, mpo)
    S = Z0.entropies(psi)
    m = Z0.measure(psi, mpo, gph=gph, gwin=6, cp_i0=2, cp_rmax=L)
    Gb = np.array(m['Gb_re']) + 1j * np.array(m['Gb_im'])
    Gf = np.array(m['Gf_re']) + 1j * np.array(m['Gf_im'])
    cp = np.array([v for r, v in m['CP']])
    return dict(E=E, var=var, S=S, Z=np.array(m['Z']), Pb=np.array(m['P_bond']), Gb=Gb, Gf=Gf, cp=cp,
                par=Z0.parity_expectation(psi), hist=hist)


def case(label, mpo, gph, h2s, h1s, L, nnn=None):
    H = V.ed_ham(L, h2s, h1s, nnn)
    worst_ed, worst_cmp = {}, {}
    for P in (0, 1):
        e_ed, psi_ed = V.ed_sector(H, L, P)
        S_ed, Z_ed, P_ed, Gb_ed, Gf_ed = V.ed_obs(psi_ed, L, 2, 6)
        cp_ed = np.array([v for r, v in V.ed_cp(psi_ed, L, 2)])
        for init in ('rand', 'prod'):
            a = run(Z0, mpo, gph, L, P, init)
            b = run(ZX, mpo, gph, L, P, init)
            ed = dict(E=abs(b['E'] - e_ed), S=np.max(np.abs(b['S'] - S_ed)), Z=np.max(np.abs(b['Z'] - Z_ed)),
                      Pb=np.max(np.abs(b['Pb'] - P_ed)), Gb=np.max(np.abs(b['Gb'] - Gb_ed)),
                      Gf=np.max(np.abs(b['Gf'] - Gf_ed)), CP=np.max(np.abs(b['cp'] - cp_ed)),
                      par=abs(b['par'] - (1 - 2 * P)), var=abs(b['var']))
            cmp = dict(E=abs(b['E'] - a['E']), S=np.max(np.abs(b['S'] - a['S'])),
                       Ehist=max(abs(x['E'] - y['E']) for x, y in zip(a['hist'], b['hist'])),
                       Z=np.max(np.abs(b['Z'] - a['Z'])), Gf=np.max(np.abs(b['Gf'] - a['Gf'])),
                       nsw=abs(len(a['hist']) - len(b['hist'])))
            for k, v in ed.items():
                worst_ed[k] = max(worst_ed.get(k, 0.0), v)
            for k, v in cmp.items():
                worst_cmp[k] = max(worst_cmp.get(k, 0.0), v)
            print(f"{label} P={P} {init}: vsED " + ' '.join(f"{k}={v:.1e}" for k, v in ed.items()) +
                  " | vsCPU " + ' '.join(f"{k}={v:.1e}" for k, v in cmp.items()), flush=True)
    return worst_ed, worst_cmp


def sweep32():
    Lc = 32
    mpo, gph, shift, info = Mo.candidate('p153', Lc, -19.97589125)
    sched = [[32, 1e-3, 1e-7, 2], [64, 1e-4, 1e-8, 2], [96, 0.0, 1e-9, 2]]
    sw = [dict(chi=c, alpha=a, res_tol=r, maxit=300, mmax=32) for c, a, r, n in sched for _ in range(n)]
    out = {}
    for name, mod in (('cpu', Z0), ('xp', ZX)):
        psi = mod.random_mps(Lc, 0, D0=8, seed=1, dtype=complex)
        t0 = time.time()
        _, hist = mod.dmrg(psi, mpo, sw)
        out[name] = hist
        print(f"L=32 {name} time {time.time() - t0:.1f}s", flush=True)
    wE, wS = 0.0, 0.0
    for x, y in zip(out['cpu'], out['xp']):
        dE = abs(x['E'] - y['E'])
        dS = abs(x['Sbar6'] - y['Sbar6']) if 'Sbar6' in x else float('nan')
        dS1 = max(abs(p - q) for p, q in zip(x['Sprof'], y['Sprof'])) if 'Sprof' in x else float('nan')
        wE = max(wE, dE)
        if dS1 == dS1:
            wS = max(wS, dS1)
        print(f"  sweep {x['sweep']} chi={x['chi']} alpha={x['alpha']} E_cpu={x['E']:.12f} dE={dE:.1e} "
              f"dSbar6={dS:.1e} dSprof={dS1:.1e} mv_cpu={x['mv_mean']:.2f} mv_xp={y['mv_mean']:.2f} "
              f"t_cpu={x['t']:.1f} t_xp={y['t']:.1f}", flush=True)
    print(f"L=32 worst dE={wE:.1e} worst dS(cut)={wS:.1e}", flush=True)


if __name__ == '__main__':
    t0 = time.time()
    print(f"backend={backend} SVD_HOST={ZX.SVD_HOST} GPU_MIN_CHI={ZX.GPU_MIN_CHI}", flush=True)
    if '--only32' in sys.argv:
        sweep32()
        sys.exit(0)
    allw = {}
    for name, gauge in (('p153', None), ('rat', None), ('sym40', 'stag'), ('sym40', None)):
        b = Mo.bond_family(name)
        hz = b['hz0']
        mpo, gph, shift, info = Mo.candidate(name, L, hz, gauge=gauge)
        h1o = [hz * Z0.ZZ_ + shift * Z0.I2] * L
        allw[f"{name}/{gauge}"] = case(f"{name}/{gauge}", mpo, gph, [info['h2']] * (L - 1), h1o, L)
    n_op = np.diag([0.0, 1.0])
    for legs in (1, 2):
        mpo, gph = Mo.xx_control(L, mu=-1.2, legs=legs)
        h2o = [-(np.kron(Z0.SP, Z0.SM) + np.kron(Z0.SM, Z0.SP))] * (L - 1) if legs == 1 else [np.zeros((4, 4))] * (L - 1)
        nnn = None if legs == 1 else [(-1.0, -1.0)] * (L - 2)
        allw[f"xx{legs}"] = case(f"xx{legs}", mpo, gph, h2o, [1.2 * n_op] * L, L, nnn=nnn)
    print('worst per case (vs ED | vs original z2dmrg.py):')
    for k, (we, wc) in allw.items():
        print(' ', k, ' '.join(f"{a}={v:.1e}" for a, v in we.items()), '|',
              ' '.join(f"{a}={v:.1e}" for a, v in wc.items()))
    if '--sweep32' in sys.argv:
        sweep32()
    print(f'time {time.time() - t0:.1f}s')
