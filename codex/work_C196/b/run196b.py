#!/usr/bin/env python3
"""C196b driver (Claude Opus 5.5 agent, 2026-11-10). Usage: python3 run196b.py <cfg.json> <outdir>

cfg keys: tag, model ('sym40', 'p153', 'xx1', 'xx2', 'tfim'), hz | mu | h, q (cell sites), p (particles per cell,
even), init ('eff' effective-model window with 2p particles, 'prod' evenly spaced product state with balanced
legs, 'rand:<seed>'), chi_init, nsw_init, chis (list), min_it, max_it, tolE (1e-10), tolS, nfull (0), res_tol (1e-9),
maxit (300), sigma (1e-2), noise (0).
Backend from the environment as in C190g: C196B_BACKEND=cupy runs chi >= C196B_GPU_MIN_CHI on the GPU,
C196B_SVD_DEV=1 does the SVDs on the device.

Per chi the window is iterated until (I1) |e_t - e_{t-2}| < tolE (twice the change of the two-iteration mean of
e, the comparison that covers a state alternating between iterations, as in C189), |Sbar_t - Sbar_{t-2}| < tolS,
and every local solve of the last two iterations has residual below 1e-7, with at least min_it iterations. The
recorded e is the mean of the last two iterations.
If max_it is reached first the chi record is written with conv.I1 = false. Each chi starts from the state of the
previous one. Outputs: <outdir>/<tag>.jsonl (one record per finished chi), <outdir>/<tag>.iter.jsonl (one line per
iteration), <outdir>/ckpt/<tag>.npz (after every iteration ckpt_every and at the end of each chi; resume is
automatic).
"""
import os
import sys
import json
import time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import z2dmrg_xp as X     # noqa: E402
import idmrg196b as B     # noqa: E402
import models as Mo       # noqa: E402

X.USE_GPU = os.environ.get('C196B_BACKEND', 'numpy') == 'cupy'
GPU_MIN_CHI = int(os.environ.get('C196B_GPU_MIN_CHI', '128'))
X.SVD_HOST = os.environ.get('C196B_SVD_DEV', '0') != '1'


def warm_alpha(cfg, t):
    """Mixer amplitude of global iteration t: cfg['warm'] = [[t_end, alpha], ...] (default none)."""
    for t_end, a in cfg.get('warm', []):
        if t < t_end:
            return a
    return 0.0


def init_window(cfg, st, log):
    n, q, p = st.n, cfg['q'], cfg['p']
    kind, _, arg = cfg['init'].partition(':')
    info = {}
    if kind == 'eff':
        b = Mo.bond_family(cfg['model'])
        par = Mo.eff_params(b, cfg['hz'])
        w_eff = complex(par['w_over_t2'])
        Ms, n0, E_eff = B.eff_window(n, 2 * p, w_eff, chis=tuple(cfg.get('eff_chis', [8, 16, 64, 64, 64])))
        info = dict(E_eff=E_eff, w_eff=[w_eff.real, w_eff.imag], N=2 * p)
    elif kind == 'u1':
        # N-definite ground state of the open window for the controls (fibdmrg without blockade)
        kw = dict(t2=cfg['u1_t2'], V2=0.0, V3=0.0, blockade=False)
        Ms, n0, E_u1 = B.eff_window(n, 2 * p, cfg['u1_w'], chis=tuple(cfg.get('eff_chis', [8, 16, 64, 64, 64])), **kw)
        info = dict(E_u1=E_u1, N=2 * p)
    elif kind == 'prod':
        N = 2 * p
        occ = sorted({int(k * n / N) + (k % 2) for k in range(N)})
        assert len(occ) == N
        Ms, n0 = B.product_window(n, occ)
        info = dict(occ=occ)
    elif kind == 'rand':
        Ms, n0 = B.random_window(n, seed=int(arg or 0))
    else:
        raise ValueError(cfg['init'])
    st.set_window(Ms, n0)
    info['dens_start'] = st.window_density()
    info['N_start'] = float(sum(info['dens_start']))
    recs = st.sweep_window(cfg.get('chi_init', 32), nsw=cfg.get('nsw_init', 2), alpha=cfg.get('alpha_init', 0.0))
    info['init_sweeps'] = [dict(E=r['E'], res=r['res_max'], disc=r['disc'], N=float(sum(r['dens']))) for r in recs]
    log(dict(init=info))
    return info


def main():
    cfg = json.load(open(sys.argv[1]))
    outdir = sys.argv[2]
    os.makedirs(os.path.join(outdir, 'ckpt'), exist_ok=True)
    tag = cfg['tag']
    fout = os.path.join(outdir, tag + '.jsonl')
    fit = os.path.join(outdir, tag + '.iter.jsonl')
    ck = os.path.join(outdir, 'ckpt', tag + '.npz')
    mpo, W, e_off = B.build_model(cfg)
    st = B.CellIDMRG(mpo, W, cfg['q'], seed=cfg.get('seed', 0), noise=cfg.get('noise', 0.0), e_off=e_off)

    def log(d):
        with open(fit, 'a') as f:
            f.write(json.dumps(d) + '\n')

    done = set()
    if os.path.exists(fout):
        for line in open(fout):
            if line.strip():
                done.add(json.loads(line)['chi'])
    meta = dict(chi=None, it=0, hist=[])
    if os.path.exists(ck):
        meta = st.load(ck)
        log(dict(resume=True, t=st.t, meta_chi=meta.get('chi'), host=os.uname().nodename))
    else:
        log(dict(start=True, cfg=cfg, host=os.uname().nodename, backend=os.environ.get('C196B_BACKEND', 'numpy')))
        init_window(cfg, st, log)
    for chi in cfg['chis']:
        if chi in done:
            continue
        xm = np
        if X.USE_GPU and chi >= GPU_MIN_CHI:
            import cupy as cp
            xm = cp
        st.to(xm)
        hist = meta['hist'] if meta.get('chi') == chi else []
        t0 = time.time()
        conv = False
        while len(hist) < cfg.get('max_it', 300):
            r = st.iterate(chi, res_tol=cfg.get('res_tol', 1e-9), maxit=cfg.get('maxit', 300),
                           sigma=cfg.get('sigma', 1e-2), nfull=cfg.get('nfull', 0), alpha=warm_alpha(cfg, st.t))
            r.pop('dens', None) if not cfg.get('log_dens', True) else None
            hist.append(dict(e=r['e_site'], S=r['Sbar'], res=r['res_max'] if r['alpha'] == 0 else 1.0))
            log(r)
            if len(hist) >= max(cfg.get('min_it', 20), 3):
                # two-iteration comparison (a state may alternate between iterations, C189): the change of the
                # two-iteration mean of e is (e_t - e_{t-2}) / 2
                dE2 = abs(hist[-1]['e'] - hist[-3]['e'])
                dS = abs(hist[-1]['S'] - hist[-3]['S'])
                if dE2 < cfg.get('tolE', 1e-10) and dS < cfg.get('tolS', 1e-6) \
                        and hist[-1]['res'] < 1e-7 and hist[-2]['res'] < 1e-7:
                    conv = True
                    break
            if len(hist) % cfg.get('ckpt_every', 5) == 0:
                st.save(ck, dict(chi=chi, it=len(hist), hist=hist))
        last = hist[-3:]
        tm0 = time.time()
        meas = {}
        try:
            meas['tm'] = st.transfer_spectrum(k=cfg.get('tm_k', 6))
        except Exception as e:     # noqa: BLE001
            meas['tm_error'] = repr(e)
        dens, Nc, vN = st.cell_number()
        rec = dict(tag=tag, chi=int(chi), cfg=cfg, e_site=0.5 * (r['e_site'] + hist[-2]['e']), e_last=r['e_site'],
                   e_prev=hist[-2]['e'], Sbar=r['Sbar'],
                   S_cell=r['S_cell'], S_mid=r['S'], disc=r['disc'], D=r['D'], k0=r['k0'], res_max=r['res_max'],
                   r_pred=r['r_pred'], fid0=r['fid0'], Ncell=Nc, varN=vN, dens=dens,
                   conv=dict(I1=conv, iters=len(hist), dE1=abs(last[-1]['e'] - last[-2]['e']),
                             dE2=abs(last[-1]['e'] - last[0]['e']), dS=abs(last[-1]['S'] - last[0]['S'])),
                   secs=time.time() - t0, meas_secs=time.time() - tm0, meas=meas, t=st.t,
                   schmidt=[float(x) for x in np.sort(st.lam_prev)[::-1][:64]])
        with open(fout, 'a') as f:
            f.write(json.dumps(rec) + '\n')
        meta = dict(chi=None, it=0, hist=[])
        st.save(ck, dict(chi=chi, it=len(hist), hist=[], done=chi))
        tm = meas.get('tm', {})
        print(f"{tag} chi={chi} it={len(hist)} I1={conv} e={rec['e_site']:.11f} Sbar={rec['Sbar']:.5f} "
              f"xi_e={tm.get('xi_even')} xi_o={tm.get('xi_odd')} N={Nc:.4f} varN={vN:.4f} disc={rec['disc']:.1e} "
              f"t={rec['secs']:.0f}s", flush=True)


if __name__ == '__main__':
    main()
