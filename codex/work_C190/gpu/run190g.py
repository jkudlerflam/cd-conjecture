#!/usr/bin/env python3
"""C190g GPU production driver (Claude Opus 5.5 agent, 2026-11-09): run190.py with z2dmrg_xp.

The sweeps run on the GPU (environment C190G_BACKEND=cupy, the default; numpy runs them on the host), sweeps
with chi below C190G_GPU_MIN_CHI on the host, the truncation SVD on the host unless C190G_SVD_DEV=1.  The output
JSON and the checkpoint have the format of run190.py, and the final measurements are the host code of z2dmrg.py.
Original docstring follows.

C190 production driver (Claude Opus 5.5 agent, 2026-11-09).

Usage: python3 run190.py <config.json> [<outdir>]

Config keys:
  tag, model ('p153', 'rat', 'sym40', 'sym80', 'xx1', 'xx2'), gauge (null or 'stag'), hz (candidates) or mu
  (controls), L, P (parity sector), init ('rand:<seed>', 'prod:<N>', 'file:<npz>', 'flip:<npz>'),
  sched: list of [chi, alpha, res_tol, nsweeps], final: {chi, res_tol, min_sweeps, max_sweeps, tolE, tolS},
  gwin (one-body window), variance (bool), maxit (Davidson iterations per local solve).
Writes <outdir>/<tag>.json after every sweep (history and, at the end, the measurements) and a checkpoint
<outdir>/ckpt/<tag>.npz after every sweep, from which a rerun resumes automatically.
"""
import sys
import os
import json
import time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
if not os.path.exists(os.path.join(HERE, 'z2dmrg.py')):
    sys.path.insert(0, os.path.dirname(HERE))
import z2dmrg as Z0
import z2dmrg_xp as Z
import models as Mo

Z.USE_GPU = os.environ.get('C190G_BACKEND', 'cupy') == 'cupy'
Z.GPU_MIN_CHI = int(os.environ.get('C190G_GPU_MIN_CHI', '0'))
Z.SVD_HOST = os.environ.get('C190G_SVD_DEV', '0') != '1'
PROFILE = os.environ.get('C190G_PROFILE', '0') == '1'


def build(cfg):
    L = cfg['L']
    m = cfg['model']
    if m in ('xx1', 'xx2'):
        mpo, gph = Mo.xx_control(L, cfg['mu'], legs=1 if m == 'xx1' else 2)
        return mpo, gph, 0.0
    mpo, gph, shift, info = Mo.candidate(m, L, cfg['hz'], gauge=cfg.get('gauge'))
    return mpo, gph, shift


def init_state(cfg, dt):
    L, P = cfg['L'], cfg['P']
    kind, _, arg = cfg['init'].partition(':')
    if kind == 'rand':
        return Z.random_mps(L, P, D0=8, seed=int(arg), dtype=dt)
    if kind == 'prod':
        N = int(arg)
        assert N % 2 == P
        occ = [int(round((k + 0.5) * L / N)) for k in range(N)]
        occ = sorted(set(min(L - 1, max(0, o)) for o in occ))
        assert len(occ) == N
        return Z.product_mps(L, occ, dtype=dt)
    if kind == 'file':
        psi = Z.MPS.load(arg)
        assert psi.P == P and psi.L == L
        return psi
    if kind == 'flip':
        psi = Z.MPS.load(arg)
        psi = Z.flip_site(psi, L // 2)
        assert psi.P == P and psi.L == L
        return psi
    raise ValueError(cfg['init'])


def main():
    cfg = json.load(open(sys.argv[1]))
    outdir = sys.argv[2] if len(sys.argv) > 2 else 'out'
    os.makedirs(os.path.join(outdir, 'ckpt'), exist_ok=True)
    tag = cfg['tag']
    fout = os.path.join(outdir, tag + '.json')
    fck = os.path.join(outdir, 'ckpt', tag + '.npz')
    L = cfg['L']
    mpo, gph, shift = build(cfg)
    dt = complex if mpo.dtype is complex else float
    maxit = cfg.get('maxit', 400)
    sweeps = []
    for chi, alpha, rt, ns in cfg['sched']:
        for _ in range(ns):
            sweeps.append(dict(chi=chi, alpha=alpha, res_tol=rt, maxit=maxit, mmax=cfg.get('mmax', 32)))
    fin = cfg['final']
    fixed = len(sweeps)
    for _ in range(fin['max_sweeps']):
        sweeps.append(dict(chi=fin['chi'], alpha=0.0, res_tol=fin['res_tol'], maxit=maxit, mmax=cfg.get('mmax', 32),
                           final=True))
    hist = []
    done = 0
    if os.path.exists(fck) and os.path.exists(fout):
        old = json.load(open(fout))
        if old.get('status') == 'done':
            print('already done', tag)
            return
        psi = Z.MPS.load(fck)
        hist = old['hist']
        done = len(hist)
        print(f'resuming {tag} after {done} sweeps', flush=True)
    else:
        psi = init_state(cfg, dt)
    t0 = time.time()
    E_off = mpo.const - L * shift
    hdr = dict(backend='cupy' if Z.USE_GPU else 'numpy', gpu_min_chi=Z.GPU_MIN_CHI, svd_host=Z.SVD_HOST,
               host=os.uname()[1], cuda_visible=os.environ.get('CUDA_VISIBLE_DEVICES'))
    if Z.USE_GPU:
        pr = Z._cp.cuda.runtime.getDeviceProperties(0)
        hdr['gpu'] = pr['name'].decode() if isinstance(pr['name'], bytes) else str(pr['name'])
    print('# c190g ' + json.dumps(hdr), flush=True)

    def stop(h):
        fh = [r for r in h if r.get('final')]
        if len(fh) < max(2, fin.get('min_sweeps', 2)):
            return False
        dE = abs(fh[-1]['E'] - fh[-2]['E'])
        dS = abs(fh[-1]['Sbar6'] - fh[-2]['Sbar6'])
        return dE < fin['tolE'] and dS < fin['tolS']

    def log(rec):
        rec['E_phys'] = rec['E'] + E_off
        if sweeps[len(hist)].get('final'):
            rec['final'] = True
        hist.append(rec)
        print(json.dumps({k: (round(v, 12) if isinstance(v, float) else v) for k, v in rec.items()
                          if k not in ('Zprof', 'Sprof')}), flush=True)
        if PROFILE:
            print('# timers ' + json.dumps({k: round(v, 3) for k, v in Z.TIMERS.items()}), flush=True)

    def ckpt(psi_, h):
        tmp = fck + '.tmp.npz'
        psi_.save(tmp, meta=dict(sweeps=len(hist)))
        os.replace(tmp, fck)
        json.dump(dict(cfg=cfg, hist=hist, status='running', shift=shift, E_off=E_off), open(fout + '.tmp', 'w'))
        os.replace(fout + '.tmp', fout)

    if done and stop(hist):
        todo = []
    else:
        todo = sweeps[done:]
    # the dmrg() log/ckpt/stop callbacks see the global history through closures
    if todo:
        psi, _ = Z.dmrg(psi, mpo, todo, log=lambda r: log(r), ckpt=ckpt,
                        stop=lambda _h: stop(hist), profile=PROFILE)
    # final measurements on the converged state
    tm = time.time()
    psi = psi.host() if isinstance(psi, Z.MPS) else psi       # final measurements: z2dmrg host code (re-exported)
    S = Z.entropies(psi)
    meas = Z.measure(psi, mpo, gph=gph, gwin=min(cfg.get('gwin', 48), L - 4), h2_list=mpo.h2, h1_list=mpo.h1_in)
    par = Z.parity_expectation(psi)
    res = dict(cfg=cfg, hist=hist, status='done', shift=shift, E_off=E_off, S=S.tolist(), meas=meas,
               parity=par, Sbar1=Z.sbar(S, L, 1), Sbar6=Z.sbar(S, L, 6), Sbar12=Z.sbar(S, L, 12),
               dims=psi.dims, t_total=time.time() - t0)
    if cfg.get('variance', True):
        tv = time.time()
        var, e1 = Z.variance(psi, mpo)
        res['variance'] = var
        res['E_mpo'] = e1
        res['E_phys'] = e1 + E_off
        res['t_var'] = time.time() - tv
    res['t_meas'] = time.time() - tm
    json.dump(res, open(fout + '.tmp', 'w'))
    os.replace(fout + '.tmp', fout)
    print('done', tag, f"E_phys={res.get('E_phys')} Sbar6={res['Sbar6']:.6f} var={res.get('variance')} "
          f"t={res['t_total']:.0f}s", flush=True)


if __name__ == '__main__':
    main()
