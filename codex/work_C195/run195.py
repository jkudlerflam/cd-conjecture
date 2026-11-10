#!/usr/bin/env python3
"""C195 production driver (Claude Opus 5.5 agent, 2026-11-10): run190g.py of C190g with three changes.

1. The solver is the C195 copy of z2dmrg_xp.py, with OFFLOAD on (environment C195_OFFLOAD=1, the default), so
   that only the active environments and MPS tensors sit on the device.
2. Stages that run to a stop rule.  Besides 'sched' (fixed stages [chi, alpha, res_tol, nsweeps]) and 'final',
   a config may have 'conv': a list of stages {chi, res_tol, min_sweeps, max_sweeps, tolE, tolS} that run between
   'sched' and 'final' and each end by the stop rule of the final stage of run190g (|dE| < tolE and
   |dSbar6| < tolS between consecutive sweeps of the stage, at least max(2, min_sweeps) sweeps).  The sweeps are
   grouped into segments, each ending with one stop-ruled stage, and dmrg() is called once per segment.  A config
   without 'conv' is one segment, the same single dmrg() call as run190g (this is what the gate compares).
   Every history record carries 'stage' (an index into the list of stages); final-stage records carry
   'final': True as in run190g.
3. The final measurements use the host MPS returned by dmrg() without another copy (memory at L = 512).

Output JSON, checkpoint format and resume behaviour are those of run190g: <outdir>/<tag>.json after every sweep,
<outdir>/ckpt/<tag>.npz after every sweep, a rerun resumes from the checkpoint and skips a finished run.
Usage: python3 run195.py <config.json> <outdir>
"""
import sys
import os
import json
import time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import z2dmrg as Z0
import z2dmrg_xp as Z
import models as Mo

Z.USE_GPU = os.environ.get('C190G_BACKEND', 'cupy') == 'cupy'
Z.GPU_MIN_CHI = int(os.environ.get('C190G_GPU_MIN_CHI', '0'))
Z.SVD_HOST = os.environ.get('C190G_SVD_DEV', '0') != '1'
Z.OFFLOAD = os.environ.get('C195_OFFLOAD', '1') == '1'
Z.OFFLOAD_TEST = os.environ.get('C195_OFFLOAD_TEST', '0') == '1'
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
    raise ValueError(cfg['init'])


def stages_of(cfg):
    """List of stages: fixed ones from 'sched', stop-ruled ones from 'conv', then 'final'."""
    st = []
    for chi, alpha, rt, ns in cfg['sched']:
        st.append(dict(kind='fix', chi=chi, alpha=alpha, res_tol=rt, n=ns))
    for c in cfg.get('conv', []):
        st.append(dict(kind='conv', alpha=0.0, n=c['max_sweeps'], **c))
    fin = cfg['final']
    st.append(dict(kind='conv', alpha=0.0, n=fin['max_sweeps'], final=True, **fin))
    return st


def stage_done(s, recs):
    """recs: the history records of stage s."""
    if s['kind'] == 'fix':
        return len(recs) >= s['n']
    if len(recs) >= s['n']:
        return True
    return rule_met(s, recs)


def rule_met(s, recs):
    if len(recs) < max(2, s.get('min_sweeps', 2)):
        return False
    dE = abs(recs[-1]['E'] - recs[-2]['E'])
    dS = abs(recs[-1]['Sbar6'] - recs[-2]['Sbar6'])
    return dE < s['tolE'] and dS < s['tolS']


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
    mmax = cfg.get('mmax', 32)
    stages = stages_of(cfg)
    hist = []
    if os.path.exists(fck) and os.path.exists(fout):
        old = json.load(open(fout))
        if old.get('status') == 'done':
            print('already done', tag)
            return
        psi = Z.MPS.load(fck)
        hist = old['hist']
        print(f'resuming {tag} after {len(hist)} sweeps', flush=True)
    else:
        psi = init_state(cfg, dt)
    t0 = time.time()
    E_off = mpo.const - L * shift
    hdr = dict(backend='cupy' if Z.USE_GPU else 'numpy', gpu_min_chi=Z.GPU_MIN_CHI, svd_host=Z.SVD_HOST,
               offload=Z.OFFLOAD, host=os.uname()[1], cuda_visible=os.environ.get('CUDA_VISIBLE_DEVICES'),
               job=os.environ.get('SLURM_JOB_ID'))
    if Z.USE_GPU:
        pr = Z._cp.cuda.runtime.getDeviceProperties(0)
        hdr['gpu'] = pr['name'].decode() if isinstance(pr['name'], bytes) else str(pr['name'])
    print('# c195 ' + json.dumps(hdr), flush=True)

    def recs_of(k):
        return [r for r in hist if r.get('stage') == k]

    def rss_gb():
        try:
            for line in open('/proc/self/status'):
                if line.startswith('VmRSS'):
                    return int(line.split()[1]) / 1e6
        except Exception:
            pass
        return None

    # segments: consecutive fixed stages closed by one stop-ruled stage
    segs, curseg = [], []
    for k, s in enumerate(stages):
        curseg.append(k)
        if s['kind'] == 'conv':
            segs.append(curseg)
            curseg = []
    # legacy histories (no 'stage' key) are not expected here
    for seg in segs:
        todo, labels = [], []
        for k in seg:
            s = stages[k]
            if stage_done(s, recs_of(k)):
                continue
            nrem = s['n'] - len(recs_of(k))
            for _ in range(nrem):
                todo.append(dict(chi=s['chi'], alpha=s['alpha'], res_tol=s['res_tol'], maxit=maxit, mmax=mmax))
                labels.append(k)
        if not todo:
            continue
        kconv = seg[-1]
        cnt = {'i': 0}

        def log(rec, labels=labels, cnt=cnt):
            k = labels[cnt['i']]
            cnt['i'] += 1
            rec['stage'] = k
            rec['E_phys'] = rec['E'] + E_off
            if stages[k].get('final'):
                rec['final'] = True
            rec['host_rss_GB'] = rss_gb()
            hist.append(rec)
            print(json.dumps({kk: (round(v, 12) if isinstance(v, float) else v) for kk, v in rec.items()
                              if kk not in ('Zprof', 'Sprof')}), flush=True)
            if PROFILE:
                print('# timers ' + json.dumps({kk: round(v, 3) for kk, v in Z.TIMERS.items()}), flush=True)

        def ckpt(psi_, h):
            tc = time.time()
            tmp = fck + '.tmp.npz'
            psi_.save(tmp, meta=dict(sweeps=len(hist)))
            os.replace(tmp, fck)
            json.dump(dict(cfg=cfg, hist=hist, status='running', shift=shift, E_off=E_off), open(fout + '.tmp', 'w'))
            os.replace(fout + '.tmp', fout)
            print(f'# ckpt {time.time() - tc:.1f}s', flush=True)
            if len(hist) == int(os.environ.get('C195_EXIT_AFTER', '-1')):   # test of the resume path only
                print('# exit for the resume test', flush=True)
                sys.exit(0)

        psi, _ = Z.dmrg(psi, mpo, todo, log=log, ckpt=ckpt, stop=lambda _h, kc=kconv: rule_met(stages[kc], recs_of(kc)),
                        profile=PROFILE)
    # final measurements on the converged state (host code of z2dmrg.py, no extra copy of the MPS)
    tm = time.time()
    if Z.USE_GPU:
        Z._cp.get_default_memory_pool().free_all_blocks()
    psi = Z0.MPS([Z.asnp(x) for x in psi.M], list(psi.n0), psi.P)
    S = Z.entropies(psi)
    print(f'# entropies {time.time() - tm:.0f}s rss {rss_gb()}', flush=True)
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
          f"t={res['t_total']:.0f}s t_meas={res['t_meas']:.0f}s", flush=True)


if __name__ == '__main__':
    main()
