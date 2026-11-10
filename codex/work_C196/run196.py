"""C196 driver (Claude Opus 5.5 agent, 2026-11-10).  python3 run196.py <cfg.json> <outdir>

cfg: tag, model (p153 | rat | sym40 | xx1 | xx2 | tfim), hz or mu or h, chis (stage list), and optional
nmin (20), win (20), tolE (1e-10), tolS (1e-5), maxit_stage ({chi: cap}), max_wall_stage (s), res_tol (1e-9),
maxit (300).

Stop rule of a chi stage (PROOF 4): at least nmin iterations, and over the last `win` iterations
max |e_t - e_{t-2}| <= tolE (per site) and max |S_t - S_{t-2}| <= tolS, where e_t = (E_t - E_{t-1}) / 2.
A stage that hits its iteration cap or wall budget ends unconverged and is flagged.

Outputs: <outdir>/<tag>.jsonl (one line per finished chi), <outdir>/logs/<tag>.iter.jsonl (every iteration),
<outdir>/ckpt/<tag>.npz (checkpoint, resumable), <outdir>/<tag>.json (header and status).
Environment: C196_BACKEND=cupy|numpy, C196_GPU_MIN_CHI (128), C196_SVD_DEV (1 on the device).
"""
import os
import sys
import time
import json
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import z2dmrg_xp as X          # noqa: E402
import idmrg196 as I           # noqa: E402

X.USE_GPU = os.environ.get('C196_BACKEND', 'numpy') == 'cupy'
X.GPU_MIN_CHI = int(os.environ.get('C196_GPU_MIN_CHI', '128'))
X.SVD_HOST = os.environ.get('C196_SVD_DEV', '0') != '1'


def backend_for(chi):
    if X.USE_GPU and chi >= X.GPU_MIN_CHI:
        return X._cp
    return np


def main(cfg_path, out):
    cfg = json.load(open(cfg_path))
    tag = cfg['tag']
    os.makedirs(os.path.join(out, 'logs'), exist_ok=True)
    os.makedirs(os.path.join(out, 'ckpt'), exist_ok=True)
    f_pts = os.path.join(out, tag + '.jsonl')
    f_it = os.path.join(out, 'logs', tag + '.iter.jsonl')
    f_ck = os.path.join(out, 'ckpt', tag + '.npz')
    f_hdr = os.path.join(out, tag + '.json')
    mpo, W, e_off = I.build_model(cfg)
    ex = I.exact_energy(cfg)
    nmin, win = cfg.get('nmin', 20), cfg.get('win', 20)
    tolE, tolS = cfg.get('tolE', 1e-10), cfg.get('tolS', 1e-5)
    caps = {int(k): v for k, v in cfg.get('maxit_stage', {}).items()}
    wall_cap = cfg.get('max_wall_stage', 6 * 3600)
    st = I.IDMRG(mpo, W, cutoff=cfg.get('cutoff', 1e-14), seed=cfg.get('seed', 0))
    stage, it_stage, hist = 0, 0, []
    if os.path.exists(f_ck):
        meta = st.load(f_ck)
        stage, it_stage, hist = meta['stage'], meta['it_stage'], meta['hist']
        print(f'# resume {tag} at stage {stage} (chi {cfg["chis"][min(stage, len(cfg["chis"]) - 1)]}), '
              f'iteration {it_stage}, t {st.t}', flush=True)
    hdr = dict(cfg=cfg, e_off=e_off, exact=ex, backend='cupy' if X.USE_GPU else 'numpy',
               gpu_min_chi=X.GPU_MIN_CHI, svd_host=X.SVD_HOST, host=os.uname()[1],
               cuda_visible=os.environ.get('CUDA_VISIBLE_DEVICES'), status='running')
    json.dump(hdr, open(f_hdr, 'w'), indent=1)
    fit = open(f_it, 'a')
    t_ck = time.time()
    while stage < len(cfg['chis']):
        chi = cfg['chis'][stage]
        st.to(backend_for(chi))
        cap = caps.get(chi, cfg.get('maxit_default', 3000))
        t_stage0 = time.time() - (hist[-1]['t_stage'] if hist and hist[-1].get('chi') == chi else 0.0)
        hist = [h for h in hist if h.get('chi') == chi]       # per-stage history only
        conv = False
        while True:
            t0 = time.time()
            rec = st.step(chi, res_tol=cfg.get('res_tol', 1e-9), maxit=cfg.get('maxit', 300),
                          alpha=I.warmup_alpha(st.nstep))
            rec['chi'] = chi
            rec['dt'] = time.time() - t0
            Eh = st.E
            rec['e'] = (Eh[-1] - Eh[-2]) / 2 + e_off if len(Eh) >= 2 else None
            it_stage += 1
            rec['it'] = it_stage
            rec['t_stage'] = time.time() - t_stage0
            hist.append(rec)
            if len(hist) >= 3 and hist[-3]['e'] is not None:
                rec['de2'] = abs(rec['e'] - hist[-3]['e'])
                rec['dS2'] = abs(rec['S'] - hist[-3]['S'])
            fit.write(json.dumps(rec) + '\n')
            if it_stage % 50 == 0:
                fit.flush()
            done = False
            if it_stage >= nmin and len(hist) >= win + 2:
                last = hist[-win:]
                if all('de2' in h and h['alpha'] == 0 for h in last):
                    deM = max(h['de2'] for h in last)
                    dSM = max(h['dS2'] for h in last)
                    if deM <= tolE and dSM <= tolS:
                        conv, done = True, True
            if it_stage >= cap or rec['t_stage'] > wall_cap:
                done = True
            if done:
                break
            if time.time() - t_ck > 600:
                st.save(f_ck, dict(stage=stage, it_stage=it_stage, hist=hist[-(win + 4):]))
                t_ck = time.time()
        fit.flush()
        # point record
        last = hist[-win:]
        t1 = time.time()
        xi = st.transfer_spectrum(k=cfg.get('tm_k', 6), m=cfg.get('tm_m', 60), tol=cfg.get('tm_tol', 1e-10))
        pt = dict(tag=tag, chi=chi, D=st.D, iters=it_stage, t_final=st.t, converged=conv,
                  e=hist[-1]['e'], e_mean2=(hist[-1]['e'] + hist[-2]['e']) / 2,
                  de2_max=max(h.get('de2', np.inf) for h in last), dS2_max=max(h.get('dS2', np.inf) for h in last),
                  res_max=max(h['res'] for h in last), r_pred_last=hist[-1]['r_pred'],
                  disc_max=max(h['disc'] for h in last), nmv_mean=float(np.mean([h['nmv'] for h in last])),
                  S_t=hist[-1]['S'], S_prev=hist[-1]['S_prev'], S=(hist[-1]['S'] + hist[-1]['S_prev']) / 2,
                  n=(hist[-1]['n1'] + hist[-1]['n2'] + hist[-2]['n1'] + hist[-2]['n2']) / 4,
                  n_sites=[hist[-1]['n1'], hist[-1]['n2']],
                  xi_even=xi['xi_even'], xi_odd=xi['xi_odd'], xi=xi['xi'], lam0=xi['lam0'], spec=xi['spec'],
                  t_tm=time.time() - t1, t_stage=hist[-1]['t_stage'])
        if ex is not None:
            pt['e_err'] = pt['e'] - ex['e']
        with open(f_pts, 'a') as f:
            f.write(json.dumps(pt) + '\n')
        print(f"{tag} chi={chi} D={st.D} it={it_stage} conv={conv} e={pt['e']:.12f} de2={pt['de2_max']:.1e} "
              f"dS2={pt['dS2_max']:.1e} res={pt['res_max']:.1e} S={pt['S']:.5f} xi_e={pt['xi_even']:.2f} "
              f"xi_o={pt['xi_odd']:.2f} lam0={pt['lam0']:.8f} n={pt['n']:.5f} disc={pt['disc_max']:.1e} "
              f"t={pt['t_stage']:.0f}s", flush=True)
        stage += 1
        it_stage = 0
        hist = hist[-(win + 4):]
        st.save(f_ck, dict(stage=stage, it_stage=0, hist=[]))
        t_ck = time.time()
    hdr['status'] = 'done'
    json.dump(hdr, open(f_hdr, 'w'), indent=1)


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
