"""C176 production driver. Usage: python3 run.py <cfg.json> [outdir]
cfg.json: a list of dicts with keys
  tag, L, N, wa, phi, [t2, V2, V3, blockade, gauge], chis (list, one entry per sweep), meas_chi (measure at every
  sweep whose chi >= meas_chi, and always at the last sweep), seed.
One jsonl line per sweep is appended to <outdir>/<tag>.jsonl, and a checkpoint is kept in <outdir>/ckpt/<tag>.pkl.
If the checkpoint exists the run resumes from it and skips the sweeps already recorded.
"""
import numpy as np, json, sys, os, time
from fibdmrg import Chain, DMRG
import analyze


def run_one(cfg, outdir):
    tag = cfg['tag']
    L, N = cfg['L'], cfg['N']
    w = cfg['wa'] * np.exp(1j * cfg['phi'])
    if abs(cfg['phi']) < 1e-12:
        w = cfg['wa']
    gauge = cfg.get('gauge', 1 if abs(cfg['phi'] - np.pi / 2) < 0.2 else 0)
    kw = {k: cfg[k] for k in ('t2', 'V2', 'V3', 'blockade') if k in cfg}
    ch = Chain(L, N, w, gauge=gauge, **kw)
    dm = DMRG(ch, seed=cfg.get('seed', 0))
    outf = os.path.join(outdir, f"{tag}.jsonl")
    ck = os.path.join(outdir, 'ckpt', f"{tag}.pkl")
    os.makedirs(os.path.dirname(ck), exist_ok=True)
    done = 0
    if os.path.exists(ck) and os.path.exists(outf):
        dm.load(ck)
        done = sum(1 for _ in open(outf))
    chis = cfg['chis']
    mchi = cfg.get('meas_chi', 10 ** 9)
    for n_sw, chi in enumerate(chis):
        if n_sw < done:
            continue
        t0 = time.time()
        E, tr = dm.sweep(chi, tol=cfg.get('tol', 1e-10))
        rec = dict(tag=tag, cfg=cfg, sweep=n_sw, chi=int(chi), E=E, trunc=tr, maxchi=int(dm.chi_bond.max()),
                   S=dm.S.tolist(), trunc_bond=dm.trunc.tolist(), secs=time.time() - t0)
        last = (n_sw == len(chis) - 1)
        if (chi >= mchi or last) and not cfg.get('nomeas', False):
            t1 = time.time()
            m = dm.measure()
            rec['meas'] = analyze.reduce_meas(m, L, N, store_G=last)
            rec['meas_secs'] = time.time() - t1
        dm.save(ck)
        with open(outf, 'a') as f:
            f.write(json.dumps(rec) + "\n")
        print(f"{tag} sweep {n_sw} chi={chi} E={E:.12f} trunc={tr:.2e} Smid={dm.S[L//2-1]:.5f} t={rec['secs']:.0f}s", flush=True)


if __name__ == "__main__":
    cfgs = json.load(open(sys.argv[1]))
    outdir = sys.argv[2] if len(sys.argv) > 2 else '.'
    if isinstance(cfgs, dict):
        cfgs = [cfgs]
    for c in cfgs:
        run_one(c, outdir)
