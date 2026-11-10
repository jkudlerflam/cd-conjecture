#!/usr/bin/env python3
"""C190 Typhon pipeline (Claude Opus 5.5 agent, 2026-11-09).

Usage: python3 pipeline190.py --node A|B --slots 8 --threads 6 [--out out] [--dry]

A simple scheduler on one node: tasks run as subprocesses of run190.py, longest first, each with a fixed number
of BLAS threads.  Field tuning (registered in RESULT.md): for each (bond, target density) three L = 128 even
runs at mu_eff = mu* - 0.03, mu*, mu* + 0.03 (mu* from C176's effective-model E(N +- 1) at the target), then a
least-squares line n_mid(mu_eff) through the three final middle-half densities fixes mu_eff for the target n_mid,
and h_z follows from the second-order map (models.hz_for_mu_eff).  The production runs of that key start then.
"""
import os
import sys
import json
import time
import argparse
import subprocess
import numpy as np
import models as Mo

TARGETS = {'n316': dict(n_mid=0.1957, mu=-2.165), 'n18': dict(n_mid=0.1326, mu=-2.215)}
SCHED = {
    64: dict(sched=[[32, 1e-3, 1e-7, 2], [64, 1e-4, 1e-8, 2], [128, 1e-5, 1e-8, 2], [192, 0.0, 1e-8, 2],
                    [256, 0.0, 1e-9, 2]],
             final=dict(chi=384, res_tol=1e-9, min_sweeps=3, max_sweeps=8, tolE=1e-9, tolS=2e-4)),
    128: dict(sched=[[32, 1e-3, 1e-7, 2], [64, 1e-4, 1e-8, 2], [128, 1e-4, 1e-8, 2], [192, 1e-5, 1e-8, 2],
                     [256, 0.0, 1e-8, 2], [384, 0.0, 1e-9, 2]],
              final=dict(chi=512, res_tol=1e-9, min_sweeps=3, max_sweeps=8, tolE=1e-9, tolS=3e-4)),
    256: dict(sched=[[32, 1e-3, 1e-7, 2], [64, 1e-4, 1e-8, 2], [128, 1e-4, 1e-8, 2], [192, 1e-5, 1e-8, 2],
                     [256, 1e-5, 1e-8, 2], [384, 0.0, 1e-8, 2], [512, 0.0, 1e-9, 2]],
              final=dict(chi=768, res_tol=1e-9, min_sweeps=2, max_sweeps=6, tolE=2e-9, tolS=5e-4)),
}
TUNE_SCHED = dict(sched=[[32, 1e-3, 1e-7, 2], [64, 1e-4, 1e-8, 2], [128, 1e-5, 1e-8, 2], [192, 0.0, 1e-8, 1]],
                  final=dict(chi=192, res_tol=1e-8, min_sweeps=2, max_sweeps=4, tolE=1e-8, tolS=1e-3))
COLLAB_HZ = {'p153': -19.97589125, 'rat': -19.9759}


def free_control_mu(L, n, legs):
    """Chemical potential halfway between the last filled and first empty open-chain level (unique GC ground
    state), and the resulting parity."""
    Lleg = L if legs == 1 else L // 2
    eps = np.sort(-2 * np.cos(np.pi * np.arange(1, Lleg + 1) / (Lleg + 1)))
    N = int(round(n * Lleg))
    mu = (eps[N - 1] + eps[N]) / 2
    Ntot = N * legs
    return float(mu), Ntot % 2, Ntot


def cand_task(key, bond, L, P, init, hz, gauge, tag=None, extra=None):
    sc = SCHED[L]
    cfg = dict(tag=tag or f"{key}_L{L}_P{P}_{init.replace(':', '')}", model=bond, gauge=gauge, hz=hz, L=L, P=P,
               init=init, sched=sc['sched'], final=sc['final'], gwin=48 if L >= 128 else 32,
               variance=(L <= 128), maxit=300)
    if extra:
        cfg.update(extra)
    return cfg


def cost(cfg):
    cx = 4 if cfg['model'] in ('p153', 'rat') else 1
    return cfg['L'] * cfg['final']['chi'] ** 3 * cx


def tasks_for_node(node):
    """Returns (tuning keys, immediate tasks, deferred task makers keyed by tuning key)."""
    tune, now, later = {}, [], {}

    def add_tune(key, bond, tgt, gauge):
        mu0 = TARGETS[tgt]['mu']
        b = Mo.bond_family(bond)
        runs = []
        for k, dmu in enumerate((-0.03, 0.0, 0.03)):
            mu = mu0 + dmu
            hz = Mo.hz_for_mu_eff(b, mu)
            cfg = dict(tag=f"tune_{key}_{k}", model=bond, gauge=gauge, hz=hz, L=128, P=0, init='rand:11',
                       sched=TUNE_SCHED['sched'], final=TUNE_SCHED['final'], gwin=32, variance=False, maxit=300,
                       mu_eff=mu)
            runs.append(cfg)
        tune[key] = dict(bond=bond, target=TARGETS[tgt]['n_mid'], runs=runs, gauge=gauge)

    def add_prod(key, bond, gauge, inits):
        def make(hz):
            out = []
            for L in (64, 128, 256):
                for P in (0, 1):
                    for init in inits(L, P):
                        out.append(cand_task(key, bond, L, P, init, hz, gauge))
            return out
        later[key] = make

    def n_init(n, L, P):
        N = int(round(n * L))
        if N % 2 != P:
            N += 1
        return f'prod:{N}'

    if node == 'A':
        # collaborator's field, phi = 1.53 unrounded: direct reproduction of the L = 128 downturn
        for L in (64, 128, 256):
            for P in (0, 1):
                inits = ['rand:1'] + ([n_init(0.22, L, P)] if L <= 128 else [])
                for init in inits:
                    now.append(cand_task('p153c', 'p153', L, P, init, COLLAB_HZ['p153'], None))
        add_tune('p153_n316', 'p153', 'n316', None)
        add_prod('p153_n316', 'p153', None,
                 lambda L, P: ['rand:1'] + ([n_init(0.1957, L, P)] if L == 128 else []))
    elif node == 'B':
        add_tune('sym40_n316', 'sym40', 'n316', 'stag')
        add_prod('sym40_n316', 'sym40', 'stag',
                 lambda L, P: ['rand:1'] + ([n_init(0.1957, L, P)] if L >= 128 else []))
        add_tune('sym40_n18', 'sym40', 'n18', 'stag')
        add_prod('sym40_n18', 'sym40', 'stag', lambda L, P: ['rand:1'])
        add_tune('sym80_n316', 'sym80', 'n316', 'stag')
        add_prod('sym80_n316', 'sym80', 'stag', lambda L, P: ['rand:1'])
        add_tune('rat_n316', 'rat', 'n316', None)
        add_prod('rat_n316', 'rat', None, lambda L, P: ['rand:1'])
        # free controls at the same sizes and chi schedule, density 0.196 (per leg for the two-leg control)
        for legs in (1, 2):
            for L in (64, 128, 256):
                mu, P, Ntot = free_control_mu(L, 0.1957, legs)
                sc = SCHED[L]
                now.append(dict(tag=f"xx{legs}_L{L}_P{P}_rand1", model=f'xx{legs}', mu=mu, L=L, P=P, init='rand:1',
                                sched=sc['sched'], final=sc['final'], gwin=48 if L >= 128 else 32, variance=False,
                                maxit=300, N_exact=Ntot))
    return tune, now, later


def finish_tuning(key, t, outdir):
    pts = []
    for cfg in t['runs']:
        r = json.load(open(os.path.join(outdir, cfg['tag'] + '.json')))
        if r.get('status') != 'done':
            return None
        dens = np.array(r['meas']['dens'])
        L = cfg['L']
        pts.append((cfg['mu_eff'], float(np.mean(dens[L // 4: 3 * L // 4])), cfg['hz']))
    mu = np.array([p[0] for p in pts])
    n = np.array([p[1] for p in pts])
    a, b = np.polyfit(mu, n, 1)
    mu_star = (t['target'] - b) / a if a > 0 else mu[1]
    mu_star = float(np.clip(mu_star, mu.min() - 0.05, mu.max() + 0.05))
    hz = Mo.hz_for_mu_eff(Mo.bond_family(t['bond']), mu_star)
    res = dict(key=key, points=pts, slope=a, intercept=b, target=t['target'], mu_eff=mu_star, hz=hz)
    json.dump(res, open(os.path.join(outdir, 'tune', key + '.json'), 'w'), indent=1)
    return hz


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--node', required=True)
    ap.add_argument('--slots', type=int, default=8)
    ap.add_argument('--threads', type=int, default=6)
    ap.add_argument('--out', default='out')
    ap.add_argument('--dry', action='store_true')
    a = ap.parse_args()
    out = a.out
    for d in (out, os.path.join(out, 'cfg'), os.path.join(out, 'logs'), os.path.join(out, 'tune')):
        os.makedirs(d, exist_ok=True)
    tune, now, later = tasks_for_node(a.node)
    queue = list(now)
    for key, t in tune.items():
        queue += t['runs']
    pending_keys = set(tune)
    if a.dry:
        for c in sorted(queue, key=cost, reverse=True):
            print(c['tag'], c['L'], c['final']['chi'], f"{cost(c):.2e}")
        for key, mk in later.items():
            for c in mk(-20.0):
                print('(after tuning)', c['tag'], f"{cost(c):.2e}")
        return
    running = {}
    env = dict(os.environ)
    for v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
        env[v] = str(a.threads)

    def done_ok(tag):
        f = os.path.join(out, tag + '.json')
        try:
            return json.load(open(f)).get('status') == 'done'
        except Exception:
            return False

    t0 = time.time()
    while queue or running or pending_keys:
        # finished processes
        for tag in list(running):
            p, cfg = running[tag]
            if p.poll() is not None:
                print(f"[{time.time() - t0:8.0f}s] finished {tag} rc={p.returncode}", flush=True)
                del running[tag]
        # tuning keys whose runs are all done
        for key in list(pending_keys):
            tags = [c['tag'] for c in tune[key]['runs']]
            if all(done_ok(tg) for tg in tags):
                hz = finish_tuning(key, tune[key], out)
                print(f"[{time.time() - t0:8.0f}s] tuned {key}: hz = {hz:.9f}", flush=True)
                queue += later[key](hz)
                pending_keys.discard(key)
            elif all(tg not in running and tg not in [c['tag'] for c in queue] for tg in tags) and \
                    not all(done_ok(tg) for tg in tags):
                print(f"[{time.time() - t0:8.0f}s] tuning {key} failed; using the predicted field", flush=True)
                b = Mo.bond_family(tune[key]['bond'])
                mu0 = tune[key]['runs'][1]['mu_eff']
                queue += later[key](Mo.hz_for_mu_eff(b, mu0))
                pending_keys.discard(key)
        # launch
        queue.sort(key=cost, reverse=True)
        tuning_first = [c for c in queue if c['tag'].startswith('tune_')]
        order = tuning_first + [c for c in queue if not c['tag'].startswith('tune_')]
        while len(running) < a.slots and order:
            cfg = order.pop(0)
            queue.remove(cfg)
            if done_ok(cfg['tag']):
                continue
            fcfg = os.path.join(out, 'cfg', cfg['tag'] + '.json')
            json.dump(cfg, open(fcfg, 'w'), indent=1)
            log = open(os.path.join(out, 'logs', cfg['tag'] + '.log'), 'a')
            p = subprocess.Popen([sys.executable, 'run190.py', fcfg, out], stdout=log, stderr=subprocess.STDOUT,
                                 env=env)
            running[cfg['tag']] = (p, cfg)
            print(f"[{time.time() - t0:8.0f}s] started {cfg['tag']} (cost {cost(cfg):.2e})", flush=True)
        time.sleep(20)
    print('pipeline done', flush=True)


if __name__ == '__main__':
    main()
