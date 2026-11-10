#!/usr/bin/env python3
"""C190g GPU pipeline (Claude Opus 5.5 agent, 2026-11-09).

Usage: python3 pipeline190g.py --out outG [--cpu-out ../c190] [--threads 4] [--dry]

One run per GPU (the GPUs of the allocation, from CUDA_VISIBLE_DEVICES), tasks taken in the priority order of
specs/C190g_gpu_port.md step 6.  Configurations are built by pipeline190.cand_task (same schedules, seeds and
initializations as the CPU pipeline), so the output JSON and checkpoints have the CPU format.  A task whose
output is already 'done' is skipped, and a rerun resumes every unfinished run from its checkpoint.
(d) Controls are run only where the CPU control (node B, outB) is not done when a GPU becomes free.
(e) The other tuned keys are added once their tune file exists in the CPU outA/tune or outB/tune (read only).
"""
import os
import sys
import json
import time
import argparse
import subprocess
import pipeline190 as P

ap = argparse.ArgumentParser()
ap.add_argument('--out', default='outG')
ap.add_argument('--cpu-out', default='/scratch/lustre/jkudlerflam/c18x/c190')
ap.add_argument('--threads', type=int, default=4)
ap.add_argument('--gpus', default=None)
ap.add_argument('--dry', action='store_true')
a = ap.parse_args()

HZ = {'p153c': P.COLLAB_HZ['p153']}


def tuned_hz(key):
    for sub in ('outA', 'outB'):
        f = os.path.join(a.cpu_out, sub, 'tune', key + '.json')
        if os.path.exists(f):
            try:
                return json.load(open(f))['hz']
            except Exception:
                pass
    return None
LATER_KEYS = [('p153_n316', 'p153', None), ('rat_n316', 'rat', None), ('sym40_n18', 'sym40', 'stag'),
              ('sym80_n316', 'sym80', 'stag')]


def n_init(n, L, Pp):
    N = int(round(n * L))
    if N % 2 != Pp:
        N += 1
    return f'prod:{N}'


def base_tasks():
    t = []
    # (a) collaborator's field, L = 256, both sectors
    for Pp in (0, 1):
        t.append(P.cand_task('p153c', 'p153', 256, Pp, 'rand:1', HZ['p153c'], None))
    # (b) L = 128 and 64, rand:1 and the product start (as in the CPU pipeline)
    for L in (128, 64):
        for init in ('rand:1', 'prod'):
            for Pp in (0, 1):
                ini = init if init != 'prod' else n_init(0.22, L, Pp)
                t.append(P.cand_task('p153c', 'p153', L, Pp, ini, HZ['p153c'], None))
    # (c) sym40_n316 at its tuned field (outB/tune/sym40_n316.json, exact value), rand:1 plus prod at L >= 128
    HZ['sym40_n316'] = tuned_hz('sym40_n316')
    assert HZ['sym40_n316'] is not None, 'sym40_n316 tune file missing'
    for L in (256, 128, 64):
        for Pp in (0, 1):
            t.append(P.cand_task('sym40_n316', 'sym40', L, Pp, 'rand:1', HZ['sym40_n316'], 'stag'))
        if L >= 128:
            for Pp in (0, 1):
                t.append(P.cand_task('sym40_n316', 'sym40', L, Pp, n_init(0.1957, L, Pp), HZ['sym40_n316'], 'stag'))
    # (d) controls, same configuration as pipeline190 node B
    for legs in (1, 2):
        for L in (256, 128, 64):
            mu, Pp, Ntot = P.free_control_mu(L, 0.1957, legs)
            sc = P.SCHED[L]
            t.append(dict(tag=f"xx{legs}_L{L}_P{Pp}_rand1", model=f'xx{legs}', mu=mu, L=L, P=Pp, init='rand:1',
                          sched=sc['sched'], final=sc['final'], gwin=48 if L >= 128 else 32, variance=False,
                          maxit=300, N_exact=Ntot, _control=True))
    return t


def tuned_tasks():
    """(e) keys whose CPU tuning has finished: same initializations as pipeline190.tasks_for_node."""
    out = []
    for key, bond, gauge in LATER_KEYS:
        hz = tuned_hz(key)
        if hz is None:
            continue
        for L in (256, 128, 64):
            for Pp in (0, 1):
                inits = ['rand:1']
                if key == 'p153_n316' and L == 128:
                    inits.append(n_init(0.1957, L, Pp))
                for ini in inits:
                    out.append(P.cand_task(key, bond, L, Pp, ini, hz, gauge))
    return out


def cpu_done(tag):
    for sub in ('outA', 'outB'):
        f = os.path.join(a.cpu_out, sub, tag + '.json')
        try:
            if json.load(open(f)).get('status') == 'done':
                return True
        except Exception:
            pass
    return False


def done_ok(tag):
    try:
        return json.load(open(os.path.join(a.out, tag + '.json'))).get('status') == 'done'
    except Exception:
        return False


def main():
    for d in (a.out, os.path.join(a.out, 'cfg'), os.path.join(a.out, 'logs')):
        os.makedirs(d, exist_ok=True)
    gpus = (a.gpus or os.environ.get('CUDA_VISIBLE_DEVICES', '0')).split(',')
    queue = base_tasks()
    seen = set(c['tag'] for c in queue)
    if a.dry:
        for c in queue + tuned_tasks():
            print(c['tag'], c['L'], c['final']['chi'], c.get('hz', c.get('mu')))
        return
    free = list(gpus)
    running = {}
    t0 = time.time()
    print(f"[{0:8.0f}s] GPUs {gpus}", flush=True)
    while queue or running:
        for tag in list(running):
            p, g = running[tag]
            if p.poll() is not None:
                print(f"[{time.time() - t0:8.0f}s] finished {tag} rc={p.returncode} on GPU {g}", flush=True)
                del running[tag]
                free.append(g)
        for c in tuned_tasks():
            if c['tag'] not in seen:
                seen.add(c['tag'])
                queue.append(c)
                print(f"[{time.time() - t0:8.0f}s] added {c['tag']} (hz = {c['hz']})", flush=True)
        while free and queue:
            cfg = queue.pop(0)
            tag = cfg['tag']
            if done_ok(tag):
                continue
            if cfg.pop('_control', False) and cpu_done(tag):
                print(f"[{time.time() - t0:8.0f}s] skip {tag} (CPU control done)", flush=True)
                continue
            g = free.pop(0)
            fcfg = os.path.join(a.out, 'cfg', tag + '.json')
            json.dump(cfg, open(fcfg, 'w'), indent=1)
            env = dict(os.environ)
            env['CUDA_VISIBLE_DEVICES'] = g
            for v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
                env[v] = str(a.threads)
            log = open(os.path.join(a.out, 'logs', tag + '.log'), 'a')
            p = subprocess.Popen([sys.executable, 'run190g.py', fcfg, a.out], stdout=log, stderr=subprocess.STDOUT,
                                 stdin=subprocess.DEVNULL, env=env)
            running[tag] = (p, g)
            print(f"[{time.time() - t0:8.0f}s] started {tag} on GPU {g}", flush=True)
        time.sleep(30)
    print('pipeline done', flush=True)


if __name__ == '__main__':
    main()
