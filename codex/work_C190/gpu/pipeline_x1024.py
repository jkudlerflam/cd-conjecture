#!/usr/bin/env python3
"""C190g: the chi = 1024 stage at L = 256 (coordinator, Claude Fable 5, 2026-11-10).

Usage: python3 pipeline_x1024.py [--src outG] [--out outX1024] [--threads 4] [--dry]

Continues finished L = 256 runs of the registered configurations from their final chi = 768 state
(<src>/ckpt/<tag>.npz, read only) with sweeps at chi = 1024.  Model, field, gauge, Davidson settings and
tolerances are those of the registered final stage (pipeline190.SCHED[256]['final']), with min_sweeps = 3
because the first sweep grows the bonds.  Output goes to a separate folder with the same tags and the format
of run190g.py, so nothing in <src> is written.  One run per GPU.  A rerun resumes an unfinished run from its
own checkpoint in <out>/ckpt and skips finished ones.  Registration: RESULT.md, section R8.  The two-leg control is continued last, so that candidate and control are
read at the same chi.
"""
import os
import sys
import json
import time
import argparse
import subprocess

TAGS = ['p153c_L256_P0_rand1', 'p153c_L256_P1_rand1', 'sym40_n316_L256_P0_rand1', 'sym40_n316_L256_P1_rand1',
        'xx2_L256_P0_rand1']
FINAL = dict(chi=1024, res_tol=1e-9, min_sweeps=3, max_sweeps=6, tolE=2e-9, tolS=5e-4)
SRC_CHI = 768

ap = argparse.ArgumentParser()
ap.add_argument('--src', default='outG')
ap.add_argument('--out', default='outX1024')
ap.add_argument('--threads', type=int, default=4)
ap.add_argument('--gpus', default=None)
ap.add_argument('--dry', action='store_true')
a = ap.parse_args()


def make_cfg(tag):
    """Configuration of the continuation, or None with a reason if the source run is not a finished chi = 768 run."""
    fsrc = os.path.join(a.src, tag + '.json')
    fck = os.path.join(a.src, 'ckpt', tag + '.npz')
    try:
        old = json.load(open(fsrc))
    except Exception as e:
        return None, f'source output unreadable ({e})'
    if old.get('status') != 'done':
        return None, 'source run not done'
    last = old['hist'][-1]
    if last['chi'] != SRC_CHI or old['cfg']['L'] != 256:
        return None, f"source run ends at chi = {last['chi']}, L = {old['cfg']['L']}"
    if not os.path.exists(fck):
        return None, 'source checkpoint missing'
    cfg = dict(old['cfg'])
    cfg['init'] = 'file:' + fck
    cfg['sched'] = []
    cfg['final'] = dict(FINAL)
    cfg['continues'] = dict(src=a.src, sweeps=len(old['hist']), chi=SRC_CHI, E_phys=last['E_phys'],
                            Sbar1=last['Sbar1'], Sbar6=last['Sbar6'], Sbar12=last['Sbar12'], trunc=last['trunc'])
    return cfg, None


def done_ok(tag):
    try:
        return json.load(open(os.path.join(a.out, tag + '.json'))).get('status') == 'done'
    except Exception:
        return False


def main():
    queue = []
    for tag in TAGS:
        cfg, why = make_cfg(tag)
        if cfg is None:
            print(f'skip {tag}: {why}', flush=True)
            continue
        queue.append(cfg)
    if a.dry:
        for c in queue:
            print(c['tag'], c['model'], c.get('gauge'), c.get('hz', c.get('mu')), c['init'], c['final'], c['continues'])
        return
    for d in (a.out, os.path.join(a.out, 'cfg'), os.path.join(a.out, 'logs')):
        os.makedirs(d, exist_ok=True)
    gpus = (a.gpus or os.environ.get('CUDA_VISIBLE_DEVICES', '0')).split(',')
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
        while free and queue:
            cfg = queue.pop(0)
            tag = cfg['tag']
            if done_ok(tag):
                print(f"[{time.time() - t0:8.0f}s] skip {tag} (done)", flush=True)
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
