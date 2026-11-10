#!/usr/bin/env python3
"""C195 production pipeline (Claude Opus 5.5 agent, 2026-11-10).

Usage: python3 pipeline195.py [--out out512] [--threads 4] [--dry]

p153c (h_z = -19.97589125, the collaborator's field, model 'p153', no gauge) at L = 512, both parity sectors,
initialization rand:1 (the registered one), one run per GPU of the allocation (run195.py).  The schedule is that
of pipeline190.SCHED[256] for the stages up to chi = 512 (same chi values, mixer amplitudes, residual targets,
two sweeps each), then chi = 768 and 1024 as stages that run to the stop rule, then the final stage.  It is
written in RESULT.md ("Schedule") before production.  A rerun resumes unfinished runs from their checkpoints and
skips finished ones.
"""
import os
import sys
import json
import time
import argparse
import subprocess

L = 512
HZ = -19.97589125
SCHED512 = dict(
    sched=[[32, 1e-3, 1e-7, 2], [64, 1e-4, 1e-8, 2], [128, 1e-4, 1e-8, 2], [192, 1e-5, 1e-8, 2],
           [256, 1e-5, 1e-8, 2], [384, 0.0, 1e-8, 2], [512, 0.0, 1e-9, 2]],
    conv=[dict(chi=768, res_tol=1e-9, min_sweeps=3, max_sweeps=6, tolE=4e-9, tolS=5e-4),
          dict(chi=1024, res_tol=1e-9, min_sweeps=3, max_sweeps=6, tolE=4e-9, tolS=5e-4)],
    final=dict(chi=1536, res_tol=1e-9, min_sweeps=3, max_sweeps=8, tolE=4e-9, tolS=5e-4))


def cfg_for(P):
    return dict(tag=f'p153c_L{L}_P{P}_rand1', model='p153', gauge=None, hz=HZ, L=L, P=P, init='rand:1',
                sched=SCHED512['sched'], conv=SCHED512['conv'], final=SCHED512['final'], gwin=48, variance=False,
                maxit=300)


ap = argparse.ArgumentParser()
ap.add_argument('--out', default='out512')
ap.add_argument('--threads', type=int, default=4)
ap.add_argument('--gpus', default=None)
ap.add_argument('--dry', action='store_true')
a = ap.parse_args()


def done_ok(tag):
    try:
        return json.load(open(os.path.join(a.out, tag + '.json'))).get('status') == 'done'
    except Exception:
        return False


def main():
    queue = [cfg_for(0), cfg_for(1)]
    if a.dry:
        for c in queue:
            print(json.dumps(c))
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
            p = subprocess.Popen([sys.executable, 'run195.py', fcfg, a.out], stdout=log, stderr=subprocess.STDOUT,
                                 stdin=subprocess.DEVNULL, env=env)
            running[tag] = (p, g)
            print(f"[{time.time() - t0:8.0f}s] started {tag} on GPU {g}", flush=True)
        time.sleep(30)
    print('pipeline done', flush=True)


if __name__ == '__main__':
    main()
