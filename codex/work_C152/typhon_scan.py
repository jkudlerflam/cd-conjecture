"""Two-stage ring-ED scan of the parity class on random field lines (Typhon, 48 cores).
Stage 1: N random directions (Jx, Jy, D, Jz) on S^3, hz in [0, 3] (31 values), c_eff from L = 12, 16.
Stage 2: L = 12, 16, 20 for the top 400 stage-1 points and every point with c_eff > 1.1.
Also runs the extra task file tasks_extra.json (if present) at L = 12, 16, 20.
Outputs: stage1.jsonl, stage2.jsonl, extra.jsonl in the working directory."""
import json, os, sys
import numpy as np
from multiprocessing import Pool
from ring_ceff import ceff

def w1(t):
    lab, *p = t
    try:
        r = ceff(tuple(p), (12, 16)); r['label'] = lab
        r['c'] = [float(x) for x in r['c']]; r['gapL'] = [float(x) for x in r['gapL']]; r['S'] = [float(x) for x in r['S']]
        return r
    except Exception as e:
        return dict(label=lab, p=p, err=str(e))
def w2(t):
    lab, *p = t
    try:
        r = ceff(tuple(p), (12, 16, 20)); r['label'] = lab
        r['c'] = [float(x) for x in r['c']]; r['gapL'] = [float(x) for x in r['gapL']]; r['S'] = [float(x) for x in r['S']]
        return r
    except Exception as e:
        return dict(label=lab, p=p, err=str(e))

if __name__ == "__main__":
    N = int(sys.argv[1]); nproc = int(sys.argv[2])
    rng = np.random.default_rng(152)
    tasks = []
    for i in range(N):
        v = rng.normal(size=4); v /= np.linalg.norm(v)
        for hz in np.linspace(0, 3, 31):
            tasks.append([f"line{i}"] + [float(x) for x in v] + [float(hz)])
    with Pool(nproc) as pool:
        if os.path.exists("tasks_extra.json"):
            ex = json.load(open("tasks_extra.json"))
            with open("extra.jsonl", "a") as f:
                for r in pool.imap_unordered(w2, ex):
                    f.write(json.dumps(r) + "\n"); f.flush()
        res = []
        with open("stage1.jsonl", "a") as f:
            for r in pool.imap_unordered(w1, tasks, chunksize=4):
                f.write(json.dumps(r) + "\n"); f.flush(); res.append(r)
        ok = [r for r in res if 'c' in r]
        ok.sort(key=lambda r: -r['c'][0])
        sel = ok[:400] + [r for r in ok[400:] if r['c'][0] > 1.1]
        t2 = [[r['label']] + r['p'] for r in sel]
        with open("stage2.jsonl", "a") as f:
            for r in pool.imap_unordered(w2, t2):
                f.write(json.dumps(r) + "\n"); f.flush()
    print("done")
