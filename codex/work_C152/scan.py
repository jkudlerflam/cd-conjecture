"""Parallel ring-ED scan. Usage: python3 scan.py tasks.json out.jsonl nproc L1,L2[,L3]
tasks.json: list of [label, Jx, Jy, D, Jz, hz]."""
import sys, json
import numpy as np
from multiprocessing import Pool
from ring_ceff import ceff

Ls = tuple(int(x) for x in sys.argv[4].split(","))
def work(t):
    lab, *p = t
    try:
        r = ceff(tuple(float(x) for x in p), Ls)
        r['label'] = lab
        r['c'] = [float(x) for x in r['c']]; r['gapL'] = [float(x) for x in r['gapL']]; r['S'] = [float(x) for x in r['S']]
        return r
    except Exception as e:
        return dict(label=lab, p=p, err=str(e))
if __name__ == "__main__":
    tasks = json.load(open(sys.argv[1]))
    with Pool(int(sys.argv[3])) as pool, open(sys.argv[2], "a") as f:
        for r in pool.imap_unordered(work, tasks, chunksize=1):
            f.write(json.dumps(r) + "\n"); f.flush()
