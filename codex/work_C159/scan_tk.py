"""C159 volume scan of the strata T and K with C155's ring estimator (Claude Opus agent, 2026-11-07).
modes: T   uniform on the unit sphere of T coordinates (Jx, Jy, Jz, Dx, Dy, Dz)
       K   uniform on the unit sphere of K coordinates (Jx, Jy, Jz, Dz, fx, fy)
       TFM ferromagnetic Heisenberg (-1,-1,-1) plus a T perturbation of radius r in [0.1, 0.9]
       KFM ferromagnetic Heisenberg plus a K perturbation (radius r in [0.1, 0.9])
       TAF antiferromagnetic Heisenberg (1,1,1) plus a T perturbation of radius r in [0.2, 1.2]
Each line: v, spiral-parent data (C154 section 2), C155 estimate at rings (12, 16), lowest levels per momentum
at L = 16 (relative to E0, times L/2pi).
usage: scan_tk.py MODE TASK NPTS [OUTDIR]"""
import sys, os, json, time
import numpy as np
import est as E
import tk

mode, task, npts = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
outdir = sys.argv[4] if len(sys.argv) > 4 else 'scan'
os.makedirs(outdir, exist_ok=True)
rng = np.random.default_rng(1000 * task + {'T': 1, 'K': 2, 'TFM': 3, 'KFM': 4, 'TAF': 5}[mode])
_rd = E.ring_data
_cache = {}


def _ring_data(h, L, r=2, nev=2):
    d = _rd(h, L, r, nev)
    _cache[L] = d
    return d


E.ring_data = _ring_data   # keep the L = 16 levels for the record (estimate() looks ring_data up at call time)
idxT = [0, 1, 2, 3, 4, 5]; idxK = [0, 1, 2, 5, 6, 7]


def sample():
    if mode == 'T':
        return tk.rand_T(rng)
    if mode == 'K':
        return tk.rand_K(rng)
    base = np.zeros(9)
    if mode in ('TFM', 'KFM'):
        base[:3] = -1.0
        lo, hi = 0.1, 0.9
    else:
        base[:3] = 1.0
        lo, hi = 0.2, 1.2
    idx = idxK if mode == 'KFM' else idxT
    d = np.zeros(9); d[idx] = rng.normal(size=6); d /= np.linalg.norm(d)
    return tk.norm(base / np.sqrt(3) + rng.uniform(lo, hi) * d)


out = open(os.path.join(outdir, f'{mode}_{task}.jsonl'), 'a')
t0 = time.time()
for i in range(npts):
    v = sample()
    J, f = E.J_f(v)
    par = tk.spiral_parent(J, f)
    h = E.bond_from_JB(J, f)
    try:
        e = E.estimate(h, Ls=(12, 16), nev=2)
    except Exception as ex:
        e = dict(error=str(ex))
    rec = dict(mode=mode, i=i, vec=v.tolist(), Delta=float(par['Delta']), theta=float(par['theta']),
               harm=float(par['harm']), hm=float(par['hm']), m=par['m'].tolist(),
               nsym=len(tk.symmetry_check(J, f)), **e)
    if 16 in _cache:
        d16 = _cache[16]
        rec['lev16'] = np.round(16 * (d16['lev'] - d16['E0']) / (2 * np.pi), 5).tolist()
    out.write(json.dumps(rec, default=float) + '\n'); out.flush()
print('done', mode, task, npts, f'{time.time() - t0:.0f}s', flush=True)
