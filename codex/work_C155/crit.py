"""C155 stage B (gap closings first) and stage A (uniform volume sample) (Claude Opus agent, 2026-11-07).

Stage B, one restart:
  B1  minimize the criticality indicator crit = |z_eff - 1| (est.estimate, rings 12 and 16) by Nelder-Mead
      from a random start (uniform on S^8 or a random point of a symmetry stratum);
  B2  at the minimum: full estimate, symmetry and known-stratum labels, and a ball test (12 random points at
      radius 0.03 and 12 at 0.08 on S^8: fraction with crit < 0.15, and their c_obj);
  B3  slide along the critical set: maximize c_obj - 2 crit from the minimum, record the end point.
Stage A: uniform sample of S^8, records (crit, c_S, c_cas, c_obj, z1, z2, ndeg) per point.
"""
import argparse, json, os, time
import numpy as np
from scipy.optimize import minimize
from est import estimate, bond_v, J_f, PAULI, I2
from search import COORDS, tangent_basis

GAPLESS = 0.15
KEYS = ('c_S', 'c_fit', 'c_cas', 'c_obj', 'z1', 'z2', 'crit', 'ndeg', 'v', 'm0', 'gap1', 'gap2')


def ev(v):
    v = np.asarray(v, float); v = v / np.linalg.norm(v)
    try:
        e = estimate(bond_v(v), Ls=(12, 16), nev=3)
    except Exception:
        return None
    return {k: e[k] for k in KEYS}


# ------------------------------------------------------------------ labels
def _charge_nullity(h, sign):
    """Dimension of {n in R^3 : [h, n.s (x) 1 + sign 1 (x) n.s] = A (x) 1 - sign 1 (x) A for some A}.
    sign = +1: uniform charge sum_j n.s_j; sign = -1: staggered charge sum_j (-1)^j n.s_j."""
    cols = []
    for a in range(3):
        q = np.kron(PAULI[a], I2) + sign * np.kron(I2, PAULI[a])
        cols.append((h @ q - q @ h).ravel())
    for a in range(3):   # telescoping freedom
        cols.append(-(np.kron(PAULI[a], I2) - sign * np.kron(I2, PAULI[a])).ravel())
    M = np.array(cols).T
    M = np.concatenate([M.real, M.imag], 0)
    s = np.linalg.svd(M, compute_uv=False)
    return int(np.sum(s < 1e-7 * max(1, s[0])))


def labels(v):
    v = np.asarray(v, float); v = v / np.linalg.norm(v)
    J, f = J_f(v)
    h = bond_v(v)
    lab = []
    if _charge_nullity(h, +1) > 0:
        lab.append('U1')
    if _charge_nullity(h, -1) > 0:
        lab.append('U1stag')
    # free fermion: common left and right null vector n of J with f parallel to n
    U, s, Vt = np.linalg.svd(J)
    if s[-1] < 1e-6:
        nl, nr = U[:, -1], Vt[-1]
        if abs(abs(nl @ nr) - 1) < 1e-6 and np.linalg.norm(f - (f @ nr) * nr) < 1e-6:
            lab.append('free')
    # discrete symmetries (principal frame), as in contrib/d2_grid/strata.py
    S = (J + J.T) / 2
    D = np.array([J[1, 2] - J[2, 1], J[2, 0] - J[0, 2], J[0, 1] - J[1, 0]]) / 2
    for name, R in [('Rx', np.diag([1, -1, -1])), ('Ry', np.diag([-1, 1, -1])), ('Rz', np.diag([-1, -1, 1]))]:
        if np.allclose(R @ S @ R.T, S, atol=1e-6) and np.allclose(R @ D, D, atol=1e-6) and np.allclose(R @ f, f, atol=1e-6):
            lab.append(name)
    for name, O in [('T', -np.eye(3)), ('Kx', np.diag([-1, 1, 1])), ('Ky', np.diag([1, -1, 1])), ('Kz', np.diag([1, 1, -1]))]:
        if np.allclose(O @ S @ O.T, S, atol=1e-6) and np.allclose(O @ D, -D, atol=1e-6) and np.allclose(O @ f, f, atol=1e-6):
            lab.append(name)
    # near-degenerate principal values (enhanced rotation symmetry when D, f allow)
    return lab


def sphere_kick(v, r, rng, coords=None):
    w = rng.normal(size=9)
    if coords is not None:
        m = np.zeros(9); m[coords] = 1; w *= m
    w -= (w @ v) * v
    w /= np.linalg.norm(w)
    u = np.cos(r) * v + np.sin(r) * w
    return u / np.linalg.norm(u)


def ball_test(v, rng, coords=None):
    out = []
    for r in (0.03, 0.08):
        for _ in range(12):
            u = sphere_kick(v, r, rng, coords)
            e = ev(u)
            if e is None:
                continue
            out.append(dict(r=r, crit=e['crit'], c_obj=e['c_obj'], c_S=e['c_S'], ndeg=e['ndeg']))
    summ = {}
    for r in (0.03, 0.08):
        rows = [o for o in out if o['r'] == r]
        summ[str(r)] = dict(n=len(rows), frac_gapless=float(np.mean([o['crit'] < GAPLESS for o in rows])) if rows else None,
                            max_c_obj=float(max(o['c_obj'] for o in rows)) if rows else None)
    return out, summ


def nm(fun, v0, coords, step, maxfev, rng):
    E = tangent_basis(v0, coords, rng)
    n = E.shape[1]
    Qr, _ = np.linalg.qr(rng.normal(size=(n, n)))
    sim = np.vstack([np.zeros(n), step * Qr])
    best = dict(f=np.inf, v=v0)

    def g(y):
        v = v0 + E @ y; v = v / np.linalg.norm(v)
        val = fun(v)
        if val < best['f']:
            best.update(f=val, v=v)
        return val
    minimize(g, np.zeros(n), method='Nelder-Mead', options=dict(initial_simplex=sim, maxfev=maxfev, xatol=1e-4, fatol=1e-4))
    return best


def restart(task, irs, rng, maxfev):
    kinds = ['rand', 'rand', 'rand', 'Z2', 'TR', 'RA', 'Z2conf', 'TRconf', 'RAconf']
    kind = kinds[(task + irs) % len(kinds)]
    st = 'NS' if kind == 'rand' else kind[:2]
    v0 = np.zeros(9); v0[COORDS[st]] = rng.normal(size=len(COORDS[st])); v0 /= np.linalg.norm(v0)
    coords = COORDS[st] if kind.endswith('conf') else COORDS['NS']

    def fcrit(v):
        e = ev(v)
        if e is None:
            return 10.0
        pen = 1.0 if e['ndeg'][-1] >= 3 else 0.0   # exclude flat-band manifolds
        return min(e['crit'], 10.0) + pen

    b1 = nm(fcrit, v0, coords, 0.25, maxfev, rng)
    vmin = b1['v']
    emin = ev(vmin)
    balls, bsum = ball_test(vmin, rng)

    def fslide(v):
        e = ev(v)
        if e is None:
            return 10.0
        if e['ndeg'][-1] >= 3:
            return 5.0
        return -(min(e['c_obj'], 4.0) - 2 * min(e['crit'], 2.0))

    b3 = nm(fslide, vmin, coords, 0.05, maxfev, rng)
    eslide = ev(b3['v'])
    return dict(task=task, restart=irs, kind=kind, v0=v0.tolist(), vmin=vmin.tolist(), crit_min=b1['f'], est_min=emin,
                labels_min=labels(vmin), ball=bsum, ball_points=balls, vslide=b3['v'].tolist(), est_slide=eslide,
                labels_slide=labels(b3['v']))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--task', type=int, required=True)
    ap.add_argument('--stage', choices=['A', 'B'], required=True)
    ap.add_argument('--hours', type=float, default=1.0)
    ap.add_argument('--npts', type=int, default=400)
    ap.add_argument('--maxfev', type=int, default=150)
    ap.add_argument('--out', default='results_crit')
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    if args.stage == 'A':
        fn = os.path.join(args.out, f'uniform_{args.task}.jsonl')
        done = sum(1 for _ in open(fn)) if os.path.exists(fn) else 0
        rng = np.random.default_rng([1551, args.task])
        pts = rng.normal(size=(args.npts, 9))
        with open(fn, 'a') as fh:
            for i in range(done, args.npts):
                v = pts[i] / np.linalg.norm(pts[i])
                e = ev(v)
                fh.write(json.dumps(dict(i=i, v=v.tolist(), est=e)) + '\n'); fh.flush()
        return
    fn = os.path.join(args.out, f'gap_{args.task}.jsonl')
    done = sum(1 for _ in open(fn)) if os.path.exists(fn) else 0
    t0 = time.time(); irs = done
    while time.time() - t0 < args.hours * 3600:
        rng = np.random.default_rng([1552, args.task, irs])
        ts = time.time()
        rec = restart(args.task, irs, rng, args.maxfev)
        rec['seconds'] = time.time() - ts
        with open(fn, 'a') as fh:
            fh.write(json.dumps(rec) + '\n')
        em, es = rec['est_min'] or {}, rec['est_slide'] or {}
        print(f"task {args.task} restart {irs} {rec['kind']} crit_min {rec['crit_min']:.3f} c_obj@min {em.get('c_obj', np.nan):.3f} "
              f"labels {rec['labels_min']} ball {rec['ball']} slide c_obj {es.get('c_obj', np.nan):.3f} crit {es.get('crit', np.nan):.3f} "
              f"{rec['seconds']:.0f}s", flush=True)
        irs += 1


if __name__ == '__main__':
    main()
