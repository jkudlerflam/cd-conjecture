"""C155 adversarial maximization of the ring central-charge estimator over B_2 (Claude Opus agent, 2026-11-07).

Each task runs Nelder-Mead restarts until its wall-time budget is spent, appending one JSON line per
restart to results/task_<t>.jsonl (checkpoint; a restarted task continues from the number of lines).

Objective modes
  min   : c_obj = min(c_S, c_cas)                      (default, both estimators must be fooled together)
  Sfit  : min(c_S, c_fit) at points with a nondegenerate L = 16 ground state (L*gap > 1e-3), else 0
Seeds
  rand  : uniform on S^8
  Z2, TR, RA : random point of the symmetry stratum (contrib/d2_grid/strata.py); 'conf' variants keep the
               whole optimization inside that stratum
  XZ, XZT : the XZ neighbourhood (h = X (x) Z, C144), seeded at distance delta in {0.1, 0.2, 0.3, 0.45} along
          YX, ZY, YX +- ZY or a random direction, in J-space coordinates; XZT stays in C144's tuned hyperplane
          dJ_zx = (2/pi) dJ_yy. No crit penalty there, since the emergent velocities are O(delta^2).
  XZB   : the Kitaev chain in a [110] field (C156): XZ plus a uniform field B (x + z)/sqrt 2, B in [0.1, 0.5],
          with a 0.05 kick in the 11 normal directions; no crit penalty.
  known : special points of B_2 (XXZ family, Heisenberg, XX at saturation, TFI and XY Ising lines, XZ point,
          staggered box boundary, U(1) + DM, multicritical intersections) plus a small random kick
The search runs in a local chart v = normalize(v0 + E y), E an orthonormal basis of the tangent space
(8 directions, or the tangent space of the stratum when confined).
"""
import argparse, json, os, time
import numpy as np
from scipy.optimize import minimize
from est import estimate, bond_v, J_f

COORDS = {'Z2': [0, 1, 2, 5, 8], 'TR': [0, 1, 2, 3, 4, 5], 'RA': [0, 1, 2, 5, 6, 7],
          'NS': list(range(9))}


def to_v(J, f):
    """Principal-frame 9-vector of a general (J, f)."""
    J = np.asarray(J, float); f = np.asarray(f, float)
    S = (J + J.T) / 2
    w, R = np.linalg.eigh(S)
    Jp = R.T @ J @ R; fp = R.T @ f
    D = np.array([Jp[1, 2] - Jp[2, 1], Jp[2, 0] - Jp[0, 2], Jp[0, 1] - Jp[1, 0]]) / 2
    v = np.concatenate([np.diag(Jp), D, fp])
    return v / np.linalg.norm(v)


def V(Jx, Jy, Jz, D=(0, 0, 0), f=(0, 0, 0)):
    v = np.array([Jx, Jy, Jz, *D, *f], float)
    return v / np.linalg.norm(v)


def known_points():
    P = []
    for d in [-0.95, -0.7, -0.5, 0, 0.5, 0.9, 1.0, 1.1]:
        P.append(('XXZ%.2f' % d, V(1, 1, d)))
    P.append(('FM_Heis', V(-1, -1, -1)))
    P.append(('XXsat', V(1, 1, 0, f=(0, 0, 2))))
    P.append(('XXZ-1_field', V(1, 1, -1, f=(0, 0, 0.5))))
    P.append(('Heis_sat', V(1, 1, 1, f=(0, 0, 4))))
    for g in [0.5, 1.0]:
        P.append(('XY_Ising_g%.1f' % g, V(-(1 + g) / 2, -(1 - g) / 2, 0, f=(0, 0, -1))))
    P.append(('TFI', V(-1, 0, 0, f=(0, 0, -1))))
    P.append(('TFI_long', V(-1, 0, 0, f=(0.2, 0, -1))))
    P.append(('XYmulticrit', V(-0.5, -0.5, 0, f=(0, 0, -1))))     # gamma -> 0 end of the Ising line
    # XZ point (compass), sigma^x sigma^z
    J = np.zeros((3, 3)); J[0, 2] = 1
    P.append(('XZ', to_v(J, np.zeros(3))))
    J2 = J.copy(); J2[1, 1] = 0.2
    P.append(('XZ+YY', to_v(J2, np.zeros(3))))
    # staggered box a(XX - YY) + uZZ + bZ near the open-region boundary u = sqrt(1 + b^2/4)
    for b in [0.5, 1.5, 3.0]:
        u = np.sqrt(1 + b * b / 4)
        P.append(('stag_b%.1f' % b, V(1, -1, 0.97 * u, f=(0, 0, b))))
        P.append(('stagneg_b%.1f' % b, V(1, -1, -0.97, f=(0, 0, b))))
    for dz in [0.5, 1.0, 2.0]:
        P.append(('XXZ_DM%.1f' % dz, V(1, 1, 0.5, D=(0, 0, dz))))
    P.append(('KitaevGamma', to_v(np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0.]]) + np.diag([0, 0, 1.]), np.zeros(3))))
    P.append(('XYZ_AT', V(1, 0.5, 0.5)))   # zero-field XYZ critical (two equal couplings)
    P.append(('XYZ_FMAT', V(1, -1, 1)))
    return P


KNOWN = known_points()


def objective_value(e, mode, penal=True):
    """Screening value. With penal, subtract 2 max(0, crit - 0.2): a CFT has gap ~ 1/L, so points whose ring
    gaps do not scale like 1/L are pushed down (this removes the level-crossing and spiral artifacts that the
    unpenalized pilot converged to, probe_pilot_best.out)."""
    if mode == 'min':
        val = e['c_obj']
    else:
        val = min(e['c_S'], e['c_fit']) if (e['Lgap'][-1] > 1e-3 and e['ndeg'][-1] < 3) else 0.0
    if not np.isfinite(val):
        val = 0.0
    val = min(val, 4.0)
    if penal:
        val -= 2 * max(0.0, min(e['crit'], 3.0) - 0.2)
    return float(val)


# ------------------------------------------------------------- XZ neighbourhood (C144, coordinator request)
# J-space coordinates w = (J flattened row-major, f), 12 entries; XZ point J0 = E_xz (h = X (x) Z).
IDX = lambda a, b: 3 * a + b
W0 = np.zeros(12); W0[IDX(0, 2)] = 1.0
E_YX = np.zeros(12); E_YX[IDX(1, 0)] = 1.0
E_ZY = np.zeros(12); E_ZY[IDX(2, 1)] = 1.0
N_T = np.zeros(12); N_T[IDX(2, 0)] = 1.0; N_T[IDX(1, 1)] = -2 / np.pi   # T: dJ_zx - (2/pi) dJ_yy = 0
N_T /= np.linalg.norm(N_T)


def w_to_v(w):
    return to_v(w[:9].reshape(3, 3), w[9:])


def complement_basis(normals):
    M = np.array(normals).T
    Q, _ = np.linalg.qr(np.concatenate([M, np.eye(12)], 1))
    return Q[:, len(normals):12]


B_XZ = complement_basis([W0])            # 11 directions
# Kitaev chain in a [110] field = XZ point with a uniform field along the bisector (x + z)/sqrt 2 (C156)
F_BIS = np.zeros(12); F_BIS[9] = F_BIS[11] = 1 / np.sqrt(2)
B_XZT = complement_basis([W0, N_T])      # 10 directions inside the tuned hyperplane


def tangent_basis(v0, coords, rng):
    idx = np.array(coords)
    M = np.zeros((9, len(idx)))
    M[idx, np.arange(len(idx))] = 1.0
    # project out v0 inside the subspace
    P = M - np.outer(v0, v0 @ M)
    Q, _ = np.linalg.qr(P)
    # drop the direction that is (numerically) parallel to v0
    U, s, _ = np.linalg.svd(P, full_matrices=False)
    E = U[:, s > 1e-8]
    return E


def run_restart(task, irs, rng, args, log_high):
    kinds = ['rand', 'rand', 'Z2', 'TR', 'RA', 'known', 'Z2conf', 'TRconf', 'RAconf', 'known',
             'XZ', 'XZT', 'XZT', 'XZ', 'XZB', 'XZB']
    kind = kinds[(task + irs) % len(kinds)]
    mode = 'Sfit' if ((task * 7 + irs) % 4 == 3) else 'min'
    penal = not kind.startswith('XZ')
    if kind.startswith('XZ'):
        Bs = B_XZT if kind == 'XZT' else B_XZ
        delta = float(rng.choice([0.1, 0.2, 0.3, 0.45]))
        pick = rng.integers(6)
        if kind == 'XZB':
            delta = float(rng.uniform(0.1, 0.5))
            kick = Bs @ rng.normal(size=Bs.shape[1]); kick *= 0.05 / np.linalg.norm(kick)
            d0 = F_BIS + kick / delta; dl = 'bisector-field'
        elif pick == 0:
            d0 = E_YX.copy(); dl = 'YX'
        elif pick == 1:
            d0 = E_ZY.copy(); dl = 'ZY'
        elif pick == 2:
            d0 = (E_YX + rng.choice([-1, 1]) * E_ZY) / np.sqrt(2); dl = 'YX+-ZY'
        else:
            d0 = Bs @ rng.normal(size=Bs.shape[1]); d0 /= np.linalg.norm(d0); dl = 'rand'
        if kind != 'XZB':
            d0 = Bs @ (Bs.T @ d0); d0 /= np.linalg.norm(d0)
        w0 = W0 + delta * d0
        E = Bs
        chart = lambda y: w_to_v(w0 + E @ y)
        step = 0.3 * delta
        label = f'{kind}:{dl}:d{delta}'
        restart_chart = lambda wbest: (lambda y: w_to_v(wbest + E @ y))
        y_to_w = lambda y: w0 + E @ y
    else:
        if kind == 'rand':
            v0 = rng.normal(size=9); coords = COORDS['NS']; label = 'rand'
        elif kind in ('Z2', 'TR', 'RA', 'Z2conf', 'TRconf', 'RAconf'):
            st = kind[:2]
            v0 = np.zeros(9); v0[COORDS[st]] = rng.normal(size=len(COORDS[st]))
            coords = COORDS[st] if kind.endswith('conf') else COORDS['NS']; label = kind
        else:
            name, vk = KNOWN[rng.integers(len(KNOWN))]
            v0 = vk + 0.05 * rng.normal(size=9); coords = COORDS['NS']; label = 'known:' + name
        v0 = v0 / np.linalg.norm(v0)
        E = tangent_basis(v0, coords, rng)
        chart = lambda y: (v0 + E @ y) / np.linalg.norm(v0 + E @ y)
        step = float(rng.choice([0.08, 0.15, 0.3]))
        y_to_w = None
    evals = []
    best = dict(val=-np.inf)
    cur = dict(chart=chart)

    def f(y):
        v = cur['chart'](y)
        try:
            e = estimate(bond_v(v), Ls=(12, 16), nev=3)
        except Exception as ex:  # ARPACK failure etc.
            return 1.0
        val = objective_value(e, mode, penal)
        evals.append(1)
        if val > best['val']:
            best.update(val=val, v=v.tolist(), est=e, y=np.array(y), stage=cur.get('stage', 0))
        raw = objective_value(e, mode, False)
        if raw >= args.high or e['c_S'] >= 1.3 and e['c_cas'] >= 1.0:
            log_high.append(dict(point=v.tolist(), val=val, raw=raw, vel=e['v'],
                                 **{k: e[k] for k in ('c_S', 'c_fit', 'c_cas', 'crit', 'm0', 'Lgap', 'ndeg', 'E0')}))
        return -val

    n = E.shape[1]
    Qr, _ = np.linalg.qr(rng.normal(size=(n, n)))
    sim = np.vstack([np.zeros(n), step * Qr])
    minimize(f, np.zeros(n), method='Nelder-Mead',
             options=dict(initial_simplex=sim, maxfev=args.maxfev, xatol=1e-3, fatol=1e-3))
    # one local restart from the best point with a smaller simplex
    if y_to_w is not None:
        wb = y_to_w(best['y'])
        cur['chart'] = lambda y: w_to_v(wb + E @ y)
        wbest = wb
    else:
        v1 = np.array(best['v'])
        E2 = tangent_basis(v1, coords, rng)
        cur['chart'] = lambda y: (v1 + E2 @ y) / np.linalg.norm(v1 + E2 @ y)
        n = E2.shape[1]
        wbest = None
    cur['stage'] = 1
    sim = np.vstack([np.zeros(n), 0.4 * step * np.eye(n)])
    minimize(f, np.zeros(n), method='Nelder-Mead',
             options=dict(initial_simplex=sim, maxfev=args.maxfev // 2, xatol=1e-3, fatol=1e-3))
    est_b = {k: (v if not isinstance(v, np.ndarray) else v.tolist()) for k, v in best['est'].items()}
    return dict(task=task, restart=irs, seed=label, mode=mode, penal=penal, step=step, nfev=len(evals),
                best_val=best['val'], best_v=best['v'], best_est=est_b,
                best_w=(wbest.tolist() if wbest is not None else None))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--task', type=int, required=True)
    ap.add_argument('--hours', type=float, default=1.0)
    ap.add_argument('--maxfev', type=int, default=300)
    ap.add_argument('--high', type=float, default=1.05)
    ap.add_argument('--out', default='results')
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    fn = os.path.join(args.out, f'task_{args.task}.jsonl')
    done = sum(1 for _ in open(fn)) if os.path.exists(fn) else 0
    t0 = time.time()
    irs = done
    while time.time() - t0 < args.hours * 3600:
        rng = np.random.default_rng([155, args.task, irs])
        highs = []
        ts = time.time()
        rec = run_restart(args.task, irs, rng, args, highs)
        rec['highs'] = sorted(highs, key=lambda d: -d['val'])[:40]
        rec['n_high'] = len(highs)
        rec['seconds'] = time.time() - ts
        with open(fn, 'a') as fh:
            fh.write(json.dumps(rec) + '\n')
        print(f"task {args.task} restart {irs} {rec['seed']} {rec['mode']} best {rec['best_val']:.3f} "
              f"nfev {rec['nfev']} highs {len(highs)} {rec['seconds']:.0f}s", flush=True)
        irs += 1


if __name__ == '__main__':
    main()
