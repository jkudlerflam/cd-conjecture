#!/usr/bin/env python3
"""Numerical maximization of A/(pi M log2 d) over bonds at small d (mechanism iv).

Parametrization.  h = sum_{mu,nu} c_{mu nu} G_mu (x) G_nu over an orthonormal Hermitian basis
{G_0 = 1/sqrt d, G_1, ...} with c_00 dropped (identity) and the one-site part restricted to
the symmetric combination G_mu (x) 1 + 1 (x) G_mu (the antisymmetric one is the gauge, which
changes neither E_0 on a ring nor M).  Options restrict to real bonds (time-reversal
symmetric), to reflection-symmetric bonds, or to a U(1)-symmetric subspace.

Objective.  E_0(L) on rings L in --Ls (full-space Lanczos, no symmetry), the Casimir fit
over all L, over the last three and the last two; the score is the smallest of the three
ratios, which penalizes bonds whose finite-size energies do not follow the Casimir form.
Nelder-Mead from random starts, then a local refinement; every evaluation is logged.

    python3 optimize_random.py --d 2 --Ls 8 10 12 14 --starts 4 --evals 400 --seed 1 --out data/opt_d2_s1.jsonl
"""
import argparse, json, time, sys
import numpy as np
from scipy import optimize
from core import Ring, hnorm, fit_cas, herm_basis, LOG2, spread


def basis(d):
    B = list(herm_basis(d))
    B = [b / np.sqrt(np.trace(b @ b).real) for b in B]
    return [np.eye(d) / np.sqrt(d)] + B


def param_shape(d, real=True, reflect=False):
    n = d * d  # basis size incl identity
    pairs = []
    for mu in range(n):
        for nu in range(n):
            if mu == 0 and nu == 0:
                continue
            if mu == 0 and nu != 0:
                continue  # one-site part handled through (mu,0)+(0,mu) symmetric combination
            if reflect and nu < mu:
                continue
            pairs.append((mu, nu))
    return pairs


def build(x, d, pairs, G, reflect=False):
    h = np.zeros((d * d, d * d), complex)
    for c, (mu, nu) in zip(x, pairs):
        if nu == 0:
            h += c * (np.kron(G[mu], G[0]) + np.kron(G[0], G[mu]))
        elif reflect and mu != nu:
            h += c * (np.kron(G[mu], G[nu]) + np.kron(G[nu], G[mu]))
        else:
            h += c * np.kron(G[mu], G[nu])
    return h


def score(h, d, Ls, nlev=1, tol=1e-8):
    """Pessimistic Casimir ratio.  Two-point estimates A(L_i, L_{i+1}) from consecutive ring
    sizes must agree for a chain in a CFT regime (they drift slowly, by marginal operators,
    for the Heisenberg chain); the score is their minimum minus half their spread, in units
    of pi M log2 d.  Fit artifacts from level crossings or incommensurate ground states give
    wildly different two-point values and are penalized."""
    E0 = []
    for L in Ls:
        w = Ring(h, d, L).lowest(k=max(nlev, 1), tol=tol)
        E0.append(float(w[0]))
    M = hnorm(h, d)
    if M < 1e-9:
        return -10.0, dict(E0=E0, M=M)
    den = np.pi * M * np.log(d) / LOG2
    two = [fit_cas(Ls[i:i + 2], E0[i:i + 2], quartic=False)[0] / den for i in range(len(Ls) - 1)]
    A = fit_cas(Ls, E0)[0]
    sc = min(two) - 0.5 * (max(two) - min(two))
    return float(sc), dict(E0=E0, M=float(M), A_fit=A, two_point=two, ratio_fit=A / den)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--d', type=int, default=2)
    ap.add_argument('--Ls', nargs='+', type=int, default=[8, 10, 12, 14])
    ap.add_argument('--starts', type=int, default=4)
    ap.add_argument('--evals', type=int, default=400)
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--complex', action='store_true')
    ap.add_argument('--reflect', action='store_true')
    ap.add_argument('--init', default='random', help='random | heis (d=2) | su (SU(d) swap) | file.json')
    ap.add_argument('--noise', type=float, default=0.3)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    d = a.d
    G = basis(d)
    pairs = param_shape(d, real=not a.complex, reflect=a.reflect)
    rng = np.random.default_rng(a.seed)
    log = open(a.out, 'a')
    nparam = len(pairs)

    def coords_of(h):
        x = []
        for (mu, nu) in pairs:
            if nu == 0:
                x.append(np.trace(h @ np.kron(G[mu], G[0])).real / 2)
            elif a.reflect and mu != nu:
                x.append(np.trace(h @ np.kron(G[mu], G[nu])).real)
            else:
                x.append(np.trace(h @ np.kron(G[mu], G[nu])).real)
        return np.array(x)

    for s in range(a.starts):
        init = a.init
        if init == 'mix':
            init = ['random', 'proj', 'proj', 'proj', 'heis' if d == 2 else 'su', 'xxz'][s % 6]
        if init == 'random':
            x0 = rng.normal(size=nparam)
        elif init == 'proj':
            # minus a random projector of rank r (Haar random subspace), plus noise
            r = int(rng.integers(1, d * d))
            Q, _ = np.linalg.qr(rng.normal(size=(d * d, r)) + (1j * rng.normal(size=(d * d, r)) if a.complex else 0))
            h0 = -(Q @ Q.conj().T)
            x0 = coords_of(h0) + a.noise * 0.3 * rng.normal(size=nparam) * np.linalg.norm(coords_of(h0)) / np.sqrt(nparam)
        elif init == 'xxz':
            if d == 2:
                sx = np.array([[0, 1], [1, 0]]) / 2; sy = np.array([[0, -1j], [1j, 0]]) / 2; sz = np.diag([.5, -.5])
                dl = rng.uniform(-1, 1)
                h0 = np.kron(sx, sx) + np.kron(sy, sy) + dl * np.kron(sz, sz)
            else:
                sz = np.diag(np.arange((d - 1) / 2, -(d + 1) / 2, -1)); sp = np.diag(np.sqrt([((d - 1) / 2) * ((d + 1) / 2) - m * (m + 1) for m in np.arange((d - 1) / 2 - 1, -(d + 1) / 2, -1)]), 1)
                sx = (sp + sp.T) / 2; sy = (sp - sp.T) / 2j
                dl = rng.uniform(-1, 1)
                h0 = (np.kron(sx, sx) + np.kron(sy, sy) + dl * np.kron(sz, sz)).real
            x0 = coords_of(h0) + a.noise * rng.normal(size=nparam) * np.linalg.norm(coords_of(h0)) / np.sqrt(nparam)
        else:
            if init == 'heis' and d == 2:
                sx = np.array([[0, 1], [1, 0]]) / 2; sy = np.array([[0, -1j], [1j, 0]]) / 2; sz = np.diag([.5, -.5])
                h0 = sum(np.kron(s_, s_) for s_ in (sx, sy, sz))
            elif init == 'su':
                h0 = np.eye(d * d).reshape(d, d, d, d).transpose(1, 0, 2, 3).reshape(d * d, d * d)
            else:
                h0 = np.array(json.load(open(init))['h_re'])
            x0 = coords_of(h0) + a.noise * rng.normal(size=nparam) * np.linalg.norm(coords_of(h0)) / np.sqrt(nparam)
        x0 = x0 / np.linalg.norm(x0)
        best = [-np.inf, None, None]
        t0 = time.time()
        neval = [0]

        def f(x):
            h = build(x, d, pairs, G, reflect=a.reflect)
            try:
                sc, info = score(h, d, a.Ls)
            except Exception as e:  # Lanczos failure
                sc, info = -10.0, {'err': str(e)}
            neval[0] += 1
            if sc > best[0]:
                best[0], best[1], best[2] = sc, x.copy(), info
                log.write(json.dumps({'start': s, 'seed': a.seed, 'eval': neval[0], 'score': sc, 'info': info,
                                      'x': x.tolist(), 'time': time.time() - t0}) + '\n'); log.flush()
            return -sc

        optimize.minimize(f, x0, method='Nelder-Mead',
                          options=dict(maxfev=a.evals, xatol=1e-4, fatol=1e-5, adaptive=True,
                                       initial_simplex=None))
        # local polish from the best point with a smaller simplex
        if best[1] is not None:
            xb = best[1]
            simplex = np.vstack([xb] + [xb + 0.05 * np.linalg.norm(xb) * e for e in np.eye(nparam)])
            optimize.minimize(f, xb, method='Nelder-Mead',
                              options=dict(maxfev=a.evals // 2, xatol=1e-5, fatol=1e-6, adaptive=True,
                                           initial_simplex=simplex))
        h = build(best[1], d, pairs, G, reflect=a.reflect)
        log.write(json.dumps({'start': s, 'seed': a.seed, 'init': init, 'final': True, 'score': best[0], 'info': best[2],
                              'x': best[1].tolist(), 'h_re': h.real.tolist(), 'h_im': h.imag.tolist(),
                              'pairs': pairs, 'reflect': a.reflect, 'Ls': a.Ls, 'nevals': neval[0],
                              'time': time.time() - t0}) + '\n'); log.flush()
        print('start', s, 'score', best[0], 'evals', neval[0], 'time', round(time.time() - t0, 1), flush=True)


if __name__ == '__main__':
    main()
