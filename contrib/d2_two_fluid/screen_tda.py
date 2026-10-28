#!/usr/bin/env python3
"""Two-body (N <= 2, Tamm-Dancoff) screen of the exact-vacuum family of twobody.py.

For each bond (t = 1, phi, V, gL, gR) find the threshold mu* at which the lowest excitation
above the exact vacuum reaches zero, then fill the isolated bands rigidly (mu = mu* + dmu)
and record the first dmu at which two separate Fermi seas exist.  For the two seas a
(lower charge q = -dE/dmu) and b, record centres k_a, k_b, charges, densities, binding
below the two-particle continuum, and the distance from 2 pi Z of every low-order vertex
momentum m k_a + n k_b.  gL = gR = 0 is the exact U(1) two-body problem.

usage: screen_tda.py SEED NPTS OUT.json [as]   (random points; 'as' = antisymmetric conversion slice)
       screen_tda.py --point phi V gL |gR| arg(gR)
"""
import sys, json, time
import numpy as np
from scipy.optimize import brentq
from twobody import isolated_levels, bond, to_pauli

nK = 720
Ks = np.linspace(-np.pi, np.pi, nK, endpoint=False)
VERTS = [(1, 0), (0, 1), (1, 1), (1, -1), (2, 0), (0, 2), (2, -1), (1, -2), (2, 1), (1, 2), (2, -2), (2, 2)]


def levels(p, mu):
    return isolated_levels(Ks, 1.0, p['phi'], p['V'], mu, p['gL'], p['gR'])


def gmin(p, mu):
    E1, E2, ed = levels(p, mu)
    return np.nanmin(np.concatenate([E1, ed]))


def threshold(p):
    return brentq(lambda mu: gmin(p, mu), -40, 40, xtol=1e-10)


def intervals(mask):
    if not mask.any():
        return []
    if mask.all():
        return [(0, nK - 1)]
    m = np.asarray(mask, int)
    starts = np.where((m - np.roll(m, 1)) == 1)[0]
    ends = np.where((m - np.roll(m, -1)) == 1)[0]
    res = []
    for s in starts:
        e = ends[ends >= s]
        res.append((s, e[0] if len(e) else ends[0]))
    return res


def wrap(x):
    return (x + np.pi) % (2 * np.pi) - np.pi


def sea_info(iv, E, q, edge):
    s, e = iv
    idx = np.arange(s, (e if e >= s else e + nK) + 1) % nK
    j = idx[np.argmin(E[idx])]
    width = len(idx)
    centre = wrap(Ks[s] + (width - 1) / 2 * 2 * np.pi / nK)
    return dict(k=float(centre), kmin=float(Ks[j]), q=float(q[j]), dens=width / nK, Emin=float(E[j]),
                bind=float(edge[j] - E[j]))


def analyze(p, dmus=np.concatenate([np.linspace(0.005, 0.1, 20), np.linspace(0.12, 2.0, 48)])):
    mus = threshold(p)
    res = dict(p={k: (v if not isinstance(v, complex) else [v.real, v.imag]) for k, v in p.items()}, mu_star=mus)
    # valleys at threshold
    E1, E2, ed = levels(p, mus)
    E1b, E2b, _ = levels(p, mus + 1e-6)
    q1 = -(E1b - E1) / 1e-6
    vals = []
    for i in range(nK):
        if np.isnan(E1[i]):
            continue
        l, r = E1[i - 1], E1[(i + 1) % nK]
        if (np.isnan(l) or E1[i] <= l) and (np.isnan(r) or E1[i] <= r):
            vals.append(dict(k=float(Ks[i]), E=float(E1[i]), q=float(q1[i]), bind=float(ed[i] - E1[i])))
    res['valleys'] = vals
    res['single_min_free'] = float(-2 - mus)
    for dm in dmus:
        mu = mus + dm
        E1, E2, ed = levels(p, mu)
        E1b, E2b, _ = levels(p, mu + 1e-6)
        q1 = -(E1b - E1) / 1e-6; q2 = -(E2b - E2) / 1e-6
        seas = [sea_info(iv, E1, q1, ed) for iv in intervals(np.nan_to_num(E1, nan=1) < 0)]
        seas += [sea_info(iv, E2, q2, ed) for iv in intervals(np.nan_to_num(E2, nan=1) < 0)]
        if len(seas) >= 2:
            seas.sort(key=lambda s: s['q'])
            a, b = seas[0], seas[-1]
            vm = {f"{m},{n}": float(abs(wrap(m * a['kmin'] + n * b['kmin']))) for m, n in VERTS}
            res['two_sea'] = dict(dmu=float(dm), seas=seas, a=a, b=b, vertex_mom=vm,
                                  dmin=min(v for k, v in vm.items() if k not in ('2,-2', '2,2')),
                                  dmin_all=min(vm.values()), Qrel=float(abs(wrap(b['kmin'] - 2 * a['kmin']))))
            break
    J, B = to_pauli(bond(1.0, p['phi'], p['V'], mus, p['gL'], p['gR']))
    res['J'] = J.tolist(); res['B'] = B.tolist()
    return res


def random_point(rng, family='all'):
    if family == 'as':          # antisymmetric conversion, gR = -gL e^{i delta}: A(0) = gL + gR small
        g = float(rng.uniform(0.2, 3)); dl = float(rng.uniform(-0.3, 0.3))
        return dict(phi=float(rng.uniform(0, np.pi)), V=float(rng.uniform(-3, 3)), gL=g,
                    gR=complex(-g * np.exp(1j * dl)))
    return dict(phi=float(rng.uniform(0, np.pi)), V=float(rng.uniform(-4, 4)), gL=float(rng.uniform(0, 3)),
                gR=complex(rng.uniform(0, 3) * np.exp(1j * rng.uniform(0, 2 * np.pi))))


if __name__ == '__main__':
    if sys.argv[1] == '--point':
        phi, V, gL, gRa, gRp = map(float, sys.argv[2:7])
        r = analyze(dict(phi=phi, V=V, gL=gL, gR=gRa * np.exp(1j * gRp)))
        print(json.dumps(r, indent=1))
        sys.exit()
    seed, npts, out = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
    family = sys.argv[4] if len(sys.argv) > 4 else 'all'
    rng = np.random.default_rng(seed)
    results = []
    t0 = time.time()
    for i in range(npts):
        p = random_point(rng, family)
        try:
            results.append(analyze(p))
        except Exception as ex:
            results.append(dict(p={k: (v if not isinstance(v, complex) else [v.real, v.imag]) for k, v in p.items()},
                                error=str(ex)))
        if i % 20 == 0:
            json.dump(results, open(out, 'w'))
    json.dump(results, open(out, 'w'))
    print('done', npts, time.time() - t0)
