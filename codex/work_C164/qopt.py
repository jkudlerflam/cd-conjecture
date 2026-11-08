"""C164 adversarial search for quantum two-valley onsets above a uniform background (Claude Opus agent, 2026-11-08).

Objective (minimized by Nelder-Mead restarts on the unit 8-sphere, or inside a sub-ensemble):
  ring ED at L1 = 12 and L2 = 14; uniform background required (unique ground state at k = 0 at both sizes);
  E(k) = trigonometric interpolant of the L2 one-particle band; A = global minimum (k_A, E_A);
  B = lowest point of E(k) with |k - k_A| > 0.9 (mod 2 pi) (a second valley if it is a local minimum);
  O = (E_B - E_A)/bw + 0.25 E_A/bw + penalties
  penalties: not uniform (+5), B not a local minimum (+0.3), minima not below the two-particle floor (+0.5),
             gap closing between L1 and L2 like 1/L (gapless rather than gapped, +1 if gap(L2)/gap(L1) < 0.8).
A point with O < 0.05 has two near-degenerate isolated valleys and a gap below 0.2 bandwidths: a two-valley
onset is one or two tunings away.
usage: qopt.py task hours kind out.jsonl  (kind as in qscreen.sample, the start is drawn from it; 'free' = S^8)
"""
import sys, json, time
import numpy as np
from scipy.optimize import minimize
import est, qband, qscreen

TWO_PI = 2 * np.pi


def cdist(a, b):
    return float(np.abs((a - b + np.pi) % TWO_PI - np.pi))


def valleys(E1, ngrid=720):
    """Interpolated band; A = global minimum; B = lowest other local minimum; barrier = the lower of the two arc
    maxima between them."""
    ks, Ek = qband.interp(E1, ngrid)
    n = len(Ek)
    loc = [i for i in range(n) if Ek[i] <= Ek[i - 1] and Ek[i] <= Ek[(i + 1) % n]]
    iA = int(np.argmin(Ek))
    others = [i for i in loc if cdist(ks[i], ks[iA]) > 0.6]
    if not others:
        return ks, Ek, iA, None, 0.0
    iB = min(others, key=lambda i: Ek[i])
    i0, i1 = sorted((iA, iB))
    arc1 = Ek[i0:i1 + 1].max(); arc2 = max(Ek[i1:].max(), Ek[:i0 + 1].max())
    return ks, Ek, iA, iB, float(min(arc1, arc2))


def iso_discrete(E1):
    """One-particle level below the discrete two-particle floor min_m1 E1(m1) + E1(m - m1) in every sector?"""
    L = len(E1)
    E2 = np.array([min(E1[m1] + E1[(m - m1) % L] for m1 in range(L)) for m in range(L)])
    return E1 < E2 - 1e-9


def tail(E1):
    """Weight of the two highest harmonics |r| >= L/2 - 1 of the sampled band, relative to the bandwidth.
    An isolated band whose particle is small compared with L has exponentially small high harmonics; a jagged
    (period-2 or level-crossing) 'band' does not."""
    L = len(E1); c = np.fft.fft(E1) / L; r = np.abs(np.fft.fftfreq(L, 1.0 / L))
    bw = max(E1.max() - E1.min(), 1e-9)
    return float(np.sum(np.abs(c[r >= L // 2 - 1])) / bw)


def bad_momentum(k, margin=0.35):
    return min(cdist(k, 0.0), cdist(k, np.pi)) < margin


def evaluate(v, L1=12, L2=14, detail=False):
    v = np.asarray(v, float); v = v / np.linalg.norm(v)
    J, f = est.J_f(v); h = est.bond_from_JB(J, f)
    b1 = qband.band(h, L1, nev=2); b2 = qband.band(h, L2, nev=2)
    info = dict(m0=(b1['m0'], b2['m0']))
    if b1['m0'] != 0 or b2['m0'] != 0:
        return 5.0, info
    E1 = b2['E1']; bw = max(E1.max() - E1.min(), 1e-6)
    ks, Ek, iA, iB, bar = valleys(E1)
    kA, EA = float(ks[iA]), float(Ek[iA])
    info.update(kA=kA, EA=EA, bw=float(bw), gap12=float(b1['E1'].min()), gap14=float(E1.min()))
    if iB is None:
        # no second valley: distance to having one, measured by the flattest far region
        far = np.array([cdist(k, kA) > 0.9 for k in ks])
        return 2.0 + float((Ek[far].min() - EA) / bw), info
    kB, EB = float(ks[iB]), float(Ek[iB])
    O = abs(EB - EA) / bw + 0.25 * max(min(EA, EB), 0) / bw
    # barrier between the valleys, at least 5 percent of the bandwidth
    brel = (bar - max(EA, EB)) / bw
    O += 4 * max(0.0, 0.05 - brel)
    # finite-size consistency: the L1 band has valleys at the same places with the same energies
    ks1, Ek1, iA1, iB1, _ = valleys(b1['E1'])
    if iB1 is None:
        O += 0.5
    else:
        pairs = sorted([(cdist(ks1[iA1], kA) + cdist(ks1[iB1], kB), abs(Ek1[iA1] - EA) + abs(Ek1[iB1] - EB)),
                        (cdist(ks1[iA1], kB) + cdist(ks1[iB1], kA), abs(Ek1[iA1] - EB) + abs(Ek1[iB1] - EA))])
        dk, dE = pairs[0]
        O += 0.5 * max(0.0, dk - 0.3) + 2 * dE / bw
    # incommensurate valleys, no pair creation from the vacuum (2kA, 2kB, kA + kB not 0 mod 2 pi)
    for k in (kA, kB):
        if bad_momentum(k):
            O += 0.5
    if cdist(kA + kB, 0.0) < 0.35:
        O += 0.5
    # isolated one-particle levels at the lattice points nearest the two valleys
    iso = iso_discrete(E1)
    mA = int(round(kA / (TWO_PI / L2))) % L2; mB = int(round(kB / (TWO_PI / L2))) % L2
    if not (iso[mA] and iso[mB]):
        O += 0.5
    g1, g2 = b1['E1'].min(), E1.min()
    if g2 / max(g1, 1e-12) < 0.8:
        O += 1.0
    t1, t2 = tail(b1['E1']), tail(E1)
    O += 5 * (max(0.0, t1 - 0.02) + max(0.0, t2 - 0.02))
    info.update(tail=(t1, t2))
    info.update(kB=kB, EB=EB, barrier=float(brel), iso=bool(iso[mA] and iso[mB]))
    if detail:
        info['E1'] = E1.tolist(); info['E1_L1'] = b1['E1'].tolist()
    return float(O), info


def run(task, hours, kind, fn, maxfev=250):
    rng = np.random.default_rng(1000 + task)
    t0 = time.time()
    with open(fn, 'a') as fo:
        r = 0
        while time.time() - t0 < hours * 3600:
            # start: draw from the ensemble until uniform
            for _ in range(200):
                v0, meta = qscreen.sample(rng, kind if kind != 'free' else 'u')
                O0, inf0 = evaluate(v0)
                if O0 < 5:
                    break
            best = [O0, v0, inf0]
            def obj(x):
                O, inf = evaluate(x)
                if O < best[0]:
                    best[:] = [O, np.asarray(x) / np.linalg.norm(x), inf]
                return O
            v0 = np.asarray(v0) / np.linalg.norm(v0)
            minimize(obj, v0, method='Nelder-Mead',
                     options=dict(maxfev=maxfev, xatol=1e-5, fatol=1e-6, initial_simplex=v0 + 0.08 * np.vstack([np.zeros(9), np.eye(9)])))
            rec = dict(task=task, restart=r, kind=kind, O0=O0, O=best[0], v=[float(x) for x in best[1]], info=best[2], meta=meta,
                       t=time.time() - t0)
            fo.write(json.dumps(rec) + '\n'); fo.flush()
            r += 1


if __name__ == '__main__':
    task, hours, kind, fn = int(sys.argv[1]), float(sys.argv[2]), sys.argv[3], sys.argv[4]
    maxfev = int(sys.argv[5]) if len(sys.argv) > 5 else 250
    run(task, hours, kind, fn, maxfev)
