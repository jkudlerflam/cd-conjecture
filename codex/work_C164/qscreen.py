"""C164 quantum screen for two-valley one-particle bands (Claude Opus agent, 2026-11-08).

For each bond: ring ED at L (momentum resolved), require a unique ground state at k = 0 (uniform background),
take the one-particle band E1(k) (lowest level above E0 in each sector), interpolate it trigonometrically and
record its local minima. A bond is flagged "two-valley" when
  (i)  the band has two local minima k_A, k_B (both lattice points or interpolated), separated by more than
       2 * 2 pi / L, with E_B - E_A < 0.5 * E_A (near-degenerate relative to the gap);
  (ii) both minima are isolated: E1(k) at the minima lies below the two-particle floor min_q E(q) + E(k - q)
       (from the interpolant), so the lowest state there is a one-particle state;
  (iii) the gap E_A is below 0.5 times the bandwidth (low-lying, so an onset is nearby in parameter space).
Ensembles (kind):
  hF    J = -1 + eps (S0 + A(D)), B = eps b, (S0, D, b) uniform on S^10 (N), b = 0 for hFT (T)
  dm    pure-DM wedge on the lattice: J = -1 + eps A(z) + eps^2 S, f = eps^2 b, (S, b) uniform on S^8 scaled by s
  u     uniform on the unit 8-sphere (N), uT with zero field (T)
  I     D = 0 (site-inversion symmetric, E(k) = E(-k) exactly), uniform; IhF the same near h_F
usage: qscreen.py seed count kind L out.jsonl
"""
import sys, json, time
import numpy as np
import est, qband

TWO_PI = 2 * np.pi


def rand_S(rng):
    S = rng.normal(size=(3, 3)); return (S + S.T) / 2


def to_v(J, f):
    Js = (J + J.T) / 2; Ja = (J - J.T) / 2
    w, R = np.linalg.eigh(Js)
    Ja2 = R.T @ Ja @ R; f2 = R.T @ f
    return np.array([w[0], w[1], w[2], Ja2[1, 2], Ja2[2, 0], Ja2[0, 1], *f2])


def A_of(D):
    A = np.zeros((3, 3))
    A[0, 1] += D[2]; A[1, 0] -= D[2]
    A[1, 2] += D[0]; A[2, 1] -= D[0]
    A[2, 0] += D[1]; A[0, 2] -= D[1]
    return A


def sample(rng, kind):
    if kind in ('hF', 'hFT'):
        eps = rng.uniform(0.3, 1.5)
        S = rand_S(rng); S -= np.trace(S) / 3 * np.eye(3)
        D = rng.normal(size=3); b = rng.normal(size=3) if kind == 'hF' else np.zeros(3)
        nrm = np.sqrt(np.sum(S * S) + D @ D + b @ b)
        J = -np.eye(3) + eps * (S + A_of(D)) / nrm; f = eps * b / nrm
        return to_v(J, f), dict(eps=float(eps))
    if kind == 'dm':
        eps = rng.uniform(0.3, 1.0); s = np.exp(rng.uniform(np.log(0.1), np.log(5.0)))
        S = rand_S(rng); b = rng.normal(size=3)
        nrm = np.sqrt(np.sum(S * S) + b @ b)
        J = -np.eye(3) + eps * A_of([0, 0, 1.0]) + eps * eps * s * S / nrm; f = eps * eps * s * b / nrm
        return to_v(J, f), dict(eps=float(eps), s=float(s))
    if kind == 'IhF':                              # inversion-symmetric (D = 0) near h_F
        eps = rng.uniform(0.3, 1.5)
        S = rand_S(rng); S -= np.trace(S) / 3 * np.eye(3); b = rng.normal(size=3)
        nrm = np.sqrt(np.sum(S * S) + b @ b)
        return to_v(-np.eye(3) + eps * S / nrm, eps * b / nrm), dict(eps=float(eps))
    v = rng.normal(size=9)
    if kind == 'I':                                # inversion-symmetric (D = 0), uniform
        v[3:6] = 0
    if kind == 'uT':
        v[6:] = 0
    return v / np.linalg.norm(v), {}


def pdist(a, b):
    d = abs((a - b) % TWO_PI)
    return min(d, TWO_PI - d)


def assess(v, L=14):
    J, f = est.J_f(v)
    h = est.bond_from_JB(J, f)
    b = qband.band(h, L, nev=2)
    lev = b['lev'] - b['E0']
    out = dict(m0=b['m0'], gap0=float(np.sort(lev.ravel())[1]))
    E1 = b['E1']
    out['gap'] = float(E1.min()); out['bw'] = float(E1.max() - E1.min())
    if b['m0'] != 0 or out['gap'] < 1e-6:
        out['flag'] = False
        return out
    ks, Ek = qband.interp(E1, ngrid=720)
    mins = sorted(qband.local_minima(ks, Ek), key=lambda t: t[1])
    out['mins'] = [(round(k, 4), round(E, 6)) for k, E in mins[:4]]
    flag = False
    if len(mins) >= 2 and out['gap'] < 0.5 * max(out['bw'], 1e-9):
        kA, EA = mins[0]
        for kB, EB in mins[1:]:
            if pdist(kA, kB) > 2 * TWO_PI / L and EB - EA < 0.5 * EA:
                E2 = qband.two_particle_floor(ks, Ek)
                iA = int(np.argmin(np.abs(ks - kA))); iB = int(np.argmin(np.abs(ks - kB)))
                iso = bool(Ek[iA] < E2[iA] and Ek[iB] < E2[iB])
                out.update(kA=kA, EA=EA, kB=kB, EB=EB, iso=iso)
                flag = iso
                break
    out['flag'] = flag
    return out


if __name__ == '__main__':
    seed, count, kind, L, fn = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3], int(sys.argv[4]), sys.argv[5]
    rng = np.random.default_rng(seed)
    t0 = time.time()
    with open(fn, 'a') as fo:
        for i in range(count):
            v, meta = sample(rng, kind)
            try:
                o = assess(v, L)
            except Exception as ex:
                o = dict(error=repr(ex))
            o.update(v=[float(x) for x in v], kind=kind, seed=seed, idx=i, **meta)
            fo.write(json.dumps(o) + '\n'); fo.flush()
    print('done', fn, f'{time.time() - t0:.0f}s')
