"""Fit the twist-relaxed pseudospin levels of the cusp ring to an N-site XXZ ring.

Pseudospin chain (C173 PROOF 4, C162): H_s = sum_l [ (J/2)(e^{i Phi/N} S+_l S-_{l+1} + h.c.) + J Delta S^z_l S^z_{l+1} ],
J > 0 (the sign of J_perp is a gauge choice on even rings up to a shift of Phi by pi N; we scan Phi over [0, 2 pi)).
Data per ring (from analyze2): D_ps = lowest |M| = 2 (S^z = +-1) level, D_4 = lowest |M| = 4, and the lowest M = 0 level
in each pseudospin momentum class kappa = 1..N-1 (kappa = 0: first excited). All relative to the relaxed ground state.
The ground state of the XXZ ring fixes its kappa = 0. Overall scale J is fitted linearly for each (Delta, Phi).
Usage: imported by fit_runs.py."""
import numpy as np
from itertools import combinations


def xxz_levels(N, Delta, Phi):
    """Return dict: (Sz, q_index) -> sorted energies, for J = 1."""
    dim = 1 << N
    states = np.arange(dim)
    out = {}
    tw = np.exp(1j * Phi / N)
    for nup in range(N + 1):
        basis = [s for s in states if bin(s).count('1') == nup]
        idx = {s: i for i, s in enumerate(basis)}
        Hm = np.zeros((len(basis), len(basis)), complex)
        for i, s in enumerate(basis):
            for l in range(N):
                a, b = (s >> l) & 1, (s >> ((l + 1) % N)) & 1
                Hm[i, i] += Delta * (0.25 if a == b else -0.25)
                if a != b:
                    t = s ^ (1 << l) ^ (1 << ((l + 1) % N))
                    # S+_l S-_{l+1}: up moves from l+1 to l  (a = 0, b = 1) gets e^{-i Phi/N}; reverse gets e^{+i Phi/N}
                    ph = tw if (a == 1 and b == 0) else np.conj(tw)
                    Hm[idx[t], i] += 0.5 * ph
        # translation operator in this sector
        T = np.zeros_like(Hm)
        for i, s in enumerate(basis):
            t = ((s << 1) | (s >> (N - 1))) & (dim - 1)
            T[idx[t], i] = 1.0
        e, v = np.linalg.eigh(Hm)
        # resolve momentum within degenerate blocks
        i = 0
        while i < len(e):
            j = i
            while j + 1 < len(e) and abs(e[j + 1] - e[i]) < 1e-9: j += 1
            V = v[:, i:j + 1]
            tm = V.conj().T @ T @ V
            ev = np.linalg.eigvals(tm)
            for x in ev:
                q = int(np.round(np.angle(x) / (2 * np.pi / N))) % N
                out.setdefault((nup - N / 2, q), []).append(e[i])
            i = j + 1
    for k in out: out[k].sort()
    return out


def features(N, Delta, Phi):
    lv = xxz_levels(N, Delta, Phi)
    Sz0 = 0.0 if N % 2 == 0 else 0.5
    cands = [(E[0], key) for key, E in lv.items() if abs(key[0]) == Sz0]
    E0, (sz0, q0) = min(cands)
    f = {}
    f['D_ps'] = min(E[0] for key, E in lv.items() if abs(key[0]) == Sz0 + 1) - E0
    if N >= 4:
        f['D_4'] = min(E[0] for key, E in lv.items() if abs(key[0]) == Sz0 + 2) - E0
    for kap in range(N):
        key = (sz0, (q0 + kap) % N)
        Es = [x for x in lv.get(key, [])]
        # also the mirror Sz -> -Sz for odd N is equivalent; for even N sz0 = 0
        if kap == 0:
            f['k0'] = Es[1] - E0 if len(Es) > 1 else None
        else:
            f[f'k{kap}'] = Es[0] - E0 if Es else None
    # global GS check: is it in Sz = Sz0?
    f['gs_is_Sz0'] = all(E[0] >= E0 - 1e-9 for E in lv.values())
    return f


def fit(N, data, Deltas=np.linspace(-1.2, 3.0, 85), Phis=np.linspace(0, 2 * np.pi, 49)[:-1], keys=None):
    """data: dict of the same keys as features() (missing allowed). Return list of (rms_rel, Delta, Phi, J)."""
    keys = keys or [k for k in data if data[k] is not None]
    y = np.array([data[k] for k in keys])
    res = []
    for D in Deltas:
        for P in Phis:
            f = features(N, D, P)
            x = np.array([f.get(k) if f.get(k) is not None else np.nan for k in keys])
            if np.any(np.isnan(x)) or not f['gs_is_Sz0']:
                continue
            J = (x @ y) / (x @ x)
            r = np.sqrt(np.mean(((J * x - y) / y) ** 2))
            res.append((r, D, P, J))
    res.sort()
    return res
