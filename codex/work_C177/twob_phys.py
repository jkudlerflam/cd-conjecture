#!/usr/bin/env python3
"""C177 (Claude Opus 5.5 agent, 2026-11-09): two-excitation binding of the full chain at the PHYSICAL chemical potential.
Answers the coordinator's C186 conflict.  twob.py/twob_fit.py measured the band and fitted V2, V3, V4 at mu_eff = -12
(particles far above the vacuum), where the energy denominators of every virtual process carry -2 mu = +12 |t2|, which
shifts the effective couplings at relative order mu/V ~ (Delta/V)^2 -- the same order as the effect being measured.
Here everything is measured at mu_eff slightly below the band bottom (default -2.35, the dressed vacuum is the ground state),
directly on full-chain rings (momentum and parity resolved, exact):
  E_vac = lowest even K=0;  E_1(k) = lowest odd at k (one excitation);  E_2(pi) = two lowest even at K = pi (vacuum excluded
  by momentum).  Direct finite-L binding B(L) = E_2(pi) - E_vac - 2 min_k [E_1(k) - E_vac].
Then the N-conserving two-body matching of twob_fit.py is repeated at this mu (band harmonics from E_1(k) at L = 12, V2, V3,
V4 fitted to the two lowest K = pi levels at L = 14..20) and solved at L = 400.
usage: twob_phys.py r[,r...] mu_eff L[,L...] out.json [Dmode]   Dmode = prod (tune.json D) or retune (|w|/|t2| = 1 at this mu)"""
import sys, json
import numpy as np
from scipy.optimize import least_squares
import scipy.sparse.linalg as spla
import est as E, band as B, twob as T

def bond_at(r, mu, Dscale):
    Tn = [d for d in json.load(open('tune.json')) if abs(d['r'] - r) < 1e-9][0]
    cp = dict(Tn['couplings'])
    cp['hz'] = cp['hz'] - (mu - Tn['target_mu']) * Tn['t2'] / Tn['dc0_dhz']
    cp['D'] = cp['D'] * Dscale
    return cp

def _sector(Tb, h, m, par, nev):
    s = Tb.sector(m); H = Tb.sector_matrix(h, m)
    pc = np.array([bin(int(x)).count('1') for x in s['R']])
    idx = np.nonzero(pc % 2 == par)[0]
    Hs = H[idx][:, idx]
    if Hs.shape[0] <= 1500:
        w, V = np.linalg.eigh(Hs.toarray()); w, V = w[:nev], V[:, :nev]
    else:
        w, V = spla.eigsh(Hs, k=nev, which='SA', tol=1e-13); o = np.argsort(w); w, V = w[o], V[:, o]
    N = (np.abs(V) ** 2 * pc[idx][:, None]).sum(0)
    return w, N


def levels(h, L, nev1=40, nev2=8):
    """Number-filtered levels: E_vac (even, K=0, lowest); E_1(k) = lowest odd state at k with <N> < 2 (one excitation;
    three-excitation states, <N> near 3, can lie lower at k far from the valleys); E_2(pi) = two lowest even K = pi
    states with <N> < 3 (two excitations).  Returns also the <N> of the chosen states."""
    Tb = E.table(L)
    w, N = _sector(Tb, h, 0, 0, 1)
    Ev = w[0]
    e1 = np.zeros(L); n1 = np.zeros(L)
    for m in range(L):
        w, N = _sector(Tb, h, m, 1, nev1)
        ok = np.nonzero(N < 2.0)[0]
        e1[m] = w[ok[0]] - Ev if len(ok) else np.nan; n1[m] = N[ok[0]] if len(ok) else np.nan
    w, N = _sector(Tb, h, L // 2, 0, nev2)
    ok = np.nonzero(N < 3.0)[0][:2]
    return Ev, e1, w[ok] - Ev, n1, N[ok]

if __name__ == '__main__':
    rs = [float(x) for x in sys.argv[1].split(',')]; mu = float(sys.argv[2]); Ls = [int(x) for x in sys.argv[3].split(',')]
    Dmode = sys.argv[5] if len(sys.argv) > 5 else 'prod'
    res = []
    for r in rs:
        Ds = 1.0
        for it in range(3 if Dmode == 'retune' else 1):
            cp = bond_at(r, mu, Ds); h = B.bond(**cp)
            Ev, e1, _, _, _ = levels(h, 12)
            f, _ = B.fit(e1, 12); t2 = -f['b2'] / 2; u = abs(f['a1'] / 2) / t2
            if Dmode == 'retune':
                Ds = Ds / u
        cp = bond_at(r, mu, Ds); h = B.bond(**cp)
        Ev, e1, _, n1, _ = levels(h, 12); f, fres = B.fit(e1, 12)
        print(f'   L=12 one-excitation <N> range {np.nanmin(n1):.3f} to {np.nanmax(n1):.3f}', flush=True); t2 = -f['b2'] / 2; u = abs(f['a1'] / 2) / t2
        tm = T.tm_from_fit(f); em = T.band_min(tm, f['c0'])
        t2LO = 1 / (2 * r)
        data = []
        for L in Ls:
            Ev, e1L, e2L, n1L, n2L = levels(h, L, nev1=6)
            BL = e2L - 2 * np.nanmin(e1L)
            bare = T.eff2(L, np.pi, tm, f['c0'], {2: -2 * t2, 3: -2 / 3 * t2}, nev=2)
            Bbare = bare - 2 * np.nanmin(e1L)
            data.append((L, e2L))
            print(f"r={r} mu={mu} L={L}: min E1/|t2|={np.nanmin(e1L) / t2:+.5f} <N2>={np.round(n2L, 3).tolist()} B_full(L)/|t2|={np.round(BL / t2, 5).tolist()}  B_effbare(L)/|t2|={np.round(Bbare / t2, 5).tolist()}", flush=True)
        def resid(p, keys):
            V = {2: -2.0 * t2, 3: -2 / 3 * t2}
            for kk, v in zip(keys, p): V[kk] = v * t2
            out = []
            for L, e2L in data:
                out += list((T.eff2(L, np.pi, tm, f['c0'], V, nev=2) - e2L) / t2)
            return out
        rec = dict(r=r, mu=mu, Dmode=Dmode, Dscale=Ds, couplings=cp, t2=t2, t2_over_LO=t2 / t2LO, u=u, fit=f, fit_resid=fres, fits={})
        for name, keys, p0 in (('A', (2, 3), [-2, -2 / 3]), ('B', (2, 3, 4), [-2, -2 / 3, 0])):
            sol = least_squares(resid, p0, args=(keys,))
            V = {2: -2.0 * t2, 3: -2 / 3 * t2}
            for kk, v in zip(keys, sol.x): V[kk] = v * t2
            Einf = T.eff2(400, np.pi, tm, f['c0'], V, nev=1)[0]
            v2, v3 = V[2] / t2, V[3] / t2
            rec['fits'][name] = dict(V={kk: V[kk] / t2 for kk in V}, rms=float(np.sqrt(np.mean(np.square(sol.fun)))),
                                     E_inf=(Einf - 2 * em) / t2, D=2 * v2 * v3 + 4 * v2 + v3 * u * u + 4 * v3 + 8, v=v2 + 2, g=-v3)
            print(f"   fit {name}: V/|t2|={ {kk: round(V[kk] / t2, 4) for kk in V} } rms={rec['fits'][name]['rms']:.1e} "
                  f"v=V2+2={v2 + 2:+.4f} g={-v3:.4f} D={rec['fits'][name]['D']:+.3f}  (E2-2emin)/|t2| at L=400: {(Einf - 2 * em) / t2:+.5f}", flush=True)
        print(f"   band at this mu: |t2|/t2_LO={t2 / t2LO:.4f} u=|w|/|t2|={u:.4f} a3/|t2|={f['a3'] / t2:+.4f} b4/|t2|={f['b4'] / t2:+.4f}", flush=True)
        res.append(rec)
    json.dump(res, open(sys.argv[4], 'w'), indent=1, default=float)
