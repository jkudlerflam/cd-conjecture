#!/usr/bin/env python3
"""C183 (Claude Opus 5.5 agent, 2026-11-09): stage T, the two-body map of the full qubit chain on the J_x + J_y = 0
stratum, over V/Delta = r and |w|.  Extends C177's twob.py / twob_fit.py (which did |w| = 1 only) and adds C175's
dilute pseudospin criterion evaluated with the fitted interactions.

For each (r, |w|):
 1. tune183.tune: D such that the exact band has the requested |w|/|t2| (phi = pi/2).
 2. Full chain at the far-below field: E_2(K = pi) - E_vac, two lowest even-parity levels, rings L in Ls (C177 twob.full2).
 3. Effective two-body model with the exact fitted band (C177 twob.eff2): fit V2, V3 (model A) and V2, V3, V4 (B)
    to the 2 x len(Ls) levels; then E_2 - 2 eps_min at L = 400 with the fitted V (negative = bound pair).
 4. C175 scat.py (bare band 2|w| cos(k - pi/2) - 2 cos 2k, |t2| = 1) with the model-A V2, V3: the like-valley
    length a_s, the unlike matrix A (Hermitian part), Delta_sigma = (a_s - A_d)/|A_o|, and the lowest two-body level
    over all K minus 2 eps_min (C175 calib.py E0 scan).  c = 2 in the dilute gas iff no bound state and
    -1 < Delta_sigma <= 1 (C175 P1, C180 claim 5).
usage: twob183.py r1,r2 wa1,wa2 L1,L2,.. out.json
"""
import sys, json
import numpy as np
import scipy.linalg as sl
from scipy.optimize import least_squares
import band as B
import twob as T
import tune183 as TU
import scat175 as SC


def E0_K(w, K, V2, V3, V4=0.0, Lr=1500):
    c = SC.rel_coeffs(K, w)
    n = Lr
    pot = np.zeros(n); pot[0] = V2; pot[1] = V3; pot[2] = V4
    ab = np.zeros((3, n), complex); ab[0, :] = c[0] + pot; ab[1, :-1] = c[-1]; ab[2, :-2] = c[-2]
    return sl.eig_banded(ab, lower=True, select='i', select_range=(0, 0), eigvals_only=True)[0]


def run(r, wa, Ls, mu=-12.0):
    Tn = TU.tune(r, wa, np.pi / 2)
    cp = TU.bond_cp(r, Tn['Jxpy'], Tn['D'], TU.hz_far(r, mu))   # C183 (after the coordinator's correction): field set by mu_eff
    h = B.bond(**cp)
    Ev, eps = B.band(h, 12)
    f, res = B.fit(eps, 12)
    t2 = -f['b2'] / 2
    tm = T.tm_from_fit(f); c0 = f['c0']; em = T.band_min(tm, c0)
    data = []
    for L in Ls:
        ev, e2 = T.full2(h, L)
        data.append((L, np.array((e2 - ev)[:2])))
    rec = dict(r=r, wa=wa, mu=mu, t2=t2, D=Tn['D'], fit=f, full=[(L, ((e - 2 * em) / t2).tolist()) for L, e in data])
    for name, keys, p0 in (('A', (2, 3), [-2, -2 / 3]), ('B', (2, 3, 4), [-2, -2 / 3, 0])):
        def model(p):
            V = {2: -2.0 * t2, 3: -2 / 3 * t2}
            for k, v in zip(keys, p):
                V[k] = v * t2
            return V
        def resid(p):
            V = model(p); out = []
            for L, ef in data:
                out += list((T.eff2(L, np.pi, tm, c0, V, nev=2) - ef) / t2)
            return out
        sol = least_squares(resid, p0)
        V = model(sol.x)
        Einf = (T.eff2(400, np.pi, tm, c0, V, nev=1)[0] - 2 * em) / t2
        rec[name] = dict(V={k: V[k] / t2 for k in V}, rms=float(np.sqrt(np.mean(np.square(sol.fun)))), Eb400=float(Einf))
    V2, V3 = rec['A']['V'][2], rec['A']['V'][3]
    w = 1j * wa
    try:
        A, ka, kb, K = SC.A_matrix(w, V2=V2, V3=V3)
        AH = (A + A.conj().T) / 2
        a_s = SC.a_like(w, 0, V2=V2, V3=V3)[0].real
        Ad = AH[0, 0].real; Ao = abs(AH[0, 1])
        rec['scat'] = dict(a_s=float(a_s), A_d=float(Ad), A_o=float(Ao), Delta_sigma=float((a_s - Ad) / Ao))
    except Exception as e:  # noqa
        rec['scat'] = dict(error=str(e))
    emb = SC.valleys(w)[0][1]
    Lr = 3000 if wa < 1 else 1500
    eK = [(float(KK), float(E0_K(w, KK, V2, V3, Lr=Lr) - 2 * emb)) for KK in np.linspace(-np.pi, np.pi, 49)]
    kmin = min(eK, key=lambda x: x[1])
    rec['minK'] = dict(K=kmin[0], E=kmin[1])
    return rec


if __name__ == '__main__':
    rs = [float(x) for x in sys.argv[1].split(',')]
    was = [float(x) for x in sys.argv[2].split(',')]
    Ls = [int(x) for x in sys.argv[3].split(',')]
    mus = [float(x) for x in sys.argv[5].split(',')] if len(sys.argv) > 5 else [-12.0]
    out = []
    for r in rs:
        for wa in was:
          for mu in mus:
            rec = run(r, wa, Ls, mu)
            out.append(rec)
            s = rec.get('scat', {})
            print(f"r={r:5.1f} |w|={wa:3.1f} mu={mu:+.1f}  A: V2={rec['A']['V'][2]:+.4f} V3={rec['A']['V'][3]:+.4f} rms={rec['A']['rms']:.1e} "
                  f"Eb400={rec['A']['Eb400']:+.5f} | B: V4={rec['B']['V'][4]:+.4f} rms={rec['B']['rms']:.1e} Eb400={rec['B']['Eb400']:+.5f} | "
                  f"Dsig={s.get('Delta_sigma', float('nan')):+.3f} a_s={s.get('a_s', float('nan')):+.3f} A_d={s.get('A_d', float('nan')):.3f} "
                  f"|A_o|={s.get('A_o', float('nan')):.3f} | minK K={rec['minK']['K']:+.3f} E={rec['minK']['E']:+.5f}", flush=True)
            json.dump(out, open(sys.argv[4], 'w'), indent=1)
