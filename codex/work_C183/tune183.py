#!/usr/bin/env python3
"""C183 (Claude Opus 5.5 agent, 2026-11-09): tune the full qubit chain near the Fibonacci cusp to a target
effective first harmonic w = |w| e^{i phi} (units of the measured |t2|) at V/Delta = r, generalizing C177 tune.py
off the J_x + J_y = 0 stratum.

Bond (band.bond): h = Jx XX + Jy YY + D (XY - YX) + Jz ZZ + (hz/2)(Z1 + 1Z) [+ Jxz (XZ+ZX) + Kxz (XZ-ZX) + (hx/2)(X1+1X)],
Delta = Jx - Jy = 1, V = 4 Jz = r.  The dressed one-particle band eps(k) = E_odd(k) - E_vac (ring ED, L = 12,
parity split, far below the band) is fitted with all harmonics to order 4 (band.fit):
    eps = c0 + a1 sin k + ac1 cos k + b2 cos 2k + ...   ;   |t2| = -b2/2,
and the effective-model first harmonic 2|w||t2| cos(k - phi) = 2|w||t2| (cos phi cos k + sin phi sin k) gives
    a1 = 2 |w| |t2| sin phi   (set by D, leading order a1 = 4 D),
    ac1 = 2 |w| |t2| cos phi  (set by J_x + J_y, leading order ac1 = 2 (J_x + J_y)).
phi and pi - phi are equivalent (k -> pi - k), so phi <= pi/2 is used with the sign of J_x + J_y fixed by the data.

Chemical potential.  e_bot(hz) = min_k eps(k) is linear in hz (C177 checked c0(hz) to 5e-4 |t2| at r = 8).
x = -e_bot/|t2| is the dressed chemical potential above the band bottom in |t2| units (C173's mu_eff = -2.2 at
|w| = 1, phi = pi/2 is x = 0.05).  hz(x) = hz0 + (-x |t2| - e_bot0) / (d e_bot / d hz).

Vacuum density n_vac: <(1 - Z)/2> in the dressed vacuum (k = 0, even sector) at the far-below field.
"""
import sys, json
import numpy as np
import scipy.sparse.linalg as spla
import est as E
import band as B


def bond_cp(r, Jxpy, D, hz, extra=None):
    cp = dict(Jx=(1.0 + Jxpy) / 2, Jy=(Jxpy - 1.0) / 2, D=float(D), Jz=r / 4.0, hz=float(hz))
    if extra:
        cp.update(extra)
    return cp


def hz_far(r, mu_eff=-12.0):
    t2 = 1.0 / (2 * r)
    return (t2 * (mu_eff + 6) - r) / 2   # C177 band.couplings: mu = t2 (mu_eff + 6), hz = (mu - 4 Jz)/2


def vac_density(h, L=12):
    T = E.table(L)
    s = T.sector(0)
    H = T.sector_matrix(h, 0)
    pc = np.array([bin(int(x)).count('1') for x in s['R']])
    idx = np.nonzero(pc % 2 == 0)[0]
    w, v = spla.eigsh(H[idx][:, idx], k=1, which='SA', tol=1e-12)
    p = np.abs(v[:, 0]) ** 2
    return float(np.sum(p * pc[idx]) / L / np.sum(p))


def measure(r, Jxpy, D, hz, L=12):
    cp = bond_cp(r, Jxpy, D, hz)
    Ev, eps = B.band(B.bond(**cp), L)
    f, res = B.fit(eps, L)
    t2 = -f['b2'] / 2
    k = np.linspace(-np.pi, np.pi, 20001)
    ef = (f['c0'] + f['a1'] * np.sin(k) + f['anti_cos1'] * np.cos(k) + f['b2'] * np.cos(2 * k) + f['a3'] * np.sin(3 * k)
          + f['b4'] * np.cos(4 * k) + f['anti_sin2'] * np.sin(2 * k) + f['anti_cos3'] * np.cos(3 * k) + f['anti_sin4'] * np.sin(4 * k))
    i = int(np.argmin(ef))
    # the two local minima (valleys) and their offset
    mins = [j for j in range(len(k)) if ef[j] < ef[j - 1] and ef[j] < ef[(j + 1) % len(k)]]
    vals = sorted((float(ef[j]), float(k[j])) for j in mins)
    return dict(cp=cp, fit=f, res=res, t2=t2, a1=f['a1'], ac1=f['anti_cos1'], e_bot=float(ef[i]), k_bot=float(k[i]),
                valleys=vals, eps=eps.tolist())


def tune(r, wa, phi, L=12, iters=4):
    t2 = 1.0 / (2 * r)
    D = wa * np.sin(phi) * t2 / 2
    Jxpy = wa * np.cos(phi) * t2
    hz = hz_far(r)
    for it in range(iters):
        m = measure(r, Jxpy, D, hz, L)
        t2 = m['t2']
        D *= abs(2 * wa * np.sin(phi) * t2 / m['a1'])  # keep D > 0 (C177/C173 convention)
        tgt = 2 * wa * np.cos(phi) * t2
        if abs(tgt) > 1e-14:
            Jxpy = Jxpy * tgt / m['ac1'] if abs(m['ac1']) > 1e-14 else Jxpy
        else:
            Jxpy = 0.0
    m = measure(r, Jxpy, D, hz, L)
    t2 = m['t2']
    wa_m = np.hypot(m['a1'], m['ac1']) / (2 * t2)
    phi_m = float(np.arctan2(abs(m['a1']), abs(m['ac1'])))
    pts = []
    for mu in (-12.0, -9.0):
        mm = measure(r, Jxpy, D, hz_far(r, mu), L)
        pts.append((hz_far(r, mu), mm['e_bot']))
    slope = (pts[1][1] - pts[0][1]) / (pts[1][0] - pts[0][0])
    nv = vac_density(B.bond(**bond_cp(r, Jxpy, D, hz)), L)
    dv = m['valleys'][1][0] - m['valleys'][0][0] if len(m['valleys']) > 1 else None
    return dict(r=r, wa=wa, phi=phi, Jxpy=float(Jxpy), D=float(D), Jz=r / 4.0, t2=float(t2), wa_meas=float(wa_m),
                phi_meas=phi_m, hz0=pts[0][0], ebot0=pts[0][1], slope=float(slope), n_vac=nv,
                valleys=m['valleys'], valley_offset_t2=(None if dv is None else dv / t2),
                harm_over_t2={kk: v / t2 for kk, v in m['fit'].items()}, fit_res=m['res'])


def hz_of_x(T, x):
    return T['hz0'] + (-x * T['t2'] - T['ebot0']) / T['slope']


def couplings(T, x, extra=None):
    cp = bond_cp(T['r'], T['Jxpy'], T['D'], hz_of_x(T, x), extra)
    return cp


if __name__ == '__main__':
    # usage: tune183.py r1,r2 wa1,wa2 phi1,phi2 out.json
    out = []
    for r in [float(v) for v in sys.argv[1].split(',')]:
        for wa in [float(v) for v in sys.argv[2].split(',')]:
            for phi in [float(v) for v in sys.argv[3].split(',')]:
                T = tune(r, wa, phi)
                out.append(T)
                print(json.dumps({k: (round(v, 6) if isinstance(v, float) else v) for k, v in T.items() if k != 'harm_over_t2'}), flush=True)
    json.dump(out, open(sys.argv[4], 'w'), indent=1)


# ---------------------------------------------------------------------------------------------------------------
# C183 addition (after the coordinator's correction, 2026-11-09): the transverse-DM bond of C186,
# Gamma (XZ - ZX) = Kxz, which breaks the Z2 parity, so the one-particle band is read without a parity split:
# eps(k) = E_1(k) - E_vac with E_vac the lowest state at k = 0 and E_1(k) the lowest state at k != 0 (second lowest at k = 0),
# valid far below the band, where the lowest state of each k != 0 sector is the dressed one-particle state.
def band_np(h, L=12):
    T = E.table(L)
    eps = np.zeros(L); Ev = None
    for m in range(L):
        H = T.sector_matrix(h, m).toarray()
        w = np.linalg.eigvalsh(H)
        if m == 0:
            Ev = w[0]; eps[0] = w[1]
        else:
            eps[m] = w[0]
    return Ev, eps - Ev


def measure_np(r, Jxpy, D, hz, extra, L=12):
    cp = bond_cp(r, Jxpy, D, hz, extra)
    Ev, eps = band_np(B.bond(**cp), L)
    f, res = B.fit(eps, L)
    k = np.linspace(-np.pi, np.pi, 20001)
    ef = (f['c0'] + f['a1'] * np.sin(k) + f['anti_cos1'] * np.cos(k) + f['b2'] * np.cos(2 * k) + f['a3'] * np.sin(3 * k)
          + f['b4'] * np.cos(4 * k) + f['anti_sin2'] * np.sin(2 * k) + f['anti_cos3'] * np.cos(3 * k) + f['anti_sin4'] * np.sin(4 * k))
    mins = [j for j in range(len(k)) if ef[j] < ef[j - 1] and ef[j] < ef[(j + 1) % len(k)]]
    vals = sorted((float(ef[j]), float(k[j])) for j in mins)
    return dict(fit=f, res=res, t2=-f['b2'] / 2, valleys=vals, e_bot=float(ef.min()), eps=eps)


def tune_dm(r, wa, Gamma, L=12, offset=0.0):
    """D from the target |w| (first harmonic magnitude), J_x + J_y by bisection so that the two valley minima differ by
    offset * |t2| (offset = 0: degenerate valleys, the effective phi = pi/2 of the full chain); Kxz = Gamma."""
    extra = dict(Kxz=float(Gamma))
    t2 = 1.0 / (2 * r)
    D = wa * t2 / 2
    hz = hz_far(r)
    def split(Jxpy, D):
        m = measure_np(r, Jxpy, D, hz, extra, L)
        v = m['valleys']
        if len(v) < 2:
            return None, m
        # signed: (energy of the valley at k > pi/2 side) - (other)
        (e1, k1), (e2, k2) = v[0], v[1]
        hi = e1 if abs(k1) > np.pi / 2 else e2
        lo = e2 if abs(k1) > np.pi / 2 else e1
        return (hi - lo) / m['t2'], m
    for it in range(3):
        m = measure_np(r, 0.0, D, hz, extra, L)
        D *= abs(wa * 2 * m['t2'] / np.hypot(m['fit']['a1'], m['fit']['anti_cos1']))
    a, b = -0.05 * t2 * 8, 0.05 * t2 * 8
    sa, _ = split(a, D); sb, _ = split(b, D)
    for it in range(50):
        c = 0.5 * (a + b); sc, m = split(c, D)
        if sc is None or sa is None:
            break
        if (sc - offset) * (sa - offset) <= 0:
            b, sb = c, sc
        else:
            a, sa = c, sc
        if abs(b - a) < 1e-12:
            break
    Jxpy = 0.5 * (a + b)
    s, m = split(Jxpy, D)
    t2m = m['t2']
    pts = []
    for mu in (-12.0, -9.0):
        mm = measure_np(r, Jxpy, D, hz_far(r, mu), extra, L)
        pts.append((hz_far(r, mu), mm['e_bot']))
    slope = (pts[1][1] - pts[0][1]) / (pts[1][0] - pts[0][0])
    # uniform Delta N = 1 element (coordinator note): measured by <X> in the dressed vacuum (full ED, ring L = 10, far
    # below the band, two-level picture <X> ~ g / E_1(0)); if nonzero, cancel it with h_x by secant iteration.
    def mx(hx):
        cp = bond_cp(r, Jxpy, D, hz, dict(extra, hx=hx))
        J = np.array([[cp['Jx'], cp['D'], cp['Kxz']], [-cp['D'], cp['Jy'], 0.0], [-cp['Kxz'], 0.0, cp['Jz']]])
        f = np.array([hx, 0.0, cp['hz']])
        import c166lib as X
        H = X.H_sparse(J, f, 10, True).toarray()
        w, v = np.linalg.eigh(H)
        psi = v[:, 0]
        X0 = X.site_op(X.PA[0], 0, 10)
        return float(np.real(np.vdot(psi, X0 @ psi)))
    m0 = mx(0.0)
    hx = 0.0
    if abs(m0) > 1e-9:
        h1 = 1e-4; m1 = mx(h1)
        for it in range(6):
            h2 = h1 - m1 * (h1 - hx) / (m1 - m0)
            hx, m0, h1, m1 = h1, m1, h2, mx(h2)
            if abs(m1) < 1e-12:
                break
        hx = h1
        extra = dict(extra, hx=float(hx))
    return dict(r=r, wa=wa, Gamma=Gamma, offset=offset, mx_vac_before=float(mx(0.0)), hx_cancel=float(hx),
                mx_vac_after=float(mx(hx)), Jxpy=float(Jxpy), D=float(D), Jz=r / 4.0, t2=float(t2m),
                wa_meas=float(np.hypot(m['fit']['a1'], m['fit']['anti_cos1']) / (2 * t2m)), valley_split_t2=s,
                valleys=m['valleys'], hz0=pts[0][0], ebot0=pts[0][1], slope=float(slope), n_vac=None, extra=extra,
                harm_over_t2={kk: v / t2m for kk, v in m['fit'].items()})
