"""C189 Part 1: Luttinger parameters K_c and K_s from C176's converged open chains (PROOF sections 2 to 5).

Runs are read from the last line of C176's jsonl files (final sweep, with 'meas'), or generated exactly for the
control ff2v (free fermions with the candidate's band, phi = pi/2, open chain, same L and N), in the same format.
"""
import numpy as np, json, os, glob
from scipy.optimize import least_squares

K0 = np.arccos(-0.25)          # d-frame valley position at |w| = 1, phi = pi/2


# ------------------------------------------------------------------ data
def last_meas(path):
    recs = [json.loads(l) for l in open(path) if l.strip()]
    r = recs[-1]
    m = None
    for x in reversed(recs):
        if 'meas' in x and 'Gf_re' in x['meas']:
            m = x['meas']; break
    return r, m


def load_run(path):
    r, m = last_meas(path)
    L, N = r['cfg']['L'], r['cfg']['N']
    G = np.array(m['Gf_re']) + 1j * np.array(m['Gf_im'])
    Gb = np.array(m['Gb_re']) + 1j * np.array(m['Gb_im'])
    D2 = {int(k.split('_')[1]): np.array(v)[:, 0] + 1j * np.array(v)[:, 1] for k, v in m['Dpair'].items() if k.startswith('2_')}
    D3 = {int(k.split('_')[1]): np.array(v)[:, 0] + 1j * np.array(v)[:, 1] for k, v in m['Dpair'].items() if k.startswith('3_')}
    return dict(tag=r['tag'], L=L, N=N, E=r['E'], trunc=r['trunc'], dens=np.array(m['dens']), G=G, Gb=Gb,
                gwin=m['gwin'], Cnn={int(k): np.array(v) for k, v in m['Cnn'].items()}, D2=D2, D3=D3,
                Wa=m['Wa_f'], Wb=m['Wb_f'], cfg=r['cfg'])


def ff_hamiltonian(L, phi=np.pi / 2, u=1.0, t2=-1.0):
    """single-particle matrix of the c-frame model without blockade or interactions (open chain):
    H = sum_j (w c+_{j+1} c_j + h.c.) + (t2 c+_{j+2} c_j + h.c.), w = u e^{i phi}"""
    w = u * np.exp(1j * phi)
    h = np.zeros((L, L), dtype=complex)
    for j in range(L - 1):
        h[j + 1, j] += w; h[j, j + 1] += np.conj(w)
    for j in range(L - 2):
        h[j + 2, j] += t2; h[j, j + 2] += np.conj(t2)
    return h


def ff_run(L, N, phi=np.pi / 2):
    """exact free-fermion open chain in the C176 measurement format (c frame)."""
    h = ff_hamiltonian(L, phi)
    e, U = np.linalg.eigh(h)
    occ = U[:, :N]
    G = np.conj(occ) @ occ.T          # G[i, j] = <c+_i c_j>
    E = float(e[:N].sum())
    dens = np.real(np.diag(G))
    refs = [L // 2 - 3, L // 2 - 2, L // 2 - 1, L // 2]
    Cnn = {i0: np.array([-abs(G[i0, j]) ** 2 for j in range(i0 + 1, L)]) for i0 in refs}
    D2 = {}
    d = 2
    for i0 in refs:
        row = []
        for j in range(i0 + d + 1, L - d):
            # <c+_{i0+d} c+_{i0} c_j c_{j+d}> = G[i0+d, j+d] G[i0, j] - G[i0+d, j] G[i0, j+d]
            row.append(G[i0 + d, j + d] * G[i0, j] - G[i0 + d, j] * G[i0, j + d])
        D2[i0] = np.array(row)
    a0, a1 = L // 4, 3 * L // 4
    Gw = G[a0:a1, a0:a1]
    ks = np.linspace(-np.pi, np.pi, 256, endpoint=False)
    x = np.arange(a1 - a0)
    ph = np.exp(1j * np.outer(ks, x))
    nk = np.real(np.einsum('ki,ij,kj->k', ph, Gw, ph.conj())) / (a1 - a0)
    Wa = float(nk[np.cos(ks) > 0].sum()); Wb = float(nk[np.cos(ks) < 0].sum())
    return dict(tag=f'ff2v_L{L}_N{N}', L=L, N=N, E=E, trunc=0.0, dens=dens, G=Gw, Gb=None, gwin=[a0, a1], Cnn=Cnn,
                D2=D2, Wa=Wa, Wb=Wb, cfg=dict(L=L, N=N, phi=phi))


def n_mid(run):
    L = run['L']
    return float(np.mean(run['dens'][L // 4: 3 * L // 4]))


# ------------------------------------------------------------------ E_c: structure factor
def Sq(run, qs, Rmax=None):
    L = run['L']
    rows = run['Cnn']
    if Rmax is None:
        Rmax = L // 2 - 4
    Cr = np.mean([np.asarray(v)[:Rmax] for v in rows.values()], axis=0)
    n0 = np.mean([run['dens'][i] for i in rows])
    C0 = n0 * (1 - n0)
    r = np.arange(1, Rmax + 1)
    return C0 + 2 * np.cos(np.outer(qs, r)) @ Cr


def K_from_Sq(run, qlo, qhi, ncomp=2, nq=40):
    """slope of S(q) on [qlo, qhi] through the origin-free line fit; S = K q / pi (ncomp = 2), K q / (2 pi) (1)."""
    qs = np.linspace(qlo, qhi, nq)
    S = Sq(run, qs)
    A = np.vstack([qs, np.ones_like(qs)]).T
    sl, b = np.linalg.lstsq(A, S, rcond=None)[0]
    return float(sl * (np.pi if ncomp == 2 else 2 * np.pi)), float(b)


# ------------------------------------------------------------------ correlation exponents
def envelope(y, w=8):
    y = np.abs(np.asarray(y))
    return np.array([y[max(0, i - w // 2): i + w // 2 + 1].max() for i in range(len(y))])


def chord(r, L):
    return (2 * L / np.pi) * np.sin(np.pi * r / (2 * L))


def slope_loglog(r, y, lo, hi, L=None, use_chord=False):
    sel = (r >= lo) & (r <= hi) & (y > 0)
    x = chord(r[sel], L) if use_chord else r[sel]
    return float(np.polyfit(np.log(x), np.log(y[sel]), 1)[0])


def G_profile(run):
    G = run['G']
    W = G.shape[0]
    i = W // 2 - 4
    Gr = np.array([np.mean([abs(G[i + d, i + d + r]) for d in range(4)]) for r in range(1, W - i - 4)])
    return np.arange(1, 1 + len(Gr)), Gr


def D_profile(run):
    Ds = [np.abs(v) for v in run['D2'].values()]
    nmin = min(len(d) for d in Ds)
    Dr = np.mean([d[:nmin] for d in Ds], axis=0)
    return np.arange(3, 3 + nmin), Dr


def eta(run, which, lo, hi, use_chord=True, w=8):
    if which == 'G':
        r, y = G_profile(run)
    else:
        r, y = D_profile(run)
    e = envelope(y, w)
    return -slope_loglog(r, e, lo, hi, run['L'], use_chord)


# ------------------------------------------------------------------ Friedel oscillations at the walls
def friedel_fit(run, xlo, xhi, ks, n0=None, fit_k=False):
    """dens(x) - n0 = sum_k A_k cos(k x' + phi_k) / d(x')^eta, x' = x + 1 (distance to the virtual site 0),
    d(x') = (2(L+1)/pi) sin(pi x'/(L+1)); both walls averaged by reflection. Returns eta and the fit."""
    L = run['L']
    dn = np.asarray(run['dens'])
    prof = 0.5 * (dn + dn[::-1])
    x = np.arange(1, L + 1, dtype=float)
    sel = (x >= xlo) & (x <= xhi)
    xs = x[sel]; ys = prof[sel]
    Lp = L + 1
    dist = (2 * Lp / np.pi) * np.sin(np.pi * xs / Lp)
    nk = len(ks)
    if n0 is None:
        n0 = n_mid(run)

    def model(p):
        et = p[0]; nb = p[1]
        out = np.full_like(xs, nb)
        for i in range(nk):
            A, B = p[2 + 2 * i], p[3 + 2 * i]
            k = p[2 + 2 * nk + i] if fit_k else ks[i]
            out += (A * np.cos(k * xs) + B * np.sin(k * xs)) / dist ** et
        return out
    p0 = [1.0, n0] + [0.05, 0.05] * nk + (list(ks) if fit_k else [])
    res = least_squares(lambda p: model(p) - ys, p0, max_nfev=20000)
    p = res.x
    rms = float(np.sqrt(np.mean(res.fun ** 2)))
    amp = [float(np.hypot(p[2 + 2 * i], p[3 + 2 * i])) for i in range(nk)]
    kout = list(p[2 + 2 * nk:]) if fit_k else list(ks)
    return dict(eta=float(p[0]), nb=float(p[1]), amp=amp, k=[float(k) for k in kout], rms=rms)


def ks_two_valley(nb):
    kF = np.pi * nb / 2
    return [2 * (K0 + kF), 2 * (K0 - kF)]


# ------------------------------------------------------------------ finite-size scaling of two-point functions
def G_sym(run):
    """|<c+_i c_j>| for pairs placed symmetrically about the centre, i = L/2 - m, j = L/2 + m, r = 2m."""
    G = run['G']; a0 = run['gwin'][0]; L = run['L']
    c = L // 2 - a0
    ms = np.arange(1, min(c, G.shape[0] - c - 1) + 1)
    return 2 * ms, np.array([abs(G[c - m, c + m]) for m in ms])


def D_row(run):
    """|D(i0, i0 + r)| averaged over the four reference sites, r = 3, 4, ..."""
    return D_profile(run)


def amp_at(r, y, r0, half=6, mode='max'):
    sel = np.abs(r - r0) <= half
    if mode == 'max':
        return float(np.max(np.abs(y[sel])))
    return float(np.sqrt(np.mean(np.abs(y[sel]) ** 2)))


def eta_fss(runs, which, x, mode='max', half=6):
    """effective exponents from pairs of sizes at fixed r/L = x: eta = -dln A / dln L."""
    Ls = sorted(runs)
    A = {}
    for L in Ls:
        run = runs[L]
        r, y = (G_sym(run) if which == 'G' else D_row(run))
        r0 = x * L
        if r0 + half > r.max() or r0 - half < r.min():
            continue
        A[L] = amp_at(r, y, r0, half, mode)
    Ls = sorted(A)
    out = {}
    for i in range(len(Ls) - 1):
        L1, L2 = Ls[i], Ls[i + 1]
        out[(L1, L2)] = -np.log(A[L2] / A[L1]) / np.log(L2 / L1)
    return out, A


def friedel_fss(runs, x, ks_fn, half=6):
    """Friedel amplitude at fixed x/L from the wall (both walls), after removing the smooth part with a
    13-site running mean; eta_F = -dln A / dln L."""
    A = {}
    for L, run in runs.items():
        dn = np.asarray(run['dens']); prof = 0.5 * (dn + dn[::-1])
        sm = np.convolve(prof, np.ones(13) / 13, mode='same')
        osc = prof - sm
        x0 = int(round(x * L))
        sel = slice(max(7, x0 - half), x0 + half + 1)
        A[L] = float(np.max(np.abs(osc[sel])))
    Ls = sorted(A)
    return {(Ls[i], Ls[i + 1]): -np.log(A[Ls[i + 1]] / A[Ls[i]]) / np.log(Ls[i + 1] / Ls[i]) for i in range(len(Ls) - 1)}, A


# ------------------------------------------------------------------ boson-frame and conformal (chord) fits
def Gb_sym(run):
    """|<b+_i b_j>| for pairs symmetric about the centre (boson frame, no string), r = 2m."""
    G = run['Gb']; a0 = run['gwin'][0]; L = run['L']
    c = L // 2 - a0
    ms = np.arange(1, min(c, G.shape[0] - c - 1) + 1)
    return 2 * ms, np.array([abs(G[c - m, c + m]) for m in ms])


def conf_fit_sym(r, y, L, lo, hi, w=8):
    """symmetric pairs: log env = a - eta log[(2L/pi) sin(pi r/2L)] + g log cos(pi r/2L); returns eta, g."""
    e = envelope(y, w)
    sel = (r >= lo) & (r <= hi) & (e > 0)
    X = np.vstack([np.ones(sel.sum()), -np.log((2 * L / np.pi) * np.sin(np.pi * r[sel] / (2 * L))),
                   np.log(np.cos(np.pi * r[sel] / (2 * L)))]).T
    c = np.linalg.lstsq(X, np.log(e[sel]), rcond=None)[0]
    return float(c[1]), float(c[2])


def conf_fit_row(r, y, L, i0, lo, hi, w=8):
    """rows from a fixed site i0 to j = i0 + r (both measured from the left wall at -1):
    log env = a - eta log d(r) + g log[d(x+y)^2 / (d(2x) d(2y))], d(u) = (2L'/pi) sin(pi u / 2L'), L' = L + 1,
    x = i0 + 1, y = i0 + 1 + r. Returns eta, g."""
    Lp = L + 1
    d = lambda u: (2 * Lp / np.pi) * np.sin(np.pi * u / (2 * Lp))
    e = envelope(y, w)
    sel = (r >= lo) & (r <= hi) & (e > 0)
    x = i0 + 1; yy = i0 + 1 + r[sel]
    X = np.vstack([np.ones(sel.sum()), -np.log(d(r[sel])), np.log(d(x + yy) ** 2 / (d(2 * x) * d(2 * yy)))]).T
    c = np.linalg.lstsq(X, np.log(e[sel]), rcond=None)[0]
    return float(c[1]), float(c[2])


def peak_points(r, y, h=4):
    """centred RMS envelope over +-h samples (unbiased to first order for a monotone power law; for an oscillating
    function it is the amplitude over sqrt 2). Points within h of either end are dropped."""
    a = np.abs(np.asarray(y)) ** 2
    idx = np.arange(h, len(a) - h)
    env = np.array([np.sqrt(np.mean(a[i - h: i + h + 1])) for i in idx])
    return r[idx], env


def fit_chord(r, y, L, lo, hi, kind='sym', i0=None, h=3, image=False):
    """exponent from peak points with the chord distance; kind 'sym' (symmetric pairs, d = (2L/pi) sin(pi r/2L))
    or 'row' (fixed i0, d(r) with L' = L + 1, optional image factor log[d(x+y)^2/(d(2x) d(2y))])."""
    rp, yp = peak_points(r, y, h)
    sel = (rp >= lo) & (rp <= hi)
    rp, yp = rp[sel], yp[sel]
    if kind == 'sym':
        X = np.vstack([np.ones(len(rp)), -np.log((2 * L / np.pi) * np.sin(np.pi * rp / (2 * L)))]).T
    else:
        Lp = L + 1
        d = lambda u: (2 * Lp / np.pi) * np.sin(np.pi * u / (2 * Lp))
        cols = [np.ones(len(rp)), -np.log(d(rp))]
        if image:
            x = i0 + 1; yy = i0 + 1 + rp
            cols.append(np.log(d(x + yy) ** 2 / (d(2 * x) * d(2 * yy))))
        X = np.vstack(cols).T
    c = np.linalg.lstsq(X, np.log(yp), rcond=None)[0]
    return float(c[1])


def D_rows(run, sep=2):
    return run['D2'] if sep == 2 else run['D3']
