"""Band structure, Fermi points, resonance margins, Fermi-point couplings and self-consistent
Hartree-Fock for a 16 x 16 two-mode bond (conventions of bonds.py).

usage: python3 -B bands.py CELL x [x ...]       (prints one JSON line per coupling)
Written by Claude (Opus 5.5) for contrib/routeA_n2_chains, 2026-10-28.
"""
import sys, json, itertools
import numpy as np
from scipy.optimize import brentq
import bonds as B


def blocks(h):
    """One-body bond h (4x4 on (j,0),(j,1),(j+1,0),(j+1,1)) -> T, M with
    H = sum_j c_j^+ M c_j + c_j^+ T c_{j+1} + h.c."""
    return h[:2, 2:], h[:2, :2] + h[2:, 2:]


def bloch(T, M, k):
    return M + T * np.exp(1j * k) + T.conj().T * np.exp(-1j * k)


def fermi_points(T, M, mu=0.0, nk=20001):
    """Roots of the bands of h(k) - mu, with band index, velocity and spinor."""
    ks = np.linspace(-np.pi, np.pi, nk)
    E = np.array([np.linalg.eigvalsh(bloch(T, M, k)) for k in ks]) - mu
    out = []
    for b in range(2):
        e = E[:, b]
        for i in np.flatnonzero(e[:-1] * e[1:] < 0):
            f = lambda k: np.linalg.eigvalsh(bloch(T, M, k))[b] - mu
            k0 = brentq(f, ks[i], ks[i + 1], xtol=1e-15)
            w, U = np.linalg.eigh(bloch(T, M, k0))
            dh = 1j * T * np.exp(1j * k0) - 1j * T.conj().T * np.exp(-1j * k0)
            v = float(np.vdot(U[:, b], dh @ U[:, b]).real)
            out.append(dict(k=float(k0), band=b, v=v, u=U[:, b]))
    return sorted(out, key=lambda r: r['k']), E + mu


def density(T, M, mu, nk=4001):
    ks = np.linspace(-np.pi, np.pi, nk, endpoint=False)
    E = np.array([np.linalg.eigvalsh(bloch(T, M, k)) for k in ks])
    return float((E < mu).sum() / nk)


def mu_for_density(T, M, nu, nk=4001):
    ks = np.linspace(-np.pi, np.pi, nk, endpoint=False)
    E = np.sort(np.array([np.linalg.eigvalsh(bloch(T, M, k)) for k in ks]).ravel())
    n = int(round(nu * nk))
    return 0.5 * (E[n - 1] + E[n])


def margins(pts):
    """delta2 = min distance between Fermi points; for every 4-fermion process
    psi^+_a psi^+_b psi_c psi_d with {a,b} != {c,d}: the momentum mismatch mod 2 pi.
    Returns the sorted list of distinct small mismatches with labels."""
    seas(pts)
    p = [r['k'] for r in pts]
    lab = [('R' if r['v'] > 0 else 'L') + str(r.get('sea', 'x')) for r in pts]
    d = lambda x: abs((x + np.pi) % (2 * np.pi) - np.pi)
    d2 = min(d(p[a] - p[b]) for a, b in itertools.combinations(range(len(p)), 2))
    pairs = list(itertools.combinations(range(len(p)), 2))       # distinct points (Pauli)
    procs = []
    for (a, b), (c, e) in itertools.combinations(pairs, 2):
        procs.append((d(p[a] + p[b] - p[c] - p[e]), f'{lab[a]}{lab[b]}<-{lab[c]}{lab[e]}'))
    procs.sort()
    return d2, procs


def seas(pts):
    """Group Fermi points into Fermi seas: in each band an occupied arc starts at a point with
    v < 0 (left edge L) and ends at the next point with v > 0 (right edge R), going up in k.
    Returns a list of (L, R) point dicts, and labels each point with its sea index."""
    out = []
    for b in (0, 1):
        q = sorted([r for r in pts if r['band'] == b], key=lambda r: r['k'])
        for i, r in enumerate(q):
            if r['v'] < 0:
                for j in range(1, len(q) + 1):
                    s_ = q[(i + j) % len(q)]
                    if s_['v'] > 0:
                        out.append((r, s_)); break
    for a, (L, R) in enumerate(out):
        L['sea'] = a; R['sea'] = a
    return out


def centres(pts):
    out = {}
    for a, (L, R) in enumerate(seas(pts)):
        kr = R['k'] if R['k'] >= L['k'] else R['k'] + 2 * np.pi
        out[a] = float(np.angle(np.exp(1j * (L['k'] + kr) / 2)))
    return out


def wedge(x, y):
    return np.array([x[i] * y[j] - x[j] * y[i] for i in range(4) for j in range(i + 1, 4)])


def amp(U, o1, o2, i1, i2):
    """C28 eq. (13): antisymmetrized continuum vertex from the bond quartic matrix U."""
    w = lambda r: np.concatenate([r['u'], np.exp(1j * r['k']) * r['u']])
    return complex(np.vdot(wedge(w(o1), w(o2)), U @ wedge(w(i1), w(i2))))


def couplings(U, pts):
    """Forward couplings G_ab = A(R_a L_b; R_a L_b) and the pair transfer p = A(R1 L1; R0 L0)
    between the two Fermi seas a = 0, 1, and C28's one-loop combination x, y
    (x < -|y|: pair transfer irrelevant at one loop)."""
    S = seas(pts)
    if len(S) != 2:
        return None
    L = {a: S[a][0] for a in (0, 1)}; R = {a: S[a][1] for a in (0, 1)}
    G = np.array([[amp(U, R[a], L[b], R[a], L[b]).real for b in (0, 1)] for a in (0, 1)])
    p = amp(U, R[1], L[1], R[0], L[0])
    v = [0.5 * (R[a]['v'] - L[a]['v']) for a in (0, 1)]
    x = -G[0, 0] / (2 * np.pi * v[0]) - G[1, 1] / (2 * np.pi * v[1]) + (G[0, 1] + G[1, 0]) / (np.pi * (v[0] + v[1]))
    y = abs(p) / np.pi * np.sqrt(1 / (2 * v[0] * v[1]) + 2 / (v[0] + v[1]) ** 2)
    gdim = np.array([[G[a, b] / (np.pi * (v[a] + v[b])) for b in (0, 1)] for a in (0, 1)])
    return dict(G=G.tolist(), pair=[p.real, p.imag], v=v, x=float(x), y=float(y), G_dimless=gdim.tolist(),
                pair_dimless=float(abs(p) / (np.pi * (v[0] + v[1]))))


def corr(T, M, mu, nk=4001):
    """rho(r)_{ab} = <c^+_{j,a} c_{j+r,b}> for r = 0, 1 at chemical potential mu."""
    ks = np.linspace(-np.pi, np.pi, nk, endpoint=False)
    r0 = np.zeros((2, 2), complex); r1 = np.zeros((2, 2), complex)
    for k in ks:
        w, Uk = np.linalg.eigh(bloch(T, M, k))
        occ = Uk[:, w < mu]
        P = occ.conj() @ occ.T                 # P[a,b] = sum u_a^* u_b
        r0 += P; r1 += P * np.exp(1j * k)
    return r0 / nk, r1 / nk


def hf_onebody(K, rho_bond):
    """Self-consistent-field one-body operator of the bond K at the Gaussian bond density matrix
    rho_bond[a,b] = <c_a^+ c_b>: F_ab = d<K>/d rho_ab, computed by Wick on the full 16 x 16 bond
    (all 3- and 4-body terms included) with a finite-difference-free formula: <K> for a Gaussian
    state is evaluated exactly and differentiated numerically (central differences)."""
    def expK(rho):
        # Gaussian state on 4 modes with <c_a^+ c_b> = rho[a,b]: number-conserving, so its density
        # matrix is exp(-sum c^+ X c)/Z with rho^T = 1/(e^X + 1); build it from rho via X.
        r = 0.5 * (rho + rho.conj().T)
        w, V = np.linalg.eigh(r.T)
        w = np.clip(w.real, 1e-12, 1 - 1e-12)
        X = V @ np.diag(np.log(1 / w - 1)) @ V.conj().T
        A = sum(X[a, b] * B.CS[a].conj().T @ B.CS[b] for a in range(4) for b in range(4))
        ev, Ev = np.linalg.eigh(A)
        rhoM = Ev @ np.diag(np.exp(-(ev - ev.min()))) @ Ev.conj().T
        rhoM /= np.trace(rhoM)
        return np.trace(rhoM @ K)
    F = np.zeros((4, 4), complex)
    eps = 1e-6
    for a in range(4):
        for b in range(4):
            for ph in (1, 1j):
                d = np.zeros((4, 4), complex); d[a, b] += eps * ph; d[b, a] += np.conj(eps * ph)
                g = (expK(rho_bond + d) - expK(rho_bond - d)) / (2 * eps)
                # <K> linear part: sum F_ab rho_ab ; derivative along d gives F_ab ph + F_ba conj(ph)
                if ph == 1:
                    re = g
                else:
                    im = g
            # Solve F_ab + F_ba (re) and i F_ab - i F_ba (im) for a != b; diagonal: re = 2 F_aa
            if a == b:
                F[a, a] = re.real / 2
            else:
                F[a, b] = 0.5 * (re - 1j * im) if False else 0.5 * (re + im / 1j)
    # enforce hermiticity
    return 0.5 * (F + F.conj().T)


def hartree_fock(K, nu=None, mu=0.0, iters=300, mix=0.3, nk=2001):
    """Self-consistent Hartree-Fock of the chain sum_j K(j,j+1), at fixed density nu per cell
    (mu adjusted) or, if nu is None, at fixed chemical potential mu (K already contains the fields).
    Returns T, M of the converged one-body bond F, the chemical potential, iterations, F (4x4)."""
    c, h, U, Hq, H4, r = B.normal_order(K)
    Fh = h.copy()
    for it in range(iters):
        Tf, Mf = blocks(Fh)
        if nu is not None:
            mu = mu_for_density(Tf, Mf, nu, nk)
        r0, r1 = corr(Tf, Mf, mu, nk)
        rb = np.block([[r0, r1], [r1.conj().T, r0]])
        Fnew = hf_onebody(K, rb)
        if np.linalg.norm(Fnew - Fh) < 1e-9 * np.linalg.norm(Fh):
            Fh = Fnew; break
        Fh = (1 - mix) * Fh + mix * Fnew
    Tf, Mf = blocks(Fh)
    if nu is not None:
        mu = mu_for_density(Tf, Mf, nu, nk)
    return Tf, Mf, mu, it, Fh


def report(cell, x, nu=None, hf=True):
    K, info = B.CELLS[cell](x)
    c, h, U, Hq, H4, r = B.normal_order(K)
    T, M = blocks(h)
    s = np.abs(T).max()
    mu0 = 0.0
    nu0 = density(T, M, mu0)
    out = dict(cell=cell, x=x, scale=float(s), nonquartic_res=r / s, quartic_norm=float(np.linalg.norm(U) / s),
               T=(T / s).round(8).tolist(), M=(M / s).round(8).tolist(), density_bare=nu0)
    pts, E = fermi_points(T / s, M / s, 0.0)
    out['bare'] = dict(points=[(round(p['k'], 6), p['band'], p['v']) for p in pts], centres=centres(pts),
                       bandwidths=[float(np.ptp(E[:, b])) for b in (0, 1)], band_ranges=[[float(E[:, b].min()), float(E[:, b].max())] for b in (0, 1)])
    if len(pts) >= 2:
        d2, procs = margins(pts)
        out['bare']['delta2'] = d2
        out['bare']['smallest_mismatches'] = [(round(m, 6), lab) for m, lab in procs[:6]]
    if len(pts) == 4:
        out['bare']['couplings'] = couplings(U / s, pts)
    if hf:
        nuh = nu0 if nu is None else nu
        Tf, Mf, muh, it, _ = hartree_fock(K / s, nuh)
        ptsh, Eh = fermi_points(Tf, Mf, muh)
        out['hf'] = dict(nu=nuh, iters=it, mu=float(muh), points=[(round(p['k'], 6), p['band'], p['v']) for p in ptsh],
                         centres=centres(ptsh))
        if len(ptsh) >= 2:
            d2, procs = margins(ptsh)
            out['hf']['delta2'] = d2
            out['hf']['smallest_mismatches'] = [(round(m, 6), lab) for m, lab in procs[:6]]
        if len(ptsh) == 4:
            out['hf']['couplings'] = couplings(U / s, ptsh)
    return out


if __name__ == '__main__':
    cell = sys.argv[1]
    for x in sys.argv[2:]:
        o = report(cell, float(x))
        print(json.dumps(o, default=lambda z: str(z)), flush=True)
