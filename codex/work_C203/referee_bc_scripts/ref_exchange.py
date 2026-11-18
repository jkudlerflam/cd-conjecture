"""Referee C203: Theorem B checks with own code.
(1) sharp V=0 bound with Haar-random W, (2) domination E_L(W) >= E_L(sc),
(3) explicit seam holonomy Theta: ring spectrum in sector m equals the union over
    eigenphases of Theta of twisted scalar spectra, (4) M >= t by SDP, (5) V <= -2t.
"""
import itertools, sys
import numpy as np
import scipy.sparse as sp
from ref_common import digits, two_site_op, ring_H, lowest, gauge_M, haar_unitary

rng = np.random.default_rng(20261118)
C_SHARP = 4 - 8 / np.pi


def bond(d0, d1, U, t, V, mu):
    d = d0 + d1
    n = np.diag([0.0] * d0 + [1.0] * d1)
    W = np.zeros((d * d, d * d), complex)
    ins = [(c, e) for c in range(d0, d) for e in range(d0)]      # H1 (x) H0
    outs = [(a, b) for a in range(d0) for b in range(d0, d)]     # H0 (x) H1
    for oi, (a, b) in enumerate(outs):
        for ii, (c, e) in enumerate(ins):
            W[a * d + b, c * d + e] = U[oi, ii]
    I = np.eye(d)
    h = -t * (W + W.conj().T) + V * np.kron(n, n) + 0.5 * mu * (np.kron(n, I) + np.kron(I, n))
    return h, W


def e0_free(t, mu):
    xs = [0.0, 1.0]
    if abs(mu) < 2 * t:
        xs.append(np.arccos(mu / (2 * t)) / np.pi)
    return min(mu * x - 2 * t / np.pi * np.sin(np.pi * x) for x in xs)


def Esc_free(L, t, mu):
    return min(mu * m - 2 * t * np.sin(np.pi * m / L) / np.sin(np.pi / L) for m in range(L + 1))


def sharp_checks():
    print("== (1),(2),(4) V=0 sharp bound, domination, M >= t ==")
    worst = 0.0
    shapes = [(1, 1), (1, 2), (2, 1), (2, 2), (1, 3), (3, 1)]
    for (d0, d1) in shapes:
        d = d0 + d1
        Lmax = {2: 12, 3: 8, 4: 7}[d]
        for trial in range(2):
            U = haar_unitary(d0 * d1, rng)
            t = [1.0, 0.63][trial]
            for mu in [-2.5, -0.9, 0.0, 0.4, 1.3, 2.6]:
                h, _ = bond(d0, d1, U, t, 0.0, mu)
                M = gauge_M(h, d)
                e0 = e0_free(t, mu)
                for L in range(2, Lmax + 1):
                    EL = lowest(ring_H(h, L, d))[0]
                    Esc = Esc_free(L, t, mu)
                    dom = EL - Esc
                    q = L * (L * e0 - EL)
                    bnd = 2 * t * L * (1 / np.sin(np.pi / L) - L / np.pi)
                    rt = q / (C_SHARP * t)
                    rM = q / (C_SHARP * M)
                    worst = max(worst, rM)
                    flag = "" if (dom > -1e-8 and q <= bnd + 1e-8 and M >= t - 1e-6) else "  <-- VIOLATION"
                    if L in (2, 3, Lmax) or flag:
                        print(f"d0,d1={d0},{d1} t={t} mu={mu:5.2f} L={L:2d} EL={EL:.10f} EL-Esc={dom:.2e} "
                              f"L(Le0-EL)={q:.6f} bnd={bnd:.6f} /((4-8/pi)t)={rt:.6f} M={M:.6f} /((4-8/pi)M)={rM:.6f}{flag}")
    print("largest L(Le0-EL)/((4-8/pi)M) over all cases:", worst)


def sector_words(L, m):
    return [s for s in itertools.product([0, 1], repeat=L) if sum(s) == m]


def fiber_index(D, s, d0):
    cls = (D >= d0).astype(int)
    return np.where((cls == np.array(s)).all(axis=1))[0]


def scalar_twisted(L, m, t, V, mu, phi):
    words = sector_words(L, m)
    pos = {w: k for k, w in enumerate(words)}
    H = np.zeros((len(words), len(words)), complex)
    for w in words:
        k = pos[w]
        H[k, k] += V * sum(w[j] * w[(j + 1) % L] for j in range(L)) + mu * m
        for j in range(L):
            jj = (j + 1) % L
            if w[j] == 1 and w[jj] == 0:
                w2 = list(w); w2[j], w2[jj] = 0, 1; w2 = tuple(w2)
                amp = -t * (np.exp(1j * phi) if jj == 0 else 1.0)  # seam bond (L-1, 0)
                H[pos[w2], k] += amp
                H[k, pos[w2]] += np.conj(amp)
    return np.linalg.eigvalsh(H)


def holonomy_check(d0, d1, L, t, V, mu):
    d = d0 + d1
    U = haar_unitary(d0 * d1, rng)
    h, W = bond(d0, d1, U, t, V, mu)
    D = digits(L, d)
    Wops = [two_site_op(W, j, (j + 1) % L, L, d, D).tocsr() for j in range(L)]
    H = ring_H(h, L, d).tocsr()
    maxdev_spec, maxdev_cons, maxdev_theta = 0.0, 0.0, 0.0
    for m in range(1, L):
        words = sector_words(L, m)
        idx = {s: fiber_index(D, s, d0) for s in words}
        sorted_w = tuple([0] * (L - m) + [1] * m)
        Us = {sorted_w: np.eye(len(idx[sorted_w]), dtype=complex)}
        queue = [sorted_w]
        while queue:
            sp_ = queue.pop(0)
            for j in range(L - 1):  # open bonds only
                if sp_[j] == 0 and sp_[j + 1] == 1:
                    s = list(sp_); s[j], s[j + 1] = 1, 0; s = tuple(s)
                    if s not in Us:
                        Wblk = Wops[j][idx[sp_]][:, idx[s]].toarray()
                        Us[s] = Us[sp_] @ Wblk
                        queue.append(s)
        assert len(Us) == len(words)
        # consistency on every open edge: U_s = U_{s'} W_j
        for s in words:
            for j in range(L - 1):
                if s[j] == 1 and s[j + 1] == 0:
                    s2 = list(s); s2[j], s2[j + 1] = 0, 1; s2 = tuple(s2)
                    Wblk = Wops[j][idx[s2]][:, idx[s]].toarray()
                    maxdev_cons = max(maxdev_cons, np.abs(Us[s] - Us[s2] @ Wblk).max())
        # seam holonomy, bond (L-1, 0)
        thetas = []
        for s in words:
            if s[L - 1] == 1 and s[0] == 0:
                s2 = list(s); s2[L - 1], s2[0] = 0, 1; s2 = tuple(s2)
                Wblk = Wops[L - 1][idx[s2]][:, idx[s]].toarray()
                thetas.append(Us[s2] @ Wblk @ np.linalg.inv(Us[s]))
        Th = thetas[0]
        for T in thetas[1:]:
            maxdev_theta = max(maxdev_theta, np.abs(T - Th).max())
        phases = np.angle(np.linalg.eigvals(Th))
        union = np.sort(np.concatenate([scalar_twisted(L, m, t, V, mu, ph) for ph in phases]))
        allidx = np.sort(np.concatenate([idx[s] for s in words]))
        spec = np.linalg.eigvalsh(H[allidx][:, allidx].toarray())
        maxdev_spec = max(maxdev_spec, np.abs(np.sort(spec) - union).max())
    return maxdev_cons, maxdev_theta, maxdev_spec


def holonomy_checks():
    print("== (3) seam holonomy: ring spectrum = union of twisted scalar spectra ==")
    for (d0, d1, L) in [(1, 2, 2), (1, 2, 3), (1, 2, 4), (1, 2, 5), (2, 1, 5), (2, 2, 2), (2, 2, 3), (2, 2, 4), (2, 2, 5),
                        (1, 3, 4), (3, 1, 5), (2, 3, 4)]:
        for (V, mu) in [(0.0, 0.3), (0.85, -0.4), (-1.3, 0.7)]:
            c, th, s = holonomy_check(d0, d1, L, 1.0, V, mu)
            print(f"d0,d1={d0},{d1} L={L} V={V} mu={mu}: path-indep dev={c:.1e}, seam-Theta s-dependence={th:.1e}, spectrum dev={s:.1e}")


def ferro_checks():
    print("== (5) V <= -2t: E_L = L min(0, V+mu) ==")
    for (d0, d1) in [(1, 2), (2, 2), (1, 3)]:
        d = d0 + d1
        U = haar_unitary(d0 * d1, rng)
        for (V, mu) in [(-2.0, 0.5), (-2.0, 3.0), (-3.1, 1.0), (-2.4, 3.0)]:
            h, _ = bond(d0, d1, U, 1.0, V, mu)
            for L in [2, 3, 5]:
                EL = lowest(ring_H(h, L, d))[0]
                print(f"d0,d1={d0},{d1} V={V} mu={mu} L={L}: EL - L*min(0,V+mu) = {EL - L * min(0, V + mu):.2e}")


def M_compare():
    print("== M(h_W) vs M(h_sc) and t ==")
    for (d0, d1) in [(1, 2), (2, 2), (1, 3), (2, 3)]:
        d = d0 + d1
        for (t, V, mu) in [(1.0, 0.0, 0.0), (1.0, 0.0, 1.7), (0.7, 1.1, -0.3), (1.0, -2.6, 0.9)]:
            U = haar_unitary(d0 * d1, rng)
            h, _ = bond(d0, d1, U, t, V, mu)
            hs, _ = bond(1, 1, np.eye(1), t, V, mu)
            print(f"d0,d1={d0},{d1} t={t} V={V} mu={mu}: M(W)={gauge_M(h, d):.6f} M(sc)={gauge_M(hs, 2):.6f}")


if __name__ == "__main__":
    which = sys.argv[1:] or ["sharp", "hol", "ferro", "M"]
    if "hol" in which: holonomy_checks()
    if "ferro" in which: ferro_checks()
    if "M" in which: M_compare()
    if "sharp" in which: sharp_checks()
