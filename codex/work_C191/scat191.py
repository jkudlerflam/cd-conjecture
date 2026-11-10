"""C191: exact two-body S-matrix of C173's effective model (3.2) at the Fermi-surface kinematics.

Written independently of C175's scat.py/smatrix.py (same physics, separate code).

Model, units |t2| = 1: hard-core particles, no two adjacent; hop amplitude a_d for a move by +d:
a_{+1} = w, a_{-1} = conj(w), a_{+-2} = -1 (the NNN hop, middle site empty by the blockade);
V(2) = V2, V(3) = V3 (defaults -2, -2/3).  Band eps(k) = sum_d a_d e^{-ikd}.

Two particles x1 < x2, r = x2 - x1 >= rmin (rmin = 2: blockade), psi = e^{iK x1} f(r),
  (h f)(r) = sum_d a_d [ e^{-iKd} f(r+d) + f(r-d) ] + V(r) f(r),  f(r) = 0 for r < rmin.
A mode f = z^r = e^{ipr} has right-particle momentum p, left-particle momentum K - p.
Relative velocity v_rel = eps'(p) - eps'(K-p) (> 0: separating, outgoing).
Ordered channels: 'ab' = left particle in valley a, right in valley b.
"""
import numpy as np


class Model:
    def __init__(self, w=1j, t2=-1.0, V2=-2.0, V3=-2/3, rmin=2, nnn=True):
        self.w, self.t2, self.V2, self.V3, self.rmin = w, t2, V2, V3, rmin
        self.a = {1: w, -1: np.conj(w)}
        if nnn:
            self.a[2] = t2
            self.a[-2] = t2
        self.D = max(abs(d) for d in self.a)

    def eps(self, k):
        return np.real(sum(ad * np.exp(-1j * k * d) for d, ad in self.a.items()))

    def deps(self, k):
        return np.real(sum(ad * (-1j * d) * np.exp(-1j * k * d) for d, ad in self.a.items()))

    def d2eps(self, k):
        return np.real(sum(ad * (-d * d) * np.exp(-1j * k * d) for d, ad in self.a.items()))

    def V(self, r):
        return {2: self.V2, 3: self.V3}.get(r, 0.0)

    def coeffs(self, K):
        """c[s]: coefficient of f(r+s) in (h f)(r), without V."""
        c = {}
        for d, ad in self.a.items():
            c[d] = c.get(d, 0) + ad * np.exp(-1j * K * d)
            c[-d] = c.get(-d, 0) + ad
        return c

    def roots(self, K, E):
        c = self.coeffs(K)
        D = self.D
        poly = [c.get(s, 0) for s in range(D, -D - 1, -1)]
        poly[D] = poly[D] - E           # the z^0 coefficient
        return np.roots(poly)

    def boundary_matrix(self, modes, K):
        """Rows r = rmin..rmin+D-1: true equation minus free equation, applied to each mode."""
        c = self.coeffs(K)
        rows = []
        for r in range(self.rmin, self.rmin + self.D):
            row = []
            for z in modes:
                val = 0j
                for s, cs in c.items():
                    if r + s < self.rmin:
                        val -= cs * z ** (r + s)
                val += self.V(r) * z ** r
                row.append(val)
            rows.append(row)
        return np.array(rows)

    def smatrix(self, K, E, tol=1e-8):
        """Flux-normalized S-matrix. Returns (S, p_in, p_out) with S[i_out, j_in]."""
        z = self.roots(K, E)
        opn = [x for x in z if abs(abs(x) - 1) < tol]
        dec = [x for x in z if abs(x) < 1 - tol]
        ps = np.angle(opn)
        v = np.array([self.deps(p) - self.deps(K - p) for p in ps])
        inc = [i for i in range(len(ps)) if v[i] < 0]
        out = [i for i in range(len(ps)) if v[i] > 0]
        modes = [opn[i] for i in out] + dec
        M = self.boundary_matrix(modes, K)
        S = np.zeros((len(out), len(inc)), complex)
        for j, ii in enumerate(inc):
            rhs = -self.boundary_matrix([opn[ii]], K)[:, 0]
            sol = np.linalg.lstsq(M, rhs, rcond=None)[0]
            S[:, j] = sol[:len(out)] * np.sqrt(np.abs(v[out]) / abs(v[ii]))
        return S, ps[inc], ps[out], v

    def valleys(self):
        ks = np.linspace(-np.pi, np.pi, 400001, endpoint=False)
        e = self.eps(ks)
        idx = np.nonzero((e < np.roll(e, 1)) & (e < np.roll(e, -1)))[0]
        out = []
        for i in idx:
            k = ks[i]
            for _ in range(50):
                k -= self.deps(k) / self.d2eps(k)
            out.append(k)
        out.sort(key=lambda k: self.eps(k))
        return out

    def fermi_points(self, kv, rho):
        """Fermi points kv + pL < kv < kv + pR of a valley at density rho: pR - pL = 2 pi rho,
        eps(kv + pR) = eps(kv + pL).  Returns (pL, pR, eF)."""
        width = 2 * np.pi * rho
        f = lambda pR: self.eps(kv + pR) - self.eps(kv + pR - width)
        lo, hi = 1e-9, width - 1e-9
        flo = f(lo)
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if np.sign(f(mid)) == np.sign(flo):
                lo = mid
            else:
                hi = mid
        pR = 0.5 * (lo + hi)
        return pR - width, pR, self.eps(kv + pR)


def unlike_amplitudes(m, ka, kb, pL, pR):
    """Physical amplitudes for the unlike pair at the Fermi-surface kinematics.
    in_ab : left a at ka+pR (moving right), right b at kb-pR (moving left)   [aR, bL]
    out_ab: left a at ka+pL (moving left),  right b at kb-pL (moving right)  [aL, bR]  (reflection, backscattering)
    out_ba: left b at kb-pR, right a at ka+pR                                  (transmission, momenta unchanged)
    Returns dict with r (in_ab->out_ab), t (in_ab->out_ba), and the same for in_ba (aL, bR)."""
    K = ka + kb
    E = m.eps(ka + pR) + m.eps(kb - pR)
    S, pin, pout, v = m.smatrix(K, E)
    def idx(arr, p):
        d = np.abs(np.angle(np.exp(1j * (arr - p))))
        i = int(np.argmin(d))
        assert d[i] < 1e-6, (arr, p)
        return i
    i_ab = idx(pin, kb - pR)       # right particle b moving left
    i_ba = idx(pin, ka + pL)       # right particle a moving left
    o_ab = idx(pout, kb - pL)      # right particle b moving right
    o_ba = idx(pout, ka + pR)      # right particle a moving right
    res = dict(E=E, K=K,
               r=S[o_ab, i_ab], t=S[o_ba, i_ab],
               rp=S[o_ba, i_ba], tp=S[o_ab, i_ba],
               unit=np.linalg.norm(S.conj().T @ S - np.eye(S.shape[1])),
               nopen=len(pin))
    return res


def like_amplitude(m, kv, pL, pR):
    """Like pair in one valley at the Fermi kinematics: in = (left at kv+pR moving right, right at kv+pL moving left),
    out = (left at kv+pL, right at kv+pR). Fermion reference S = -1 (f ~ sin). Returns S and delta_s with S = -e^{2i delta_s}."""
    K = 2 * kv + pL + pR
    E = m.eps(kv + pR) + m.eps(kv + pL)
    S, pin, pout, v = m.smatrix(K, E)
    def idx(arr, p):
        d = np.abs(np.angle(np.exp(1j * (arr - p))))
        i = int(np.argmin(d)); assert d[i] < 1e-6, (arr, p); return i
    i0 = idx(pin, kv + pL)
    o0 = idx(pout, kv + pR)
    s = S[o0, i0]
    return dict(S=s, delta=0.5 * np.angle(-s), nopen=len(pin), unit=np.linalg.norm(S.conj().T @ S - np.eye(S.shape[1])), K=K, E=E)
