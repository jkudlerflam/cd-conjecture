"""C191: one valley alone in the asymptotic Bethe ansatz (ABA), with the exact like-pair phase of model (3.2).

Bethe-Yang equations for N particles of valley a on a ring of L sites (momenta p_j measured from k_a):
   p_j L = 2 pi I_j + sum_{l != j} Phi(p_j, p_l),   Phi(p_j, p_l) = -tau(p_j, p_l) if v_j > v_l, +tau(p_l, p_j) otherwise,
tau = arg(-S) the like-pair phase relative to free fermions (hard rods: tau = -2 q a, Phi = +a (p_j - p_l), exact).
E = sum eps(k_a + p_j).  K and u from E(N +- 2) and the current J (shift of all I_j).
This treats the intra-valley interaction to all orders in the phase (exact for hard rods and for integrable kernels).
"""
import sys
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scat191 import Model, like_amplitude


def tau_grid(m, kv, pmax, npts=161):
    ps = np.linspace(-pmax, pmax, npts)
    T = np.zeros((npts, npts))
    for i, p1 in enumerate(ps):
        for j, p2 in enumerate(ps):
            if i == j:
                continue
            v1, v2 = m.deps(kv + p1), m.deps(kv + p2)
            if v1 <= v2:          # need left particle faster: in = (left at p1 moving faster)
                continue
            K = 2 * kv + p1 + p2
            E = m.eps(kv + p1) + m.eps(kv + p2)
            S, pin, pout, v = m.smatrix(K, E)
            # incoming mode: right particle at p2 (slower), outgoing: right particle at p1
            di = np.abs(np.angle(np.exp(1j * (pin - (kv + p2)))))
            do = np.abs(np.angle(np.exp(1j * (pout - (kv + p1)))))
            ii, oo = int(np.argmin(di)), int(np.argmin(do))
            T[i, j] = np.angle(-S[oo, ii])
            T[j, i] = T[i, j]
    # fill the diagonal by continuity (tau -> 0 at equal momenta for the like pair)
    return ps, T


class ABA:
    def __init__(self, m, kv, pmax, npts=161):
        self.m, self.kv = m, kv
        ps, T = tau_grid(m, kv, pmax, npts)
        self.interp = RegularGridInterpolator((ps, ps), T)

    def Phi(self, p):
        P1, P2 = np.meshgrid(p, p, indexing='ij')
        tau = self.interp(np.stack([P1, P2], -1))
        v = self.m.deps(self.kv + p)
        faster = v[:, None] > v[None, :]
        Ph = np.where(faster, -tau, tau)
        np.fill_diagonal(Ph, 0.0)
        return Ph

    def solve(self, L, I, p0=None, tol=1e-12):
        I = np.asarray(I, float)
        p = 2 * np.pi * I / L if p0 is None else p0.copy()
        for it in range(500):
            pn = (2 * np.pi * I + self.Phi(p).sum(1)) / L
            if np.max(np.abs(pn - p)) < tol:
                p = pn
                break
            p = 0.5 * p + 0.5 * pn
        return p, float(np.sum(self.m.eps(self.kv + p)))

    def ground(self, L, N, J=0):
        """Minimize over the centre of the quantum-number block (the valley is not inversion-symmetric)."""
        best = None
        base = np.arange(N) - (N - 1) / 2.0
        for c in np.arange(-3, 3.01, 1.0):
            I = base + c + J / 2.0
            p, E = self.solve(L, I)
            if best is None or E < best[1]:
                best = (p, E, c)
        return best

    def KU(self, L, N):
        E0 = self.ground(L, N)
        Ep = self.ground(L, N + 2)
        Em = self.ground(L, N - 2)
        d2N = (Ep[1] + Em[1] - 2 * E0[1]) / 4.0
        u_over_K = L / np.pi * d2N
        # current: shift all I by +-1 around the optimal centre
        c = E0[2]
        base = np.arange(N) - (N - 1) / 2.0 + c
        Ej = [self.solve(L, base + s)[1] for s in (-1, 0, 1)]
        d2J = (Ej[0] + Ej[2] - 2 * Ej[1])      # E = (pi/2L) uK J^2 with J = 2 s
        uK = L / (4 * np.pi) * d2J
        return np.sqrt(uK / u_over_K), np.sqrt(uK * u_over_K)


if __name__ == "__main__":
    # check on hard rods: K = (1 - rho)^2 exactly
    h = Model(w=1.0, nnn=False, V2=0, V3=0)
    kh = h.valleys()[0]
    abah = ABA(h, kh, 1.2, npts=81)
    for rho in (0.1, 0.2):
        L = 400; N = int(rho * L)
        K, u = abah.KU(L, N)
        print(f"hard rods rho={rho}: ABA K={K:.5f} exact (1-rho)^2={(1-rho)**2:.5f}")
    m = Model(w=1j)
    ka = m.valleys()[0]
    aba = ABA(m, ka, 0.75, npts=121)
    for n in [float(eval(a)) for a in sys.argv[1:]] or [1/16, 3/32, 1/8, 5/32, 3/16]:
        rho = n / 2
        for L in (256, 512):
            N = int(round(rho * L))
            K, u = aba.KU(L, N)
            print(f"n={n:.5f} rho={rho:.5f} L={L} N={N}: one-valley ABA K_par={K:.4f} u_par={u:.4f}")
