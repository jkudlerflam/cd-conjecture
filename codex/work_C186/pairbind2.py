"""C186: binding of two dressed quasiparticles in the mixed cusp model (PROOF 3.1) on N-capped rings, tuned on the ring.

Ring of L sites, particle number <= Nmax (ring_mixed.Ring). E1(q) = lowest level with <N> < 1.7 at grid momentum q.
Tuning (on the ring): phi such that the two valley minima of E1 are degenerate, then mu such that the lower one is 0
(the dilute onset of single quasiparticles). Then B = T(K*) - E_low(K*) with K* = q_a + q_b (grid minima) and
T(K) = min_q E1(q) + E1(K - q). Bound pairs: E_low < T, B > 0. Calibrated at r = 0 against the exact two-body
binding energies (u = 2: ring L = 48 contains the valley momenta).
"""
import numpy as np, json, sys, time
from scipy.optimize import brentq
from ring_mixed import Ring, model


class Tuner:
    def __init__(self, L, Nmax, u, v, rmode):
        self.ring = Ring(L, Nmax); self.L = L; self.u = u; self.v = v
        self.r = 0.0 if rmode == "zero" else np.sqrt(v)
        s = np.arcsin(u / 4)
        g = lambda k: int(round((k % (2 * np.pi)) * L / (2 * np.pi))) % L
        self.ia0, self.ib0 = g(-s), g(np.pi + s)

    def E1(self, i, phi, mu):
        e, Nx = self.ring.lowest(model(self.u, phi, self.v, mu, r=self.r), 2 * np.pi * i / self.L, k=5)
        for ei, ni in zip(e, Nx):
            if 0.5 < ni < 1.7:
                return ei
        return np.inf

    def valley(self, i0, phi, mu):
        vals = [(self.E1((i0 + d) % self.L, phi, mu), (i0 + d) % self.L) for d in (-3, -2, -1, 0, 1, 2, 3)]
        return min(vals)

    def tune(self, phi0=np.pi / 2, span=0.15):
        u = self.u
        mu = -2 - u ** 2 / 4
        phi = phi0
        for it in range(3):
            f = lambda ph: self.valley(self.ia0, ph, mu)[0] - self.valley(self.ib0, ph, mu)[0]
            if self.r == 0:
                phi = phi0
            else:
                fa, fb = f(phi0 - span), f(phi0 + span)
                if fa * fb > 0:
                    raise RuntimeError(f"phi window: {fa} {fb}")
                phi = brentq(f, phi0 - span, phi0 + span, xtol=1e-6)
            g = lambda m_: min(self.valley(self.ia0, phi, m_)[0], self.valley(self.ib0, phi, m_)[0])
            hi = mu + 0.1
            while not g(hi) < 0:
                hi -= 0.05
            lo = hi - 0.2
            while not g(lo) > 0:
                lo -= 0.2
            mu = brentq(g, lo, hi, xtol=1e-7)
        self.phi, self.mu = phi, mu
        (ea, ia), (eb, ib) = self.valley(self.ia0, phi, mu), self.valley(self.ib0, phi, mu)
        return phi, mu, ea, ia, eb, ib


def run(L, Nmax, u, v, rmode):
    t0 = time.time()
    tu = Tuner(L, Nmax, u, v, rmode)
    phi, mu, ea, ia, eb, ib = tu.tune()
    E1 = np.array([tu.E1(i, phi, mu) for i in range(L)])
    iK = (ia + ib) % L
    T = min(E1[j] + E1[(iK - j) % L] for j in range(L))
    e, Nx = tu.ring.lowest(model(u, phi, v, mu, r=tu.r), 2 * np.pi * iK / L, k=4)
    out = dict(L=L, Nmax=Nmax, u=u, v=v, rmode=rmode, r=tu.r, phi=phi, mu=mu, qa=2 * np.pi * ia / L,
               qb=2 * np.pi * ib / L, ea=ea, eb=eb, K=2 * np.pi * iK / L, T=float(T), E_low=float(e[0]),
               N_low=float(Nx[0]), B=float(T - e[0]), levels=[float(x) for x in e], Ns=[float(x) for x in Nx],
               sec=round(time.time() - t0, 1))
    print(json.dumps(out), flush=True)
    return out


if __name__ == "__main__":
    u = float(sys.argv[1]); rmode = sys.argv[2]
    vs = [float(x) for x in sys.argv[3].split(',')]
    L = int(sys.argv[4]); Nmaxs = [int(x) for x in sys.argv[5].split(',')]
    for v in vs:
        for Nmax in Nmaxs:
            run(L, Nmax, u, v, rmode)
