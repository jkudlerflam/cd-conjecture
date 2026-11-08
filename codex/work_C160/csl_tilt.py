"""C160 (Claude Opus agent, 2026-11-08). Classical chiral soliton lattice (CSL) in the K stratum near h_F
and its out-of-plane (tilt) instability (PROOF section 3.4).

Continuum functional (length unit 1/eps, energy unit eps^2 per site, from the lattice bond
v = (-1, -1, -1 + kappa eps^2/2, 0, 0, -eps, -h eps^2, 0, 0) in est.py coordinates):
    E[n] = int dx [ |n'|^2/2 - (n x n')_z + (kappa/2) n_z^2 - h n_x ].
Planar ansatz n = (cos t, sin t, 0): E = int [ t'^2/2 - t' - h cos t ] (sine-Gordon with twist, CSL for h < pi^2/16).
Out-of-plane Hessian on the planar CSL: H_z = -d^2 + W0(x) + kappa, W0 = -t'^2 + 2 t' + h cos t.
Outputs, for several h:
  - CSL period T(h);
  - kappa_c(h) = -lambda_0(-d^2 + W0) (periodic ground state) and the Bloch band minimum position (must be k = 0);
  - the symplectic pairing  Omega = int t' psi_0 dx / (||t'|| ||psi_0||)  (must be nonzero, and positive);
  - full two-component minimization across kappa_c: tilt amplitude, period and energy against the uniform canted state.
"""
import numpy as np
from scipy.optimize import minimize, minimize_scalar


def csl_planar(h, T, N=256):
    """minimize the planar functional on a ring of length T with winding 2 pi."""
    dx = T / N
    x = np.arange(N) * dx
    def f(dt):
        t = 2 * np.pi * x / T + dt
        tp = np.roll(t, -1) - t; tp[-1] += 2 * np.pi
        tp /= dx
        E = np.sum(0.5 * tp ** 2 - tp - h * np.cos(t)) * dx
        g = np.zeros(N)
        # d/dt_i of 0.5 tp_i^2 - tp_i (tp_i = (t_{i+1}-t_i)/dx) and of tp_{i-1}
        gi = (tp - 1)
        g += -gi + np.roll(gi, 1)
        g += h * np.sin(t) * dx
        return E, g
    r = minimize(f, np.zeros(N), jac=True, method='L-BFGS-B', options=dict(maxiter=5000, gtol=1e-10))
    t = 2 * np.pi * x / T + r.x
    return r.fun / T, t, dx


def csl(h, N=256):
    res = minimize_scalar(lambda T: csl_planar(h, T, N)[0], bounds=(2.0, 60.0), method='bounded',
                          options=dict(xatol=1e-5))
    T = res.x
    e, t, dx = csl_planar(h, T, N)
    return T, e, t, dx


def hill(W, dx, k=0.0):
    """Bloch Hamiltonian of -d^2 + W on a ring of N points with Bloch phase k*T."""
    N = len(W)
    H = np.diag(2 / dx ** 2 + W).astype(complex)
    ph = np.exp(1j * k * N * dx)
    for i in range(N - 1):
        H[i, i + 1] = H[i + 1, i] = -1 / dx ** 2
    H[N - 1, 0] = -ph / dx ** 2; H[0, N - 1] = -np.conj(ph) / dx ** 2
    return np.linalg.eigh(H)


def full_energy(h, kappa, T, N=128, init=None, tilt0=0.3):
    """two-component minimization on a ring of length T, in-plane winding 2 pi."""
    dx = T / N
    x = np.arange(N) * dx
    def unpack(p):
        th = np.pi / 2 - p[:N]            # p[:N] = latitude psi (n_z = sin psi)
        ph = 2 * np.pi * x / T + p[N:]
        return th, ph
    def f(p):
        th, ph = unpack(p)
        n = np.stack([np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)], 1)
        n1 = np.roll(n, -1, 0)
        d = (n1 - n) / dx
        cz = n[:, 0] * n1[:, 1] - n[:, 1] * n1[:, 0]
        E = np.sum(0.5 * np.sum(d ** 2, 1) - cz / dx + 0.5 * kappa * n[:, 2] ** 2 - h * n[:, 0]) * dx
        # gradient w.r.t. n
        G = np.zeros_like(n)
        G += (-(d) + np.roll(d, 1, 0)) / dx * dx      # from |d|^2/2 terms
        G[:, 0] += -n1[:, 1] / dx * dx; G[:, 1] += n1[:, 0] / dx * dx
        n0 = np.roll(n, 1, 0)
        G[:, 0] += n0[:, 1] / dx * dx; G[:, 1] += -n0[:, 0] / dx * dx
        G[:, 2] += kappa * n[:, 2] * dx
        G[:, 0] += -h * dx
        dth = np.stack([np.cos(th) * np.cos(ph), np.cos(th) * np.sin(ph), -np.sin(th)], 1)
        dph = np.stack([-np.sin(th) * np.sin(ph), np.sin(th) * np.cos(ph), 0 * th], 1)
        gpsi = -np.sum(G * dth, 1)
        gph = np.sum(G * dph, 1)
        return E, np.r_[gpsi, gph]
    p0 = np.r_[tilt0 * np.ones(N), np.zeros(N)] if init is None else init
    r = minimize(f, p0, jac=True, method='L-BFGS-B', options=dict(maxiter=20000, gtol=1e-11))
    th, ph = unpack(r.x)
    return r.fun / T, np.cos(th), r.x


def uniform_min(h, kappa):
    c = np.linspace(-1, 1, 200001)          # c = n_x, in the xz plane (n_y = 0 is optimal)
    s2 = 1 - c ** 2
    return np.min(-h * c + 0.5 * kappa * s2)


if __name__ == '__main__':
    print('h      T_CSL     e_CSL     kappa_c   argmin_k(band)  pairing   (full-min checks below)')
    rows = []
    for h in (0.05, 0.15, 0.3, 0.45, 0.55):
        T, e, t, dx = csl(h)
        tp = (np.roll(t, -1) - t); tp[-1] += 2 * np.pi; tp /= dx
        tp = 0.5 * (tp + np.roll(tp, 1))
        W0 = -tp ** 2 + 2 * tp + h * np.cos(t)
        w, v = hill(W0, dx)
        lam0 = w[0]; psi0 = v[:, 0].real; psi0 *= np.sign(psi0.sum())
        ks = np.linspace(0, np.pi / T, 9)
        band = [hill(W0, dx, k)[0][0] for k in ks]
        kmin = ks[int(np.argmin(band))]
        pair = np.dot(tp, psi0) / np.linalg.norm(tp) / np.linalg.norm(psi0)
        rows.append((h, T, e, -lam0))
        print(f'{h:.2f}  {T:8.4f}  {e:+.5f}  {-lam0:+.5f}   {kmin:.4f}          {pair:.4f}   min psi0 = {psi0.min():.2e} (nodeless if > 0)')
    print()
    print('full two-component minimization across kappa_c (period optimized), compared with the uniform state')
    for h, T0, e0, kc in rows[1:4]:
        for dk in (+0.15, +0.05, -0.05, -0.15):
            kappa = kc + dk
            best = None
            for T in np.linspace(0.8 * T0, 1.6 * T0, 9):
                ee, nz, p = full_energy(h, kappa, T)
                if best is None or ee < best[0]:
                    best = (ee, T, np.abs(nz).max(), nz.mean())
            eu = uniform_min(h, kappa)
            print(f'h={h:.2f} kappa={kappa:+.4f} (kappa-kappa_c={dk:+.2f}): texture e={best[0]:+.5f} T={best[1]:.3f} '
                  f'max|n_z|={best[2]:.3f} <n_z>={best[3]:+.3f} ; uniform e={eu:+.5f} -> '
                  f'{"texture" if best[0] < eu - 1e-6 else "uniform"}')
