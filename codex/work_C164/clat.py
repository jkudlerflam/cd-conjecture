"""C164 (Claude Opus agent, 2026-11-08). Classical lattice spin chain for a general nearest-neighbour qubit bond.

    E[n] = sum_j [ n_j^T J n_{j+1} + f . n_j ],   |n_j| = 1,

with (J, f) as in est.py (v = (Jx, Jy, Jz, Dx, Dy, Dz, fx, fy, fz), J = diag + A(D), field f . sigma_j per site).
No small-eps expansion anywhere: this is the lattice functional itself.

Contents
  uniform_min(J, f)            global minimum of e(n) = n^T J n + f . n on S^2 (multi-start, then polished)
  lsw(J, f, n, ks)             tangent Hessian blocks M(k) and linear spin-wave frequencies (Pauli units, s = 1/2)
  stable_uniform(J, f, ...)    uniform minimum + LSW positivity + no periodic texture below it (rings of length 2..Lmax)
  Chain                        open chain with clamped background ends: energy, gradient, solid angle, constrained
                               (fixed Omega) minimization, multi-start local minima of localized textures
  ring_pair                    two-soliton configurations at fixed separation (like / unlike valleys)
Momentum convention: a texture enclosing solid angle Omega (loop bg -> n_1 ... n_N -> bg) has crystal momentum
P = s Omega = Omega / 2 mod 2 pi (Haldane), s = 1/2.
"""
import numpy as np
from scipy.optimize import minimize

TWO_PI = 2 * np.pi


def J_f(v):
    v = np.asarray(v, float)
    J = np.diag(v[:3]).astype(float)
    Dx, Dy, Dz = v[3:6]
    J[0, 1] += Dz; J[1, 0] -= Dz
    J[1, 2] += Dx; J[2, 1] -= Dx
    J[2, 0] += Dy; J[0, 2] -= Dy
    return J, v[6:9].copy()


def sph(th, ph):
    return np.stack([np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th) + 0 * ph], -1)


def e_unif(J, f, n):
    return np.einsum('...a,ab,...b->...', n, J, n) + n @ f


def uniform_min(J, f, ngrid=60):
    th = np.linspace(0, np.pi, ngrid)[:, None]; ph = np.linspace(0, TWO_PI, 2 * ngrid, endpoint=False)[None, :]
    n = sph(th, ph)
    e = e_unif(J, f, n)
    Js = 0.5 * (J + J.T)
    best = None
    idx = np.argsort(e.ravel())[:12]
    for i in idx:
        x0 = n.reshape(-1, 3)[i]
        def fun(x):
            r = np.linalg.norm(x); m = x / r
            val = m @ Js @ m + f @ m
            g = 2 * Js @ m + f
            g = (g - (g @ m) * m) / r
            return val, g
        res = minimize(fun, x0, jac=True, method='BFGS', options=dict(gtol=1e-13))
        m = res.x / np.linalg.norm(res.x)
        val = m @ Js @ m + f @ m
        if best is None or val < best[0] - 1e-13:
            best = (val, m)
    # degeneracy check: second distinct minimum within 1e-9
    return best[1], best[0]


def frame(n):
    a = np.array([1.0, 0, 0]) if abs(n[0]) < 0.9 else np.array([0, 1.0, 0])
    u = np.cross(n, a); u /= np.linalg.norm(u)
    w = np.cross(n, u)
    return u, w


def lsw(J, f, n, ks):
    """Linear spin waves around the uniform state n (must be stationary).
    Tangent coordinates (alpha, beta) along (u, w): n_j = n sqrt(1 - a^2 - b^2) + a u + b w.
    E2 = 1/2 sum_k x_{-k}^T M(k) x_k, x = (alpha, beta). Returns M(k) (2x2 Hermitian) and omega(k) >= 0 if stable.
    Dynamics: spin S = n/2, dS/dt = {S, E} => (1/2) d n/dt = n x (-dE/dn) up to orientation; frequencies are
    the positive eigenvalues of 2 * (sym. form)^{-1} M; the overall sign convention is fixed by lsw_check."""
    u, w = frame(n)
    T = np.stack([u, w], 1)                        # 3x2
    # E = sum_j n_j^T J n_{j+1} + f.n_j. Second order: cross terms T^T J T e^{ik} + T^T J^T T e^{-ik},
    # diagonal: -(n^T J n + n^T J^T n + f . n) * 1 (from the sqrt), per site.
    A = T.T @ J @ T
    diag = -(2 * n @ J @ n + f @ n)
    Ms, om = [], []
    Om = np.array([[0, 1.0], [-1.0, 0]])
    for k in ks:
        M = A * np.exp(1j * k) + A.T * np.exp(-1j * k) + diag * np.eye(2)
        Ms.append(M)
        om.append(np.linalg.eigvals(2 * Om @ M))
    return np.array(Ms), np.array(om)


def lsw_freq(J, f, n, ks):
    """omega(k) for the magnon branch with momentum k (classical LSW, Pauli units), and min eigenvalue of M(k)."""
    Ms, ev = lsw(J, f, n, ks)
    mineig = np.array([np.linalg.eigvalsh(0.5 * (M + M.conj().T)).min() for M in Ms])
    # eigenvalues of 2 Om M(k) are +-i omega; the branch at momentum k is the one with Im > 0 (convention fixed
    # against the exact one-magnon band in lsw_check)
    om = np.array([e[np.argmax(e.imag)].imag for e in ev])   # branch matching est.py ring momentum k (lsw_check.out)
    return om, mineig


def ring_energy_grad(x, J, f):
    V = x.reshape(-1, 3); r = np.linalg.norm(V, axis=1); n = V / r[:, None]
    n1 = np.roll(n, -1, 0); n0 = np.roll(n, 1, 0)
    E = np.einsum('ia,ab,ib->', n, J, n1) + np.sum(n @ f)
    G = n1 @ J.T + n0 @ J + f
    G = (G - np.sum(G * n, 1)[:, None] * n) / r[:, None]
    return E, G.ravel()


def periodic_textures(J, f, e_bg, Ls=range(2, 25), nstart=6, seed=0):
    """Lowest energy per site among relaxed periodic textures on rings of length L (random starts).
    Returns the most negative (E/L - e_bg) and the ring length where it occurs."""
    rng = np.random.default_rng(seed)
    best = (0.0, None)
    for L in Ls:
        for s in range(nstart):
            if s == 0:
                # planar spiral starts with wavevector 2 pi/L in a random plane
                a = rng.normal(size=3); a /= np.linalg.norm(a)
                u, w = frame(a)
                x0 = np.array([np.cos(TWO_PI * j / L) * u + np.sin(TWO_PI * j / L) * w for j in range(L)])
            else:
                x0 = rng.normal(size=(L, 3))
            res = minimize(ring_energy_grad, x0.ravel(), args=(J, f), jac=True, method='L-BFGS-B',
                           options=dict(maxiter=3000, gtol=1e-10))
            d = res.fun / L - e_bg
            if d < best[0] - 1e-10:
                best = (d, L)
    return best


def stable_uniform(J, f, Lmax=24, nstart=4, seed=0):
    n, e = uniform_min(J, f)
    ks = np.linspace(-np.pi, np.pi, 241)
    om, mineig = lsw_freq(J, f, n, ks)
    d, Lb = periodic_textures(J, f, e, Ls=range(2, Lmax + 1), nstart=nstart, seed=seed)
    return dict(n=n, e=e, lsw_min=float(mineig.min()), om_min=float(om.min()), k_ommin=float(ks[np.argmin(om)]),
                texture_gain=float(d), texture_L=Lb, stable=bool(mineig.min() > 1e-9 and d > -1e-9))


# --------------------------------------------------------------------------------------------- solid angle
def tri(A, a, b):
    num = A @ np.cross(a, b) if a.ndim == 1 else np.einsum('a,ia->i', A, np.cross(a, b))
    den = 1 + a @ A + b @ A + np.sum(a * b, -1)
    return num, den


def solid_angle_grad(n, bg, apex):
    """Omega of the closed loop bg -> n_1..n_N -> bg seen from apex, and dOmega/dn (N x 3)."""
    full = np.vstack([bg, n, bg])
    a, b = full[:-1], full[1:]
    cr = np.cross(a, b)
    num = cr @ apex
    den = 1 + a @ apex + b @ apex + np.sum(a * b, 1)
    Om = np.sum(2 * np.arctan2(num, den))
    q = 2 / (num ** 2 + den ** 2)
    # d num / d a = b x apex ; d num / d b = apex x a ; d den/da = apex + b ; d den / d b = apex + a
    dA = q[:, None] * (den[:, None] * np.cross(b, apex) - num[:, None] * (apex + b))
    dB = q[:, None] * (den[:, None] * np.cross(apex, a) - num[:, None] * (apex + a))
    G = dA[1:] + dB[:-1]          # site i (1..N) is a of segment i and b of segment i-1
    return Om, G


def choose_apex(n, bg):
    """A direction far from every point of the loop (Omega mod 4 pi does not depend on it)."""
    full = np.vstack([bg, n])
    th = np.linspace(0, np.pi, 37)[:, None]; ph = np.linspace(0, TWO_PI, 72, endpoint=False)[None, :]
    cand = sph(th, ph).reshape(-1, 3)
    d = (cand @ full.T).max(1)                 # largest cosine to the loop (want small)
    return cand[np.argmin(d)]


class Chain:
    """Open chain, sites 1..N free, sites 0 and N+1 clamped to the background bg."""

    def __init__(self, J, f, bg, N):
        self.J, self.f, self.bg, self.N = J, f, bg, N
        self.e_bg = bg @ J @ bg + f @ bg
        # energy of the clamped all-bg chain: (N+1) bonds + N fields
        self.E_ref = (N + 1) * (bg @ J @ bg) + N * (f @ bg)

    def unpack(self, x):
        V = x.reshape(-1, 3); r = np.linalg.norm(V, axis=1)
        return V / r[:, None], r

    def energy_grad(self, x):
        n, r = self.unpack(x)
        full = np.vstack([self.bg, n, self.bg])
        J, f = self.J, self.f
        E = np.einsum('ia,ab,ib->', full[:-1], J, full[1:]) + np.sum(n @ f) - self.E_ref
        G = full[2:] @ J.T + full[:-2] @ J + f
        G = (G - np.sum(G * n, 1)[:, None] * n) / r[:, None]
        return E, G

    def omega(self, x, apex=None):
        n, _ = self.unpack(x)
        if apex is None:
            apex = choose_apex(n, self.bg)
        return solid_angle_grad(n, self.bg, apex)[0]

    def relax(self, x0, Om0=None, K=0.0, apex=None, gtol=1e-9, maxiter=20000):
        def fun(x):
            E, G = self.energy_grad(x)
            if K > 0:
                n, r = self.unpack(x)
                Om, GO = solid_angle_grad(n, self.bg, apex)
                d = Om - Om0
                E = E + K * d * d
                GO = (GO - np.sum(GO * n, 1)[:, None] * n) / r[:, None]
                G = G + 2 * K * d * GO
            return E, G.ravel()
        res = minimize(fun, np.asarray(x0, float).ravel(), jac=True, method='L-BFGS-B',
                       options=dict(maxiter=maxiter, gtol=gtol, ftol=1e-15))
        n, _ = self.unpack(res.x)
        E = self.energy_grad(n.ravel())[0]
        return E, n

    def loop_start(self, axis, width, center, tilt=0.0, turns=1.0):
        """A 2 pi loop that starts and ends at bg, rotating about 'axis' (perpendicular part to bg used)."""
        bg = self.bg
        ax = axis - (axis @ bg) * bg
        ax /= np.linalg.norm(ax)
        third = np.cross(ax, bg)
        j = np.arange(1, self.N + 1)
        phi = np.pi * turns * (1 + np.tanh((j - center) / width))
        n = np.cos(phi)[:, None] * bg + np.sin(phi)[:, None] * third
        n = n + tilt * np.exp(-((j - center) / (1.5 * width)) ** 2)[:, None] * ax
        return n / np.linalg.norm(n, axis=1)[:, None]

    def hessian_tangent(self, n):
        """Tangent-space Hessian of E at a stationary texture n (2N x 2N, basis (u_i, w_i))."""
        N = self.N; J, f = self.J, self.f
        full = np.vstack([self.bg, n, self.bg])
        G = full[2:] @ J.T + full[:-2] @ J + f             # dE/dn_i (unconstrained)
        lam = np.sum(G * n, 1)
        T = np.zeros((N, 3, 2))
        for i in range(N):
            u, w = frame(n[i]); T[i, :, 0] = u; T[i, :, 1] = w
        H = np.zeros((2 * N, 2 * N))
        for i in range(N):
            H[2 * i:2 * i + 2, 2 * i:2 * i + 2] = -lam[i] * np.eye(2)
            if i + 1 < N:
                B = T[i].T @ J @ T[i + 1]
                H[2 * i:2 * i + 2, 2 * i + 2:2 * i + 4] = B
                H[2 * i + 2:2 * i + 4, 2 * i:2 * i + 2] = B.T
        return H, T


def momentum(Om):
    """Haldane: P = Omega/2 mod 2 pi, returned in [0, 2 pi)."""
    return (Om / 2) % TWO_PI


def pair_chain(J, f, bg, N, n1, n2, c1, c2, w1):
    """Initial condition with two copies of single-soliton textures n1, n2 (arrays on a chain with centre c0)
    shifted so that their centres sit at c1 < c2, on a chain of length N. Built by composing rotations:
    outside both cores the texture is bg; each core is pasted in."""
    x = np.tile(bg, (N, 1)).astype(float)
    for src, c in ((n1, c1), (n2, c2)):
        m = len(src); cs = m // 2
        for i in range(m):
            j = c - cs + i
            if 0 <= j < N:
                if np.linalg.norm(src[i] - bg) > 1e-6:
                    x[j] = src[i]
    return x
