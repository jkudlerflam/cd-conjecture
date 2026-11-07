"""Core library for the crossed-channel Fock analysis (contrib/cd_crossed_fock).

Conventions follow codex/work_C125 (PROOF section 2.1):
  column = periodic history s = (s_0, ..., s_{2m-1}) of one site, slice spacing beta/(2m);
  K_e(s,t) = prod_{j even} g(s_{j+1}, t_{j+1}; s_j, t_j),  K_o the same with j odd,
  g = exp(-delta h'), delta = beta/m, h' the sublattice-rotated bond,
  Z_m(L,beta) = Tr (K_e K_o)^{L/2} for even L.
K_e = G^{(x)m} on the pairs (0,1),(2,3),..., with G the realigned gate,
G[(s_j,s_{j+1}),(t_j,t_{j+1})] = g[(s_{j+1},t_{j+1}),(s_j,t_j)];  K_o is the same on (1,2),...,(2m-1,0).
U_sym = K_e^{1/2} K_o K_e^{1/2} >= 0 has eigenvalues u_j = exp(-2 eps_j), so Z = sum_j exp(-L eps_j).
Crossed momentum: rotation by two slices (one Trotter step), eigenvalue exp(2 pi i n/m), p = 2 pi n/beta.
Charge: Q(s) = sum_{j even} (s_j - s_{j+1}) for d = 2 XXZ (winding of the Sz world lines), conserved.
Only numpy and scipy are used.
"""
import numpy as np
import scipy.sparse as sp
from scipy.linalg import expm, eigh

sx = np.array([[0., 1.], [1., 0.]]) / 2
sy = np.array([[0., -1j], [1j, 0.]]) / 2
sz = np.diag([.5, -.5])


def bond_xxz(delta):
    """Original bond h = SxSx + SySy + Delta SzSz (real)."""
    return (np.kron(sx, sx) + np.kron(sy, sy) + delta * np.kron(sz, sz)).real


def bond_xxz_rot(delta):
    """Bond after a pi rotation about y on the second site: h' = -SxSx + SySy - Delta SzSz.
    -h' = Sx (x) conj Sx + Sy (x) conj Sy + Delta Sz (x) conj Sz, reflection positive for Delta >= 0."""
    return (-np.kron(sx, sx) + np.kron(sy, sy) - delta * np.kron(sz, sz)).real


def bond_positive(Jx=0., Jy=0., Jz=0., gx=0., gz=0.):
    """General d=2 bond in the reflection-positive frame (spin units S = sigma/2):
    h' = -Jx SxSx + Jy SySy - Jz SzSz - (gx Sx + gz Sz)(x)1/2 - 1(x)(gx Sx + gz Sz)/2,
    -h' = sum_a J_a S^a (x) conj(S^a) + B(x)1 + 1(x)conj(B), B = (gx Sx + gz Sz)/2 real, J_a >= 0.
    XXZ(Delta) at zero field in the original frame is (Jx,Jy,Jz) = (1,1,Delta); a uniform gz here is a
    staggered field there; (Jx, 0, 0, 0, gz) with Sx -> sigma/2 is the transverse-field Ising chain."""
    I = np.eye(2)
    B = gx * sx + gz * sz
    h = -Jx * np.kron(sx, sx) + Jy * np.kron(sy, sy) - Jz * np.kron(sz, sz) - 0.5 * (np.kron(B, I) + np.kron(I, B))
    return h.real


def gauge_M(h, d=2, tries=8, seed=0):
    """M = min over scalar and o(x)1 - 1(x)o of (lambda_max - lambda_min)/2, by Nelder-Mead over Hermitian o."""
    from scipy.optimize import minimize
    rng = np.random.default_rng(seed)
    basis = []
    for a in range(d):
        for b in range(d):
            E = np.zeros((d, d), complex)
            if a == b:
                E[a, a] = 1
            elif a < b:
                E[a, b] = E[b, a] = 1
            else:
                E[a, b] = 1j; E[b, a] = -1j
            basis.append(E)
    I = np.eye(d)
    def spread(x):
        o = sum(c * E for c, E in zip(x, basis))
        ev = np.linalg.eigvalsh(h + np.kron(o, I) - np.kron(I, o))
        return (ev[-1] - ev[0]) / 2
    best = spread(np.zeros(len(basis)))
    for t in range(tries):
        x0 = rng.normal(size=len(basis)) * (0.3 if t else 0)
        r = minimize(spread, x0, method='Nelder-Mead', options=dict(xatol=1e-10, fatol=1e-12, maxiter=20000))
        best = min(best, r.fun)
    return best


def half_spread(h):
    """M for XXZ-type bonds: gauge-optimized half spread. For the bonds used here
    (SU(2) or U(1)xZ2 symmetric) the optimal one-site gauge term vanishes and the scalar shift
    centres the spectrum, so M = (lambda_max - lambda_min)/2 (as in C125's checks)."""
    ev = np.linalg.eigvalsh(h)
    return (ev[-1] - ev[0]) / 2


def realign(g, d):
    g4 = g.reshape(d, d, d, d)          # g4[a,b,c,e] = <a b| g |c e>, (a,b) at slice j+1
    G4 = np.einsum('qupr->pqru', g4)    # G4[s_j,s_{j+1},t_j,t_{j+1}] = g4[s_{j+1},t_{j+1},s_j,t_j]
    return G4.reshape(d * d, d * d)


def psd_sqrt(G):
    w, v = eigh((G + G.T) / 2)
    if w.min() < -1e-12 * max(1., abs(w).max()):
        raise ValueError('realigned gate not PSD, min eig %g' % w.min())
    w = np.clip(w, 0, None)
    return (v * np.sqrt(w)) @ v.T, w


class Column:
    """Column (history) space of 2m slices of a d-level site, with integer basis index
    s = sum_j s_j d^(2m-1-j) (C order, axis j = slice j)."""

    def __init__(self, d, m):
        self.d, self.m, self.n = d, m, 2 * m
        self.N = d ** (2 * m)
        idx = np.arange(self.N, dtype=np.int64)
        self.digits = np.array([(idx // d ** (2 * m - 1 - j)) % d for j in range(2 * m)])  # (2m, N)
        # permutation for rolling by one slice: (P psi)[s] = psi[s shifted]
        # roll1[s] = index of the string (s_1, s_2, ..., s_{2m-1}, s_0)
        pw = d ** (2 * m - 1 - np.arange(2 * m))
        self.roll1 = (np.roll(self.digits, -1, axis=0) * pw[:, None]).sum(0)
        self.roll2 = (np.roll(self.digits, -2, axis=0) * pw[:, None]).sum(0)

    def charge(self, mode='u1'):
        """Conserved column charge (d=2).
        'u1': winding Q = sum_{j even}(s_j - s_{j+1}) in Sz units (positive-frame XXZ, which conserves
              the staggered magnetization, so each gate obeys s_j - s_{j+1} = t_j - t_{j+1} at even j);
        'z2': parity of the number of flips at even slices (bonds with only pair flips, e.g. XYZ, TFIM);
        'none': no charge."""
        if mode == 'none':
            return np.zeros(self.N, int)
        if mode == 'z2':
            return ((self.digits[0::2] != self.digits[1::2]).sum(0) % 2).astype(int)
        sZ = 0.5 - self.digits
        return np.rint(sZ[0::2].sum(0) - sZ[1::2].sum(0)).astype(int)


def kron_power(A, m):
    out = sp.csr_matrix(A)
    Asp = sp.csr_matrix(A)
    for _ in range(m - 1):
        out = sp.kron(out, Asp, format='csr')
    return out


class Transfer:
    """Checkerboard transfer for bond h (already rotated if needed), inverse temperature beta,
    Trotter number m (2m slices)."""

    def __init__(self, h, d, beta, m, col=None):
        self.d, self.m, self.beta = d, m, beta
        self.delta = beta / m
        self.g = expm(-self.delta * h)
        self.G = realign(self.g, d)
        self.Gh, self.Gw = psd_sqrt(self.G)
        self.col = col if col is not None else Column(d, m)
        c = self.col
        Ke_half = kron_power(self.Gh, m)          # acts on pairs (0,1),(2,3),...
        Ke = kron_power(self.G, m)
        # K_o = P^T K_e P with (P psi)[s] = psi[roll1[s]]: pairs (1,2),...,(2m-1,0)
        P = sp.csr_matrix((np.ones(c.N), (np.arange(c.N), c.roll1)), shape=(c.N, c.N))
        self.P = P
        self.Ke_half = Ke_half.tocsr()
        self.Ke = Ke.tocsr()
        self.Ko = (P.T @ Ke @ P).tocsr()
        self.Ko_half = (P.T @ Ke_half @ P).tocsr()

    def U_matvec(self, x):
        return self.Ke_half @ (self.Ko @ (self.Ke_half @ x))

    def sector(self, idx):
        """restrict the three factors to a set of basis indices (a conserved-charge sector)"""
        Kh = self.Ke_half[idx][:, idx]
        Ko = self.Ko[idx][:, idx]
        return Kh.tocsr(), Ko.tocsr()


def momentum_basis(col, idx, n):
    """Orbit basis of rotation by two slices restricted to basis indices idx (closed under rotation).
    Returns sparse matrix B (len(idx) x D) with orthonormal columns |r,n>,
    |r,n> = p_r^{-1/2} sum_{k<p_r} exp(-2 pi i n k/m) |R^k r>, R = rotation by two slices."""
    m = col.m
    pos = -np.ones(col.N, dtype=np.int64)
    pos[idx] = np.arange(len(idx))
    seen = np.zeros(len(idx), bool)
    rows, cols, vals = [], [], []
    D = 0
    r2 = col.roll2
    for a in range(len(idx)):
        if seen[a]:
            continue
        orbit = [idx[a]]
        s = r2[idx[a]]
        while s != idx[a]:
            orbit.append(s)
            s = r2[s]
        p = len(orbit)
        for s in orbit:
            seen[pos[s]] = True
        if (n * p) % m != 0:
            continue
        ph = np.exp(-2j * np.pi * n * np.arange(p) / m) / np.sqrt(p)
        rows.extend(pos[np.array(orbit)])
        cols.extend([D] * p)
        vals.extend(ph)
        D += 1
    return sp.csr_matrix((vals, (rows, cols)), shape=(len(idx), D))


def sector_spectrum(T, idx, n, chunk=400, vectors=False):
    """Dense diagonalization of U_sym in the sector (charge set idx, momentum n).
    Returns eigenvalues u (descending) and optionally eigenvectors in the sector basis plus B."""
    Kh, Ko = T.sector(idx)
    B = momentum_basis(T.col, idx, n)
    D = B.shape[1]
    if D == 0:
        return (np.zeros(0), None, B) if vectors else np.zeros(0)
    Mt = np.zeros((D, D), complex)
    BH = B.conj().T.tocsr()
    for c0 in range(0, D, chunk):
        X = B[:, c0:c0 + chunk].toarray()
        Y = Kh @ (Ko @ (Kh @ X))
        Mt[:, c0:c0 + chunk] = BH @ Y
    Mt = (Mt + Mt.conj().T) / 2
    if vectors:
        u, V = eigh(Mt)
        o = np.argsort(-u)
        return u[o], V[:, o], B
    u = np.linalg.eigvalsh(Mt)
    return np.sort(u)[::-1]


# ---------------- direct channel -----------------

def ring_bond_ops(h, d, L):
    """List of dense two-site operators on an L-ring for bond h, as sparse matrices."""
    ops = []
    I = sp.identity
    for x in range(L):
        y = (x + 1) % L
        if y == x + 1:
            op = sp.kron(sp.kron(I(d ** x), sp.csr_matrix(h)), I(d ** (L - x - 2)))
        else:  # bond (L-1, 0): permute
            # build via explicit basis action
            N = d ** L
            idx = np.arange(N)
            dig = np.array([(idx // d ** (L - 1 - j)) % d for j in range(L)])
            hr = h.reshape(d, d, d, d)
            rows, cols, vals = [], [], []
            for a in range(d):
                for b in range(d):
                    for c in range(d):
                        for e in range(d):
                            v = hr[a, b, c, e]   # <a b|h|c e> on sites (L-1, 0)
                            if v == 0:
                                continue
                            mask = (dig[L - 1] == c) & (dig[0] == e)
                            src = idx[mask]
                            tgt = src + (a - c) * 1 + (b - e) * d ** (L - 1)
                            rows.append(tgt); cols.append(src); vals.append(np.full(len(src), v))
            op = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(N, N))
        ops.append(op.tocsr())
    return ops


def sz_sectors(L):
    N = 2 ** L
    idx = np.arange(N)
    pop = np.array([bin(i).count('1') for i in range(N)])
    return [idx[pop == k] for k in range(L + 1)]


def direct_logZ(h, L, beta, m=None, sectors=None):
    """log Z for the even L-ring (d=2, Sz-conserving bond h).
    m=None: exact log Tr exp(-beta H_L). Otherwise the checkerboard Trotter trace
    Tr[(exp(-delta H_e) exp(-delta H_o))^m], delta = beta/m, H_e = bonds (x,x+1) with x even."""
    from scipy.special import logsumexp
    ops = ring_bond_ops(h, 2, L)
    He = sum(ops[x] for x in range(0, L, 2))
    Ho = sum(ops[x] for x in range(1, L, 2))
    sectors = sectors if sectors is not None else sz_sectors(L)
    betas = np.atleast_1d(beta)
    out = []
    logs = [[] for _ in betas]
    for idx in sectors:
        A = He[idx][:, idx].toarray(); B = Ho[idx][:, idx].toarray()
        if m is None:
            ev = np.linalg.eigvalsh(A + B)
            for i, b in enumerate(betas):
                logs[i].append(logsumexp(-b * ev))
        else:
            wa, va = np.linalg.eigh(A); wb, vb = np.linalg.eigh(B)
            for i, b in enumerate(betas):
                dl = b / m
                Eh = (va * np.exp(-dl * wa / 2)) @ va.T
                Eo = (vb * np.exp(-dl * wb)) @ vb.T
                Vs = Eh @ Eo @ Eh
                ev = np.linalg.eigvalsh((Vs + Vs.T) / 2)
                ev = np.clip(ev, 1e-300, None)
                logs[i].append(logsumexp(m * np.log(ev)))
    res = np.array([logsumexp(l) for l in logs])
    return res if np.ndim(beta) else res[0]


def direct_logZ_dense(h, L, beta, m=None):
    """generic (no symmetry) version of direct_logZ for small even L"""
    from scipy.special import logsumexp
    ops = ring_bond_ops(h, 2, L)
    He = sum(ops[x] for x in range(0, L, 2)).toarray()
    Ho = sum(ops[x] for x in range(1, L, 2)).toarray()
    if m is None:
        return logsumexp(-beta * np.linalg.eigvalsh(He + Ho))
    wa, va = np.linalg.eigh(He); wb, vb = np.linalg.eigh(Ho); dl = beta / m
    Eh = (va * np.exp(-dl * wa / 2)) @ va.T; Eo = (vb * np.exp(-dl * wb)) @ vb.T
    ev = np.clip(np.linalg.eigvalsh(Eh @ Eo @ Eh), 1e-300, None)
    return logsumexp(m * np.log(ev))
