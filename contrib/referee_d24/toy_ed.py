"""Exact diagonalization of a small analogue of the vacancy completion.

Toy model (local dimension d = 1 + 2b).  A site is either a vacancy or an
occupied letter (c, s) with crystal colour c in Z_b and a spin-1/2 s.
The bond is

    h = -t sum_a (|0,a><a,0| + h.c.)
        + J P_occ [ K (1 - delta_{c', c+1 mod b})
                    + (Jxy/2)(S+S- + S-S+) + Jz SzSz
                    + (eps_c + eps_c')/2 ] P_occ ,

so, exactly as in the d=24 claim, vacancy hops are species blind and the
whole internal bond (penalty, exchange and its half on-site terms) acts only
on physically adjacent occupied pairs.  Colours never change, so the cyclic
colour word in squeezed coordinates is conserved: for b >= 2 the internal
chain is a frozen period-b crystal (N must be a multiple of b for a legal
word), with a critical XXZ payload (c_pay = 1) riding on it.  The colour
energies eps_c couple the charge to the period-b pattern, which produces the
fractional marker harmonics of the d=24 argument.  b = 1 is the plain
vacancy + XXZ chain (a t-J-type model with c = 2), used as calibration.

Basis: (occupied set, rotation r of a fixed necklace w, spin string of the
occupied sites in physical order starting at site 0).  A particle hopping
across the cut rotates the internal word; no fermionic sign appears because
the physical vacancy swap is bosonic.  Every basis identity is checked
against a brute-force product-basis construction for small L by
`selftest()`.

numpy and scipy only.
"""
import itertools
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla


def necklace_period(w):
    n = len(w)
    for p in range(1, n + 1):
        if n % p == 0 and all(w[(k + p) % n] == w[k] for k in range(n)):
            return p
    return n


class Toy:
    def __init__(self, L, N, w, n_up, b):
        self.L, self.N, self.w, self.b = L, N, list(w), b
        self.P = necklace_period(w)
        self.occs = np.array([sum(1 << x for x in c)
                              for c in itertools.combinations(range(L), N)], dtype=np.int64)
        self.spins = np.array([sum(1 << k for k in c)
                               for c in itertools.combinations(range(N), n_up)], dtype=np.int64)
        self.occ_rank = -np.ones(1 << L, dtype=np.int64)
        self.occ_rank[self.occs] = np.arange(len(self.occs))
        self.spin_rank = -np.ones(1 << N, dtype=np.int64)
        self.spin_rank[self.spins] = np.arange(len(self.spins))
        self.no, self.ns = len(self.occs), len(self.spins)
        self.dim = self.no * self.P * self.ns
        # positions of particles for each occupation
        bits = ((self.occs[:, None] >> np.arange(L)[None, :]) & 1).astype(bool)
        self.pos = np.array([np.nonzero(r)[0] for r in bits])            # (no, N)
        self.bits = bits
        # particle index at each site, -1 if empty
        self.at = -np.ones((self.no, L), dtype=np.int64)
        for k in range(N):
            self.at[np.arange(self.no), self.pos[:, k]] = k
        g = np.arange(self.dim)
        self.oi = g // (self.P * self.ns)
        self.r = (g // self.ns) % self.P
        self.si = g % self.ns

    def index(self, oi, r, si):
        return (oi * self.P + r) * self.ns + si

    def colour(self, r, k):
        w = np.array(self.w)
        return w[(k + r) % self.N]

    @staticmethod
    def rot_right(s, N):
        """Spin string after the last particle becomes the first."""
        last = (s >> (N - 1)) & 1
        return ((s << 1) & ((1 << N) - 1)) | last

    def hamiltonian(self, t, J, K, Jxy, Jz, eps, mu=0.0, theta=0.0, eta=None):
        """eta[c]: extra energy of an adjacent legal pair with colours (c, c+1)."""
        L, N, b = self.L, self.N, self.b
        rows, cols, vals = [], [], []
        diag = np.zeros(self.dim)
        occ = self.occs[self.oi]
        spin = self.spins[self.si]
        eps = np.asarray(eps, float)
        for x in range(L):
            y = (x + 1) % L
            ox = (occ >> x) & 1
            oy = (occ >> y) & 1
            # vacancy hop: particle x -> y (to the right), amplitude -t e^{i theta}
            m = (ox == 1) & (oy == 0)
            if m.any():
                g = np.nonzero(m)[0]
                no = occ[g] ^ (1 << x) ^ (1 << y)
                r2 = self.r[g].copy()
                s2 = spin[g].copy()
                if y == 0:      # crossed the cut: last particle becomes first
                    s2 = self.rot_right(s2, N)
                    r2 = (r2 - 1) % self.P
                tgt = self.index(self.occ_rank[no], r2, self.spin_rank[s2])
                amp = -t * np.exp(1j * theta)
                rows += [tgt, g]
                cols += [g, tgt]
                vals += [np.full(len(g), amp), np.full(len(g), np.conj(amp))]
            # internal bond on adjacent occupied pair (x, y)
            m = (ox == 1) & (oy == 1)
            if m.any():
                g = np.nonzero(m)[0]
                kx = self.at[self.oi[g], x]
                ky = self.at[self.oi[g], y]
                cx = self.colour(self.r[g], kx)
                cy = self.colour(self.r[g], ky)
                sx = (spin[g] >> kx) & 1
                sy = (spin[g] >> ky) & 1
                d = Jz * (sx - 0.5) * (sy - 0.5) + 0.5 * (eps[cx] + eps[cy])
                d = d + K * (cy != (cx + 1) % b)
                if eta is not None:
                    d = d + np.asarray(eta, float)[cx] * (cy == (cx + 1) % b)
                diag[g] += J * d
                f = sx != sy
                if f.any():
                    gg = g[f]
                    s2 = spin[gg] ^ (1 << kx[f]) ^ (1 << ky[f])
                    tgt = self.index(self.oi[gg], self.r[gg], self.spin_rank[s2])
                    rows.append(tgt)
                    cols.append(gg)
                    vals.append(np.full(len(gg), 0.5 * J * Jxy))
        diag -= mu * N
        rows.append(np.arange(self.dim))
        cols.append(np.arange(self.dim))
        vals.append(diag.astype(complex))
        H = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                          shape=(self.dim, self.dim))
        return H

    def translation(self):
        """Permutation T: every particle moves one site to the right."""
        L, N = self.L, self.N
        occ = self.occs[self.oi]
        spin = self.spins[self.si]
        wrap = ((occ >> (L - 1)) & 1) == 1
        no = ((occ << 1) & ((1 << L) - 1)) | (occ >> (L - 1))
        s2 = np.where(wrap, self.rot_right(spin, N), spin)
        r2 = np.where(wrap, (self.r - 1) % self.P, self.r)
        return self.index(self.occ_rank[no], r2, self.spin_rank[s2])

    def site_ops(self):
        """Diagonal n_x and S^z_x arrays, shape (L, dim)."""
        occ = self.occs[self.oi]
        spin = self.spins[self.si]
        n = np.array([(occ >> x) & 1 for x in range(self.L)], float)
        sz = np.zeros((self.L, self.dim))
        for x in range(self.L):
            k = self.at[self.oi, x]
            ok = k >= 0
            sz[x, ok] = ((spin[ok] >> k[ok]) & 1) - 0.5
        return n, sz

    def local_letter(self):
        """Letter index at each site: 0 vacancy, 1 + 2c + s otherwise. (dim, L)."""
        out = np.zeros((self.dim, self.L), dtype=np.int64)
        spin = self.spins[self.si]
        for x in range(self.L):
            k = self.at[self.oi, x]
            ok = k >= 0
            c = self.colour(self.r[ok], k[ok])
            s = (spin[ok] >> k[ok]) & 1
            out[ok, x] = 1 + 2 * c + s
        return out

    def entropy(self, psi, ell):
        """Von Neumann entropy of sites [0, ell) for state psi."""
        lett = self.local_letter()
        d = 1 + 2 * self.b
        left = np.zeros(self.dim, dtype=np.int64)
        right = np.zeros(self.dim, dtype=np.int64)
        for x in range(ell):
            left = left * d + lett[:, x]
        for x in range(ell, self.L):
            right = right * d + lett[:, x]
        ul, li = np.unique(left, return_inverse=True)
        ur, ri = np.unique(right, return_inverse=True)
        nz = np.abs(psi) > 0
        M = sp.csr_matrix((psi[nz], (li[nz], ri[nz])), shape=(len(ul), len(ur)))
        # block structure by left quantum numbers keeps this small
        rho = (M @ M.conj().T).tocsr()
        ncomp, lab = sp.csgraph.connected_components(rho != 0, directed=False)
        S = 0.0
        for cblock in range(ncomp):
            idx = np.nonzero(lab == cblock)[0]
            blk = rho[idx][:, idx].toarray()
            ev = np.linalg.eigvalsh(blk)
            ev = ev[ev > 1e-14]
            S -= float(np.sum(ev * np.log(ev)))
        return S


def lowest(H, k=6):
    if H.shape[0] <= 400:
        ev, vec = np.linalg.eigh(H.toarray())
        return ev[:k], vec[:, :k]
    ev, vec = spla.eigsh(H, k=k, which='SA', tol=1e-11, ncv=max(2 * k + 10, 30))
    o = np.argsort(ev)
    return ev[o], vec[:, o]


def momenta(Tperm, vecs, L):
    """Momenta of (possibly degenerate) eigenvectors: diagonalize T in the block."""
    Tv = np.zeros_like(vecs)
    Tv[Tperm] = vecs
    Tm = vecs.conj().T @ Tv
    ev = np.linalg.eigvals(Tm)
    return np.sort(np.round(np.angle(ev) * L / (2 * np.pi), 6))


def lanczos_pole(H, v0, steps=120, wthresh=0.05):
    """Lowest pole of the spectral function of v0 with weight above wthresh."""
    nrm = np.linalg.norm(v0)
    if nrm < 1e-12:
        return None, None
    v = v0 / nrm
    vprev = np.zeros_like(v)
    a, bb = [], []
    beta = 0.0
    V = [v]
    for _ in range(steps):
        w = H @ v - beta * vprev
        al = np.vdot(v, w).real
        w = w - al * v
        for u in V:      # full reorthogonalization
            w -= np.vdot(u, w) * u
        beta = np.linalg.norm(w)
        a.append(al)
        if beta < 1e-10:
            break
        bb.append(beta)
        vprev, v = v, w / beta
        V.append(v)
    T = np.diag(a) + np.diag(bb[:len(a) - 1], 1) + np.diag(bb[:len(a) - 1], -1)
    ev, U = np.linalg.eigh(T)
    wts = np.abs(U[0]) ** 2
    for e, wt in zip(ev, wts):
        if wt > wthresh:
            return e, wt
    return None, None


def selftest():
    """Compare with a brute-force product-basis construction (L=4..6; b=1,2,3)."""
    errs = []
    for (L, N, b, n_up) in [(5, 2, 2, 1), (5, 4, 2, 2), (5, 3, 2, 1), (6, 3, 1, 1), (6, 4, 1, 2), (4, 3, 3, 1)]:
        w = [k % b for k in range(N)]
        toy = Toy(L, N, w, n_up, b)
        prm = dict(t=1.0, J=0.7, K=3.0, Jxy=1.1, Jz=0.4, eps=np.linspace(-0.5, 0.5, b),
                   eta=np.linspace(0.3, -0.2, b))
        H = toy.hamiltonian(**prm).toarray()
        d = 1 + 2 * b
        # brute force on the full product space
        hb = np.zeros((d * d, d * d), complex)
        for a in range(1, d):
            hb[0 * d + a, a * d + 0] -= prm['t']
            hb[a * d + 0, 0 * d + a] -= prm['t']
        for a1 in range(1, d):
            for a2 in range(1, d):
                c1, s1 = divmod(a1 - 1, 2)
                c2, s2 = divmod(a2 - 1, 2)
                i = a1 * d + a2
                hb[i, i] += prm['J'] * (prm['Jz'] * (s1 - .5) * (s2 - .5)
                                        + .5 * (prm['eps'][c1] + prm['eps'][c2])
                                        + prm['K'] * (c2 != (c1 + 1) % b)
                                        + prm['eta'][c1] * (c2 == (c1 + 1) % b))
                if s1 != s2:
                    j = (1 + 2 * c1 + s2) * d + (1 + 2 * c2 + s1)
                    hb[j, i] += .5 * prm['J'] * prm['Jxy']
        Hf = np.zeros((d ** L, d ** L), complex)
        for x in range(L):
            y = (x + 1) % L
            # embed hb on (x, y)
            for idx in range(d ** L):
                digs = [(idx // d ** (L - 1 - z)) % d for z in range(L)]
                i = digs[x] * d + digs[y]
                col = hb[:, i]
                for j in np.nonzero(col)[0]:
                    dd = list(digs)
                    dd[x], dd[y] = divmod(j, d)
                    jdx = sum(v * d ** (L - 1 - z) for z, v in enumerate(dd))
                    Hf[jdx, idx] += col[j]
        # map toy basis into product basis
        lett = toy.local_letter()
        emb = np.zeros((d ** L, toy.dim))
        for g in range(toy.dim):
            jdx = sum(int(v) * d ** (L - 1 - z) for z, v in enumerate(lett[g]))
            emb[jdx, g] = 1
        Hp = emb.T @ Hf @ emb
        leak = np.abs((np.eye(d ** L) - emb @ emb.T) @ Hf @ emb).max()
        errs.append(dict(L=L, N=N, b=b, dim=toy.dim, matrix_error=float(np.abs(Hp - H).max()),
                         leakage=float(leak)))
        Tp = toy.translation()
        assert len(np.unique(Tp)) == toy.dim
    return errs


if __name__ == '__main__':
    import json
    print(json.dumps(selftest(), indent=2))
