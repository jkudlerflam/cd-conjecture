"""Exact-diagonalization checks of the modular-route identities on interacting chains
(Heisenberg, XXZ Delta = 0.5, SU(3) swap chain).

For each chain:
  1. open intervals of ell sites with the quadratic weight w_x = (x+1)(ell-1-x)/ell on the ell-1 bonds:
     E_0(K_w), the defect e_0 W - E_0(K_w), its slope in log ell (CFT: A/(12 pi)),
     and the layer-cake lower bound sum_k Delta_k E_0^open(n_k) (an exact inequality).
  2. open chains: f_b(n) = E_0^open(n bonds) - n e_0 (CFT: f_b - pi A/(24 n), f_b < 0).
  3. rings of L sites, reduced state rho_I of the ground state on ell = L/2 consecutive sites:
     S(rho_I), the top Schmidt vector chi and its K_w energy, and the relative entropy
     S(rho_I || tau_beta) with tau_beta = exp(-beta K_w)/Z for beta = 2 pi n/v, n = 1/2, 1, 2
     (CFT slope in log ell: beta A/(12 pi) + pi B/(3 beta) - C/3, which vanishes at beta = 2 pi/v
     for a single velocity). On a ring e_0 is replaced by E_0(L)/L.
Charge-sector restricted sparse ED (S^z for d = 2, colour counts for d = 3).
"""
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spla, itertools, sys
from scipy.linalg import eigh

sx = np.array([[0, 1], [1, 0]]) / 2; sy = np.array([[0, -1j], [1j, 0]]) / 2; sz = np.diag([.5, -.5])
def xxz(D): return (np.kron(sx, sx) + np.kron(sy, sy) + D * np.kron(sz, sz)).real
def swap3():
    h = np.zeros((9, 9))
    for a in range(3):
        for b in range(3): h[3 * a + b, 3 * b + a] = 1
    return h

def weights(ell):
    x = np.arange(ell - 1)
    return (x + 1) * (ell - 1 - x) / ell

class Sector:
    """Basis of d^L product states with prescribed colour counts, and sparse bond operators."""
    def __init__(self, d, L, counts):
        self.d, self.L = d, L
        states = []
        for s in itertools.product(range(d), repeat=L):
            if all(s.count(c) == counts[c] for c in range(d - 1)):
                states.append(s)
        self.S = np.array(states, dtype=np.int8)
        self.N = len(states)
        pw = d ** np.arange(L - 1, -1, -1)
        self.pw = pw
        self.codes = self.S @ pw
        order = np.argsort(self.codes)
        self.sorted_codes = self.codes[order]; self.order = order
    def index(self, codes):
        pos = np.searchsorted(self.sorted_codes, codes)
        return self.order[pos]
    def bond_op(self, h, i, j):
        d, S = self.d, self.S
        rows, cols, vals = [], [], []
        hb = h.reshape(d, d, d, d)   # hb[c,dd,a,b] = <c dd| h |a b>
        for a in range(d):
            for b in range(d):
                mask = (S[:, i] == a) & (S[:, j] == b)
                idx = np.nonzero(mask)[0]
                if len(idx) == 0: continue
                for c in range(d):
                    for dd in range(d):
                        v = hb[c, dd, a, b]
                        if abs(v) < 1e-14: continue
                        T = S[idx].copy(); T[:, i] = c; T[:, j] = dd
                        rows.append(self.index(T @ self.pw)); cols.append(idx); vals.append(np.full(len(idx), v))
        rows = np.concatenate(rows); cols = np.concatenate(cols); vals = np.concatenate(vals)
        return sp.csr_matrix((vals, (rows, cols)), shape=(self.N, self.N))
    def full_vector(self, v):
        out = np.zeros(self.d ** self.L); out[self.codes] = v; return out

def ground(Hs):
    if Hs.shape[0] <= 400:
        w, v = np.linalg.eigh(Hs.toarray()); return w[0], v[:, 0]
    w, v = spla.eigsh(Hs, k=1, which='SA', tol=1e-11); return w[0], v[:, 0]

def weighted_open(d, h, counts, wts):
    L = len(wts) + 1; sec = Sector(d, L, counts)
    H = sum(wts[x] * sec.bond_op(h, x, x + 1) for x in range(L - 1))
    return ground(H)[0]

def ring(d, h, counts, L):
    sec = Sector(d, L, counts)
    H = sum(sec.bond_op(h, x, (x + 1) % L) for x in range(L))
    E, v = ground(H)
    return E, sec.full_vector(v)

def dense_Kw(d, h, ell):
    wts = weights(ell); K = np.zeros((d ** ell, d ** ell))
    for x in range(ell - 1):
        K += wts[x] * np.kron(np.kron(np.eye(d ** x), h), np.eye(d ** (ell - x - 2)))
    return K, wts

def analyse(name, d, h, counts_of, A, v, c, e0_exact=None, ells_open=(), ns_bdy=(), rings=()):
    B, C = c / v, c
    print(f"\n==== {name}: d={d}, A={A:.5f}, v={v:.5f}, c={c}")
    # rings first, to get e_0 if not supplied
    ring_data = {}
    for L in rings:
        E, psi = ring(d, h, counts_of(L), L); ring_data[L] = (E, psi)
        print(f"  ring L={L:2d}  E0/L={E/L:.8f}  (E0/L + pi A/(6 L^2) = {E/L + np.pi*A/(6*L*L):.8f})")
    if e0_exact is None:
        Ls = sorted(ring_data); e0 = np.mean([ring_data[L][0] / L + np.pi * A / (6 * L * L) for L in Ls[-2:]])
        print(f"  e_0 estimated from the two largest rings: {e0:.8f}")
    else:
        e0 = e0_exact; print(f"  e_0 exact = {e0:.8f}")
    # 2. open-chain boundary energy
    for n in ns_bdy:
        E = weighted_open(d, h, counts_of(n + 1), np.ones(n))
        fb = E - n * e0
        print(f"  open n={n:2d} bonds: f_b(n)={fb:.6f}   f_b(n)+pi A/(24 n)={fb + np.pi*A/(24*n):.6f}")
    # 1. weighted intervals
    print(f"  weighted intervals, CFT slope A/(12 pi) = {A/(12*np.pi):.6f}")
    prev = None; Eopen = {}
    for ell in ells_open:
        wts = weights(ell); W = wts.sum()
        E0 = weighted_open(d, h, counts_of(ell), wts)
        defect = e0 * W - E0
        # layer cake
        lc = 0.0; prev_t = 0.0
        for t in np.unique(wts):
            n_k = int((wts >= t).sum())
            if n_k not in Eopen: Eopen[n_k] = weighted_open(d, h, counts_of(n_k + 1), np.ones(n_k))
            lc += (t - prev_t) * Eopen[n_k]; prev_t = t
        line = f"  ell={ell:2d}  defect={defect:.6f}  layer-cake: E0-lc={E0-lc:.5f} (must be >= 0)"
        if prev is not None:
            line += f"   slope vs ell={prev[0]}: {(defect-prev[1])/np.log(ell/prev[0]):.6f}"
        print(line); prev = (ell, defect)
    # 3. rings reduced to half
    for L in sorted(ring_data):
        E, psi = ring_data[L]; ell = L // 2; eL = E / L
        rho = psi.reshape(d ** ell, -1); rho = rho @ rho.T
        p, U = np.linalg.eigh(rho); p = np.clip(p, 1e-300, 1)
        S = -(p * np.log(p)).sum()
        chi = U[:, -1]
        K, wts = dense_Kw(d, h, ell); W = wts.sum()
        ek, Uk = eigh(K)
        E0K = ek[0]
        Erho = np.trace(rho @ K); Echi = chi @ K @ chi
        out = f"  ring L={L:2d} -> ell={ell}: S(rho_I)={S:.5f}  -log p_max={-np.log(p[-1]):.5f} (CFT S/2={S/2:.5f})  defect(eL)={eL*W-E0K:.5f}  chi-defect={eL*W-Echi:.5f}  rho-defect={eL*W-Erho:.5f}"
        rel = []
        for nn in (0.5, 1.0, 2.0):
            beta = 2 * np.pi * nn / v
            logZ = np.logaddexp.reduce(-beta * ek)
            rel.append((nn, -S + beta * Erho + logZ))
        out += "   S(rho||tau_beta) at beta v/(2pi)=1/2,1,2: " + ", ".join(f"{r:.5f}" for _, r in rel)
        print(out)
    print("  CFT slopes of S(rho||tau) in log ell at beta v/(2pi)=n: ", [round(2*np.pi*nn/v*A/(12*np.pi)+np.pi*B/(3*2*np.pi*nn/v)-C/3, 5) for nn in (0.5, 1.0, 2.0)])

if __name__ == '__main__':
    which = sys.argv[1] if len(sys.argv) > 1 else 'all'
    half = lambda L: (L // 2, L - L // 2)
    third = lambda L: (L // 3 + (L % 3 > 0), L // 3 + (L % 3 > 1), L // 3)   # most balanced colour counts
    if which in ('all', 'heis'):
        analyse("Heisenberg", 2, xxz(1.0), half, np.pi / 2, np.pi / 2, 1, e0_exact=0.25 - np.log(2),
                ells_open=range(4, 19), ns_bdy=(5, 9, 13, 17), rings=(10, 12, 14, 16, 18))
    if which in ('all', 'xxz'):
        analyse("XXZ Delta=0.5", 2, xxz(0.5), half, 3 * np.sqrt(3) / 4, 3 * np.sqrt(3) / 4, 1, e0_exact=-3 / 8,
                ells_open=range(4, 19), ns_bdy=(5, 9, 13, 17), rings=(10, 12, 14, 16, 18))
    if which in ('all', 'su3'):
        analyse("SU(3) swap", 3, swap3(), third, 4 * np.pi / 3, 2 * np.pi / 3, 2, e0_exact=None,
                ells_open=range(3, 14), ns_bdy=(5, 8, 11), rings=(9, 12, 15))
