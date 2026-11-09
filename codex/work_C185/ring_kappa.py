"""C185: ring ED of the transverse-DM cusp effective model with its dN = 1 term (no N conservation), to see the
valley content of the ground state against phi. Model (units |t2| = 1), blockaded ring of L sites:
  H = sum_j [ w b+_{j+1} b_j + t2 b+_{j+2} b_j + h.c. + V2 n_j n_{j+2} + V3 n_j n_{j+3} - mu n_j ]
      + kappa sum_s [ b+_s (n_{s+2} - n_{s-2}) + h.c. ],  w = u e^{i phi}, t2 = -1, V2 = -2 + v, V3 = -2/3.
Particles never pass (blockade), so hard-core bosons; n(k) from <b+_i b_j>. Valley weights: W_a = sum of n(k) over
cos k > 0, W_b over cos k < 0 (valleys at k_a = -asin(u/4), k_b = pi + asin(u/4) when phi = pi/2).
usage: python3 ring_kappa.py L u v kappa mu phi1,phi2,...
"""
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spla, sys

def basis(L):
    out = []
    for s in range(1 << L):
        if s & ((s << 1) | (s >> (L - 1))) & ((1 << L) - 1):
            continue
        out.append(s)
    return np.array(out, dtype=np.int64)

def build(L, u, phi, v, kappa, mu, B, idx):
    w = u * np.exp(1j * phi); t2 = -1.0; V2 = -2 + v; V3 = -2 / 3
    rows, cols, vals = [], [], []
    occ = ((B[:, None] >> np.arange(L)[None, :]) & 1).astype(np.int8)
    diag = -mu * occ.sum(1).astype(float)
    for j in range(L):
        diag += V2 * occ[:, j] * occ[:, (j + 2) % L] + V3 * occ[:, j] * occ[:, (j + 3) % L]
    rows.append(np.arange(len(B))); cols.append(np.arange(len(B))); vals.append(diag.astype(complex))
    def add(src_mask_ok, newB, amp):
        sel = np.nonzero(src_mask_ok)[0]
        tgt = newB[sel]
        ok = np.isin(tgt, B)
        sel = sel[ok]; tgt = tgt[ok]
        rows.append(idx(tgt)); cols.append(sel); vals.append(np.full(len(sel), amp, complex))
    for j in range(L):
        for dist, amp in ((1, w), (2, t2)):
            k = (j + dist) % L
            # j -> k with amplitude amp (b+_k b_j), and k -> j with conj
            m = (occ[:, j] == 1) & (occ[:, k] == 0)
            add(m, B ^ (1 << j) ^ (1 << k), amp)
            m = (occ[:, k] == 1) & (occ[:, j] == 0)
            add(m, B ^ (1 << j) ^ (1 << k), np.conj(amp))
        if kappa != 0:
            for c, sg in ((2, 1.0), (-2, -1.0)):
                q = (j + c) % L
                m = (occ[:, j] == 0) & (occ[:, q] == 1)
                add(m, B ^ (1 << j), sg * kappa)            # creation at j (blockade checked by membership)
                m = (occ[:, j] == 1) & (occ[:, q] == 1)
                add(m, B ^ (1 << j), sg * np.conj(kappa))   # annihilation at j
    H = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(len(B), len(B)))
    return H, occ

def nk(L, psi, B, idx, occ):
    G = np.zeros((L, L), complex)
    for i in range(L):
        G[i, i] = np.vdot(psi, occ[:, i] * psi)
        for j in range(L):
            if i == j: continue
            m = (occ[:, j] == 1) & (occ[:, i] == 0)
            sel = np.nonzero(m)[0]
            nb = B[sel] ^ (1 << i) ^ (1 << j)
            ok = np.isin(nb, B); sel = sel[ok]; nb = nb[ok]
            G[i, j] = np.vdot(psi[idx(nb)], psi[sel])   # <b+_i b_j>
    ks = 2 * np.pi * np.arange(L) / L
    n = np.array([np.real(np.sum(G * np.exp(1j * k * (np.arange(L)[:, None] - np.arange(L)[None, :])))) / L for k in ks])
    return ks, n, np.real(np.trace(G))

if __name__ == "__main__":
    L, u, v, kappa, mu = int(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5])
    phis = [float(x) for x in sys.argv[6].split(',')]
    B = basis(L); order = np.argsort(B); Bs = B[order]
    idx = lambda t: np.searchsorted(Bs, t)
    B = Bs
    print(f"# L={L} u={u} v={v} kappa={kappa} mu={mu} dim={len(B)}")
    print("# phi   E0   <N>/L   W_a   W_b   W_a/W_b   k_peak_a  k_peak_b   gap(E1-E0)")
    for phi in phis:
        H, occ = build(L, u, phi, v, kappa, mu, B, idx)
        e, V = spla.eigsh(H, k=2, which='SA', tol=1e-10)
        o = np.argsort(e); e = e[o]; psi = V[:, o[0]]
        ks, n, N = nk(L, psi, B, idx, occ)
        a = np.cos(ks) > 1e-9; b = np.cos(ks) < -1e-9
        Wa, Wb = n[a].sum(), n[b].sum()
        ka = ks[a][np.argmax(n[a])]; kb = ks[b][np.argmax(n[b])]
        print(f"{phi:.4f} {e[0]:.6f} {N/L:.4f} {Wa:.4f} {Wb:.4f} {Wa/max(Wb,1e-12):.3f} {ka:.3f} {kb:.3f} {e[1]-e[0]:.5f}", flush=True)
