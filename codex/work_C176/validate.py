"""C176 validation of fibdmrg.py against C173's ED (eff_ed.py, copied unmodified) at L <= 24."""
import numpy as np, sys, time
from fibdmrg import Chain, DMRG, ground
from eff_ed import build, entropies, basis_fixed
import scipy.sparse.linalg as spla

def ed(L, N, w, **kw):
    b, H = build(L, N, w, **kw)
    if H.shape[0] < 500:
        e, v = np.linalg.eigh(H.toarray()); return b, e[0], v[:, 0], e[1]
    e, v = spla.eigsh(H, k=2, which='SA', tol=1e-12); o = np.argsort(e)
    return b, e[o[0]], v[:, o[0]], e[o[1]]

def ed_obs(b, psi, L):
    dens = np.array([np.sum(np.abs(psi) ** 2 * ((b >> j) & 1)) for j in range(L)])
    return dens

if __name__ != "__main__": cases = []
else: cases = [
    # (L, N, |w|, phi, kw for Chain, kw for ED, gauge)
    (12, 3, 1.0, np.pi / 2, {}, {}, 1),
    (12, 3, 1.0, np.pi / 2, {}, {}, 0),
    (16, 4, 1.0, 1.3, {}, {}, 0),
    (18, 3, 1.0, 1.52, {}, {}, 1),
    (20, 4, 1.0, np.pi / 2, {}, {}, 1),
    (24, 4, 1.0, np.pi / 2, {}, {}, 1),
    (24, 5, 0.7, 0.4, {}, {}, 0),
    (22, 9, 1.0, 1.3, {}, {}, 0),            # dense, near close packing
    (20, 5, 1.0, 0.0, dict(t2=0.0, V2=0.0, V3=0.0), dict(t2=0.0, V2=0.0, V3=0.0), 0),   # c=1 control
    (16, 5, 0.0, 0.0, dict(V2=0.0, V3=0.0, blockade=False), dict(V2=0.0, V3=0.0, blockade=False), 0),  # c=2 control
]
if __name__ == "__main__":
    worst = 0
    import builtins
    for (L, N, wa, phi, ckw, ekw, g) in cases:
        w = wa * np.exp(1j * phi)
        t0 = time.time()
        b, e0, psi, e1 = ed(L, N, w, **ekw)
        Sed = entropies(b, psi, L)
        dens_ed = ed_obs(b, psi, L)
        dm, hist = ground(L, N, w, chis=[16, 32, 64, 128, 128], ch_kw=dict(gauge=g, **ckw), tol=1e-12)
        E = hist[-1][1]
        m = dm.measure(gwin=(0, L), refs=[L // 2])
        dS = np.max(np.abs(dm.S - Sed))
        dd = np.max(np.abs(np.array(m['dens']) - dens_ed))
        # one-body density matrix (fermion, c frame) vs ED
        from ring_ed import popcount
        G = np.array(m['Gf_re']) + 1j * np.array(m['Gf_im'])
        idx = {int(s): k for k, s in enumerate(b)}
        Ged = np.zeros((L, L), complex)
        for i in range(L):
            for j in range(L):
                if i == j:
                    Ged[i, i] = dens_ed[i]; continue
                occ = np.nonzero(((b >> j) & 1) & (1 - ((b >> i) & 1)))[0]
                if len(occ) == 0: continue
                ts = (b[occ] & ~(1 << j)) | (1 << i)
                pos = np.searchsorted(b, ts); pos = np.minimum(pos, len(b) - 1); good = b[pos] == ts
                lo, hi = min(i, j), max(i, j)
                mask = ((1 << hi) - 1) & ~((1 << (lo + 1)) - 1)
                sg = (-1.0) ** popcount(b[occ] & mask)
                Ged[i, j] = np.sum(np.conj(psi[pos[good]]) * psi[occ[good]] * sg[good])
        dG = np.max(np.abs(G - Ged))
        # P_j = c+_{j+1}c_j + h.c. expectation vs ED
        Ped = np.array([2 * np.real(Ged[j + 1, j]) for j in range(L - 1)])
        dP = np.max(np.abs(np.array(m['P_exp']) - Ped))
        gap = e1 - e0
        print(f"L={L} N={N} |w|={wa} phi={phi:.3f} gauge={g} {ckw}: E_dmrg-E_ed={E-e0:+.2e} gap={gap:.3e} max|dS|={dS:.1e} "
              f"max|dn|={dd:.1e} max|dG|={dG:.1e} max|dP|={dP:.1e} nn_adj={m['nn_adjacent_max']:.1e} trunc={hist[-1][2]:.1e} ({time.time()-t0:.0f}s)", flush=True)
        worst = max(worst, abs(E - e0))
    print("worst |dE| =", worst)
