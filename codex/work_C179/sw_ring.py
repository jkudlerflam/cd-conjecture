"""Exact Schrieffer-Wolff (des Cloizeaux) effective Hamiltonian of the full qubit RING (translation invariant, no edges),
L = 12, in the spin basis. Tests the bulk selection rules:
 - separation-2 pair creation from the vacuum, <b,b+2|H|0>: zero on Jx + Jy = 0 to all orders (U-odd uniform term);
   nonzero, = -2 Delta Re(w)/V at leading order, off it;
 - separation-3 pair creation <b,b+3|H|0> = +2 Delta^3/(3V^2) at leading order;
 - the largest Delta N = 4 element, predicted O(Delta^6/V^5)."""
import numpy as np
from sw_exact import op, X, Y, Z
L = 12
def ring_heff(Delta, V, w, mu):
    Jz = V / 4; hz = (mu - 4 * Jz) / 2; Jx = (w.real + Delta) / 2; Jy = (w.real - Delta) / 2; D = w.imag / 2
    H = np.zeros((2 ** L, 2 ** L), complex)
    for j in range(L):
        k = (j + 1) % L
        H += Jx * op(L, {j: X, k: X}) + Jy * op(L, {j: Y, k: Y}) + Jz * op(L, {j: Z, k: Z}) + D * (op(L, {j: X, k: Y}) - op(L, {j: Y, k: X})) + hz * op(L, {j: Z})
    import scipy.linalg as sl
    e, v = sl.eigh(H, driver='evr')
    occ = lambda s: [j for j in range(L) if (s >> (L - 1 - j)) & 1]
    P = [s for s in range(2 ** L) if not any(((s >> (L - 1 - j)) & 1) and ((s >> (L - 1 - (j + 1) % L)) & 1) for j in range(L))]
    m = len(P); A = v[P, :m]; u, sv, vh = np.linalg.svd(A); U = u @ vh
    return U @ np.diag(e[:m]) @ U.conj().T, {s: i for i, s in enumerate(P)}, occ, e[m] - e[m - 1]
idx = lambda sites: sum(1 << (L - 1 - j) for j in sites)
for wt in (1j, 0.6 + 0.8j):
    for V in (20.0, 40.0, 80.0):
        t2 = 1 / (2 * V); w = wt * t2; mu = -2.2 * t2 + 3 / V
        He, pos, occ, gap = ring_heff(1.0, V, w, mu)
        p2 = [He[pos[idx([b, (b + 2) % L])], pos[0]] for b in range(L)]
        p3 = [He[pos[idx([b, (b + 3) % L])], pos[0]] for b in range(L)]
        S = list(pos.keys()); N = np.array([len(occ(s)) for s in S])
        m4 = np.abs(He)[np.abs(N[:, None] - N[None, :]) == 4].max()
        print(f"w/t2={wt} V={V:3.0f}: sep-2 <b,b+2|H|0>: max {np.max(np.abs(p2)):.2e} (pred {abs(2*w.real/V):.2e}); "
              f"sep-3 mean {np.mean(p3).real:.4e} (pred {2/(3*V**2):.4e}), spread {np.ptp(np.abs(p3)):.1e}; max|dN=4| = {m4:.2e} (V^-5 = {V**-5:.2e})", flush=True)
