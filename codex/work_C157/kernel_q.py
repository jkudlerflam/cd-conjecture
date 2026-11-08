"""Second-order YY kernel K_jl = int dtau <da_j(tau) da_l(0)>_A and its Fourier transform, for
(i) the critical A ring (XZ point) and (ii) a dimerized (massive) A ring with hoppings 1 + d (-1)^j.

Same free-fermion formula as C144 kernel.py, generalized to an arbitrary quadratic A:
    K_jl = 2 sum_{mu,nu} W^j_{mu nu} conj(W^l_{mu nu}) / (lam_mu + lam_nu),  W^j = s_j (v_j v_{j+1}^T - v_{j+1} v_j^T),
where da_j = a_j - <a_j> and a_j = i g_j g_{j+1} with the hopping sign of bond j absorbed into s_j = sign.
Outputs: K(r) for r = 0..8 (both sublattices), the Fourier transform Kt(q) = sum_r K(r) e^{iqr} (sublattice average),
its minimum (positivity check), the location of its maximum, the RZFA-matched ratio g = 2 beta^2 K(2), and the sign
of K(2) needed for the repulsive (c = 3/2) side of RZFA, which is K(2) < 0.
"""
import numpy as np


def chain(L, hop, bc=-1):
    A = np.zeros((L, L))
    for j in range(L):
        k = (j + 1) % L
        s = bc if j == L - 1 else 1.0
        A[j, k] += hop[j] * s
        A[k, j] -= hop[j] * s
    return A


def kernel_rows(L, hop, bc=-1, rows=(0, 1)):
    A = chain(L, hop, bc)
    w, U = np.linalg.eigh(2j * A)
    pos = w > 1e-12
    lam = w[pos]; V = U[:, pos]
    den = lam[:, None] + lam[None, :]

    def W(j):
        k = (j + 1) % L
        s = bc if j == L - 1 else 1.0
        return s * (np.outer(V[j], V[k]) - np.outer(V[k], V[j]))

    out = {}
    Ws = [W(l) for l in range(L)]
    for j in rows:
        out[j] = np.array([2 * np.sum(Ws[j] * np.conj(Ws[(j + r) % L]) / den) for r in range(L)]).real
    return out, -0.5 * lam.sum()


if __name__ == "__main__":
    L = 400
    q = np.linspace(0, np.pi, 9)
    for d in (0.0, 0.05, 0.2, 0.5):
        hop = 1.0 + d * (-1.0) ** np.arange(L)
        K, E0 = kernel_rows(L, hop)
        Kav = 0.5 * (K[0] + K[1])
        r = np.arange(L); rr = np.minimum(r, L - r)
        Kt = lambda qq: np.sum(Kav * np.cos(qq * rr))
        Ktq = np.array([Kt(x) for x in q])
        # full matrix positivity check on a smaller ring
        Ls = 60
        hs = 1.0 + d * (-1.0) ** np.arange(Ls)
        Ks, _ = kernel_rows(Ls, hs, rows=range(Ls))
        Kmat = np.array([np.roll(Ks[j], j) for j in range(Ls)])
        Kmat = 0.5 * (Kmat + Kmat.T)
        print(f"d = {d}:  E0/L = {E0/L:.6f}")
        print("   K(r), even j :", np.round(K[0][:9], 6).tolist())
        print("   K(r), odd j  :", np.round(K[1][:9], 6).tolist())
        print("   Kt(q) at q = k pi/8, k = 0..8:", np.round(Ktq, 5).tolist())
        print(f"   min eigenvalue of K (ring L = {Ls}): {np.linalg.eigvalsh(Kmat).min():.3e};  "
              f"K(2) avg = {Kav[2]:+.6f}  -> RZFA g/beta^2 = 2 K(2) = {2*Kav[2]:+.5f} (attractive if > 0)")
        tail = [round(float(Kav[x] * (-1) ** x * 2 * np.pi * x), 4) for x in (2, 4, 8, 16, 32)]
        print("   K(r)(-1)^r 2 pi r at r = 2,4,8,16,32:", tail)
