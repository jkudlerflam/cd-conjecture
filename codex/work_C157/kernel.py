"""Second-order kernel K_{jl} = <0| da_j Q (H_A - E0)^{-1} Q da_l |0> of the critical A Majorana chain.

H_A = sum_j a_j, a_j = i g_j g_{j+1} (ring of L Majoranas, BC sign bc on the closing bond).
Free-fermion formula: with h = 2iA, positive eigenpairs (lam_mu, v_mu),
    G_mn(tau) = <g_m(tau) g_n> = 2 sum_mu exp(-lam_mu tau) v_m,mu conj(v_n,mu),
    K_jl = int_0^inf dtau [G_jl G_{j+1,l+1} - G_{j,l+1} G_{j+1,l}]
         = 2 sum_{mu,nu} W^j_{mu nu} conj(W^l_{mu nu}) / (lam_mu + lam_nu),
    W^j_{mu nu} = s_j (v_{j,mu} v_{j+1,nu} - v_{j+1,mu} v_{j,nu})   (s_j = bc on the closing bond).
Checks: (i) against ED of the full XZ qubit ring (L = 8, 10), (ii) reality and symmetry,
(iii) the analytic tail K(r) -> (-1)^r / (2 pi r) derived in PROOF.md section 3.1.
"""
import sys
import numpy as np
from xzlib import chain_A, term, ring_operator


def kernel_row(L, bc=-1):
    A = chain_A(L, bc)
    h = 2j * A
    w, U = np.linalg.eigh(h)
    pos = w > 1e-12
    lam = w[pos]; V = U[:, pos]
    E0 = -0.5 * lam.sum()
    den = lam[:, None] + lam[None, :]

    def W(j):
        k = (j + 1) % L
        s = bc if j == L - 1 else 1.0
        return s * (np.outer(V[j], V[k]) - np.outer(V[k], V[j]))

    W0 = W(0)
    K = np.array([2 * np.sum(W0 * np.conj(W(l)) / den) for l in range(L)])
    return K, E0


def ed_kernel(L):
    """K_{0,l} from the full XZ qubit ring by dense ED (manifold = all states within 1e-8 of E0)."""
    h = term("xz")
    N = 1 << L
    op = ring_operator(h, L)
    H = np.zeros((N, N), complex)
    for x in range(N):
        e = np.zeros(N, complex); e[x] = 1
        H[:, x] = op.matvec(e)
    w, U = np.linalg.eigh(H)
    E0 = w[0]
    man = w < E0 + 1e-8
    psi = U[:, 0]
    # a_l operators as matrices
    ops = []
    for l in range(L):
        hl = np.zeros((N, N), complex)
        # build X_l Z_{l+1}
        x = np.arange(N)
        a = L - 1 - l; b = L - 1 - ((l + 1) % L)
        zb = 1 - 2 * ((x >> b) & 1)
        y = x ^ (1 << a)
        hl[y, x] = zb
        ops.append(hl)
    Ex = w[~man] - E0
    Ut = U[:, ~man]
    K = []
    v0 = Ut.conj().T @ (ops[0] @ psi)
    for l in range(L):
        vl = Ut.conj().T @ (ops[l] @ psi)
        K.append(np.sum(np.conj(v0) * vl / Ex))
    return np.array(K), E0, man.sum()


if __name__ == "__main__":
    for L in (8, 10):
        Ked, E0ed, D = ed_kernel(L)
        out = {}
        for bc in (-1, 1):
            K, E0 = kernel_row(L, bc)
            out[bc] = (K, E0)
        bc = min(out, key=lambda b: out[b][1])
        print(f"L={L}: ED E0={E0ed:.10f} deg={D}; free E0(bc=-1)={out[-1][1]:.10f} E0(bc=+1)={out[1][1]:.10f}")
        print("   ED K(0,l)  :", np.round(Ked.real, 8).tolist(), " max|Im|", np.abs(Ked.imag).max())
        print("   free K(0,l):", np.round(out[bc][0].real, 8).tolist(), " max|Im|", np.abs(out[bc][0].imag).max())
        print("   max diff", np.abs(Ked - out[bc][0]).max())
    print("Tail check, K(r) * 2 pi r * (-1)^r  (should -> 1):")
    for L in (200, 800):
        K, E0 = kernel_row(L, -1)
        rs = [2, 3, 4, 5, 8, 16, 32, 64]
        print(f"  L={L}: E0/L={E0/L:.8f} (-2/pi={-2/np.pi:.8f}); max|Im K|={np.abs(K.imag).max():.1e}; "
              + " ".join(f"r={r}:{(K[r].real*2*np.pi*r*(-1)**r):.5f}" for r in rs))
        Kpi = np.sum(K.real * (-1) ** np.arange(L))
        print(f"        K~(pi) = {Kpi:.6f};  (1/pi) log L = {np.log(L)/np.pi:.6f};  difference {Kpi - np.log(L)/np.pi:.6f}")
    K, _ = kernel_row(400, -1)
    np.save("kernel_L400.npy", K.real)
    print("K(r), r=0..8 at L=400:", np.round(K.real[:9], 8).tolist())
