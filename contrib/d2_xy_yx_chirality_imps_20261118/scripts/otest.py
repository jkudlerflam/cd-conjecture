"""Independent checks of ofd.py: (1) target, open chain, L=12, sparse ED of the ORIGINAL-frame H
a XY + b YX + ZZ - h Z (built by hand, no TeNPy) vs finite DMRG in the rotated frame, both parity sectors;
(2) XX x XX control, open chain, exact free-fermion energy vs finite DMRG."""
import numpy as np, scipy.sparse as sp
from scipy.sparse.linalg import eigsh
from tenpy.algorithms import dmrg
import ims, ofd

X = sp.csr_matrix([[0, 1], [1, 0]], dtype=complex); Y = sp.csr_matrix([[0, -1j], [1j, 0]]); Zm = sp.csr_matrix([[1, 0], [0, -1]], dtype=complex)


def op(o, i, L):
    m = sp.identity(1, format="csr", dtype=complex)
    for j in range(L):
        m = sp.kron(m, o if j == i else sp.identity(2, dtype=complex), format="csr")
    return m


def ed(L, d, h):
    a, b = (1 + d) / 2, (1 - d) / 2
    H = sp.csr_matrix((2 ** L, 2 ** L), dtype=complex)
    for i in range(L - 1):
        H = H + a * op(X, i, L) @ op(Y, i + 1, L) + b * op(Y, i, L) @ op(X, i + 1, L) + op(Zm, i, L) @ op(Zm, i + 1, L)
    for i in range(L):
        H = H - h * op(Zm, i, L)
    ndn = np.array([bin(s).count("1") for s in range(2 ** L)])   # bit=1 <-> Z=-1 (kron order)
    out = {}
    for p in (0, 1):
        idx = np.where(ndn % 2 == p)[0]
        Hs = H[idx][:, idx]
        out[p] = eigsh(Hs, k=1, which="SA")[0][0].real
    return out


def dm(P, L, par, chi=64):
    M = ofd.fmodel(P, L); psi = ofd.init_state(M, P, L, par, 1)
    eng = dmrg.TwoSiteDMRGEngine(psi, M, {"trunc_params": {"chi_max": chi, "svd_min": 1e-14}, "mixer": True,
                                          "max_sweeps": 60, "min_sweeps": 20, "max_E_err": 1e-14})
    return eng.run()[0]


def ff(L, t, h):
    """open XX chain -t sum (s+s- + h.c.) - h sum Z: E0 over both fermion parities (JW, Z = 2n-1)."""
    k = np.pi * np.arange(1, L + 1) / (L + 1)
    e = np.sort(-2 * t * np.cos(k) - 2 * h)
    pre = np.concatenate([[0], np.cumsum(e)]) + h * L
    return {n % 2: min(pre[m] for m in range(L + 1) if m % 2 == n % 2) for n in (0, 1)}


if __name__ == "__main__":
    for d, h in ((0.1875, 1.783), (0.125, 1.7405)):
        L = 12
        E = ed(L, d, h)
        P = ims.params({"model": "target", "d": str(d), "h": str(h)})
        print("target L=%d d=%g h=%g  ED even/odd-#down %.12f %.12f | DMRG par0 %.12f par1 %.12f" %
              (L, d, h, E[0], E[1], dm(P, L, 0), dm(P, L, 1)), flush=True)
    L = 32
    P = ims.params({"model": "ctrl", "tA": "1", "hA": "0.3", "tB": "1", "J": "0", "g": "0.3"})
    A = ff(L, 1, 0.3); B = ff(L, 1, 0.3)
    print("ctrl XXxXX L=%d  exact (A even) %.12f (A odd) %.12f | DMRG par0 %.12f par1 %.12f" %
          (L, A[0] + min(B.values()), A[1] + min(B.values()), dm(P, L, 0, 128), dm(P, L, 1, 128)), flush=True)
