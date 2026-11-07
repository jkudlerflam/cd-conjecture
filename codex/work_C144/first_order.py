"""First-order projection P V P of every tangent direction at XZ onto the degenerate manifold (ring L).
Fits P V P = c0 + c1 P (sum_j b_j) P and prints c0, c1 and the residual norm. Predictions (PROOF 1.2):
zx: c1 = 1; yy: c1 = <a>_L = E0/L; xz: c0 = E0, c1 = 0; all others (xx, zz, By, Bx, Bz, xy, yx, yz, zy): P V P = 0."""
import numpy as np
from xzlib import term, ring_operator

L = 10
N = 1 << L
def dense(h):
    op = ring_operator(h, L); H = np.zeros((N, N), complex)
    for x in range(N):
        e = np.zeros(N, complex); e[x] = 1; H[:, x] = op.matvec(e)
    return H
H0 = dense(term("xz"))
w, U = np.linalg.eigh(H0)
man = w < w[0] + 1e-8
P = U[:, man]
E0 = w[0]
print(f"L={L}: manifold dim {man.sum()}, E0 = {E0:.10f}, <a> = E0/L = {E0/L:.10f}, gap above = {w[man.sum()]-E0:.6f}")
Bsum = P.conj().T @ dense(term("zx")) @ P
D = P.shape[1]
for name in ["xz", "zx", "yy", "xx", "zz", "xy", "yx", "yz", "zy", "Bx", "By", "Bz"]:
    V = P.conj().T @ dense(term(name)) @ P
    # least squares V = c0 I + c1 Bsum
    A = np.stack([np.eye(D).ravel(), Bsum.ravel()], 1)
    c, *_ = np.linalg.lstsq(A, V.ravel(), rcond=None)
    res = np.linalg.norm(V - c[0] * np.eye(D) - c[1] * Bsum)
    print(f"{name:3s}: c0 = {c[0].real:+.10f}  c1 = {c[1].real:+.10f}  residual = {res:.1e}  ||PVP|| = {np.linalg.norm(V):.3e}")
