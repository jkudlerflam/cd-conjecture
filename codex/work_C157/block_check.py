"""Check of the exact rewriting of the bisector-field chain (PROOF section 4):
XZ frame H = sum X_j Z_{j+1} + B sum (X_j + Z_j)/sqrt2 on a ring of L = 2n qubits.
Hadamard on even sites -> compass frame; blocks (2i-1, 2i): tau^x = X1 X2, tau^z = Z1, nu^z = Z1 Z2, nu^x = X2;
sigma^z_i = tau^z_i nu^z_i tau^z_{i+1}, sigma^x_i = nu^x_i, G_i = nu^x_{i-1} tau^x_i nu^x_i, gx_i = tau^z_i;
claimed: H = sum_i [ sigma^z_i + G_{i+1} sigma^x_i sigma^x_{i+1} + (B/sqrt2)( sigma^x_i (1 + G_{i+1}) + gx_{i+1} (1 + sigma^z_i) ) ].
Checks the Pauli algebra of (sigma_i, g_i) and the operator identity, on L = 8 and 10."""
import numpy as np
from functools import reduce
X = np.array([[0, 1], [1, 0]], complex); Z = np.diag([1., -1]).astype(complex); I = np.eye(2); Y = 1j * X @ Z
Hd = np.array([[1, 1], [1, -1]]) / np.sqrt(2)
for L in (8, 10):
    n = L // 2
    def op(d):
        return reduce(np.kron, [d.get(k, I) for k in range(L)])
    B = 0.37
    Hq = sum(op({j: X, (j + 1) % L: Z}) for j in range(L)) + (B / np.sqrt(2)) * sum(op({j: X}) + op({j: Z}) for j in range(L))
    U = reduce(np.kron, [Hd if k % 2 == 1 else I for k in range(L)])   # Hadamard on even sites (1-based) = odd 0-based
    Hc = U @ Hq @ U.conj().T
    q1 = lambda i: 2 * i; q2 = lambda i: 2 * i + 1
    taux = [op({q1(i): X, q2(i): X}) for i in range(n)]; tauz = [op({q1(i): Z}) for i in range(n)]
    nuz = [op({q1(i): Z, q2(i): Z}) for i in range(n)]; nux = [op({q2(i): X}) for i in range(n)]
    sz = [tauz[i] @ nuz[i] @ tauz[(i + 1) % n] for i in range(n)]; sx = nux
    G = [nux[(i - 1) % n] @ taux[i] @ nux[i] for i in range(n)]; gx = tauz
    Hb = sum(sz[i] + G[(i + 1) % n] @ sx[i] @ sx[(i + 1) % n] + (B / np.sqrt(2)) * (sx[i] @ (np.eye(2 ** L) + G[(i + 1) % n])
             + gx[(i + 1) % n] @ (np.eye(2 ** L) + sz[i])) for i in range(n))
    comm = lambda a, b: np.abs(a @ b - b @ a).max(); acomm = lambda a, b: np.abs(a @ b + b @ a).max()
    errs = [acomm(sz[0], sx[0]), acomm(G[0], gx[0]), comm(sz[0], sx[1]), comm(sz[0], G[0]), comm(sz[0], G[1]), comm(sx[0], G[0]),
            comm(gx[0], sz[0]), comm(gx[1], sz[0]), comm(gx[0], G[1]), comm(sx[0], gx[0])]
    print(f"L = {L}: algebra max error {max(errs):.1e};  ||H_compass - H_block|| = {np.abs(Hc - Hb).max():.1e}")
