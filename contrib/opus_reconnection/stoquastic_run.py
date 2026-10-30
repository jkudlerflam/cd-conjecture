"""Stoquastic Trotter network: the exact cross terms of coarse seam records.

Checkerboard Trotterization with K (even) slices of a stoquastic qubit bond
h' = -(Sx Sx + Sy Sy) + Delta Sz Sz (the even-ring image of antiferromagnetic
XXZ under the alternating Z rotation). Column histories x in {0,1}^K; the
spatial transfer matrix B = (A C)^{L/2} is entrywise nonnegative, and
Z_L = Tr B, Z_2L = Tr B^2 for the Trotterized model (checked against the
dense Trotter product).

With records c(x) = x at N equally spaced slices and coherent record weights
a_s = sum_{c(x)=s} B_xx (C121 eq. 20):
  sum_s a_s^2 = D + X,   X = sum_{x!=y, c(x)=c(y)} B_xx B_yy,
  Z_2L        = D + Y,   Y = sum_{x!=y} B_xy B_yx,   D = sum_x B_xx^2.
Also checks the exact Jensen representation
  R_L = -log E_{mu x mu} exp(-DF),  DF(x,y) = log[B_xx B_yy/(B_xy B_yx)],
so R_L <= E_{mu x mu} DF, with mu_x = B_xx/Z_L.
Output: stoquastic_results.json."""
import json
from itertools import product
from pathlib import Path
import numpy as np
from scipy.linalg import expm

ROOT = Path(__file__).resolve().parent

def weights(delta, dt):
    sx = np.array([[0, 1], [1, 0]]) / 2; sy = np.array([[0, -1j], [1j, 0]]) / 2
    sz = np.diag([.5, -.5])
    h = (-(np.kron(sx, sx) + np.kron(sy, sy)) + delta * np.kron(sz, sz)).real
    w = expm(-dt * h)
    assert w.min() >= -1e-15
    return np.maximum(w, 0.0)         # w[(a',b'),(a,b)] = <a'b'|e^{-dt h}|ab>

def transfer(K, w):
    cols = np.array(list(product(range(2), repeat=K)))   # cols[x, t]
    nx = len(cols)
    A = np.ones((nx, nx)); C = np.ones((nx, nx))
    for t in range(K):
        u = (t + 1) % K
        idx_from = 2 * cols[:, t][:, None] + cols[:, t][None, :]
        idx_to = 2 * cols[:, u][:, None] + cols[:, u][None, :]
        fac = w[idx_to, idx_from]
        if t % 2 == 0: A *= fac
        else: C *= fac
    return A, C, cols

def dense_trotter(L, delta, dt, K):
    from zeno_run import _embed
    sx = np.array([[0, 1], [1, 0]]) / 2; sy = np.array([[0, -1j], [1j, 0]]) / 2
    sz = np.diag([.5, -.5])
    h = (-(np.kron(sx, sx) + np.kron(sy, sy)) + delta * np.kron(sz, sz)).real
    # site 0 most significant in _embed; bonds (i,i+1)
    He = sum(_embed(h.reshape(2, 2, 2, 2), i, i + 1, L, 2) for i in range(0, L, 2))
    Ho = sum(_embed(h.reshape(2, 2, 2, 2), i, (i + 1) % L, L, 2) for i in range(1, L, 2))
    Ue, Uo = expm(-dt * He), expm(-dt * Ho)
    P = np.eye(2 ** L)
    for t in range(K):
        P = (Ue if t % 2 == 0 else Uo) @ P
    return np.trace(P)

def main():
    rows = []
    K = 12
    for delta in [0.0, 1.0]:
        for beta in [2.0, 6.0]:
            dt = beta / K
            w = weights(delta, dt)
            A, Cm, cols = transfer(K, w)
            T = A @ Cm
            for L in [4, 6]:
                B = np.linalg.matrix_power(T, L // 2)
                ZL = np.trace(B); Z2 = np.sum(B * B.T)
                if L == 4:
                    zd = dense_trotter(L, delta, dt, K)
                    assert abs(ZL - zd) / zd < 1e-10, (ZL, zd)
                    zd2 = dense_trotter(2 * L, delta, dt, K)
                    assert abs(Z2 - zd2) / zd2 < 1e-10, (Z2, zd2)
                R = 2 * np.log(ZL) - np.log(Z2)
                diag = np.diag(B)
                mu = diag / ZL
                with np.errstate(divide='ignore'):
                    lB = np.log(B)
                ok = np.isfinite(lB) & np.isfinite(lB.T)
                DF = (np.log(diag)[:, None] + np.log(diag)[None, :] - lB - lB.T)
                mm = mu[:, None] * mu[None, :]
                jensen_exact = -np.log(np.sum(mm[ok] * np.exp(-DF[ok])))
                meanDF = np.sum(mm[ok] * DF[ok]) + (np.inf if np.any(mm[~ok] > 0) else 0.0)
                assert abs(jensen_exact - R) < 1e-9 and R <= meanDF + 1e-12
                D = np.sum(diag ** 2); Y = Z2 - D
                for N in [1, 2, 3, 4, 6, 12]:
                    rec = cols[:, ::K // N][:, :N] @ (2 ** np.arange(N))
                    same = rec[:, None] == rec[None, :]
                    off = same & ~np.eye(len(rec), dtype=bool)
                    X = np.sum(np.outer(diag, diag)[off])
                    Yc = np.sum((B * B.T)[off])
                    a = np.bincount(rec, weights=diag, minlength=2 ** N)
                    assert abs(a @ a - (D + X)) / (D + X) < 1e-10
                    row = dict(Delta=delta, beta=beta, K=K, L=L, N=N, R=float(R),
                               mean_DF=float(meanDF),
                               log_coh_ratio=float(np.log((D + X) / Z2)),
                               X_over_Yc=(float(X / Yc) if Yc > 0 else None), X_over_Y=float(X / Y),
                               D_over_Z2L=float(D / Z2))
                    rows.append(row)
                    print(f"Delta={delta} beta={beta} L={L} N={N} R={R:.4f} <DF>={meanDF:.4f} "
                          f"log(sum a^2/Z2L)={row['log_coh_ratio']:.4f} X/Yc={row["X_over_Yc"]} "
                          f"X/Y={row['X_over_Y']:.4f} D/Z2L={row['D_over_Z2L']:.2e}")
    (ROOT / 'stoquastic_results.json').write_text(json.dumps(rows, indent=1) + '\n')

if __name__ == '__main__':
    main()
