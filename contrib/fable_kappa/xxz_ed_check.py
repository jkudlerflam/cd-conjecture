"""Exact diagonalization check of the gradient identity and of the closed form for D_2 on
the XXZ line, including the approach to the ferromagnetic point Delta -> -1.

Rings of even length L in the S^z_tot = 0 sector (Lanczos).  The infinite-chain two-site
marginal is exact (Yang-Yang e_0 and its Delta-derivative, see xxz_closed_form.py), so
D_2(L) = L^2 ||rho_2^inf - rho_2^L||_1 needs no fit.  Printed for each L:
    L^2 (g_zz^inf - g_zz^L)       -> (pi/6) A'(Delta)          [gradient identity, S^zS^z direction]
    L^2 (2g_xx^inf - 2g_xx^L)     -> (pi/6) (A - Delta A')     [Euler / bond pairing]
    L^2 (e_0 - E_0/L) 6/pi        -> A = v(Delta)
    D_2(L)                        -> closed form
Usage: python3 xxz_ed_check.py [Lmax]   (default 20; L=24 needs about 3 GB and a few minutes)
"""
import sys, numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spla
from xxz_closed_form import v_of, vprime_of, e0_of, gzz_inf, gxx_inf, D2_closed

Lmax = int(sys.argv[1]) if len(sys.argv) > 1 else 20

def sz0_basis(L):
    # all L-bit integers with L/2 bits set, sorted
    from itertools import combinations
    states = []
    for comb in combinations(range(L), L//2):
        s = 0
        for b in comb: s |= 1 << b
        states.append(s)
    return np.array(sorted(states), dtype=np.int64)

def build_H(L, Delta, states):
    n = len(states)
    idx = np.arange(n)
    diag = np.zeros(n)
    rows = []; cols = []; vals = []
    for i in range(L):
        j = (i+1) % L
        bi = (states >> i) & 1; bj = (states >> j) & 1
        same = bi == bj
        diag += np.where(same, Delta/4, -Delta/4)
        flip = states[~same] ^ ((1 << i) | (1 << j))
        k = np.searchsorted(states, flip)
        rows.append(idx[~same]); cols.append(k); vals.append(np.full(k.shape, 0.5))
    rows = np.concatenate(rows); cols = np.concatenate(cols); vals = np.concatenate(vals)
    H = sp.coo_matrix((vals, (rows, cols)), shape=(n, n)).tocsr()
    H = H + sp.diags(diag)
    return H

def two_site(L, states, psi):
    # bond (0,1): g_zz = <S^z_0 S^z_1>, g_xx = <S^x_0 S^x_1>  (= <S^y S^y>)
    b0 = states & 1; b1 = (states >> 1) & 1
    same = b0 == b1
    p = np.abs(psi)**2
    gzz = 0.25*(p[same].sum() - p[~same].sum())
    flip = states[~same] ^ 3
    k = np.searchsorted(states, flip)
    # <S^+_0 S^-_1 + h.c.>/2 = sum over differing configs psi_s psi_{s'} /2 ; S^xS^x = that /2
    gxx = 0.5*np.sum(psi[~same]*psi[k])/2
    return gzz, gxx

sx = np.array([[0,1],[1,0]])/2; sy = np.array([[0,-1j],[1j,0]])/2; sz = np.diag([.5,-.5])
SxSx = np.kron(sx,sx).real; SySy = np.kron(sy,sy).real; SzSz = np.kron(sz,sz)
def rho2(gxx, gzz):
    return 0.25*np.eye(4) + 4*gxx*(SxSx + SySy) + 4*gzz*SzSz

for Delta in [0.5, 0.0, -0.5, -0.8, -0.9, -0.95]:
    A = v_of(Delta); Ap = vprime_of(Delta); e0 = e0_of(Delta)
    gzzI = gzz_inf(Delta); gxxI = gxx_inf(Delta)
    D2c, a, b = D2_closed(Delta)
    rinf = rho2(gxxI, gzzI)
    print(f"== XXZ Delta = {Delta}:  A = v = {A:.6f}, A' = {Ap:.6f}, e_0 = {e0:.8f}, g_zz^inf = {gzzI:.7f}, g_xx^inf = {gxxI:.7f}")
    print(f"   targets: L^2 dg_zz -> (pi/6)A' = {a:.5f};  L^2 d(2g_xx) -> (pi/6)(A - Delta A') = {b:.5f};  D_2 closed = {D2c:.5f}")
    print("    L     E0/L        L^2(e0-E0/L)6/pi   L^2 dg_zz   L^2 d(2g_xx)   D_2(L)")
    for L in range(8, Lmax+1, 2):
        states = sz0_basis(L)
        H = build_H(L, Delta, states)
        w, v = spla.eigsh(H, k=1, which='SA', tol=1e-12)
        psi = v[:,0]; E0 = w[0]
        gzz, gxx = two_site(L, states, psi)
        dr = rinf - rho2(gxx, gzz)
        tn = np.sum(np.abs(np.linalg.eigvalsh(dr)))
        print(f"  {L:3d}  {E0/L:12.8f}  {L**2*(e0-E0/L)*6/np.pi:14.5f}  {L**2*(gzzI-gzz):10.5f}  {L**2*2*(gxxI-gxx):11.5f}  {L**2*tn:10.5f}")
    print()
