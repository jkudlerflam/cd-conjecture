#!/usr/bin/env python3
"""C170: checks of C165 PROOF 1.2.
(i) open chain: spectrum of the parent (a arbitrary, Dz, hz) equals that of the AF XXZ chain mu.mu_perp + Dz mu^z mu^z + hz mu^z;
(ii) ring with a L in 2 pi Z, L even: in the mu frame (sigma_j = R_z(-a j + pi j) mu_j) the parent commutes with T and S^z separately,
     while the full bond with random harmonics commutes only with T exp(i alpha S^z_tot) for one alpha = +-(a - pi).  This is the
     exact operator identity behind the selection rule P + Q (a - pi) = 0 mod 2 pi used in PROOF 4.2."""
import numpy as np, scipy.sparse as sp
from qed import ham, op, P
from lr2 import parent
rng = np.random.default_rng(7)
# (i)
L = 8; a = 0.77; Dz = 1.3; hz = 0.4
Hs = ham(parent(a, Dz), [0, 0, hz], L)[0].toarray()
Hm = ham(np.diag([1, 1, Dz]), [0, 0, hz], L)[0].toarray()
print('(i) max spectral difference, open L=8:', np.max(np.abs(np.linalg.eigvalsh(Hs) - np.linalg.eigvalsh(Hm))))
# (ii)
L = 8; a = 2 * np.pi * 3 / 8
dJ = 0.2 * rng.standard_normal((3, 3)); dJ[2, 2] = 0; dJ[:2, :2] -= np.trace(dJ[:2, :2]) / 2 * np.eye(2)  # no charge-0 in-plane part
dJ[0, 1] = dJ[1, 0] = 0.1   # keep a symmetric charge-2 part, drop charge-0 DM_z
dB = np.array([0.15, -0.1, 0.0])
H0 = ham(parent(a, 1.0), [0, 0, 0.0], L, pbc=True)[0]
H1 = ham(parent(a, 1.0) + dJ, dB, L, pbc=True)[0]
# unitary U = prod_j exp(-i phi_j Z_j / 2), rotating sigma_j -> R_z(phi_j) sigma_j ; we need mu_j = R_z(a j - pi j) sigma_j
dim = 2 ** L
Zs = [op(P[2], j, L).diagonal().real for j in range(L)]
phis = [(a - np.pi) * j for j in range(L)]
u = np.exp(-0.5j * sum(phi * z for phi, z in zip(phis, Zs)))
U = sp.diags(u)
# translation T: site j -> j+1 (cyclic) on basis states
idx = np.arange(dim)
bits = (idx[:, None] >> (L - 1 - np.arange(L))[None, :]) & 1
tb = np.roll(bits, 1, axis=1)
tidx = (tb * (1 << (L - 1 - np.arange(L)))[None, :]).sum(1)
T = sp.csr_matrix((np.ones(dim), (tidx, idx)), shape=(dim, dim))
Sz = sp.diags(0.5 * sum(Zs))
def comm(A, B): return abs((A @ B - B @ A)).max()
for name, H in (('parent', H0), ('parent+harmonics', H1)):
    Hm = U @ H @ U.conj().T
    print(f'(ii) {name}: [H_mu, T] = {comm(Hm, T):.2e}, [H_mu, Sz] = {comm(Hm, Sz):.2e}')
    for s in (+1, -1):
        G = T @ sp.diags(np.exp(1j * s * (a - np.pi) * Sz.diagonal()))
        print(f'      [H_mu, T exp({"+" if s>0 else "-"}i(a-pi)Sz)] = {comm(Hm, G):.2e}')
