"""D_2 for n decoupled XX chains stacked into one chain with d = 2^n.

Bond h = sum_a lambda_a h_a^XX (h_a^XX acting on the a-th spin-1/2 factor of the site), so the
sectors are n free fermions with velocities v_a = lambda_a, A = sum_a lambda_a, c = n.
The ground state is a product and does not depend on the lambda_a at all, so D_2 is
velocity-blind: slowing a sector leaves its contribution to rho_2^inf - rho_2^L unchanged.
This script computes D_2(n) = L^2 || (x)_a rho_a^inf - (x)_a rho_a^L ||_1 exactly from the
two-mode Gaussian marginals (L = 2000), and compares with n D_2^XX and with log_2 d = n.
Half filling throughout (symmetric Fermi seas: antiperiodic fermions for N even, periodic for N odd).
"""
import numpy as np
from scipy.linalg import expm
from functools import reduce

sm = np.array([[0,1],[0,0]],dtype=float); Z = np.diag([1.,-1.]); I2 = np.eye(2)
c1 = np.kron(sm, I2); c2 = np.kron(Z, sm); cs=[c1,c2]

def rho_from_G(G):
    w,U = np.linalg.eigh(G); w = np.clip(w,1e-15,1-1e-15)
    K = U @ np.diag(np.log((1-w)/w)) @ U.conj().T
    Kop = sum(K[i,j]*cs[i].conj().T@cs[j] for i in range(2) for j in range(2))
    r = expm(-Kop); return (r/np.trace(r)).real

def G_ring(L):
    N=L//2
    ks = 2*np.pi*np.arange(L)/L if N%2==1 else 2*np.pi*(np.arange(L)+0.5)/L
    e=-np.cos(ks); occ=np.argsort(e)[:N]
    G=np.zeros((2,2),dtype=complex)
    for a in range(2):
        for b in range(2): G[a,b]=np.sum(np.exp(1j*ks[occ]*(a-b)))/L
    return G.real

L = 2000
rinf = rho_from_G(np.array([[0.5,1/np.pi],[1/np.pi,0.5]]))
rL = rho_from_G(G_ring(L))
X1 = rinf - rL
D2xx = L**2*np.sum(np.abs(np.linalg.eigvalsh(X1)))
print(f"single XX chain: D_2 = {D2xx:.6f}  (exact (pi+2)/3 = {(np.pi+2)/3:.6f})")
print()
print(" n   d=2^n   D_2(n)      D_2/n     D_2/(n D_2^XX)   (pi^2/6) log2 d")
# the product marginal lives on (C^2 x C^2)^{x n} with site ordering (site1 flavours..., site2 flavours...);
# a reordering of tensor factors is unitary and does not change the trace norm, so we can take
# the plain Kronecker product over flavours of the two-site marginals.
for n in range(1, 6):
    Rinf = reduce(np.kron, [rinf]*n); RL = reduce(np.kron, [rL]*n)
    D2 = L**2*np.sum(np.abs(np.linalg.eigvalsh(Rinf - RL)))
    print(f"{n:2d}  {2**n:5d}  {D2:10.5f}  {D2/n:9.5f}  {D2/(n*D2xx):12.5f}   {np.pi**2/6*n:12.5f}")
print()
print("Leading-order decomposition: Rinf - RL = sum_a (rinf^{(x)(a-1)} x X1 x rinf^{(x)(n-a)}) + O(L^-4).")
print("Trace norm of that sum (should match D_2(n) to O(1/L^2)):")
for n in range(1, 6):
    S = sum(reduce(np.kron, [rinf]*a + [X1] + [rinf]*(n-1-a)) for a in range(n))
    print(f"  n={n}: {L**2*np.sum(np.abs(np.linalg.eigvalsh(S))):.5f}")
