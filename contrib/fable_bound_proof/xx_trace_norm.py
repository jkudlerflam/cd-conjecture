"""XX chain (free fermions): two-site reduced density matrix of the ring ground
state versus the infinite chain, and the coefficient D_2 = lim L^2 ||rho_2^L - rho_2^inf||_1.
h = S^x S^x + S^y S^y, so t = 1/2, e_0 = -1/pi, v = 1, c = 1, M = 1/2, A = 1.
Conjecture-level target: D_2 <= (pi^2/6) log_2 d = 1.6449 ; lower bound from A: (pi/6) A/M = pi/3 = 1.0472.
"""
import numpy as np
from scipy.linalg import logm, expm

sm = np.array([[0,1],[0,0]],dtype=float)  # sigma^- (annihilates up=0 ->?); define c = |0><1| on basis (empty, occupied)
Z = np.diag([1.,-1.])
I2 = np.eye(2)
c1 = np.kron(sm, I2); c2 = np.kron(Z, sm)
cs=[c1,c2]

def rho_from_G(G):
    # Gaussian number-conserving state with <c_i^dag c_j> = G_ij on 2 modes
    w,U = np.linalg.eigh(G)
    w = np.clip(w,1e-14,1-1e-14)
    K = U @ np.diag(np.log((1-w)/w)) @ U.conj().T
    Kop = sum(K[i,j]*cs[i].conj().T@cs[j] for i in range(2) for j in range(2))
    r = expm(-Kop); return r/np.trace(r)

def G_ring(L):
    N=L//2
    if N%2==1: ks=2*np.pi*np.arange(L)/L            # periodic
    else: ks=2*np.pi*(np.arange(L)+0.5)/L           # antiperiodic
    e=-np.cos(ks); occ=np.argsort(e)[:N]
    G=np.zeros((2,2),dtype=complex)
    for a in range(2):
        for b in range(2):
            G[a,b]=np.sum(np.exp(1j*ks[occ]*(a-b)))/L
    return G, 2*0.5*np.sum(e[occ])   # energy with t=1/2: eps_k=-cos k

Ginf=np.array([[0.5,1/np.pi],[1/np.pi,0.5]])
rinf=rho_from_G(Ginf)
hXX=np.zeros((4,4)); # S^xS^x+S^yS^y = (1/2)(S+S- + S-S+) in spin basis, up to ordering = fermion hop
hXX = 0.5*(c1.conj().T@c2 + c2.conj().T@c1)  # with t=1/2 the bond is -(1/2)(c1^dag c2+h.c.) ; sign irrelevant for check via unitary
hXX = -hXX
print("bond eigenvalues", np.round(np.linalg.eigvalsh(hXX),6))
e0=-1/np.pi
print(" L   L^2*(e0 - E0/L)*6/pi   L^2*Tr[h drho]*6/pi   L^2||drho||_1   M*L^2||drho||_1*6/pi")
for L in [8,10,12,16,20,24,32,48,64,100,200,400,1000,2000]:
    G,E0=G_ring(L)
    r=rho_from_G(G)
    dr=rinf-r
    tn=np.sum(np.abs(np.linalg.eigvalsh(dr)))
    print(f"{L:5d} {L**2*(e0-E0/L)*6/np.pi:12.6f} {L**2*np.trace(hXX@dr).real*6/np.pi:14.6f} {L**2*tn:14.6f} {0.5*L**2*tn*6/np.pi:14.6f}")
