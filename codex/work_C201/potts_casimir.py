"""Potts ring: M of the bond, e_0 from the XXZ Bethe integral, and L(L e_0 - E_L)/M.
Claude Fable 5.2, work_C201."""
import numpy as np, itertools
from scipy.integrate import quad
import scipy.sparse as sp, scipy.sparse.linalg as spl

def eps0(gamma, s):
    """ground energy per site of H_X=(1/2)sum[XX+YY+Delta(ZZ-1)], Delta = s*cos(gamma), gamma in (0,pi/2]."""
    g = gamma if s > 0 else np.pi-gamma      # Delta = cos(g)
    D = np.cos(g)
    if abs(g) < 1e-12:
        e = 0.25-np.log(2)
    else:
        f = lambda x: 1.0/(2*np.cosh(np.pi*x)*(np.cosh(2*g*x)-np.cos(g)))
        I = quad(f, -40, 40, limit=400)[0]
        e = D/4-np.sin(g)**2*I
    return 2*e-D/2

def potts_bond(Q, s):
    E = np.ones((Q, Q))/np.sqrt(Q); I = np.eye(Q)
    F = np.zeros((Q*Q, Q*Q))
    for a in range(Q): F[a*Q+a, a*Q+a] = np.sqrt(Q)
    return -s*(F+(np.kron(E, I)+np.kron(I, E))/2)

def ring_E(Q, L, s):
    h = sp.csr_matrix(potts_bond(Q, s))
    H = sp.csr_matrix((Q**L, Q**L))
    # build via permutation of site order: bond (j,j+1)
    for j in range(L):
        if j < L-1:
            H = H+sp.kron(sp.kron(sp.identity(Q**j), h), sp.identity(Q**(L-j-2)), format='csr')
        else:
            # bond (L-1,0): build explicitly
            dim = Q**L
            idx = np.arange(dim)
            first = idx//Q**(L-1); last = idx % Q; mid = (idx//Q) % Q**(L-2) if L > 2 else np.zeros(dim, int)
            hd = potts_bond(Q, s)
            rows, cols, vals = [], [], []
            for a in range(Q):
                for b in range(Q):
                    for a2 in range(Q):
                        for b2 in range(Q):
                            v = hd[a2*Q+b2, a*Q+b]       # (last=a, first=b) -> (a2,b2)
                            if v == 0: continue
                            sel = idx[(last == a) & (first == b)]
                            new = b2*Q**(L-1)+((sel//Q) % Q**(L-2))*Q+a2 if L > 2 else np.full(len(sel), b2*Q+a2)
                            rows.extend(new); cols.extend(sel); vals.extend([v]*len(sel))
            H = H+sp.csr_matrix((vals, (rows, cols)), shape=(dim, dim))
    if Q**L <= 1500:
        return np.linalg.eigvalsh(H.toarray())[0]
    return spl.eigsh(H, k=1, which='SA', tol=1e-12)[0][0]

for Q, Lmax in ((3, 11), (4, 9), (2, 14)):
    gam = np.arccos(np.sqrt(Q)/2)
    for s in (+1, -1):
        hb = potts_bond(Q, s); ev = np.linalg.eigvalsh(hb); M = (ev[-1]-ev[0])/2
        e0 = 2*eps0(gam, s)
        print(f"Q={Q} s={s:+d}: bond spectrum min={ev[0]:.6f} max={ev[-1]:.6f} M={M:.6f} e0 per Potts site={e0:.10f}")
        for L in range(2, Lmax+1):
            E = ring_E(Q, L, s)
            print(f"   L={L:2d} E_L={E:.8f}  L e0 - E_L={L*e0-E:+.6f}  L(L e0-E_L)/M={L*(L*e0-E)/M:+.5f}")
