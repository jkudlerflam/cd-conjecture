"""Scalar (two-token resolvent) formula for the lambda^2 g^2 logical Bloch blocks of the 9-letter bond,
in arbitrary precision.  Units: lambda^2 alpha^4 g^2 = 1.
Token graph per species s: express vertices X_0..X_n (index j), data vertices D_0..D_n (index n+1+j).
"""
import mpmath as mp

def Hs(n, a, D, phi, td, psi, sw=1):
    N = 2*n+2; H = mp.zeros(N, N)
    for j in range(n+1):
        H[j,j] = a + D; H[n+1+j, n+1+j] = a
    for j in range(n):
        H[j+1, j] += -mp.expjpi(phi/mp.pi) if phi else mp.mpf(-1); H[j, j+1] += -mp.expjpi(-phi/mp.pi) if phi else mp.mpf(-1)
        H[n+2+j, n+1+j] += -1; H[n+1+j, n+2+j] += -1
    H[0, n+1] += -sw; H[n+1, 0] += -sw                       # switch D_0 <-> X_0
    if td:
        H[n, 2*n+1] += -td*mp.expj(psi); H[2*n+1, n] += -td*mp.expj(-psi)   # closing link D_n -> X_n
    return H

class Res:
    """R_st B for the Sylvester operator X -> H_s X + X H_t^T"""
    def __init__(self, Hs_, Ht_):
        self.es, self.Us = mp.eigh(Hs_); self.et, self.Ut = mp.eigh(Ht_)
        self.UsH = self.Us.H; self.Utc = self.Ut.conjugate(); self.UtT = self.Ut.T
        N = Hs_.rows
        self.den = mp.matrix(N, N)
        for i in range(N):
            for j in range(N): self.den[i,j] = 1/(self.es[i] + self.et[j])
    def __call__(self, B):
        C = self.UsH * B * self.Utc
        N = C.rows
        for i in range(N):
            for j in range(N): C[i,j] *= self.den[i,j]
        return self.Us * C * self.UtT

def blocks(n, a=4, D=30, phiA=None, phiB=0, td=0, psiA=None, psiB=None, sw=1, signs=(-1,1,1)):
    """returns T (left<-right), muL, muR as mp matrices, in units lambda^2 alpha^4 g^2."""
    if phiA is None: phiA = mp.pi/(2*(n+1))
    if psiA is None: psiA = mp.pi
    if psiB is None: psiB = mp.pi/2
    H = {'A': Hs(n, a, D, phiA, td, psiA, sw), 'B': Hs(n, a, D, phiB, td, psiB, sw)}
    R = {(s,t): Res(H[s], H[t]) for s in 'AB' for t in 'AB'}
    N = 2*n+2; srcL, srcR = n, n+1                      # X_n and D_0
    B0 = mp.zeros(N, N); B0[srcL, srcR] = 1
    X = {st: R[st](B0) for st in [('A','B'),('B','A')]}
    Dv = lambda j: n+1+j
    def col(M, c):
        Y = mp.zeros(N, N)
        for l in range(N): Y[l,c] = M[l,c]
        return Y
    def row(M, r):
        Y = mp.zeros(N, N)
        for c in range(N): Y[r,c] = M[r,c]
        return Y
    def ip(A_, B_):
        return mp.fsum(mp.conj(A_[i,j])*B_[i,j] for i in range(N) for j in range(N))
    T = mp.zeros(n, n); muL = mp.zeros(n, n); muR = mp.zeros(n, n)
    ph = mp.mpc(0,1)                                    # e^{i pi/4} * conj(e^{-i pi/4})
    # inter-cell: left emits at j (B->A), right absorbs at j' (A->B); start (B,A), end (A,B)
    rightfirst = [R[('B','B')](col(X[('B','A')], Dv(jp))) for jp in range(n)]
    leftfirst  = [R[('A','A')](row(X[('B','A')], Dv(j))) for j in range(n)]
    for j in range(n):
        for jp in range(n):
            t1 = ip(row(X[('A','B')], Dv(j)), row(rightfirst[jp], Dv(j)))
            t2 = ip(col(X[('A','B')], Dv(jp)), col(leftfirst[j], Dv(jp)))
            T[j,jp] = -signs[0]*ph*(t1 + t2)
    # intra-cell, left cell: (i) start (B,A): emit j then absorb j'; (ii) start (A,B): absorb j' then emit j
    emitL = [R[('A','A')](row(X[('B','A')], Dv(j))) for j in range(n)]
    absL  = [R[('B','B')](row(X[('A','B')], Dv(jp))) for jp in range(n)]
    emitR = [R[('A','A')](col(X[('A','B')], Dv(j))) for j in range(n)]
    absR  = [R[('B','B')](col(X[('B','A')], Dv(jp))) for jp in range(n)]
    for j in range(n):
        for jp in range(n):
            a1 = ip(row(X[('B','A')], Dv(jp)), row(emitL[j], Dv(jp)))
            a2 = ip(row(X[('A','B')], Dv(j)), row(absL[jp], Dv(j)))
            muL[j,jp] = -signs[1]*(a2 - a1)
            b1 = ip(col(X[('A','B')], Dv(jp)), col(emitR[j], Dv(jp)))
            b2 = ip(col(X[('B','A')], Dv(j)), col(absR[jp], Dv(j)))
            muR[j,jp] = -signs[2]*(b2 - b1)
    return T, muL, muR

if __name__ == '__main__':
    import numpy as np, sys
    from twocell9 import analyse, one_particle
    mp.mp.dps = 40
    for n in (2, 3):
        a, D, td = 2.6, 3.0, 0.07
        p = dict(lam=1.0, g=1.0, Delta=a-2*np.cos(np.pi/(n+2)), D=D, nu=0.0, f=0.0, phiA=np.pi/(2*(n+1)), td=td)
        r = analyse(n, p); h, vacE = one_particle(r['h2'], n, r['bitlist'])
        al4 = (2/(n+2))**2*np.sin(np.pi/(n+2))**4
        T, muL, muR = blocks(n, a=mp.mpf(a), D=mp.mpf(D), td=mp.mpf(td))
        c = lambda M: np.array(M.tolist(), dtype=complex)*al4
        for name, bf, sc in [('T', h[:n,n:], c(T)), ('muL', h[:n,:n], c(muL)), ('muR', h[n:,n:], c(muR))]:
            ratio = bf/np.where(np.abs(sc)>1e-30, sc, 1)
            print(n, name, 'max|bf|=%.3e'%np.abs(bf).max(), 'max|bf-sc|=%.2e'%np.abs(bf-sc).max(), 'max|bf+sc|=%.2e'%np.abs(bf+sc).max())
            if np.abs(bf-sc).max() > 1e-9*np.abs(bf).max() and n==2: print('   ratio bf/sc:\n', np.round(ratio,4))
