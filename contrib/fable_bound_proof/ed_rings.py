"""Exact diagonalization of small rings: Casimir coefficient A from E_0(L), and the two-site
ring-vs-line distinguishability D_2(L) = L^2 ||rho_2^inf - rho_2^L||_1.
XXZ at Delta = 0.5 and 1 (d = 2, L up to 18) and the SU(3) swap chain (d = 3, L = 6, 9, 12).
rho_2^inf is estimated by an elementwise fit rho_2^L = rho_inf + X/L^2 + Y/L^4 over the largest L
(at Delta = 1 marginal logarithms make the fit rough; the exact answer there is known).
"""
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spla, itertools, sys

def ring_gs(d, hbond, L):
    dim=d**L
    H=sp.csr_matrix((dim,dim))
    hb=sp.csr_matrix(hbond)
    for i in range(L):
        j=(i+1)%L
        if j>i:
            op=sp.kron(sp.kron(sp.identity(d**i),hb),sp.identity(d**(L-i-2)),format='csr')
        else:
            # bond (L-1,0): permute so that sites 0 and L-1 are adjacent: use swap-based construction
            # build via reshaping: op acts on indices (s0, s_{L-1})
            idx=np.arange(dim).reshape([d]*L)
            perm=np.moveaxis(idx,[L-1,0],[0,1]).reshape(-1)   # new order: s_{L-1}, s_0, rest
            Pm=sp.csr_matrix((np.ones(dim),(np.arange(dim),perm)),shape=(dim,dim))
            op0=sp.kron(hb,sp.identity(d**(L-2)),format='csr')
            op=Pm.T@op0@Pm
        H=H+op
    w,v=spla.eigsh(H,k=1,which='SA',tol=1e-11)
    return w[0],v[:,0]

def rho2(psi,d,L):
    t=psi.reshape(d,d,d**(L-2))
    return np.einsum('abk,cdk->abcd',t,t.conj()).reshape(d*d,d*d)

sx=np.array([[0,1],[1,0]])/2; sy=np.array([[0,-1j],[1j,0]])/2; sz=np.diag([.5,-.5])
def xxz(D): return (np.kron(sx,sx)+np.kron(sy,sy)+D*np.kron(sz,sz)).real

results={}
for name,d,h,Ls,e0,Aexact in [
    ("XXZ Delta=0.5",2,xxz(0.5),[8,10,12,14,16,18], None, None),
    ("Heisenberg",2,xxz(1.0),[8,10,12,14,16,18], 0.25-np.log(2), np.pi/2),
    ("SU(3) swap",3,None,[6,9,12], None, 4*np.pi/3)]:
    if h is None:
        h=np.zeros((9,9))
        for a in range(3):
            for b in range(3): h[3*a+b,3*b+a]=1
    print("==",name)
    data=[]
    for L in Ls:
        E,psi=ring_gs(d,h,L); r=rho2(psi,d,L)
        data.append((L,E,r)); print(f"  L={L:3d} E0/L={E/L:.10f}")
    # fit E0/L = e0 - a/L^2 (+ b/L^4)
    Ls_=np.array([x[0] for x in data],float); eL=np.array([x[1]/x[0] for x in data])
    X=np.vstack([np.ones_like(Ls_),Ls_**-2,Ls_**-4]).T
    coef,*_=np.linalg.lstsq(X[-4:],eL[-4:],rcond=None)
    e0fit=coef[0]; Afit=-coef[1]*6/np.pi
    print(f"  fit: e0={e0fit:.8f}  A=-6*coef/pi={Afit:.5f}" + (f"  (exact e0={e0:.8f})" if e0 else "") + (f"  (exact A={Aexact:.5f})" if Aexact else ""))
    # rho_2^inf by elementwise fit over last 4 sizes
    R=np.array([x[2] for x in data])
    cf=np.linalg.lstsq(X[-4:],R[-4:].reshape(len(data[-4:]),-1),rcond=None)[0]
    rinf=cf[0].reshape(r.shape)
    print("  L   D_2(L)=L^2||rho_inf-rho_L||_1   L^2 Tr[h(rho_inf-rho_L)] 6/pi")
    for L,E,r in data:
        dr=rinf-r; tn=np.sum(np.abs(np.linalg.eigvalsh(dr)))
        print(f"  {L:3d}  {L**2*tn:10.5f}  {L**2*np.trace(h@dr).real*6/np.pi:10.5f}")
    # also the model-free pairwise estimate: ||rho_L - rho_{L'}||_1 * 1/(1/L^2-1/L'^2)
    print("  pairwise (L,L'):  ||rho_L-rho_L'||_1/(L^-2-L'^-2)")
    for (L,E,r),(L2,E2,r2) in zip(data[:-1],data[1:]):
        tn=np.sum(np.abs(np.linalg.eigvalsh(r-r2))); print(f"   ({L},{L2}) {tn/(L**-2-L2**-2):10.5f}")
