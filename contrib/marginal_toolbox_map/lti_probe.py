import numpy as np, cvxpy as cp, scipy.sparse as sp, scipy.sparse.linalg as sla, sys
X=np.array([[0,1],[1,0]],complex);Y=np.array([[0,-1j],[1j,0]]);Z=np.diag([1.,-1]).astype(complex);I2=np.eye(2)
def ptrace_first(n):  # returns function for Tr_1 and Tr_n as cvxpy expr
    pass
def lti(h,n):
    D=2**n; r=cp.Variable((D,D),hermitian=True)
    cons=[r>>0, cp.real(cp.trace(r))==1]
    left=cp.partial_trace(r,[2,2**(n-1)],axis=0); right=cp.partial_trace(r,[2**(n-1),2],axis=1)
    cons.append(left==right)
    H=np.kron(h,np.eye(2**(n-2)))
    p=cp.Problem(cp.Minimize(cp.real(cp.trace(H@r))),cons); p.solve(solver='CLARABEL')
    return p.value
def ring(h,L):
    hs=sp.csr_matrix(h); E=None; H=sp.csr_matrix((2**L,2**L),dtype=complex)
    for x in range(L-1):
        H=H+sp.kron(sp.kron(sp.identity(2**x),hs),sp.identity(2**(L-x-2)))
    # wrap bond (L-1,0): permute
    h4=h.reshape(2,2,2,2)
    hw=np.zeros((2**L,2**L),complex) if L<=10 else None
    # build wrap via operator basis decomposition
    P=[I2,X,Y,Z]
    for a in range(4):
        for b in range(4):
            c=np.trace(np.kron(P[a],P[b]).conj().T@h)/4
            if abs(c)<1e-14: continue
            op=sp.kron(sp.kron(sp.csr_matrix(P[b]),sp.identity(2**(L-2))),sp.csr_matrix(P[a]))
            H=H+c*op
    if L<=8: return np.linalg.eigvalsh(H.toarray()).min()/L
    return sla.eigsh(H,k=1,which='SA')[0][0].real/L
rng=np.random.default_rng(1)
def rand_bond(complex_=True):
    A=rng.normal(size=(4,4))+(1j*rng.normal(size=(4,4)) if complex_ else 0); return (A+A.conj().T)/2
psi=np.array([0,1,1j,0])/np.sqrt(2)
bonds={'Heis':(np.kron(X,X)+np.kron(Y,Y)+np.kron(Z,Z))/4,'XX':(np.kron(X,X)+np.kron(Y,Y))/4,
 'TFI1':-np.kron(Z,Z)-0.5*(np.kron(X,I2)+np.kron(I2,X)),'chiral':-np.outer(psi,psi.conj()),
 'rnd1':rand_bond(),'rnd2':rand_bond(),'rnd3':rand_bond(False)}
nmax=int(sys.argv[1]) if len(sys.argv)>1 else 6; Lmax=int(sys.argv[2]) if len(sys.argv)>2 else 16
for name,h in bonds.items():
    rings={L:ring(h,L) for L in range(2,Lmax+1)}
    out=[]
    for n in range(2,nmax+1):
        v=lti(h,n); m=min(rings[L] for L in range(n,Lmax+1)); Lst=min(range(n,Lmax+1),key=lambda L:rings[L])
        out.append(f"n={n}: LTI={v:.6f} minRing(L>=n)={m:.6f}@L={Lst} diff={m-v:.2e}")
    print(name,'| E_L/L:',' '.join(f"{L}:{rings[L]:.5f}" for L in range(2,Lmax+1)))
    for o in out: print('   ',o)
    sys.stdout.flush()
