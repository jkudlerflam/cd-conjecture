"""Closed-book C9 checks. Only generated files in this directory are used.

Abstract two-port cells, not a nearest-neighbor realization inside a cell.
"""
import json
import time
import numpy as np
from scipy.linalg import eigh
from scipy.sparse.linalg import LinearOperator, eigsh

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.diag([1., -1.]).astype(complex)

def kron_all(xs):
    a = np.ones((1, 1), complex)
    for x in xs:
        a = np.kron(a, x)
    return a

def clifford(r):
    m = (r + 1) // 2
    gammas = []
    for k in range(m):
        gammas += [kron_all([Z]*k + [X] + [I2]*(m-k-1)),
                   kron_all([Z]*k + [Y] + [I2]*(m-k-1))]
    parity = kron_all([Z]*m)
    return gammas[:r], parity, gammas

def cell(r=2, s=0.1):
    """q=2, r=2 or 3, two physical ports; unused Pauli directions map to zero."""
    gs, parity, allgs = clifford(r)
    etas = [1j*g@parity for g in gs]
    n = len(parity)
    fs = [X, Y, Z][:r]
    bs = [np.kron(f, I2) for f in fs] + [np.kron(I2, f) for f in fs]
    ts = gs + etas
    Q = 4
    J = np.eye(Q*n, dtype=complex)
    for b, t in zip(bs, ts):
        J += s*np.kron(b.T, t)
    J /= Q
    vals, vecs = eigh(J)
    assert min(vals) > -1e-12
    # J_(i,alpha),(j,beta) = Phi(|i><j|)_(alpha,beta).
    # Kraus K_mu maps target -> input and Phi(A)=sum K_mu^dagger A K_mu.
    ks = [np.sqrt(v)*vecs[:, k].conj().reshape(Q, n)
          for k, v in enumerate(vals) if v > 1e-13]
    G = np.stack(ks, axis=1).reshape(Q*len(ks), n)
    D = len(G)
    left = [np.kron(np.kron(f, I2), np.eye(len(ks))) for f in fs]
    right = [np.kron(np.kron(I2, f), np.eye(len(ks))) for f in fs]
    P = G@G.conj().T
    H = np.eye(D)-P
    return dict(r=r, s=s, n=n, D=D, G=G, P=P, H=H, left=left,
                right=right, gammas=gs, etas=etas, parity=parity,
                all_gammas=allgs, choi_eigenvalues=vals)

def norm(a):
    return float(np.linalg.norm(a, ord=2))

def two_cell_check(c):
    G, H, n, D = [c[k] for k in ('G', 'H', 'n', 'D')]
    W = np.kron(G, G)
    def V_mat(x):
        return sum(a@x@b.T for a, b in zip(c['right'], c['left']))
    def V_vec(v):
        return V_mat(v.reshape(D, D)).reshape(-1)
    VW = np.stack([V_vec(W[:, k]) for k in range(n*n)], axis=1)
    V0 = W.conj().T@VW
    expected = c['s']**2*sum(np.kron(a, b) for a,b in zip(c['etas'], c['gammas']))
    # H0 inverse on excited complement is P Q + Q P + (Q Q)/2.
    p, q = c['P'], H
    def R_vec(v):
        x = v.reshape(D, D)
        return (p@x@q.T+q@x@p.T+0.5*q@x@q.T).reshape(-1)
    RVW = np.stack([R_vec(VW[:, k]) for k in range(n*n)], axis=1)
    H2 = -VW.conj().T@RVW
    def op(lam):
        def mv(v):
            x = v.reshape(D, D)
            return (H@x+x@H.T+lam*V_mat(x)).reshape(-1)
        return LinearOperator((D*D,D*D), matvec=mv, dtype=np.complex128)
    spectra=[]
    for lam in [0.04, 0.02, 0.01, 0.005]:
        operator = op(lam)
        raw, vectors = eigsh(operator, k=n*n+1, which='SA', tol=2e-12, maxiter=5000)
        eig_residual = max(np.linalg.norm(operator@vectors[:,i]-raw[i]*vectors[:,i]) for i in range(len(raw)))
        exact = np.sort(raw)
        approx = eigh(lam*V0+lam*lam*H2, eigvals_only=True)
        residual = float(np.max(np.abs(exact[:n*n]-approx)))
        spectra.append(dict(lam=lam, eigenpair_residual=float(eig_residual), eigenvalues=exact[:n*n].tolist(),
                            next_eigenvalue=float(exact[n*n]),
                            second_order_error=residual,
                            error_over_lam_cubed=residual/lam**3))
    result = dict(first_order_error=norm(V0-expected),
                  H2_scalar_removed_norm=norm(H2-np.trace(H2)/(n*n)*np.eye(n*n)),
                  H2_logical_parity_commutator=norm(H2@np.kron(c['parity'],c['parity'])-np.kron(c['parity'],c['parity'])@H2),
                  spectra=spectra)
    if c['r']==2:
        s=c['s']; T=expected/s**2
        analytic=-(c['r']*(1+2*s*s)/2)*np.eye(n*n)+1.5*s**4*T@T
        result['analytic_H2_error']=norm(H2-analytic)
        result['H2_eigenvalues']=eigh(H2,eigvals_only=True).tolist()
    return result

def ladder_check(r, cells=4):
    gs, p, allgs=clifford(r)
    n=len(p)
    es=[1j*g@p for g in gs]
    dim=n**cells
    T=np.zeros((dim,dim),complex)
    globals_=[]
    for j in range(cells):
        globals_.append([kron_all([p]*j+[g]+[np.eye(n)]*(cells-j-1)) for g in allgs])
    JW=np.zeros_like(T)
    for j in range(cells-1):
        for a in range(r):
            T+=kron_all([np.eye(n)]*j+[es[a],gs[a]]+[np.eye(n)]*(cells-j-2))
            JW+=1j*globals_[j][a]@globals_[j+1][a]
    ev=eigh(T,eigvals_only=True)
    out=dict(r=r,cells=cells,dim=dim,jordan_wigner_error=norm(T-JW),
             ground_multiplicity=int(np.sum(abs(ev-ev[0])<1e-9)),
             ground_energy=float(ev[0]))
    if r%2:
        out['unused_majorana_commutator_norms']=[norm(T@g[-1]-g[-1]@T) for g in globals_]
    return out

def main():
    start=time.time()
    out={}
    for r in (2,3):
        c=cell(r,0.1)
        G=c['G']; n=c['n']
        errors=[]
        for A,B in zip(c['left']+c['right'],c['gammas']+c['etas']):
            errors.append(norm(G.conj().T@A@G-c['s']*B))
        cross=max(norm(G.conj().T@a@b@G) for a in c['left'] for b in c['right'])
        out[f'cell_r{r}']=dict(N=n,D=c['D'],s=c['s'],
            choi_min=float(min(c['choi_eigenvalues'])),
            isometry_error=norm(G.conj().T@G-np.eye(n)),
            compression_error=max(errors),cross_port_product_compression=cross,
            zero_eigenvalues=n,gap=1.0,
            two_cells=two_cell_check(c))
        print(f'Completed r={r}', flush=True)
    out['ladder_r2']=ladder_check(2)
    out['ladder_r3']=ladder_check(3)
    out['runtime_seconds']=time.time()-start
    with open('checks_results.json','w') as f:
        json.dump(out,f,indent=2)
    print(json.dumps(out,indent=2))

if __name__=='__main__':
    main()
