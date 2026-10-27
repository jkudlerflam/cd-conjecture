"""Two-port spectral-weight synthesis of an arbitrary rank-one logical bond.

This is an abstract D=40 cell, not a local cell implementation.
"""
import json
import time
import numpy as np
from scipy.linalg import eigh, null_space
from scipy.sparse.linalg import LinearOperator, eigsh
from checks import I2, X, Y, Z, norm


def main():
    start=time.time()
    q=2; n=2; r=3; k=3; m=2*r*n; D=40; aux=D//(q*q)
    es=np.array([1.,4.,16.]); p=np.ones(k)/k
    C=1/(es[:,None]+es[None,:])
    z=null_space(np.array([np.ones(k),C@p])).reshape(k)
    z/=np.sqrt(z@C@z)
    eps=0.5*np.min(p/np.abs(z))
    G=np.zeros((D,n),complex)
    for a in range(q):
        for b in range(q):
            for i in range(n):
                G[(a*q+b)*aux+((a*q+b)*n+i),i]=1/q
    left=[np.kron(np.kron(f,I2),np.eye(aux)) for f in (X,Y,Z)]
    right=[np.kron(np.kron(I2,f),np.eye(aux)) for f in (X,Y,Z)]
    S=np.hstack([a@G for a in left+right])
    weights=[]
    for e in range(k):
        M=p[e]*np.eye(m,dtype=complex)
        M[:n,:n]-=eps*z[e]*Z
        M[r*n:(r+1)*n,r*n:(r+1)*n]+=eps*z[e]*X
        weights.append(M)
    roots=[]
    for M in weights:
        vs,us=eigh(M)
        roots.append((us*np.sqrt(vs))@us.conj().T)
    W=np.vstack(roots)
    Wfull=np.hstack([W,null_space(W.conj().T)])
    Uexc=np.hstack([S,null_space(np.hstack([G,S]).conj().T)])
    J=Uexc[:,:k*m]@Wfull.conj().T
    projectors=[J[:,e*m:(e+1)*m]@J[:,e*m:(e+1)*m].conj().T for e in range(k)]
    remainder=np.eye(D)-G@G.conj().T-J@J.conj().T
    H=sum(e*P for e,P in zip(es,projectors))+32*remainder
    H=(H+H.conj().T)/2
    E,U=eigh(H)
    den=E[:,None]+E[None,:]
    inv=np.zeros_like(den)
    inv[den>1e-8]=1/den[den>1e-8]
    def v_action(v):
        x=v.reshape(D,D)
        return sum(a@x@b.T for a,b in zip(right,left)).reshape(-1)
    def r_action(v):
        x=U.conj().T@v.reshape(D,D)@U.conj()
        return (U@(inv*x)@U.T).reshape(-1)
    GG=np.kron(G,G)
    VG=np.stack([v_action(GG[:,i]) for i in range(n*n)],axis=1)
    H1=GG.conj().T@VG
    H2=-VG.conj().T@np.stack([r_action(VG[:,i]) for i in range(n*n)],axis=1)
    a=-r*(p@C@p)
    intended=a*np.eye(n*n)+eps**2*np.kron(X,Z)
    spectra=[]
    for lam in [.04,.02,.01,.005]:
        def mv(v):
            x=v.reshape(D,D)
            return (H@x+x@H.T).reshape(-1)+lam*v_action(v)
        op=LinearOperator((D*D,D*D),matvec=mv,dtype=np.complex128)
        raw,vectors=eigsh(op,k=n*n+1,which='SA',tol=1e-12,maxiter=5000)
        eig_residual=max(np.linalg.norm(op@vectors[:,i]-raw[i]*vectors[:,i]) for i in range(len(raw)))
        exact=np.sort(raw)
        pred=eigh(lam*H1+lam**2*H2,eigvals_only=True)
        err=float(max(abs(exact[:n*n]-pred)))
        spectra.append(dict(lam=lam,eigenpair_residual=float(eig_residual),eigenvalues=exact[:n*n].tolist(),
                            next_eigenvalue=float(exact[n*n]),second_order_error=err,
                            error_over_lam_cubed=err/lam**3,
                            exact_bandwidth=float(exact[n*n-1]-exact[0]),
                            predicted_bandwidth=2*lam**2*eps**2))
    out=dict(N=n,D=D,port_q=q,levels=es.tolist(),eps=float(eps),
             exchange_coefficient=float(eps**2),scalar_coefficient=float(a),
             min_weight_eigenvalue=float(min(min(eigh(M,eigvals_only=True)) for M in weights)),
             G_isometry_error=norm(G.conj().T@G-np.eye(n)),
             syndrome_isometry_error=norm(S.conj().T@S-np.eye(m)),
             weight_sum_error=norm(sum(weights)-np.eye(m)),
             weight_realization_error=max(norm(S.conj().T@P@S-M) for P,M in zip(projectors,weights)),
             H_ground_error=norm(H@G),first_excited_energy=float(E[n]),
             first_order_norm=norm(H1),second_order_synthesis_error=norm(H2-intended),
             spectra=spectra,runtime_seconds=time.time()-start)
    with open('spectral_cell_results.json','w') as f:
        json.dump(out,f,indent=2)
    print(json.dumps(out,indent=2))

if __name__=='__main__':
    main()
