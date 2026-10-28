"""Checks for resumed C8: exact-encoding leakage and flat-parent second order.
Only numpy/scipy, no external data. Typical runtime a few seconds with one BLAS thread.
"""
import json
import numpy as np
from scipy.linalg import eigh, svd
from c8_checks import clifford, tensor


def odd_basis(q):
    """Hermitian off-diagonal blocks under diag(+1,...,+1,-1,...,-1)."""
    out=[]
    for a in range(q//2):
        for b in range(q//2,q):
            f=np.zeros((q,q),complex)
            f[a,b]=f[b,a]=np.sqrt(q/2)
            out.append(f)
            f=np.zeros((q,q),complex)
            f[a,b]=-1j*np.sqrt(q/2)
            f[b,a]=1j*np.sqrt(q/2)
            out.append(f)
    return out


def joint_code(q,m):
    fs=odd_basis(q)[:m]
    assert len(fs)==m and m%2==0
    gs,rs,p=clifford(m//2)
    D=len(p); k=q*q
    alpha=1/(4*m*np.sqrt(q))
    J=np.eye(k*D,dtype=complex)
    for f,g,r in zip(fs,gs,rs):
        J+=alpha*(np.kron(np.kron(f.T,np.eye(q)),g)
                  +np.kron(np.kron(np.eye(q),f.T),r))
    J/=k
    ev,u=eigh(J)
    assert ev.min()>0
    S=(u*np.sqrt(ev)[None,:])@u.conj().T
    V=S.reshape(k,D,k*D).transpose(0,2,1).conj().reshape(k*k*D,D)
    vt=V.reshape(k,k*D,D)
    K_L=[np.einsum('ij,jeb->ieb',np.kron(f,np.eye(q)),vt).reshape(-1,D) for f in fs]
    K_R=[np.einsum('ij,jeb->ieb',np.kron(np.eye(q),f),vt).reshape(-1,D) for f in fs]
    return fs,gs,rs,p,alpha,V,K_L,K_R


def run_case(q,m):
    fs,gs,rs,p,alpha,V,kl,kr=joint_code(q,m)
    D=len(p)
    all_k=kl+kr; all_o=gs+rs
    err_iso=np.linalg.norm(V.conj().T@V-np.eye(D))
    err_comp=max(np.linalg.norm(V.conj().T@k-alpha*o) for k,o in zip(all_k,all_o))
    err_mom=0.
    for a,ka in enumerate(all_k):
        for b,kb in enumerate(all_k):
            same_side=(a<m)==(b<m)
            target=np.eye(D) if same_side and a==b else np.zeros((D,D))
            err_mom=max(err_mom,np.linalg.norm(ka.conj().T@kb-target))
    # Leakage identity for physically commuting opposite endpoints.
    la=kl[0]-V@(V.conj().T@kl[0]); ra=kr[0]-V@(V.conj().T@kr[0])
    left=alpha*gs[0]; right=alpha*rs[0]
    comm=right@left-left@right
    leak_identity=la.conj().T@ra-ra.conj().T@la
    err_leak=np.linalg.norm(comm-leak_identity)
    # Exact finite-dimensional reduction of the needed excited transition vectors.
    syn=np.concatenate([k-V@(V.conj().T@k) for k in all_k],axis=1)
    U,s,_=svd(syn,full_matrices=False)
    U=U[:,s>1e-11]
    frame=np.concatenate([V,U],axis=1)
    local_dim=frame.shape[1]
    G=np.zeros((local_dim,D),complex); G[:D]=np.eye(D)
    l=[frame.conj().T@k for k in kl]
    r=[frame.conj().T@k for k in kr]
    W1=sum(tensor([a,b,G]) for a,b in zip(r,l))
    W2=sum(tensor([G,a,b]) for a,b in zip(r,l))
    W=W1+W2
    exc=(np.arange(local_dim)>=D).astype(float)
    denom=(exc[:,None,None]+exc[None,:,None]+exc[None,None,:]).reshape(-1)
    inv=np.zeros_like(denom); inv[denom>0]=1/denom[denom>0]
    actual=-W.conj().T@(inv[:,None]*W)
    eye=np.eye(D)
    t1=sum(tensor([a,b,eye]) for a,b in zip(rs,gs))
    t2=sum(tensor([eye,a,b]) for a,b in zip(rs,gs))
    scalar=-2*m*(.5+alpha**2)*np.eye(D**3)
    same=1.5*alpha**4*(t1@t1+t2@t2)
    cross=alpha**4*(t1@t2+t2@t1)
    predicted=scalar+same+cross
    actual_cross=-(W1.conj().T@(inv[:,None]*W2)+W2.conj().T@(inv[:,None]*W1))
    G3=tensor([G,G,G])
    err_first=np.linalg.norm(G3.conj().T@W-alpha**2*(t1+t2))
    err_second=np.linalg.norm(actual-predicted)
    err_cross=np.linalg.norm(actual_cross-cross)
    assert max(err_iso,err_comp,err_mom,err_leak,err_first,err_second,err_cross)<1e-10
    return dict(q=q,m=m,logical_dimension=D,physical_block_dimension=len(V),
                active_local_dimension=local_dim,alpha=alpha,
                isometry_error=float(err_iso),compression_error=float(err_comp),
                second_moment_error=float(err_mom),leakage_identity_error=float(err_leak),
                leakage_product=float(np.linalg.norm(la,2)*np.linalg.norm(ra,2)),
                leakage_product_lower_bound=alpha**2,
                first_order_error=float(err_first),second_order_error=float(err_second),
                adjacent_bond_error=float(err_cross),
                adjacent_term_norm=float(np.linalg.norm(cross,2)),
                mass_kernel_coefficient=1.5+2*np.cos(np.pi))


def main():
    result={'status':'Finite-dimensional consistency checks, not CFT numerics',
            'cases':[run_case(2,2),run_case(3,4)]}
    with open('resume_checks.json','w') as f:
        json.dump(result,f,indent=2)
    print(json.dumps(result,indent=2))




def spectral_parent_check(q=3,m=4):
    """Two-energy parent with parity and a reversed mass-channel sign."""
    from scipy.linalg import null_space
    fs,gs,rs,p,alpha,V,kl,kr=joint_code(q,m)
    D=len(p); dim=len(V); a2=alpha**2
    Os=gs+rs; ks=kl+kr
    syn=np.concatenate([k-V@(V.conj().T@k) for k in ks],axis=1)
    gram=syn.conj().T@syn
    ns=len(gram)
    ev,U=eigh(gram)
    ginv=(U/np.sqrt(ev)[None,:])@U.conj().T
    source=syn@ginv
    UL=np.concatenate(gs,axis=1); UR=np.concatenate(rs,axis=1)
    energies=[1.,3.]
    ps=[.5,.5]; ss=[-a2/2,-a2/2]; ts=[0.,-a2]
    weights=[]
    for pe,se,te in zip(ps,ss,ts):
        w=pe*np.eye(ns,dtype=complex)
        w[:m*D,:m*D]+=se*UL.conj().T@UL
        w[m*D:,m*D:]+=se*UR.conj().T@UR
        w[:m*D,m*D:]+=te*UL.conj().T@UR
        w[m*D:,:m*D]+=te*UR.conj().T@UL
        weights.append(w)
    roots=[]
    for w in weights:
        e,u=eigh(w); assert e.min()>0
        roots.append((u*np.sqrt(e)[None,:])@u.conj().T)
    raw=np.vstack(roots)
    W=raw@ginv
    # Explicit microscopic grading on ports and the Choi environment.
    z=np.diag([1]*(q//2)+[-1]*(q-q//2))
    omega=tensor([z,z,z,z,p])
    omega_diag=np.real(np.diag(omega))
    source_parity=np.tile(-np.diag(p).real,2*m)
    dilation_parity=np.tile(source_parity,2)
    eigenvectors=np.zeros((dim,2*ns),complex)
    # Complete the isometries separately in each parity sector.
    for sign in [-1,1]:
        phys=np.flatnonzero(omega_diag==sign)
        src=np.flatnonzero(source_parity==sign)
        rows=np.flatnonzero(dilation_parity==sign)
        code=np.flatnonzero(np.diag(p).real==sign)
        sframe=source[np.ix_(phys,src)]
        gf=V[np.ix_(phys,code)]
        complement=null_space(np.concatenate([gf,sframe],axis=1).conj().T)
        F=np.concatenate([sframe,complement[:,:len(rows)-len(src)]],axis=1)
        wp=W[np.ix_(rows,src)]
        unitary=np.concatenate([wp,null_space(wp.conj().T)],axis=1)
        eigenvectors[np.ix_(phys,rows)]=F@unitary.conj().T
    excited_energies=np.repeat(energies,ns)
    P=V@V.conj().T
    H0=(eigenvectors*excited_energies[None,:])@eigenvectors.conj().T
    H0+=3*(np.eye(dim)-P-eigenvectors@eigenvectors.conj().T)
    err_kernel=np.linalg.norm(H0@V)
    err_parity=np.linalg.norm(H0@omega-omega@H0)
    err_weights=max(np.linalg.norm(syn.conj().T@eigenvectors[:,e*ns:(e+1)*ns]
                                   @eigenvectors[:,e*ns:(e+1)*ns].conj().T@syn-weights[e])
                    for e in range(2))
    # Use physical transition amplitudes in the explicit excited eigenbasis.
    frame=np.concatenate([V,eigenvectors],axis=1)
    KK=[frame.conj().T@k for k in ks]
    L=KK[:m]; R=KK[m:]
    Wbond=sum(np.kron(a,b) for a,b in zip(R,L))
    ee=np.concatenate([np.zeros(D),excited_energies])
    den=(ee[:,None]+ee[None,:]).reshape(-1)
    inv=np.zeros_like(den); inv[den>0]=1/den[den>0]
    h2=-Wbond.conj().T@(inv[:,None]*Wbond)
    T=sum(np.kron(a,b) for a,b in zip(rs,gs))
    scalar=-m*(7/24+3*a2/4)
    A=25/24*alpha**4
    B=1/3*alpha**4
    target=scalar*np.eye(D**2)+A*T@T
    err_bond=np.linalg.norm(h2-target)
    local_inv=np.zeros_like(ee); local_inv[D:]=1/ee[D:]
    adj=np.zeros((D**3,D**3),complex)
    for a in range(m):
        for b in range(m):
            middle=L[a].conj().T@(local_inv[:,None]*R[b])
            adj-=a2*tensor([rs[a],middle+middle.conj().T,gs[b]])
    T1=np.kron(T,np.eye(D)); T2=np.kron(np.eye(D),T)
    err_adj=np.linalg.norm(adj-B*(T1@T2+T2@T1))
    assert max(err_kernel,err_parity,err_weights,err_bond,err_adj)<1e-9
    np.savez_compressed('spectral_parent_q3_m4.npz',H0=H0,V=V,
                        port_operators=np.array(fs),grading=omega,
                        weights=np.array(weights),energies=energies,
                        second_order_bond=h2,second_order_adjacent=adj)
    return dict(q=q,m=m,dimension=dim,alpha=alpha,
                minimum_weight_eigenvalue=float(min(eigh(w,eigvals_only=True).min() for w in weights)),
                gap=float(eigh(H0,eigvals_only=True)[D]),
                kernel_error=float(err_kernel),parity_error=float(err_parity),
                weight_error=float(err_weights),bond_error=float(err_bond),
                adjacent_error=float(err_adj),
                same_bond_coefficient_over_alpha4=A/alpha**4,
                adjacent_coefficient_over_alpha4=B/alpha**4,
                mass_kernel_over_alpha4=(A-2*B)/alpha**4)


if __name__=='__main__':
    main()
    spectral=spectral_parent_check()
    with open('spectral_parent_checks.json','w') as f:
        json.dump(spectral,f,indent=2)
    print(json.dumps(spectral,indent=2))
