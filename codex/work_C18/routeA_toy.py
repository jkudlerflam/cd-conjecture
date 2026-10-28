"""Literal d=64, n=2 port-idle sweep cell with a nontrivial second-order bond.
The full legal cell has dimension 1920. Four 120-dimensional Bell blocks
exhaust it, with multiplicities 1,3,3,9. This is a two-state diagonal
logical interaction, not a high-channel compiler or a CFT calculation.
"""
import json
from pathlib import Path
import numpy as np
from scipy.linalg import eigh
from scipy.sparse import csr_matrix, eye, kron as skron, block_diag
from scipy.sparse.linalg import eigsh
from clock_checks import neighbors,C,selection

I=np.eye(2); Z=np.diag([1.,-1.]); P1=np.diag([0.,1.])
SWAP=np.array([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1.]])


def main():
    n=2; c=C(n); start=('PL','I','Z','PR')
    words=[start]; prev=None
    while True:
        cand=[z for z in neighbors(words[-1]) if z!=prev]
        if not cand: break
        assert len(cand)==1
        prev=words[-1]; words.append(cand[0])
    assert len(words)==c
    Ws=[np.eye(4)]; swap_times=[]
    for t,(a,b) in enumerate(zip(words[:-1],words[1:])):
        outward=a[1].startswith('D') and b[1].startswith('A')
        u=SWAP if outward else np.eye(4)
        if outward: swap_times.append(t)
        Ws.append(u@Ws[-1])
    assert len(swap_times)==2**n
    phi=np.sqrt(2/(c+1))*np.sin(np.pi*np.arange(1,c+1)/(c+1))
    E0=-2*np.cos(np.pi/(c+1))
    # Gauged block includes the hopping clock and the right c-bit input pin.
    h=-np.diag(np.ones(c-1),1)-np.diag(np.ones(c-1),-1)-E0*np.eye(c)
    base=np.kron(h,np.eye(4))
    initial_pin=words.index(('PL','L','I','PR'))
    base[4*initial_pin:4*(initial_pin+1),4*initial_pin:4*(initial_pin+1)]+=np.kron(I,P1)
    A=I+.25*Z
    vL=np.kron(A,I); vR=np.kron(I,A)
    G0=np.column_stack([np.kron(I[:,i],I[:,0]) for i in range(2)])
    VG=np.kron(phi[:,None],G0)
    blocks={}; spectral=[]
    for left,right in ((0,0),(1,0),(0,1),(1,1)):
        H=base.copy()
        for t,w in enumerate(Ws):
            H[4*t:4*(t+1),4*t:4*(t+1)]+=w.T@(left*vL+right*vR)@w
        ev,U=eigh(H)
        multiplicity=(3 if left else 1)*(3 if right else 1)
        spectral.extend(ev.tolist()*multiplicity)
        blocks[left,right]=(H,ev,U)
    spectral=np.sort(spectral)
    assert len(spectral)==1920
    assert max(abs(spectral[:2]))<1e-12 and spectral[2]>0
    H0=blocks[0,0][0]
    residual=float(np.linalg.norm(H0@VG))
    assert residual<1e-12
    er,Ur=blocks[0,1][1:]; el,Ul=blocks[1,0][1:]
    Br=Ur.T@VG; Bl=Ul.T@VG
    T=np.einsum('ik,jl->ijkl',Br,Bl).reshape((4*c)**2,4)
    den=(er[:,None]+el[None,:]).reshape(-1)
    K2=-3*T.T@((1/den)[:,None]*T)
    # Each of the three Pauli syndrome labels gives the same contribution.
    zz=np.kron(Z,Z); zi=np.kron(Z,I); iz=np.kron(I,Z)
    coeff=[float(np.trace(K2@op)/4) for op in (np.eye(4),zi,iz,zz)]
    assert abs(coeff[-1])>1e-5
    assert np.linalg.norm(K2-np.diag(np.diag(K2)))<1e-12
    # The bare port map vanishes identically because the two Bell pairs
    # remain fixed in every ground-history slice, independent of c bits.
    # Independent finite-lambda spectra in the invariant pair-syndrome
    # subspace. Four Bell labels (a,a) times two clocks gives dimension 3600.
    X=np.array([[0.,1.],[1.,0.]])
    Y=np.array([[0.,-1j],[1j,0.]])
    bell=np.array([1.,0.,0.,1.])/np.sqrt(2)
    SB=np.column_stack([bell]+[np.kron(f,I)@bell for f in (X,Y,Z)])
    Os=[SB.conj().T@np.kron(f,I)@SB for f in (X,Y,Z)]
    vpair=sum(np.kron(o,o) for o in Os)
    ids=[0,5,10,15]
    assert np.linalg.norm(vpair[np.ix_([x for x in range(16) if x not in ids],ids)])<1e-12
    v4=vpair[np.ix_(ids,ids)].real
    finite=[]
    for lam in (.004,.002,.001):
        shifts=[]; residuals=[]
        for ii in range(2):
            for jj in range(2):
                ri=np.array([2*ii+4*t for t in range(c)])
                lj=np.array([2*jj+4*t for t in range(c)])
                hr=blocks[0,1][0][np.ix_(ri,ri)]
                hl=blocks[1,0][0][np.ix_(lj,lj)]
                free=skron(csr_matrix(h),eye(c))+skron(eye(c),csr_matrix(h))
                excited=skron(csr_matrix(hr),eye(c))+skron(eye(c),csr_matrix(hl))
                Hpair=block_diag([free,excited,excited,excited],format='csr')
                Hpair+=lam*skron(csr_matrix(v4),eye(c*c),format='csr')
                ee,vv=eigsh(Hpair,k=1,which='SA',tol=1e-11,v0=np.ones(4*c*c))
                shifts.append(float(ee[0]))
                residuals.append(float(np.linalg.norm(Hpair@vv[:,0]-ee[0]*vv[:,0])))
        error=float(np.max(np.abs(np.array(shifts)-lam**2*np.diag(K2))))
        finite.append(dict(lam=lam,dimension=4*c*c,energies=shifts,
                           second_order_error=error,error_over_lam3=error/lam**3,
                           eigenpair_residual=max(residuals)))
    delta,gamma,en,clockgap=selection(n)
    out=dict(d0=64,period_b=4,bulk_n=2,clock_vertices=c,
             legal_cell_dimension=1920,block_dimension=4*c,ground_dimension=2,
             full_legal_gap=float(spectral[2]),ground_residual=residual,
             outward_swap_edges=swap_times,initialization_time=initial_pin,
             second_order_diagonal=np.diag(K2).tolist(),
             second_order_coefficients_I_ZI_IZ_ZZ=coeff,
             first_order='zero analytically from both Bell pairs',
             selector_delta=delta,selector_mu=2-delta,
             length_excess=gamma,bare_clock_gap=clockgap,
             finite_lambda=finite,
             scope='Explicit small cell, not the unbounded-target compiler.')
    Path('routeA_toy_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))

if __name__=='__main__': main()
