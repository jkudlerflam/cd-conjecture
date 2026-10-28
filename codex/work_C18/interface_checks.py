"""Weighted hopping histories and a positive Bell-check obstruction."""
import json
import time
from pathlib import Path
import numpy as np
from scipy.linalg import eigh

I=np.eye(2); X=np.array([[0.,1.],[1.,0.]])
Y=np.array([[0.,-1j],[1j,0.]])
Z=np.diag([1.,-1.]); PAULI=[X,Y,Z]
bell=np.array([1.,0.,0.,1.])/np.sqrt(2)
Q=np.eye(4)-np.outer(bell,bell)


def kron(*xs):
    out=np.array([[1.]])
    for x in xs: out=np.kron(out,x)
    return out


def history(C):
    phi=np.sqrt(2/(C+1))*np.sin(np.arange(1,C+1)*np.pi/(C+1))
    h=-np.diag(np.ones(C-1),1)-np.diag(np.ones(C-1),-1)
    e=-2*np.cos(np.pi/(C+1))
    return phi,h-e*np.eye(C)


def no_go():
    # Physical order pL,sL,logical,sR,pR, with two independent ports.
    c=9; phi,h=history(c); d=32; N=2
    G=np.column_stack([np.kron(np.kron(bell,I[:,i]),bell) for i in range(N)])
    qL=kron(Q,I,I,I); qR=kron(I,I,I,Q)
    H=kron(h,np.eye(d)).astype(complex)
    vL=0.7+0.13*np.arange(c); vR=1.2+0.08*np.arange(c)[::-1]
    for t in range(c):
        H[t*d:(t+1)*d,t*d:(t+1)*d]+=vL[t]*qL+vR[t]*qR
    VG=np.kron(phi[:,None],G)
    ev,U=eigh(H)
    assert np.max(np.abs(ev[:N]))<1e-12 and ev[N]>0
    pL=[kron(np.eye(c),f,I,I,I,I) for f in PAULI]
    pR=[kron(np.eye(c),I,I,I,I,f) for f in PAULI]
    BsL=[U.conj().T@f@VG for f in pL]
    BsR=[U.conj().T@f@VG for f in pR]
    first=max(np.linalg.norm(VG.conj().T@f@VG) for f in pL+pR)
    resolvents=[]
    for zz in (0.,0.3,1.):
        inv=np.zeros_like(ev); mask=ev>1e-10; inv[mask]=1/(ev[mask]+zz)
        err=0.
        for Bs in (BsL,BsR):
            for i in range(3):
                for j in range(3):
                    R=Bs[i].conj().T@(inv[:,None]*Bs[j])
                    target=np.eye(N)*np.trace(R)/N if i==j else np.zeros((N,N))
                    err=max(err,float(np.linalg.norm(R-target)))
        resolvents.append(dict(z=zz,nonscalar_error=err))
    T=sum(np.einsum('ik,jl->ijkl',br,bl).reshape(len(ev)**2,N*N)
          for br,bl in zip(BsR,BsL))
    den=(ev[:,None]+ev[None,:]).reshape(-1)
    inv=np.zeros_like(den); inv[den>1e-10]=1/den[den>1e-10]
    K2=-T.conj().T@(inv[:,None]*T)
    scalar=np.trace(K2).real/(N*N)
    err=float(np.linalg.norm(K2-scalar*np.eye(N*N)))
    assert first<1e-12 and err<1e-11
    assert max(x['nonscalar_error'] for x in resolvents)<1e-11
    return dict(clock_vertices=c,cell_dimension=c*d,logical_dimension=N,
                cell_gap=float(ev[N]),first_order_norm=first,
                resolvents=resolvents,second_order_scalar=float(scalar),
                second_order_nonscalar_error=err)


def swap_gate():
    return np.array([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1.]])


def nonff():
    # Three clock slices, port p idle, SWAP(s,l), then SWAP(s,l).
    C=3; phi,h=history(C); d=8
    S=kron(I,swap_gate()); Ws=[np.eye(d),S,np.eye(d)]
    Qps=kron(Q,I)
    hbar=sum(phi[t]**2*(w.T@Qps@w) for t,w in enumerate(Ws))
    ev,V=eigh(hbar); G=V[:,:2]
    target_port=float(np.linalg.norm(G.conj().T@kron(X,I,I)@G,2))
    assert target_port>0.1
    out=[]
    for J in (10.,100.,1000.):
        HH=J*kron(h,np.eye(d)).astype(complex)
        # Work in the gauged frame; source ports are unchanged since W is idle.
        for t,w in enumerate(Ws):
            HH[t*d:(t+1)*d,t*d:(t+1)*d]+=w.T@Qps@w
        ee,VV=eigh(HH); GG=VV[:,:2]
        portnorm=float(np.linalg.norm(GG.conj().T@kron(np.eye(C),X,I,I)@GG,2))
        out.append(dict(J=J,energy_shift=float(ee[0]),ground_splitting=float(ee[1]-ee[0]),
                        port_compression_norm=portnorm))
    return dict(clock_vertices=C,cell_dimension=C*d,
                averaged_ground_energy=float(ev[0]),
                averaged_port_compression_norm=target_port,finite_J=out)


def weighted_resolvent():
    # Exact code on p,s,l. Richer positive term Qps tensor A_t,
    # not a bare Bell check. Logical syndrome dependence evades the no-go.
    C=5; phi,h=history(C); d=8
    A=[np.array([[1.+.1*t,.13], [.13,1.7-.07*t]]) for t in range(C)]
    Vblocks=[kron(Q,a) for a in A]
    Abar=sum(phi[t]**2*v for t,v in enumerate(Vblocks))
    G=np.column_stack([np.kron(bell,I[:,i]) for i in range(2)])
    VG=np.kron(phi[:,None],G)
    source=kron(np.eye(C),X,I,I)@VG
    source0=kron(X,I,I)@G
    ev0,U0=eigh(Abar); inv0=np.zeros(d); mask=ev0>1e-10
    inv0[mask]=1/(ev0[mask]+.2)
    target=source0.conj().T@((U0*inv0)@U0.conj().T)@source0
    Delta=float(ev0[2]); M=max(np.linalg.norm(v,2) for v in Vblocks)
    entries=[]
    for J in (100.,1000.,10000.):
        H=J*kron(h,np.eye(d))
        for t,v in enumerate(Vblocks): H[t*d:(t+1)*d,t*d:(t+1)*d]+=v
        ev,U=eigh(H); inv=np.zeros_like(ev); mask=ev>1e-9
        inv[mask]=1/(ev[mask]+.2)
        R=source.conj().T@((U*inv)@U.conj().T)@source
        Gamma=J*4*np.sin(3*np.pi/(2*(C+1)))*np.sin(np.pi/(2*(C+1)))
        delta=M*M/(Gamma-M+.2)
        bound=delta/((Delta+.2)*(Delta+.2-delta))
        error=float(np.linalg.norm(R-target,2))
        assert error<bound and np.linalg.norm(H@VG)<1e-10
        entries.append(dict(J=J,resolvent_error=error,bound=bound,
                            ground_residual=float(np.linalg.norm(H@VG))))
    return dict(clock_vertices=C,cell_dimension=C*d,
                nonuniform_weights=(phi**2).tolist(),
                target_resolvent_nonscalar_norm=float(np.linalg.norm(target-np.trace(target)/2*np.eye(2))),
                finite_J=entries)

if __name__=='__main__':
    start=time.monotonic()
    out=dict(positive_bell_no_go=no_go(),nonfrustration_free_example=nonff(),
             richer_syndrome_terms=weighted_resolvent())
    out['seconds']=time.monotonic()-start
    Path('interface_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
