"""Exact local massive-source tails in an onsite-dressed XX chain."""
import json
from pathlib import Path
import numpy as np
from scipy.linalg import eigh,expm,norm
from scipy.integrate import quad
from local_dressing_checks import embed

I=np.eye(2);X=np.array([[0.,1],[1,0]]);Y=np.array([[0,-1j],[1j,0]]);Z=np.diag([1.,-1])


def direct_check(L,J=.2,Delta=1.,theta=.1):
    U=expm(.5j*theta*np.kron(Z,Y))
    uall=np.array([[1.]])
    for _ in range(L):uall=np.kron(uall,U)
    Xa=np.kron(I,X);Za=np.kron(I,Z);Zs=np.kron(Z,I)
    err=norm(U@Xa@U.conj().T-np.cos(theta)*Xa-np.sin(theta)*Zs@Za,2)
    assert err<1e-14
    Xs=np.kron(X,I);Ys=np.kron(Y,I)
    h=-J/2*(np.kron(Xs,Xs)+np.kron(Ys,Ys))
    A=Delta/2*(np.eye(4)-Za)
    base=np.zeros((4**L,4**L),complex)
    literal=np.zeros_like(base)
    def arbitrary(ops):
        out=np.array([[1.]])
        for j in range(L):out=np.kron(out,ops.get(j,np.eye(4)))
        return out
    for j in range(L):
        base+=embed(A,j,1,L,4)
        literal+=embed(U.conj().T@A@U,j,1,L,4)
        for O in (Xs,Ys):
            base-=J/2*arbitrary({j:O,(j+1)%L:O})
            op=U.conj().T@O@U
            literal-=J/2*arbitrary({j:op,(j+1)%L:op})
    transformed=uall.conj().T@base@uall
    herr=float(norm(literal-transformed,2));assert herr<1e-13
    charge=sum(embed(Zs,j,1,L,4) for j in range(L))
    qerr=float(norm(literal@charge-charge@literal,2));assert qerr<1e-13
    e,v=eigh(literal);ground=v[:,0]
    obs=embed(Xa,0,1,L,4)
    coeff=abs(v.conj().T@obs@ground)**2
    ebase,vbase=eigh(base);gbase=vbase[:,0]
    zs=embed(Zs,0,1,L,4)
    zcoeff=abs(vbase.conj().T@zs@gbase)**2
    # Choose even L so XX has a unique half-filled vacuum here.
    rows=[]
    for tau in (.1,1.,10.):
        val=float(coeff@np.exp(-tau*(e-e[0]))-abs(np.vdot(ground,obs@ground))**2)
        zval=float(zcoeff@np.exp(-tau*(ebase-ebase[0]))-abs(np.vdot(gbase,zs@gbase))**2)
        expected=np.cos(theta)**2*np.exp(-Delta*tau)+np.sin(theta)**2*zval
        assert abs(val-expected)<1e-12
        rows.append(dict(time=tau,direct=val,predicted=float(expected),error=abs(val-expected)))
    source_matrix=v.conj().T@obs@v
    local_v=Delta*abs(np.sin(theta/2))+2*J
    windows=[]
    for width in (J,2*J,.5*Delta):
        selected=np.flatnonzero(e<=e[0]+width+1e-12)
        compression=float(norm(source_matrix[np.ix_(selected,selected)],2))
        weight=float(np.sum(coeff[selected]))
        bound=(width+2*local_v)/Delta
        assert compression<=bound+1e-12 and weight<=bound**2+1e-12
        windows.append(dict(width=width,states=len(selected),source_norm=compression,
                            norm_bound=bound,ground_weight=weight,weight_bound=bound**2))
    return dict(length=L,dimension=4**L,J=J,theta=theta,
                hamiltonian_error=herr,charge_error=qerr,
                observable_identity_error=float(err),correlations=rows,
                true_spectral_windows=windows)


def thermodynamic(J=.2,Delta=1.,theta=.1):
    rows=[]
    for tau in (1.,10.,30.,100.,300.,1000.):
        a=quad(lambda x:np.exp(-2*J*tau*np.sin(x)),0,np.pi/2,
               epsabs=1e-14,epsrel=1e-12)[0]/np.pi
        slow=4*np.sin(theta)**2*a*a
        fast=np.cos(theta)**2*np.exp(-Delta*tau)
        scaled=slow*(np.pi*J*tau)**2/np.sin(theta)**2
        rows.append(dict(time=tau,slow=slow,fast=fast,total=slow+fast,
                         scaled_power_tail=scaled))
    assert abs(rows[-1]['scaled_power_tail']-1)<2e-5
    return dict(J=J,Delta=Delta,theta=theta,rows=rows)


def main():
    out=dict(status='exact soluble source test; not a no-go for local elimination',
             finite=[direct_check(4),direct_check(4,J=.02,theta=.02)],thermodynamic=thermodynamic())
    Path('massive_source_tails_checks.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':main()
