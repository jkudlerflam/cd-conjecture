"""An exactly solvable d=4 test of first-order criticality vs fixed coupling."""
import json
from pathlib import Path
import numpy as np
from scipy.linalg import eigvalsh


def bloch(k,eps,Delta=1.,g=1.):
    k=np.asarray(k);A=np.zeros(k.shape+(4,4),complex)
    t=eps*(1-np.exp(-1j*k));u=g*eps
    A[...,0,1]=t;A[...,1,0]=-t.conj()
    A[...,0,2]=u;A[...,2,0]=-u
    A[...,1,3]=-u;A[...,3,1]=u
    A[...,2,3]=Delta;A[...,3,2]=-Delta
    return 1j*A


def bands(k,eps,Delta=1.,g=1.):
    w=np.sin(np.asarray(k)/2)**2;u=g*eps
    tr=Delta**2+2*u*u+4*eps*eps*w
    det2=u**4+4*Delta*eps*(Delta*eps+u*u)*w
    high=np.sqrt((tr+np.sqrt(np.maximum(0,tr*tr-4*det2)))/2)
    low=np.sqrt(det2)/high
    return low,high


def energy(L,eps,g=1.):
    k=2*np.pi*(np.arange(L)+.5)/L
    lo,hi=bands(k,eps,g=g)
    return -.5*np.sum((lo+hi-1).astype(np.longdouble),dtype=np.longdouble)


def entropy_data(eps,g,lengths,grid=32768):
    k=2*np.pi*(np.arange(grid)+.5)/grid
    ev,U=np.linalg.eigh(bloch(k,eps,g=g))
    sign=(U*np.sign(ev)[:,None,:])@np.swapaxes(U.conj(),-1,-2)
    # C=i Gamma=-sign(i A). Half-shifted Fourier grid.
    coeff=np.fft.ifft(-sign,axis=0)
    rows=[]
    for ell in lengths:
        dist=np.arange(ell)[:,None]-np.arange(ell)[None,:]
        block=coeff[dist%grid]*np.exp(1j*np.pi*dist/grid)[...,None,None]
        C=block.transpose(0,2,1,3).reshape(4*ell,4*ell)
        assert np.max(np.abs(C-C.conj().T))<1e-12
        nu=eigvalsh(C,overwrite_a=True,check_finite=False)
        p=np.clip((1+nu)/2,0,1);mask=(p>0)&(p<1)
        h=-p[mask]*np.log(p[mask])-(1-p[mask])*np.log1p(-p[mask])
        rows.append(dict(length=ell,entropy=float(h.sum()/2)))
    for a,b in zip(rows,rows[1:]):
        a['c_entropy_apparent']=3*(b['entropy']-a['entropy'])/np.log(b['length']/a['length'])
    return rows


def spin_ground(L,eps):
    I=np.eye(2);X=np.array([[0,1],[1,0]]);Y=np.array([[0,-1j],[1j,0]]);Z=np.diag([1,-1])
    XI=np.kron(X,I);XZ=np.kron(X,Z)
    onsite=(np.eye(4)-np.kron(I,Z))/2-eps*np.kron(Z,I)/2+eps*(np.kron(Y,X)+np.kron(X,Y))/2
    H=np.zeros((4**L,4**L),complex)
    def term(ops):
        z=np.array([[1.]])
        for j in range(L):z=np.kron(z,ops.get(j,np.eye(4)))
        return z
    for j in range(L):
        H+=term({j:onsite})-eps/2*term({j:XZ,(j+1)%L:XI})
    return float(eigvalsh(H,subset_by_index=(0,0),check_finite=False)[0])


def main():
    tests=[]
    for eps in (.25,.1,.02,.005):
        k=np.linspace(0,2*np.pi,4097)
        ev=np.linalg.eigvalsh(bloch(k,eps))
        low,high=bands(k,eps)
        err=float(np.max(np.abs(ev[:,2:]-np.array([low,high]).T)))
        gap=2*eps*eps/(np.sqrt(1+4*eps*eps)+1)
        assert err<1e-12 and abs(low.min()-gap)<1e-13
        tests.append(dict(epsilon=eps,bulk_gap=float(gap),gap_over_epsilon2=float(gap/eps**2),
                          band_error=err,gap_grid_points=len(k)))
    dense=[]
    for L in (3,4,5):
        exact=spin_ground(L,.1);free=float(energy(L,.1))
        err=abs(exact-free);assert err<1e-11
        dense.append(dict(length=L,dimension=4**L,spin_ground=exact,free_ground=free,error=err))
    estimators=[]
    for eps,g in ((.02,0),(.02,1),(.005,1)):
        lengths=[8,16,32,64,128,256]
        ent=entropy_data(eps,g,lengths)
        cas=[]
        for L in (8,16,32,64,128,256,512):
            c=float(4*np.longdouble(L)*(energy(2*L,eps,g)-2*energy(L,eps,g))/(np.pi*eps))
            cas.append(dict(length=L,c_casimir_apparent=c))
        if g==0:
            assert abs(cas[-1]['c_casimir_apparent']-.5)<1e-5
            assert abs(ent[-2]['c_entropy_apparent']-.5)<3e-4
        estimators.append(dict(epsilon=eps,mixing=g,entropy_grid=32768,entropy=ent,casimir=cas))
    out=dict(status='soluble crossover control, not the C8 microscopic candidate',
             band_checks=tests,dense_spin_checks=dense,estimators=estimators)
    Path('critical_crossover_checks.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':main()
