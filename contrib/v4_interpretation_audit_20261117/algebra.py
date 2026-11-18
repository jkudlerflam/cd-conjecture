#!/usr/bin/env python3
"""Small exact algebra checks, without eigenstate optimization."""
import time
import numpy as np
from threadpoolctl import threadpool_limits
from common import HERE,write


def main():
    start=time.monotonic();L=8;D=1<<L;states=np.arange(D)
    bits=(states[:,None]>>np.arange(L))&1;n=bits.sum(axis=1)
    z=1-2*bits;H=np.zeros((D,D),complex);rotated=H.copy();B=H.copy()
    for j in range(L):
        k=(j+1)%L;dest=states^(1<<j)^(1<<k)
        xy=1j*z[:,k];yx=1j*z[:,j];yy=-z[:,j]*z[:,k]
        H[dest,states]+=.5*(1-yy)+(xy-yx)/16
        rotated[dest,states]+=(9*xy+7*yx)/16
        B[dest,states]+=(bits[:,j]==bits[:,k])
    diagonal=(z*np.roll(z,-1,axis=1)).sum(axis=1)-1.75*z.sum(axis=1)
    H[states,states]+=diagonal;rotated[states,states]+=diagonal
    R=np.exp(-1j*(bits@(np.pi/2+np.pi*np.arange(L))))
    U=np.exp(-1j*np.pi/8*z.sum(axis=1))
    parity=(-1.)**n
    T=np.zeros((D,D));T[((states<<1)&(D-1))|(states>>(L-1)),states]=1
    errors=dict(antiunitary_H=float(np.max(abs(R[:,None]*H.conj()*R.conj()[None,:]-H))),
        antiunitary_T=float(np.max(abs(R[:,None]*T*R.conj()[None,:]-parity[:,None]*T))),
        uniform_rotation=float(np.max(abs(U[:,None]*H*U.conj()[None,:]-rotated))),
        parity_commutator=float(np.max(abs(H*(parity[:,None]-parity[None,:])))))
    number_commutator=float(np.max(abs(H*(n[:,None]-n[None,:]))))
    adjacent=(bits*np.roll(bits,-1,axis=1)).sum(axis=1)
    p=np.flatnonzero(adjacent==0);q=np.flatnonzero(adjacent>0)
    exact=-(B[np.ix_(p,q)]/(4*adjacent[q])[None,:])@B[np.ix_(q,p)]
    expected=np.zeros_like(exact);lookup={int(s):i for i,s in enumerate(p)}
    for col,s in enumerate(p):
        expected[col,col]=-L/4+3*n[s]/4-sum(bits[s]*np.roll(bits[s],-2))/4-sum(bits[s]*np.roll(bits[s],-3))/12
        for j in range(L):
            k=(j+2)%L
            if bits[s,j]!=bits[s,k]:
                dest=int(s^(1<<j)^(1<<k))
                if dest in lookup:expected[lookup[dest],col]-=1/8
    errors['leading_projected_second_order']=float(np.max(abs(exact-expected)))
    assert max(errors.values())<1e-12 and number_commutator>0
    write(HERE/'algebra.json',dict(L=L,errors=errors,number_commutator_max=number_commutator,
        projected_dimension=len(p),seconds=time.monotonic()-start,
        scope='Periodic spin algebra and leading second-order pair perturbation about H0=4 sum nn. Not a finite-V error bound.'))
    print(errors,'seconds',time.monotonic()-start)


if __name__=='__main__':
    with threadpool_limits(limits=1):main()
