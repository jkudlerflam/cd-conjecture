"""Nested Schur checks of alpha jets and first residue-changing witnesses.

Uses the sibling's generic Schur recursion, not its displayed sixth formula.
The independent leading.py matrix is the comparison. Seventh coefficients
are valid here before the first O(epsilon^2) effective metric correction;
no eighth-order static coefficient is interpreted as a Hermitian Hamiltonian.
"""
import sys
from pathlib import Path
import numpy as np
from scipy import sparse

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'contrib/qubit_cancellation_classification_20261116'))
import algebra as A


def series(L,q,alpha=0j,beta=0j,order=7,tau=0j,mu6=0.):
    d=1<<L;I=sparse.eye(d,format='csr',dtype=complex)
    zero=lambda:sparse.csr_matrix((d,d),dtype=complex)
    u6=50 if q==1 else 332/3;t6=-24 if q==1 else -80;c6=2 if q==1 else 18
    H=[A.chain(L,A.bond(v=1)),A.chain(L,A.bond(a=1,b=-1)),
       A.chain(L,A.bond(u=2,t=-1,p=q)),A.chain(L,A.bond(p=alpha)),
       A.chain(L,A.bond(u=-(4+4*q+2*q*q),t=2+2*q,p=beta))+L*q*q*I,
       A.chain(L,A.bond(u=-(4+4*q)*np.real(alpha),t=2*np.real(alpha)))+2*L*q*np.real(alpha)*I,
       A.chain(L,A.bond(u=u6-mu6-2*abs(alpha)**2-(4+4*q)*np.real(beta),t=t6+2*np.real(beta)+tau))
            +L*(abs(alpha)**2+2*q*np.real(beta)-c6)*I]
    H += [zero() for _ in range(order-6)]
    P=A.projectors(L);S=A.schur_series(H,P)
    keep=np.flatnonzero(abs(S[2].diagonal())<1e-9);K=A.schur_series(S[2:],keep);P=P[keep]
    keep=np.flatnonzero(abs(K[2].diagonal())<1e-9)
    return P[keep],A.schur_series(K[2:],keep),H
