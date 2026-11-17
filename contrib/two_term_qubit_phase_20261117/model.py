"""Literal two-term chain in an exactly equivalent real compass basis.

The physical H is sum(2 X_j Y_(j+1)+g Z_j Z_(j+1)). Rotating odd
sites by pi/2 around Z gives -2 XX on even bonds, +2 YY on odd bonds,
and g ZZ everywhere. Periodic length must be even. Two exact parities
permit real blocks of dimension 2**(L-2); no excitation sector is discarded.
"""
import numpy as np
from scipy import sparse


class Block:
    def __init__(self,L,parity,flip):
        assert L%2==0 and parity in (-1,1) and flip in (-1,1)
        self.L=L;self.parity=parity;self.flip=flip;self.mask=(1<<L)-1
        self.states=np.array([s for s in range(1<<(L-1)) if (-1)**s.bit_count()==parity],dtype=np.int64)
        self.index=np.full(1<<L,-1,dtype=np.int64)
        self.index[self.states]=np.arange(len(self.states))
        self.index[self.mask^self.states]=np.arange(len(self.states))
        self.z=1-2*((self.states[:,None]>>np.arange(L))&1)

    def pair(self,j,kind):
        """Projected Pauli XX, YY or ZZ on oriented adjacent sites."""
        k=(j+1)%self.L;s=self.states
        if kind=='zz':return sparse.diags(self.z[:,j]*self.z[:,k],format='csr',dtype=float)
        dest=s^(1<<j)^(1<<k)
        phase=np.where(dest<(1<<(self.L-1)),1.,self.flip)
        if kind=='yy':phase=phase*(-self.z[:,j]*self.z[:,k])
        elif kind!='xx':raise ValueError(kind)
        return sparse.coo_matrix((phase,(self.index[dest],np.arange(len(s)))),shape=(len(s),)*2).tocsr()

    def operators(self,periodic=True):
        bonds=range(self.L if periodic else self.L-1)
        self.A=[self.pair(j,'xx')*(-1) if j%2==0 else self.pair(j,'yy') for j in bonds]
        self.B=[self.pair(j,'yy') if j%2==0 else -self.pair(j,'xx') for j in bonds]
        self.D=[self.pair(j,'zz') for j in bonds]
        return sum(self.A)*2,sum(self.D)

    def translation(self,steps=1):
        """Physical one-site translation, or ordinary two-site shift.

        One-site translation stays in a fixed flip block only at parity +1.
        At parity -1 it maps flip to -flip, giving exact paired spectra.
        """
        s=self.states;L=self.L
        dest=((s<<steps)&self.mask)|(s>>(L-steps))
        target_flip=self.flip if steps%2==0 else self.parity*self.flip
        phase=np.ones(len(s),complex)
        if steps%2:
            odd=sum((s>>j)&1 for j in range(1,L,2))
            even=sum((s>>j)&1 for j in range(0,L,2))
            phase=(1j)**(even-odd)
        phase*=np.where(dest<(1<<(L-1)),1,target_flip)
        return sparse.coo_matrix((phase,(self.index[dest],np.arange(len(s)))),shape=(len(s),)*2).tocsr()

    def expand(self,v):
        out=np.zeros(1<<self.L,dtype=np.result_type(v,complex))
        out[self.states]=v/np.sqrt(2);out[self.mask^self.states]=self.flip*v/np.sqrt(2)
        return out


def pauli_full(L,sites):
    """Independent bit assembly of a Pauli word in the physical basis."""
    s=np.arange(1<<L,dtype=np.int64);dest=s.copy();phase=np.ones(len(s),complex)
    for j,p in sites.items():
        z=1-2*((s>>j)&1)
        if p=='z':phase*=z
        elif p=='x':dest^=1<<j
        elif p=='y':dest^=1<<j;phase*=1j*z
        else:raise ValueError(p)
    return sparse.coo_matrix((phase,(dest,s)),shape=(len(s),)*2).tocsr()


def physical_full(L,g,periodic=True):
    return sum(2*pauli_full(L,{j:'x',(j+1)%L:'y'})+g*pauli_full(L,{j:'z',(j+1)%L:'z'})
               for j in range(L if periodic else L-1))


def compass_phases(L):
    s=np.arange(1<<L,dtype=np.int64)
    return (1j)**sum((s>>j)&1 for j in range(1,L,2))
