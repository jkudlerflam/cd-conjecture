"""Independent occupation-space implementation of the two leading models.

No imports from the sibling perturbation calculation. Bit j is occupation j.
Use compatible rings L divisible by m when interpreting residue sequences.
"""
import numpy as np
from scipy import sparse


def parameters(q):
    if q==1: return dict(m=4,t=1.,forbidden=(1,2,3),g=1.,delta=1.)
    if q==3: return dict(m=3,t=128/3,forbidden=(1,2,4),g=16/3,delta=2/3)
    raise ValueError(q)


def bits(s,L): return tuple(j for j in range(L) if s>>j&1)


def valid(s,L,q):
    occ=bits(s,L)
    return not any((s>>((j+r)%L))&1 for j in occ for r in parameters(q)['forbidden'])


def basis(L,q,N=None):
    m=parameters(q)['m']; out=[]
    def grow(xs,left):
        if left==0:
            s=sum(1<<j for j in xs)
            if valid(s,L,q):out.append(s)
            return
        start=xs[-1]+m if xs else 0
        for j in range(start,L-m*(left-1)):
            if q==3 and xs and j-xs[-1]==4:continue
            grow(xs+[j],left-1)
    for n in (range(L//m+1) if N is None else [N]):grow([],n)
    return np.array(sorted(out),dtype=np.int64 if L<63 else object)


def dr(s,L,r):return sum((s>>((j+r)%L))&1 for j in bits(s,L))


def g8(s,L):
    return sum(((s>>((j+8)%L))&1)*(1-((s>>((j+3)%L))&1))*(1-((s>>((j+5)%L))&1)) for j in bits(s,L))


def residue_word(s,L,m):
    w=tuple(j%m for j in bits(s,L))
    return min(w[i:]+w[:i] for i in range(len(w))) if w else ()


def matrix(L,q,states=None,alpha=0j,shift=True,tau=0j,source=0j):
    """H6 plus the checked alpha jet; shift adds 2t N.

    alpha's facilitated creation coefficient is g*conj(alpha). tau is the
    freely adjustable nearest-neighbour hop, source the uniform creation.
    Do not restrict to fixed N when alpha or source is nonzero.
    """
    p=parameters(q);m,t=p['m'],p['t']
    states=basis(L,q) if states is None else states
    pos={int(s):i for i,s in enumerate(states)};rr=[];cc=[];vv=[]
    for col,ss in enumerate(states):
        s=int(ss);N=s.bit_count()
        diag=-2*t*dr(s,L,m)+p['delta']*abs(alpha)**2*dr(s,L,m)
        diag+=(10/3*dr(s,L,5)+14/3*dr(s,L,6)) if q==1 else 8*g8(s,L)
        if shift:diag+=2*t*N
        rr.append(col);cc.append(col);vv.append(diag)
        for j in range(L):
            for r,amp in ((m,-t),(1,tau)):
                k=(j+r)%L
                if (s>>j)&1 and not (s>>k)&1:
                    dest=s^(1<<j)^(1<<k)
                    if dest in pos:
                        row=pos[dest];rr.extend((row,col));cc.extend((col,row));vv.extend((amp,np.conj(amp)))
            if not (s>>j)&1:
                dest=s|(1<<j)
                if dest in pos:
                    amp=p['g']*np.conj(alpha)*(((s>>((j-m)%L))&1)-((s>>((j+m)%L))&1))+source
                    if amp:
                        row=pos[dest];rr.extend((row,col));cc.extend((col,row));vv.extend((amp,np.conj(amp)))
    H=sparse.coo_matrix((vv,(rr,cc)),shape=(len(states),len(states)),dtype=complex).tocsr()
    H.eliminate_zeros()
    return states,H


def laplacian_and_penalty(L,q,states,H):
    """Graph degree from actual transitions, penalty independently from words."""
    diag=H.diagonal().real
    off=H-sparse.diags(diag);degree=-np.asarray(off.sum(axis=1)).ravel().real
    lap=off+sparse.diags(degree)
    if q==1:potential=np.array([16/3*dr(int(s),L,5)+20/3*dr(int(s),L,6)+2*dr(int(s),L,7) for s in states])
    else:
        t=parameters(q)['t']
        potential=np.array([2*t*(dr(int(s),L,5)+dr(int(s),L,7))+8*g8(int(s),L) for s in states])
    return lap,potential,diag-degree-potential
