"""A literal reversible work machine on the 15-letter sweeping clock.
Checks a permutation completion, local transport, and a boundary readout.
The readout is one global parity, not a high-channel Clifford family.
"""
from itertools import product
import json
import numpy as np
from sweeping_clock_checks import path


def complete_permutation(size,partial,grades=None):
    if grades is None: grades=np.zeros(size,dtype=int)
    assert len(set(partial.values()))==len(partial)
    assert all(grades[a]==grades[b] for a,b in partial.items())
    perm=np.full(size,-1,dtype=int)
    for a,b in partial.items():perm[a]=b
    for g in np.unique(grades):
        dom=[a for a in range(size) if grades[a]==g and perm[a]<0]
        used=set(partial.values())
        ran=[b for b in range(size) if grades[b]==g and b not in used]
        assert len(dom)==len(ran)
        for a,b in zip(dom,ran):perm[a]=b
    assert sorted(perm)==list(range(size))
    return perm


STATES=list(product((0,1),(-1,0,1),(-1,0,1)))
INDEX={s:i for i,s in enumerate(STATES)}
partial={}
for i,(x,l,r) in enumerate(STATES):
    if (l>=0)+(r>=0)==1:
        out=(x,l^x if l>=0 else -1,r^x if r>=0 else -1)
        partial[i]=INDEX[out]
COLLISION=complete_permutation(len(STATES),partial,np.array([s[0] for s in STATES]))


def F(state):
    x,l,r=state
    return STATES[COLLISION[INDEX[(x,r,l)]]]


def evolve(n,x):
    words,steps=path(n)
    work=[(bit,-1,-1) for bit in x]
    snapshots=[tuple(work)]
    for label,j in steps:
        if label=='activate':
            # The local unitary exchanges empty and R0, fixing other states.
            bit,l,r=work[j]
            if l==-1 and r==-1:work[j]=(bit,-1,0)
            elif l==-1 and r==0:work[j]=(bit,-1,-1)
        elif label=='out':
            bit,l,r=work[j];bit2,l2,r2=work[j+1]
            work[j]=(bit,l,l2);work[j+1]=(bit2,r,r2)
        elif label=='back':work[j+1]=F(work[j+1])
        elif label=='finish':work[0]=F(work[0])
        snapshots.append(tuple(work))
    return words,snapshots


def check(n):
    inputs=list(product((0,1),repeat=n));C=len(path(n)[0]);D=2**n
    compression=np.zeros((D,D))
    active_counts=[]
    for a,x in enumerate(inputs):
        words,snapshots=evolve(n,x)
        count=0;value=0
        for w,work in zip(words,snapshots):
            # On site 1, project onto an L-lane head and read its accumulator Z.
            bit,l,r=work[0]
            if l>=0 and r==-1:
                count+=1;value+=1-2*l
                assert l==sum(x)%2
            assert tuple(z[0] for z in work)==x
        compression[a,a]=value/C
        active_counts.append(count)
    assert len(set(active_counts))==1 and active_counts[0]>0
    alpha=active_counts[0]/C
    expected=np.diag([alpha*(-1)**sum(x) for x in inputs])
    residual=float(np.max(abs(compression-expected)))
    assert residual<1e-14
    # The off-diagonal accumulator operator cannot change any input data bit.
    # Since physical data strings are retained, all its compression entries vanish.
    xcomp=np.zeros((D,D))
    allsnaps=[evolve(n,x)[1] for x in inputs]
    for a in range(D):
        for b in range(D):
            for wa,wb in zip(allsnaps[a],allsnaps[b]):
                x,l,r=wb[0]
                if l<0 or r!=-1:continue
                changed=((x,1-l,r),)+wb[1:]
                xcomp[a,b]+=int(changed==wa)/C
    assert np.max(abs(xcomp))==0
    return dict(n=n,clock_vertices=C,logical_dimension=D,work_site_dimension=len(STATES),active_vertices=active_counts[0],parity_scale=alpha,parity_residual=residual,offdiagonal_readout_norm=float(np.linalg.norm(xcomp)))


if __name__=='__main__':
    out={'collision_partial_states':len(partial),'collision_total_states':len(STATES),'checks':[check(n) for n in range(2,6)]}
    with open('sweep_machine_checks.json','w') as f:json.dump(out,f,indent=2)
    print(json.dumps(out,indent=2))
