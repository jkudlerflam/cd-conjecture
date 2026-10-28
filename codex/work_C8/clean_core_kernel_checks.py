"""Apply the actual full-space local Hamiltonian to sparse initialized histories."""
from collections import defaultdict
from itertools import product
import json
from clean_core_bond import CoreBond,Site
from clean_readout_checks import clock_stream


def check(bits):
    q=4;H=CoreBond(q);M=H.M;c=M.initialize(q,bits);n=len(c.roles);support=[]
    for clocks,lab,j in clock_stream(n):
        if lab=='out':M.bond(c,j,q)
        elif lab=='back':M.onsite(c,j+1,q)
        elif lab=='finish':M.onsite(c,0,q)
        word=[]
        for i in range(n):
            head=H.head_encode(c.control,c.p,c.x,c.phase) if i==c.pos else 0
            word.append(H.encode(Site(clocks[i],c.roles[i],*c.tape[i],head)))
        support.append(tuple(word))
    assert len(set(support))==len(support)
    residual=defaultdict(float);entries=0
    # Cache columns only within this small code support; no full Hilbert matrix.
    bulk={};left={};right={}
    for word in support:
        for j in range(n-1):
            pair=word[j:j+2]
            if pair not in bulk:bulk[pair]=H.bulk(*pair)
            for dest,v in bulk[pair].items():
                residual[word[:j]+dest+word[j+2:]]+=v;entries+=1
        if word[0] not in left:left[word[0]]=H.boundary(word[0],'L')
        for dest,v in left[word[0]].items():residual[(dest,)+word[1:]]+=v;entries+=1
        if word[-1] not in right:right[word[-1]]=H.boundary(word[-1],'R')
        for dest,v in right[word[-1]].items():residual[word[:-1]+(dest,)]+=v;entries+=1
    err=max(map(abs,residual.values()))
    assert err==0
    return dict(input=list(bits),history_support=len(support),applied_matrix_entries=entries,max_unnormalized_residual=err,
                distinct_bulk_columns=len(bulk))

if __name__=='__main__':
    rows=[check(bits) for bits in product(range(2),repeat=2)]
    open('clean_core_kernel_checks.json','w').write(json.dumps(rows,indent=2)+'\n')
    print(json.dumps(rows,indent=2))
