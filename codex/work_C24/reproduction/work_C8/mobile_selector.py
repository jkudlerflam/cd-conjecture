"""First-order mobile-marker selector checks, not a full finite-u CFT study."""
import itertools
import json
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import eigsh


def overlap(b):
    i=np.arange(1,b)
    return float(np.sum(np.sqrt(2/b)*np.sin(np.pi*i/b)
                        *np.sqrt(2/(b+1))*np.sin(np.pi*i/(b+1))))


def transition_checks():
    # Alphabet #,L,X,R as arrows between 0=L and 1=R.
    arrows=[(1,0),(0,0),(0,1),(1,1)]
    legal=lambda a,b: arrows[a][1]==arrows[b][0]
    # Both moves preserve all three affected legality indicators.
    moves=[((1,2),(2,3)),((0,1),(3,0))]
    failures=[]
    for old,new in moves:
        for a,d in itertools.product(range(4),repeat=2):
            x=[a,*old,d]; y=[a,*new,d]
            if [legal(*x[j:j+2]) for j in range(3)]!=[legal(*y[j:j+2]) for j in range(3)]:
                failures.append((x,y))
    assert not failures
    return {'local_legality_failures':len(failures)}


def rod_graph(b,M,r):
    L=b*M+r
    states=set()
    for slots in itertools.combinations(range(M),r):
        gaps=np.full(M,b,dtype=int); gaps[list(slots)]+=1
        positions=np.r_[0,np.cumsum(gaps)[:-1]]
        for shift in range(L):
            states.add(tuple(sorted((positions+shift)%L)))
    states=sorted(states); index={s:i for i,s in enumerate(states)}
    rows=[]; cols=[]
    for n,state in enumerate(states):
        for a,x in enumerate(state):
            left=(x-state[a-1])%L
            right=(state[(a+1)%M]-x)%L
            if left==b and right==b+1:
                target=list(state); target[a]=(x+1)%L; target=tuple(sorted(target))
                j=index[target]
                rows.extend([n,j]); cols.extend([j,n])
    H=coo_matrix((-np.ones(len(rows)),(rows,cols)),shape=(len(states),len(states))).tocsr()
    assert (H-H.T).nnz==0
    components=connected_components(H,directed=False,return_labels=False)
    vals=np.sort(eigsh(H,k=min(8,len(states)-2),which='SA',tol=1e-10,
                       v0=np.random.default_rng(713).normal(size=len(states)),return_eigenvectors=False))
    prediction=-2*np.sin(np.pi*r/M)/np.sin(np.pi/M)
    assert abs(vals[0]-prediction)<1e-8
    gap=float(vals[1]-vals[0])
    return dict(b=b,M=M,r=r,L=L,dimension=len(states),components=int(components),
                energy=float(vals[0]),free_gap_chain_energy=float(prediction),
                gap=gap,L_gap=L*gap,L2_gap=L*L*gap)


def main():
    out={'local_checks':transition_checks(),'overlaps':{str(b):overlap(b) for b in [5,8,16]},
         'dense_fluid':[rod_graph(5,M,M//2) for M in [6,8,10,12]],
         'one_defect':[rod_graph(5,M,1) for M in [6,8,10,12]]}
    with open('mobile_selector_checks.json','w') as f: json.dump(out,f,indent=2)
    print(json.dumps(out,indent=2))


if __name__=='__main__': main()
