"""C88 pin conservation and single-letter deletion witnesses.
All enumeration is literal, and no mismatch cutoff is used in the
finite-ring conservation check. The separate witness graph is legal.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import itertools,json
import numpy as np
from scipy.sparse.linalg import eigsh
import bond13 as b
from check_frozen import matrix

def pin(w):
    m=[a for a in w if a[0]!='V']
    return sum(a==('L',1) and b.marker(m[(i+1)%len(m)]) for i,a in enumerate(m))

def transitions(h,d):
    return {i:[(int(j),h[j,i]) for j in np.flatnonzero(h[:,i]) if j!=i] for i in range(d*d)}

def conserve():
    L=4;d=len(b.LET);h=b.bond(b.params(K=37.,field=.021,musel=.93),lam=.17)
    tr=transitions(h,d);count=0
    for w in itertools.product(range(d),repeat=L):
        npin=pin([b.LET[a] for a in w])
        for i,a in enumerate(w):
            j=(i+1)%L
            for q,v in tr[a*d+w[j]]:
                z=list(w);z[i],z[j]=divmod(q,d)
                assert pin([b.LET[a] for a in z])==npin
                count+=1
    return dict(L=L,basis_words=d**L,directed_bond_transitions=count)

def witnesses():
    seed=[('L',1),('R',0),('e',0),('L',0),('R',0),('e',0),('V',0)]
    d=len(b.LET);L=len(seed)
    tr=transitions(b.bond(b.params(n=2,K=0.),lam=.17),d)
    seed=tuple(b.IX[a] for a in seed);states=[seed];index={seed:0};found={};edges=0
    for w in states:
        assert pin([b.LET[a] for a in w])==0
        for i,a in enumerate(w):
            j=(i+1)%L
            for q,v in tr[a*d+w[j]]:
                z=list(w);z[i],z[j]=divmod(q,d);z=tuple(z)
                if any(not b.allowed(b.LET[a],b.LET[z[(k+1)%L]]) for k,a in enumerate(z)):continue
                edges+=1
                if z not in index:index[z]=len(states);states.append(z)
                for x in set(z)-set(w):
                    if x not in found:found[x]=dict(input=[b.LET[a] for a in w],output=[b.LET[a] for a in z],bond=i,amplitude=[v.real,v.imag])
    assert set(found)==set(range(d)),sorted(set(range(d))-set(found))
    return dict(legal_component_dimension=len(states),directed_bond_transitions=edges,
                witnesses={str(b.LET[x]):found[x] for x in sorted(found)})

def frozen_pins():
    out=[];markers=[1,3]
    choices=[[i for i,a in enumerate(b.MAT) if b.marker(a)==(j in markers)] for j in range(4)]
    states=list(itertools.product(*choices))
    pins=np.array([pin([b.MAT[a] for a in w]) for w in states])
    for lam in (0.,.003):
        H=matrix(markers,b.params(n=1,K=100.,musel=1.1,g=.03,td=.001),lam)
        rr,cc=H.nonzero();assert np.all(pins[rr]==pins[cc])
        energies=[];sizes=[];res=[]
        for q in (0,1,2):
            ids=np.flatnonzero(pins==q);A=H[ids][:,ids]
            e,v=eigsh(A,k=1,which='SA',tol=1e-12,v0=np.ones(len(ids)))
            energies.append(float(e[0]));sizes.append(len(ids));res.append(float(np.linalg.norm(A@v[:,0]-e[0]*v[:,0])))
        assert energies[0]<min(energies[1:])
        out.append(dict(lambda_value=lam,dimensions=sizes,lowest_by_pin_count=energies,max_residual=max(res)))
    return out

if __name__=='__main__':
    out=dict(tag='exact_enumeration_and_numerical_frozen_ED',conservation=conserve(),deletion=witnesses(),frozen_pin_selection=frozen_pins())
    json.dump(out,open('sector_checks.json','w'),indent=2)
    print(json.dumps({**out,'deletion':{k:v for k,v in out['deletion'].items() if k!='witnesses'}}),flush=True)
