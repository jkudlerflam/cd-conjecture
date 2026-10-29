"""Bond audits and modest Hubbard ED; numpy/scipy only."""
import itertools,json
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh
from ed_checks import entropies

def hubbard(L,Nu,Nd,U):
    masks=lambda N:[sum(1<<i for i in c) for c in itertools.combinations(range(L),N)]
    states=[tuple(((u>>i)&1)+2*((v>>i)&1) for i in range(L)) for u in masks(Nu) for v in masks(Nd)]
    # Interleaved Fock ordering, a_i then b_i. The same tensor bond on the seam
    # gives the parity-twisted fermionic ring, as in the other tensor-chain tests.
    index={s:k for k,s in enumerate(states)}; rows=[];cols=[];vals=[]
    for k,s in enumerate(states):
        rows.append(k);cols.append(k);vals.append(U*sum(a==3 for a in s))
        for i in range(L):
            j=(i+1)%L
            for flavor in [0,1]:
                bit=1<<flavor
                if bool(s[i]&bit)==bool(s[j]&bit):continue
                ss=list(s);ss[i]^=bit;ss[j]^=bit
                # Local two-site matrix, including i=L-1,j=0 in oriented order.
                sign=(-1)**(((s[i]>>1)&1) if flavor==0 else (s[j]&1))
                rows.append(index[tuple(ss)]);cols.append(k);vals.append(-sign)
    H=coo_matrix((vals,(rows,cols)),shape=(len(states),len(states))).tocsr()
    assert (H-H.T).nnz==0
    w,v=eigsh(H,k=3,which='SA',tol=1e-11,v0=np.random.default_rng(76).normal(size=len(states)))
    ix=np.argsort(w);w=w[ix];v=v[:,ix]
    return dict(L=L,Nu=Nu,Nd=Nd,U=U,dimension=len(states),energies=w.tolist(),entropy=entropies(v[:,0],states,4,[(0,0),(1,0),(0,1),(1,1)]),residual=float(np.linalg.norm(H@v[:,0]-w[0]*v[:,0])))

def osp(n,m):
    d=n+2*m; A=np.zeros((d,d));A[:n,:n]=np.eye(n)
    for j in range(m):A[n+2*j,n+2*j+1]=1;A[n+2*j+1,n+2*j]=-1
    p=[0]*n+[1]*(2*m);P=np.zeros((d*d,d*d))
    for a in range(d):
        for b in range(d):P[b*d+a,a*d+b]=(-1)**(p[a]*p[b])
    E=np.outer(A.ravel(),np.linalg.inv(A).ravel());z=n-2*m
    h=P+2/(2-z)*E
    return {'n':n,'m':m,'d':d,'E_squared_minus_zE':float(np.linalg.norm(E@E-z*E)),'h_nonhermiticity_Frobenius':float(np.linalg.norm(h-h.T)),'expected_nonhermiticity':float(8*np.sqrt(n*m)/abs(2-z))}

out={'osp':[osp(n,m) for n,m in [(1,1),(2,1),(3,1),(1,2),(4,1)] if n-2*m!=2], 'hubbard':[]}
for U in [0.,4.]:
    for L in [4,8,12]:
        row=hubbard(L,L//4,L//4,U);out['hubbard'].append(row);print(json.dumps(row),flush=True)
with open('extra_results.json','w') as f:json.dump(out,f,indent=2)
print(json.dumps(out['osp']),flush=True)
