"""Sp(4) invariant-sign audit and finite-size Casimir check."""
import itertools,json
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh
D=4;P=np.zeros((16,16));v=np.zeros(16)
for a in range(D):
    for b in range(D):P[b*D+a,a*D+b]=1
for a in range(D):v[a*D+3-a]=1 if a<2 else -1
K=np.outer(v,v)

def bonds(h,L,sector=False):
    states=[s for s in itertools.product(range(4),repeat=L) if not sector or (s.count(0)==s.count(3) and s.count(1)==s.count(2))]
    ix={s:k for k,s in enumerate(states)};out=[]
    for j in range(L):
        rows=[];cols=[];vals=[];k=(j+1)%L
        for src,s in enumerate(states):
            ab=s[j]*4+s[k]
            for cd in np.flatnonzero(h[:,ab]):
                ss=list(s);ss[j]=cd//4;ss[k]=cd%4
                rows.append(ix[tuple(ss)]);cols.append(src);vals.append(h[cd,ab])
        out.append(coo_matrix((vals,(rows,cols)),shape=(len(states),len(states))).tocsr())
    return out

out={'sign_audit':[],'finite_size':[]}
for sign in [-1,1]:
    hs=[x.toarray() for x in bonds(P+sign*K/3,4)]
    H=sum(hs);Q=sum(hs[j]@hs[(j+1)%4]-hs[(j+1)%4]@hs[j] for j in range(4))
    out['sign_audit'].append({'K_sign':sign,'L':4,'commutator_norm':float(np.linalg.norm(H@Q-Q@H))})
ebulk=-2*(2*np.log(2)/3+np.pi/(9*np.sqrt(3))-1/3)
vslow=np.pi/3
for L in [4,6,8]:
    H=sum(bonds(P+K/3,L,True))
    w,V=eigsh(H,k=2,which='SA',tol=1e-12,v0=np.random.default_rng(76).normal(size=H.shape[0]))
    k=np.argmin(w);E=float(w[k]);res=np.linalg.norm(H@V[:,k]-E*V[:,k])
    out['finite_size'].append({'L':L,'dimension':H.shape[0],'E0':E,'e_bulk_reference':ebulk,'v_slow_reference':vslow,'velocity_weighted_effective_c':-6*L*(E-L*ebulk)/(np.pi*vslow),'residual':float(res)})
assert out['sign_audit'][1]['commutator_norm']<1e-11
with open('symplectic_results.json','w') as f:json.dump(out,f,indent=2)
print(json.dumps(out,indent=2))
