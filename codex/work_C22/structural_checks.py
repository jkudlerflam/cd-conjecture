import itertools
import json
import numpy as np
from scipy.linalg import eigh
from spatial_bus import walker_cell,h2,kron_all,one_bus,data_ops


def grammar():
    bulk=['L','R']+['X'+a for a in 'v0+-']
    marks=['Mempty']+['M'+a for a in 'v0+-']
    arrows={'L':('L','L'),'R':('R','R'),'Mempty':('R','L')}
    arrows.update({'X'+a:('L','R') for a in 'v0+-'})
    arrows.update({'M'+a:('L','L') for a in 'v0+-'})
    ok=lambda a,b:arrows[a][1]==arrows[b][0]
    counts=[]
    for n in range(1,7):
        count=0
        for word in itertools.product(bulk,repeat=n):
            for mark in marks:
                w=('Mempty',)+word+(mark,)
                if all(ok(a,b) for a,b in zip(w,w[1:])):
                    count+=1
                    assert sum(x.startswith('X') or (x.startswith('M') and x!='Mempty') for x in w[1:])==1
        assert count==4*(n+1)
        counts.append([n,count])
    letters=bulk+marks
    rules=[]
    for a in 'v0+-':
        rules.extend([(('L','X'+a),('X'+a,'R')),
                      (('X'+a,'Mempty'),('L','M'+a))])
    checked=0
    for before,after in rules:
        for left,right in itertools.product(letters,repeat=2):
            w=(left,*before,right);v=(left,*after,right)
            assert [ok(a,b) for a,b in zip(w,w[1:])]==[ok(a,b) for a,b in zip(v,v[1:])]
            checked+=1
    return dict(legal_words=counts,context_rule_checks=checked,local_dimension=17)


def direct_second_order(n,M):
    H,G,ports,a=walker_cell(n,g=.2,phases=(.41,.31,-.07))
    e,U=eigh(H)
    e[abs(e)<1e-10]=0
    Ge=U.conj().T@G
    L,R=[[U.conj().T@F@G for F in p] for p in ports]
    W=[]
    for j in range(M-1):
        bond=np.kron(R[0],L[0])+1j*np.kron(R[1],L[2])-1j*np.kron(R[2],L[1])
        W.append(np.kron(np.kron(kron_all([Ge]*j),bond),kron_all([Ge]*(M-j-2))))
    energies=np.array([0.])
    for _ in range(M):energies=(energies[:,None]+e[None,:]).reshape(-1)
    inv=np.zeros_like(energies);inv[energies>1e-10]=1/energies[energies>1e-10]
    total=sum(W)
    K=-total.conj().T@(inv[:,None]*total)
    cross=max([float(np.linalg.norm(W[i].conj().T@(inv[:,None]*W[j])))
               for i in range(len(W)) for j in range(i)] or [0.])
    Kb,_=h2(n,g=.2,phases=(.41,.31,-.07))
    pred=sum(np.kron(np.kron(np.eye(2**(n*j)),Kb),np.eye(2**(n*(M-j-2)))) for j in range(M-1))*a**4
    assert np.linalg.norm(K-pred,2)<1e-12 and cross<1e-12
    return dict(n=n,cells=M,dimension=len(H)**M,error=float(np.linalg.norm(K-pred,2)),cross=cross)


def gauge_check(n):
    A,_=one_bus(n,g=.2,phases=(.41,.31,-.07),gauged=False)
    B,_=one_bus(n,g=.2,phases=(.41,.31,-.07),gauged=True)
    _,zs,_,_=data_ops(n);N=2**n
    from scipy.linalg import block_diag
    Us=[];u=np.eye(N)
    for j in range(n+1):
        Us.extend([u.copy(),np.eye(N),np.eye(N)])
        if j<n:u=u@zs[j]
    U=block_diag(*Us)
    return float(np.linalg.norm(U.conj().T@A@U-B))


if __name__=='__main__':
    result={'grammar':grammar(),'gauge_errors':[gauge_check(n) for n in (1,2,3)],'H2':[]}
    print(json.dumps(result),flush=True)
    for n,M in ((1,2),(1,3),(2,2),(2,3)):
        row=direct_second_order(n,M);result['H2'].append(row);print(json.dumps(row),flush=True)
    with open('structural_results.json','w') as f:json.dump(result,f,indent=2)
