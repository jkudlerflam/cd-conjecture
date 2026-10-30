"""C123 full-spectrum, symmetry-block ED with no spectral truncation."""
from pathlib import Path
from itertools import product
from collections import Counter
import json
import time
import numpy as np
from scipy.linalg import eigh
from scipy.special import logsumexp

ROOT=Path(__file__).resolve().parent

def sectors(n,d,delta=1.):
    powers=d**np.arange(n,dtype=np.int64)
    basis=np.arange(d**n,dtype=np.int64)
    digits=(basis[:,None]//powers[None,:])%d
    counts=np.stack([(digits==c).sum(axis=1) for c in range(d)],axis=1)
    groups={}
    for s,c in enumerate(counts): groups.setdefault(tuple(c),[]).append(s)
    ans=[]
    for key,states in groups.items():
        if d==2:
            if key[0]<key[1]: continue
            mult=1 if key[0]==key[1] else 2
        else:
            if tuple(sorted(key,reverse=True))!=key: continue
            mult=6
            for value in Counter(key).values():
                mult//= [1,1,2,6][value]
        states=np.array(states)
        lookup=np.full(d**n,-1,dtype=int); lookup[states]=np.arange(len(states))
        h=np.zeros((len(states),len(states)))
        dd=digits[states]
        for i in range(n):
            j=(i+1)%n
            x,y=dd[:,i],dd[:,j]
            changed=states+(y-x)*powers[i]+(x-y)*powers[j]
            col=np.arange(len(states))
            if d==3: h[lookup[changed],col]+=1.
            else:
                h[col,col]+=delta*(.5-x)*(.5-y)
                keep=x!=y
                h[lookup[changed[keep]],col[keep]]+=.5
        ev,vec=eigh(h,check_finite=False,driver='evd')
        ans.append(dict(ev=ev,vec=vec,mult=mult,digits=dd,counts=key))
    assert sum(len(s['ev'])*s['mult'] for s in ans)==d**n
    ground=min(float(s['ev'][0]) for s in ans)
    return dict(sectors=ans,E0=ground,n=n,d=d,
                maxblock=max(len(s['ev']) for s in ans))

def channel(spec,beta,N,sites):
    total=0.
    for s in spec['sectors']:
        ev,vec=s['ev'],s['vec']
        e=(vec*np.exp(-beta*(ev-spec['E0'])/(2*N)))@vec.T
        labels=sum(spec['d']**j*s['digits'][:,site] for j,site in enumerate(sites))
        mask=labels[:,None]==labels[None,:]
        rho=e.copy()
        for j in range(N):
            rho*=mask
            if j<N-1: rho=e@rho@e
        total+=s['mult']*np.sum(rho*e.T)
    return float(total)

def histories(spec,beta,N):
    # Reconstruct omitted charge/color permutations for the small numerator ring.
    # Fixed-record weights require the actual color permutation in each sector,
    # rather than multiplication by the total symmetry multiplicity.
    n,d=spec['n'],spec['d']
    allvals=np.zeros(d**N)
    labels=list(product(range(d),repeat=N))
    for s in spec['sectors']:
        from itertools import permutations
        color_perms=list(set(permutations(range(d))))
        transformed={}
        for perm in color_perms:
            key=tuple(s['counts'][perm.index(c)] for c in range(d))
            transformed.setdefault(key,perm)
        for perm in transformed.values():
            ev,vec=s['ev'],s['vec']
            f=(vec*np.exp(-beta*(ev-spec['E0'])/(4*N)))@vec.T
            e=f@f
            colors=np.array(perm)[s['digits'][:,0]]
            for k,history in enumerate(labels):
                w=f.copy()
                for j,label in enumerate(history):
                    w=(colors==label)[:,None]*w
                    w=(f if j==N-1 else e)@w
                allvals[k]+=np.sum(w*w)
    return allvals

def logz(spec,beta):
    return float(logsumexp(np.concatenate([
        -beta*(s['ev']-spec['E0'])+np.log(s['mult']) for s in spec['sectors']])))

def one_case(small,big,a,N,M,model):
    L=small['n']; d=small['d']; beta=N*a*L/M
    q=histories(small,beta,N)
    ql=channel(small,beta,N,[0])
    qd=channel(big,beta,N,[0,L])
    err=abs(q.sum()-ql)/ql
    assert err<2e-10,err
    defect=big['E0']-2*small['E0']
    G=beta*defect+np.log(q@q)-np.log(qd)
    R=beta*defect+2*logz(small,beta)-logz(big,beta)
    loss=logz(small,beta)-np.log(ql)
    lossd=logz(big,beta)-np.log(qd)
    assert -2e-9<=loss<=N*np.log(d)+2e-9
    assert -2e-9<=lossd<=2*N*np.log(d)+2e-9
    assert R-3*N*np.log(d)-2e-9<=G<=R+2*N*np.log(d)+2e-9
    return dict(model=model,L=L,N=N,a=a,G=float(G),G_per_N=float(G/N),
                R=float(R),loss=float(loss),loss_D=float(lossd),
                history_error=float(err),delta=float(defect),maxblock=big['maxblock'])

def main():
    start=time.time(); rows=[]
    for model,d,delta,M,lengths in [
        ('XX',2,0.,.5,[4,6]),('XXZ_m09',2,-.9,.5,[4,6]),
        ('XXZ_m05',2,-.5,.5,[4,6]),('XXZ_p05',2,.5,.5,[4,6]),
        ('Heisenberg',2,1.,.5,[4,6]),('SU3',3,1.,1.,[3,4])]:
        for L in lengths:
            small,big=sectors(L,d,delta),sectors(2*L,d,delta)
            for a in [.25,1.,2.,4.]:
                for N in [1,2,4,8]:
                    row=one_case(small,big,a,N,M,model)
                    if model=='XX':
                        from checks import xx_case
                        gg=xx_case(L,a,N)
                        row['gaussian_error']=abs(row['G']-gg['G'])
                        assert row['gaussian_error']<2e-9
                    rows.append(row)
            (ROOT/'ed_results.json').write_text(json.dumps(dict(ed=rows,seconds=time.time()-start),indent=2)+'\n')
            print(model,L,'block',big['maxblock'],'seconds',time.time()-start,flush=True)
    print('completed',len(rows),'cases',time.time()-start,flush=True)

if __name__=='__main__':main()
