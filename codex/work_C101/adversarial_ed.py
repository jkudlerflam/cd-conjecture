"""Spin-orbital ring screen. Small mode L=4,6,8. Large mode only inside SLURM.
A two-size Casimir fit is a diagnostic, not evidence of conformality.
"""
import argparse,itertools,json,os
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh

def masks(L):return np.array([sum(1<<i for i in c) for c in itertools.combinations(range(L),L//2)],dtype=np.int64)
def operators(L):
    ms=masks(L);m=len(ms);ix={int(s):i for i,s in enumerate(ms)}
    dim=m*m;cols=np.arange(dim);out=[]
    for typ in range(3):
        rows=[]
        for j in range(L):
            k=(j+1)%L
            flip=(((ms>>j)^(ms>>k))&1)*((1<<j)|(1<<k))
            si=np.array([ix[int(s)] for s in ms^flip])
            if typ==0: rr=(si[:,None]*m+np.arange(m)[None,:]).ravel()
            elif typ==1:rr=(np.arange(m)[:,None]*m+si[None,:]).ravel()
            else:rr=(si[:,None]*m+si[None,:]).ravel()
            rows.append(rr)
        out.append(coo_matrix((np.ones(dim*L),(np.concatenate(rows),np.tile(cols,L))),shape=(dim,dim)).tocsr())
    return out

def solve(H):
    w,V=eigsh(H,k=2,which='SA',tol=2e-11,v0=np.random.default_rng(101).normal(size=H.shape[0]))
    i=np.argmin(w);E=float(w[i]);res=float(np.linalg.norm(H@V[:,i]-E*V[:,i]))
    return dict(E0=E,gap_same_sector=float(max(w)-min(w)),residual=res,dimension=H.shape[0])

def main():
    p=argparse.ArgumentParser();p.add_argument('--cluster',action='store_true');args=p.parse_args()
    if args.cluster and not os.environ.get('SLURM_JOB_ID'):raise SystemExit('Large mode requires a SLURM job; login-node execution refused.')
    lengths=[4,6,8,12] if args.cluster else [4,6,8]
    # Interpolate two Heisenberg legs to SU(4), including unequal leg speeds and opposite couplings.
    pars=[(1-c,1-c,c) for c in [0,.2,.4,.6,.8,1]]
    pars += [(a,b,c) for a,b,c in [(1,.2,.2),(1,.2,1),(1,0,1),(0,1,1),(-.2,-.2,1),(-.5,-.5,1),(.2,.8,-.2),(1,1,-.5)]]
    results=[]
    for a,b,c in pars:
        es=[a*s+b*t+c*s*t for s,t in itertools.product([-1,1],repeat=2)]
        results.append(dict(a=a,b=b,c=c,M=(max(es)-min(es))/2,rings=[]))
    for L in lengths:
        ops=operators(L)
        for r in results:
            H=r['a']*ops[0]+r['b']*ops[1]+r['c']*ops[2]
            rr=solve(H);rr['L']=L;r['rings'].append(rr)
    for r in results:
        e={x['L']:x['E0'] for x in r['rings']}
        for L1,L2 in [(4,8)]+([(8,12)] if args.cluster else []):
            A=-6/np.pi*(e[L1]/L1-e[L2]/L2)/(L1**-2-L2**-2)
            r[f'A_fit_{L1}_{L2}']=A;r[f'Q_fit_{L1}_{L2}']=A/(2*np.pi*r['M'])
    fname='adversarial_cluster.json' if args.cluster else 'adversarial_small.json'
    with open(fname,'w') as f:json.dump(results,f,indent=2)
    print('a b c M Q(4,8), finite-size diagnostic only')
    for r in results:print(r['a'],r['b'],r['c'],r['M'],r['Q_fit_4_8'])
    print('Largest residual',max(x['residual'] for r in results for x in r['rings']))
if __name__=='__main__':main()
