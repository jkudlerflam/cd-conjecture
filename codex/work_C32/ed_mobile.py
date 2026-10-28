"""ED with active lengths and ALL legal core excitations, fast-wrapper limit.

This is not ED of the full 2736-state microscopic bond. It tests the
controlled fast-wrapper carrier Hamiltonian defined in RESULT.md.
Usage: python3 -B ed_mobile.py n cells. Each case is saved immediately.
"""
import itertools
import json
import sys
import time
import signal
from pathlib import Path
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh
import mobile as m


def run(n,M):
    start=time.monotonic()
    c=m.core(n,g=.45,mass=.5,express=.2,loop=.08)
    cs={p:m.charge(c['k'],p,.16) for p in (1,2)}
    hc={p:m.flat_model(c,ch) for p,ch in cs.items()}
    measures={p:m.measure(cc) for p,cc in hc.items()}
    bonds={(p,q):m.bond(measures[p],measures[q]) for p in cs for q in cs}
    pats=[tuple(2 if j==a else 1 for j in range(M)) for a in range(M)]
    ip={p:i for i,p in enumerate(pats)}
    d=c['d']; N=2**n
    nbrH={p:[[(int(j),hc[p]['H'][j,i]) for j in np.flatnonzero(abs(hc[p]['H'][:,i])>1e-14)] for i in range(d)] for p in cs}
    nbrF={side:{f:[[(int(j),F[j,i]) for j in np.flatnonzero(F[:,i])] for i in range(d)] for f,F in fs.items()} for side,fs in c['F'].items()}
    out={'status':'numerical','model':'fast-wrapper limit, all legal core excitations',
         'n':n,'cells':M,'core_dimension':d,'unrestricted_dimension':M*d**M,
         'length_patterns':[list(x) for x in pats],'sectors':[]}
    path=Path(f'ed_n{n}_M{M}.json')
    if path.exists():
        saved=json.loads(path.read_text())
        out['sectors']=saved.get('sectors',[])
    for q in (0,1):
        old=next((a for a in out['sectors'] if a['charge']==q),None)
        if old is not None and len(old['runs'])==3: continue
        configs=[x for x in itertools.product(range(d),repeat=M)
                 if sum(c['Q'][i] for i in x)==q and sum(c['ex'][i] for i in x)%2==0]
        idx={x:i for i,x in enumerate(configs)}; D=len(configs); DD=D*M
        code=[x for x in itertools.product(range(N),repeat=M) if sum(int(i).bit_count() for i in x)==q]
        nq=len(code); nqtot=nq*M
        G=np.zeros((DD,nqtot),complex)
        rows0=[]; cols0=[]; vals0=[]; rowsv=[]; colsv=[]; valsv=[]
        rowsu=[]; colsu=[]; valsu=[]
        for a,p in enumerate(pats):
            for ix,x in enumerate(configs):
                src=a*D+ix
                for z,bits in enumerate(code):
                    amp=np.prod([c['G'][x[j],bits[j]] for j in range(M)])
                    G[src,a*nq+z]=amp
                for site in range(M):
                    for j,v in nbrH[p[site]][x[site]]:
                        xx=x[:site]+(j,)+x[site+1:]
                        rows0.append(a*D+idx[xx]); cols0.append(src); vals0.append(v)
                for site in range(M-1):
                    for fr,fl,v in m.TERMS:
                        for j,vr in nbrF['R'][fr][x[site]]:
                            for k,vl in nbrF['L'][fl][x[site+1]]:
                                xx=x[:site]+(j,k)+x[site+2:]
                                rowsv.append(a*D+idx[xx]); colsv.append(src); valsv.append(v*vr*vl)
                    if p[site]!=p[site+1]:
                        pp=p[:site]+(p[site+1],p[site])+p[site+2:]
                        rowsu.append(ip[pp]*D+ix); colsu.append(src); valsu.append(-.2)
        def sparse(rr,cc,vv):return sp.coo_matrix((vv,(rr,cc)),shape=(DD,DD)).tocsr()
        H0=sparse(rows0,cols0,vals0); V=sparse(rowsv,colsv,valsv); U=sparse(rowsu,colsu,valsu)
        pred=np.zeros((nqtot,nqtot),complex)
        for a,p in enumerate(pats):
            K=np.zeros((N**M,N**M),complex)
            for site in range(M-1):
                K+=np.kron(np.kron(np.eye(N**site),bonds[p[site],p[site+1]]),np.eye(N**(M-site-2)))
            inds=[sum(bits[j]*N**(M-1-j) for j in range(M)) for bits in code]
            pred[a*nq:(a+1)*nq,a*nq:(a+1)*nq]=K[np.ix_(inds,inds)]
        pred+=G.conj().T@(U@G)
        rec={'charge':q,'dimension':DD,'code_dimension':nqtot,
             'code_residual':float(np.linalg.norm(H0@G)),
             'first_order_residual':float(np.linalg.norm(G.conj().T@(V@G))),
             'runs':[]}
        if old is not None: rec=old
        else: out['sectors'].append(rec)
        for lam in (.12,.06,.03):
            if any(abs(a['lambda']-lam)<1e-12 for a in rec['runs']): continue
            H=H0+lam*V+lam**2*U
            # Bound each individual diagonalization to minutes.
            def timeout(signum,frame): raise TimeoutError('ED exceeded 180 seconds')
            signal.signal(signal.SIGALRM,timeout); signal.alarm(180)
            if DD<300:
                e,W=np.linalg.eigh(H.toarray()); e=e[:nqtot+1]; W=W[:,:nqtot+1]
            else:
                e,W=eigsh(H,k=nqtot+1,which='SA',tol=2e-12,
                          ncv=min(DD-1,6*(nqtot+1)+20),maxiter=2000)
                order=np.argsort(e); e=e[order]; W=W[:,order]
            signal.alarm(0)
            resid=float(np.linalg.norm(H@W-W*e,axis=0).max())
            A=G.conj().T@W[:,:nqtot]; u,sv,vh=np.linalg.svd(A); rot=u@vh
            eff=(rot*e[:nqtot])@rot.conj().T
            err=float(np.linalg.norm(eff-lam**2*pred,2))
            rr={'lambda':lam,'error':err,'error_over_lambda4':err/lam**4,
                'max_eigenpair_residual':resid,'gap_above_retained_band':float(e[nqtot]-e[nqtot-1]),
                'minimum_code_overlap_singular_value':float(min(sv))}
            rec['runs'].append(rr); print(n,M,q,DD,rr,flush=True)
            path.write_text(json.dumps(out,indent=2)+'\n')
    out['seconds']=time.monotonic()-start
    Path(f'ed_n{n}_M{M}.json').write_text(json.dumps(out,indent=2)+'\n')
    print('elapsed',out['seconds'],flush=True)

if __name__=='__main__':run(int(sys.argv[1]),int(sys.argv[2]))
