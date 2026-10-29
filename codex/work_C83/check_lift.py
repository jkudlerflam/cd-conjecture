"""C83 full literal-bond checks, invariant closures and leakage tests."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import itertools,json
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh
import bond14 as b

def closure(seed,p,lam,colored):
    ls=b.LET if colored else b.MAT+[('V',0)];d=len(ls)
    ix={a:i for i,a in enumerate(ls)}
    h=b.bond(p,lam=lam,colored=colored)
    trans={i:[(int(j),h[j,i]) for j in np.flatnonzero(h[:,i])] for i in range(d*d)}
    seed=tuple(ix[a] for a in seed);states=[seed];index={seed:0}
    rows=[];cols=[];vals=[];pen=[];bad=[];cur=0;L=len(seed)
    while cur<len(states):
        w=states[cur]
        for i,a in enumerate(w):
            j=(i+1)%L
            for q,v in trans[a*d+w[j]]:
                z=list(w);z[i],z[j]=divmod(q,d);z=tuple(z)
                if z not in index:index[z]=len(states);states.append(z)
                rows.append(index[z]);cols.append(cur);vals.append(v)
        pen.append(sum((colored or (a<12 and w[(i+1)%L]<12))
                       and not b.allowed(ls[a],ls[w[(i+1)%L]]) for i,a in enumerate(w)))
        mat=[a for a in w if a<12]
        bad.append(sum(not b.allowed(ls[a],ls[mat[(i+1)%len(mat)]]) for i,a in enumerate(mat)))
        cur+=1
        if len(states)>25000:raise RuntimeError('small ED cap')
    H=sp.coo_matrix((vals,(rows,cols)),shape=(cur,cur)).tocsc()
    assert np.max(abs((H-H.getH()).data),initial=0)<1e-13
    if colored:assert np.all(np.array(bad)<=pen)
    return H,np.array(pen),np.array(bad),[tuple(ls[a] for a in w) for w in states]

def run(H,pen,bad,K):
    A=H+sp.diags(K*pen.astype(float),format='csc')
    sigma=-float(np.max(np.asarray(abs(H).sum(axis=1))))-1
    es,vs=eigsh(A,k=1,sigma=sigma,which='LM',tol=1e-12,ncv=35,v0=np.ones(len(pen)))
    e,v=es[0],vs[:,0]
    rec=dict(K=K,dimension=len(pen),energy=float(e),residual=float(np.linalg.norm(A@v-e*v)),
             illegal_material_weight=float(np.linalg.norm(v[bad>0])**2),
             extended_bad_expectation=float(np.vdot(v,pen*v).real))
    assert rec['residual']<1e-8
    return rec

def main():
    import bond_shift
    out=dict(parent_checks=[],commutators={},runs=[])
    for n in (2,13):
        p=b.params(n,K=37.,musel=.93,field=.021,flux=.37,theta=-.21)
        h=b.material_bond(p,.17);ref=bond_shift.bond(p,12,lam=.17)
        err=float(np.max(abs(h-ref)));assert err<1e-13
        assert all(b.allowed(a,c)==bond_shift.allowed(a,c,12) for a in b.MAT for c in b.MAT)
        out['parent_checks'].append(dict(n=n,entry_error=err))
    h=b.bond(b.params(K=37.,musel=.93,field=.021),lam=.17)
    assert np.max(abs(h-h.conj().T))==0
    for name,q in [('N',[int(a[0]!='V') for a in b.LET]),
                   ('M',[int(b.marker(a)) for a in b.LET]),('Q',[b.charge(a) for a in b.LET])]:
        qb=np.add.outer(q,q).ravel()
        err=float(np.max(abs((qb[:,None]-qb[None,:])*h)))
        assert err==0;out['commutators'][name]=err
    count=0
    for w in itertools.product(b.LET,repeat=4):
        mat=[a for a in w if a[0]!='V']
        if not mat:continue
        nb=sum(not b.allowed(a,w[(i+1)%4]) for i,a in enumerate(w))
        bm=sum(not b.allowed(a,mat[(i+1)%len(mat)]) for i,a in enumerate(mat))
        assert bm<=nb;count+=1
    out['four_site_grammar_words']=count
    cases=[('two_cells_zero_bits',[('L',0),('R',0),('R',0),('e',0)]*2+[('V',0)],0.),
           ('full_one_vacancy',[('L',1),('R',0),('e',0),('V',0)],.1),
           ('full_two_vacancies',[('L',1),('R',0),('e',0),('V',0),('V',0)],.1)]
    for name,seed,lam in cases:
        mats={}
        for colored in (False,True):
            p=b.params(n=2,K=0.)
            ss=[a if a[0]!='V' else ('V',0) for a in seed]
            H,pen,bad,states=closure(ss,p,lam,colored)
            mats[colored]=(H,pen,bad,states)
            for K in (10.,100.,1000.,10000.):
                rec=dict(case=name,colored=colored,**run(H,pen,bad,K))
                out['runs'].append(rec);print(json.dumps(rec),flush=True)
        A,pa,ba,sa=mats[True];B,pb,bb,sb=mats[False]
        ia=np.flatnonzero(pa==0);ib=np.flatnonzero(bb==0)
        lookup={tuple(a if a[0]!='V' else ('V',0) for a in sa[i]):i for i in ia}
        assert len(lookup)==len(ia)==len(ib)
        order=[lookup[sb[i]] for i in ib]
        diff=A[order][:,order]-B[ib][:,ib]
        assert np.max(abs(diff.data),initial=0)<1e-13
        out.setdefault('projections',[]).append(dict(case=name,legal_dimension=len(ia),entry_error=0.))
    with open('lift_checks.json','w') as f:json.dump(out,f,indent=2)
if __name__=='__main__':main()
