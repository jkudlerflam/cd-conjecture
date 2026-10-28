"""Rational Sp(2n) periodic ground energies in the zero Cartan sector."""
import time,json
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh

def ground(d,L,zero_cartan=True):
    t0=time.time();n=d//2;powers=d**np.arange(L,dtype=np.int64)
    ids=np.arange(d**L,dtype=np.int64)
    charges=np.zeros((n,len(ids)),dtype=np.int8)
    for power in powers:
        digit=(ids//power)%d
        for a in range(n):charges[a]+=(digit==a).astype(np.int8)-(digit==a+n).astype(np.int8)
    if zero_cartan: ids=ids[np.all(charges==0,axis=0)]
    del charges
    digits=np.array([(ids//p)%d for p in powers],dtype=np.int8)
    rows=[];cols=[];vals=[];ix=np.arange(len(ids))
    for i in range(L):
        j=(i+1)%L;a=digits[i];b=digits[j]
        dest=ids+(b.astype(np.int64)-a)*powers[i]+(a.astype(np.int64)-b)*powers[j]
        rows.append(np.searchsorted(ids,dest));cols.append(ix);vals.append(np.ones(len(ids)))
        mask=(b==(a+n)%d);src=ix[mask];aa=a[mask];bb=b[mask]
        ss=np.where(aa<n,1.,-1.)
        for c in range(d):
            e=(c+n)%d;sgn=1 if c<n else -1
            dest=ids[mask]+(c-aa.astype(np.int64))*powers[i]+(e-bb.astype(np.int64))*powers[j]
            rr=np.searchsorted(ids,dest)
            assert np.array_equal(ids[rr],dest)
            rows.append(rr);cols.append(src);vals.append(ss*sgn/(n+1))
    mat=coo_matrix((np.concatenate(vals),(np.concatenate(rows),np.concatenate(cols))),shape=(len(ids),len(ids))).tocsr()
    assert (mat-mat.T).nnz==0
    ev,vec=eigsh(mat,k=1,which='SA',tol=2e-11,v0=np.random.default_rng(35).normal(size=len(ids)))
    residual=np.linalg.norm(mat@vec[:,0]-ev[0]*vec[:,0])
    row={'d':d,'L':L,'zero_cartan':zero_cartan,'sector_dimension':len(ids),'E0_sector':float(ev[0]),'residual':float(residual),'seconds':time.time()-t0}
    print(json.dumps(row),flush=True);return row

out={'assumption':'Zero-Cartan sector contains the global ground state. No independent proof for these rings supplied.','energies':[],'fits':[]}
for d,sizes in [(4,[4,6,8,10]),(6,[4,6,8]),(8,[4,6])]:
    data=[ground(d,L) for L in sizes];out['energies']+=data
    for subset in [data, data[-3:], data[-2:]]:
        if len(subset)<2:continue
        L=np.array([r['L'] for r in subset],dtype=float);y=np.array([r['E0_sector'] for r in subset])/L
        # Two-coefficient fit even at 3+ sizes; a 3-coefficient fit separately.
        X=np.array([np.ones(len(L)),1/L**2]).T
        z=np.linalg.lstsq(X,y,rcond=None)[0];A=-6*z[1]/np.pi
        out['fits'].append({'d':d,'sizes':L.tolist(),'form':'e+a/L^2','A':float(A),'q':float(A/(np.pi*np.log2(d)))})
    if len(data)>=3:
        subset=data[-3:];L=np.array([r['L'] for r in subset],dtype=float);y=np.array([r['E0_sector'] for r in subset])/L
        z=np.linalg.solve(np.array([np.ones(len(L)),1/L**2,1/L**4]).T,y);A=-6*z[1]/np.pi
        out['fits'].append({'d':d,'sizes':L.tolist(),'form':'e+a/L^2+b/L^4','A':float(A),'q':float(A/(np.pi*np.log2(d)))})
    with open('sp_ed.json','w') as f:json.dump(out,f,indent=2)
out['full_space_checks']=[ground(d,L,False) for d,L in [(4,4),(4,6),(6,4),(6,6),(8,4)]]
for row in out['full_space_checks']:
    ref=next(x for x in out['energies'] if x['d']==row['d'] and x['L']==row['L'])
    assert abs(row['E0_sector']-ref['E0_sector'])<1e-9
out['fit_verdict']='REJECT all mixed-length fits as Casimir estimates: strong size-class drift; apparent violations are not evidence.'
a,b=[x for x in out['energies'] if x['d']==4 and x['L'] in (4,8)]
A=-6/np.pi*(a['E0_sector']/4-b['E0_sector']/8)/(1/16-1/64)
out['same_length_class_diagnostic']={'d':4,'sizes':[4,8],'A':A,'q':A/(2*np.pi),'status':'Two small sizes only; not a thermodynamic estimate.'}
with open('sp_ed.json','w') as f:json.dump(out,f,indent=2)
print(json.dumps(out['same_length_class_diagnostic'],indent=2))
