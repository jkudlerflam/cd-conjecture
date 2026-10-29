"""Small sparse ED checks. numpy/scipy only. OMP/BLAS threads should be 1."""
import json, itertools, math
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh
from scipy.special import digamma

def basis_counts(counts):
    L=sum(counts); out=[]
    def rec(prefix, left):
        if len(prefix)==L: out.append(tuple(prefix)); return
        for a,n in enumerate(left):
            if n:
                nxt=list(left); nxt[a]-=1; rec(prefix+[a],nxt)
    rec([],list(counts)); return out

def chain(counts,kind='swap',periodic=True,ts=None,J=0.):
    states=basis_counts(counts); index={s:i for i,s in enumerate(states)}
    L=sum(counts); d=len(counts); rows=[]; cols=[]; vals=[]
    for k,s in enumerate(states):
        diag=0.
        for i in range(L if periodic else L-1):
            j=(i+1)%L; a,b=s[i],s[j]; amp=0.
            if kind=='swap':
                if a==b: diag+=1.
                else: amp=1.
            elif kind=='tj':
                if a==0 and b==0: pass
                elif a==0 or b==0: amp=-ts[max(a,b)-1]
                elif a!=b: amp=J/2; diag-=J/2
            if amp:
                ss=list(s); ss[i],ss[j]=b,a
                rows.append(index[tuple(ss)]); cols.append(k); vals.append(amp)
        rows.append(k); cols.append(k); vals.append(diag)
    H=coo_matrix((vals,(rows,cols)),shape=(len(states),len(states))).tocsr()
    assert (H-H.T).nnz==0
    return H,states

def entropies(v,states,d,charge_map=None):
    L=len(states[0]); ans={}
    for ell in range(1,L//2+1):
        # Compress by the left charge sector before taking singular values.
        groups={}
        for z,s in zip(v,states):
            key=tuple(s[:ell].count(a) for a in range(d)) if charge_map is None else tuple(sum(charge_map[a][b] for a in s[:ell]) for b in range(len(charge_map[0])))
            groups.setdefault(key,[]).append((s[:ell],s[ell:],z))
        ent=0.
        for entries in groups.values():
            left={s:i for i,s in enumerate(sorted(set(x[0] for x in entries)))}
            right={s:i for i,s in enumerate(sorted(set(x[1] for x in entries)))}
            M=np.zeros((len(left),len(right)))
            for a,b,z in entries: M[left[a],right[b]]=z
            p=np.linalg.svd(M,compute_uv=False)**2; p=p[p>1e-15]
            ent-=np.dot(p,np.log(p))
        ans[ell]=float(ent)
    return ans

def main():
    rng=np.random.default_rng(76); out={'periodic':[]}
    for kind,d,Ls in [('swap',3,[6,9,12]),('swap',4,[4,8]),('tj',3,[6,9,12])]:
        for L in Ls:
            H,states=chain([L//d]*d,kind,ts=[1.]*(d-1),J=2.)
            w,v=eigsh(H,k=3,which='SA',tol=2e-11,v0=rng.normal(size=H.shape[0]))
            ix=np.argsort(w); w=w[ix]; v=v[:,ix]
            ent=entropies(v[:,0],states,d)
            ls=np.array([ell for ell in ent if ell>=2]); x=np.log(L/np.pi*np.sin(np.pi*ls/L))
            cfit=float(3*np.polyfit(x,[ent[int(ell)] for ell in ls],1)[0]) if len(ls)>=2 else None
            row={'model':kind,'d':d,'L':L,'dimension':H.shape[0],'energies':w.tolist(),'residual':float(np.linalg.norm(H@v[:,0]-w[0]*v[:,0])),'entropy':ent,'c_interval_fit_lmin2':cfit,'ground_degenerate':bool(w[1]-w[0]<1e-8)}
            if row['ground_degenerate']:
                row['c_interval_fit_lmin2']=None
                row['entropy_warning']='Arbitrary vector in a degenerate ground space; do not use for central charge.'
            if kind=='swap':
                bulk=1+2/d*(digamma(1/d)+np.euler_gamma); vel=2*np.pi/d
                row.update(bulk_energy=float(bulk),velocity=vel,c_casimir=float(-6*L*(w[0]-L*bulk)/(np.pi*vel)))
            out['periodic'].append(row)
            print(json.dumps(row),flush=True)
    # Pure hopping at fixed spin composition. Open boundaries make order conservation exact.
    H,states=chain([4,2,2],'tj',False,[1.,1.],0.)
    w=np.linalg.eigvalsh(H.toarray())
    spinless=np.array(sorted(sum(-2*np.cos(np.pi*np.array(ks)/9)) for ks in itertools.combinations(range(1,9),4)))
    err=float(np.max(np.abs(w-np.repeat(spinless,6))))
    from scipy.sparse.csgraph import connected_components
    k,labels=connected_components(H,directed=False)
    H2,_=chain([4,2,2],'tj',False,[1.,0.7],0.)
    k2,_=connected_components(H2,directed=False)
    out['frozen_order']={'L':8,'counts':[4,2,2],'dimension':len(states),'components_equal_t':int(k),'component_sizes':np.bincount(labels).tolist(),'components_unequal_t':int(k2),'full_spectrum_spinless_sixfold_error':err,'ground_multiplicity':int(np.sum(abs(w-w[0])<1e-9))}
    print(json.dumps(out['frozen_order']),flush=True)
    with open('ed_results.json','w') as f: json.dump(out,f,indent=2)
if __name__=='__main__': main()
