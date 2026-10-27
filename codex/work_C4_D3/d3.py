#!/usr/bin/env python3
"""Closed-book finite-size ED campaign. Only NumPy/SciPy; no external datasets.
One-site translations, periodic NN bonds, additive-charge and momentum sectors.
Run `python3 d3.py --help`. Set BLAS thread counts before launching.
"""
import os
for _key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ.setdefault(_key,'1')
import argparse, json, math, time
from pathlib import Path
import numpy as np
from scipy import sparse, linalg
from scipy.sparse.linalg import eigsh

VERSION = 1

def spin(s):
    m=np.arange(s,-s-1,-1); z=np.diag(m)
    plus=np.zeros((len(m),len(m)))
    for j in range(1,len(m)):
        plus[j-1,j]=np.sqrt(s*(s+1)-m[j]*(m[j]+1))
    return (plus+plus.T)/2,(plus-plus.T)/(2j),z

def projectors(s):
    ss=spin(s); d=int(2*s+1)
    x=sum(np.kron(a,a) for a in ss).real
    vals=[(j*(j+1)-2*s*(s+1))/2 for j in range(d)]
    ps=[]
    for j,v in enumerate(vals):
        p=np.eye(d*d)
        for k,w in enumerate(vals):
            if j!=k: p=p@(x-w*np.eye(d*d))/(v-w)
        p[np.abs(p)<1e-13]=0
        ps.append(p)
    return x,ps

def swap(d):
    return np.eye(d*d).reshape(d,d,d,d).transpose(1,0,2,3).reshape(d*d,d*d)

def local_components(family):
    if family in ('uls3','uls4'):
        d=int(family[-1]); return d,[swap(d)],np.eye(d,dtype=int)[:,:-1]
    if family in ('tb1','spin1'):
        x,p=projectors(1); sx,sy,sz=spin(1)
        if family=='tb1': hs=[2*p[1]+3*p[2]]
        else: hs=[x,x@x,np.kron(sz,sz),.5*(np.kron(sz@sz,np.eye(3))+np.kron(np.eye(3),sz@sz))]
        return 3,hs,np.array([[1],[0],[-1]])
    if family in ('tb32','su2','sp4'):
        x,p=projectors(1.5)
        if family=='tb32': hs=[sum(2*sum(1/k for k in range(1,j+1))*p[j] for j in range(4))]
        elif family=='su2': hs=p
        else:
            pp=swap(4); hs=[(np.eye(16)-pp)/2-p[0],(np.eye(16)+pp)/2]
        return 4,hs,np.array([[3],[1],[-1],[-3]])
    if family=='orbital':
        sx,sy,sz=spin(.5); ident=np.eye(2)
        sl=[np.kron(a,ident) for a in (sx,sy,sz)]
        tl=[np.kron(ident,a) for a in (sx,sy,sz)]
        a=sum(np.kron(b,b) for b in sl).real; b=sum(np.kron(c,c) for c in tl).real
        return 4,[a@b,a,b],np.array([[1,1],[1,-1],[-1,1],[-1,-1]])
    raise ValueError(family)

def target_charge(family,L,charges):
    if family.startswith('uls'):
        d=len(charges); n=np.full(d,L//d); n[:L%d]+=1; return n@charges
    if family=='orbital': return np.full(2,L%2,dtype=int)
    if family in ('spin1','tb1'): return np.array([0])
    return np.array([L%2])

def enumerate_basis(d,L,charges,target):
    # Meet-in-the-middle enumeration visits d^(L/2) half-words, not d^L.
    ell=L//2; other=L-ell
    def half(n):
        code=np.arange(d**n,dtype=np.int64);tmp=code.copy()
        q=np.zeros((len(code),charges.shape[1]),dtype=np.int16)
        for i in range(n):q+=charges[tmp%d];tmp//=d
        groups={}
        uq,inv=np.unique(q,axis=0,return_inverse=True)
        order=np.argsort(inv,kind='stable');counts=np.bincount(inv);starts=np.r_[0,np.cumsum(counts)]
        for i,qi in enumerate(uq):groups[tuple(qi)]=code[order[starts[i]:starts[i+1]]]
        return groups
    left=half(ell);right=half(other);out=[]
    for q,ls in left.items():
        rs=right.get(tuple(np.array(target)-q))
        if rs is not None:out.append((ls[:,None]+d**ell*rs[None,:]).ravel())
    return np.sort(np.concatenate(out)) if out else np.empty(0,dtype=np.int64)

class Sector:
    def __init__(self,family,L,target=None,cartan=False,cache_blocks=True):
        self.family=family; self.L=L; self.cache_blocks=cache_blocks
        self.d,self.hs,self.charges=local_components(family)
        if cartan:
            if family!='sp4': raise ValueError('Cartan option is for Sp4')
            self.charges=np.array([[1,0],[0,1],[0,-1],[-1,0]])
            if target is None: target=[L%2,0]
        self.target=target_charge(family,L,self.charges) if target is None else np.array(target,dtype=int)
        self.basis=enumerate_basis(self.d,L,self.charges,self.target)
        b=self.basis; rot=b.copy(); rep=b.copy(); ashift=np.zeros(len(b),dtype=np.int16)
        period=np.zeros(len(b),dtype=np.int16)
        for a in range(1,L+1):
            rot=rot//self.d+(rot%self.d)*self.d**(L-1)
            ix=rot<rep; rep[ix]=rot[ix]; ashift[ix]=a
            ix=(rot==b)&(period==0); period[ix]=a
        self.reps=b[rep==b]; self.rid=np.searchsorted(self.reps,rep)
        self.period=period; self.shift=(period-ashift)%period
        self.rperiod=period[rep==b]
        self._blocks={}
    def blocks(self,k):
        k=k%self.L
        if k in self._blocks: return self._blocks[k]
        L=self.L; d=self.d
        allowed=(k*self.rperiod)%L==0
        ids=np.flatnonzero(allowed); reps=self.reps[ids]
        inverse=np.full(len(self.reps),-1,dtype=np.int64); inverse[ids]=np.arange(len(ids))
        n=len(ids); ans=[]
        if n==0:
            self._blocks[k]=(ids,[sparse.csr_matrix((0,0)) for h in self.hs])
            return self._blocks[k]
        for h in self.hs:
            rows=[]; cols=[]; data=[]
            for site in range(L):
                j=(site+1)%L; pi=d**site; pj=d**j
                ai=(reps//pi)%d; aj=(reps//pj)%d
                pair=ai*d+aj
                for inp in range(d*d):
                    ci=np.flatnonzero(pair==inp)
                    if not len(ci): continue
                    a,b=divmod(inp,d)
                    for op in np.flatnonzero(np.abs(h[:,inp])>1e-12):
                        c,e=divmod(int(op),d)
                        new=reps[ci]+(c-a)*pi+(e-b)*pj
                        bi=np.searchsorted(self.basis,new)
                        if np.any(bi>=len(self.basis)) or np.any(self.basis[bi]!=new):
                            raise ValueError('Local h violates requested charge')
                        ri=inverse[self.rid[bi]]; use=ri>=0
                        v=h[op,inp]*np.sqrt(self.rperiod[ids[ci[use]]]/self.period[bi[use]])
                        if k: v=v*np.exp(2j*np.pi*k*self.shift[bi[use]]/L)
                        rows.append(ri[use]); cols.append(ci[use]); data.append(v)
            H=sparse.coo_matrix((np.concatenate(data),(np.concatenate(rows),np.concatenate(cols))),shape=(n,n)).tocsr() if data else sparse.csr_matrix((n,n))
            H.sum_duplicates(); H.eliminate_zeros()
            err=H-H.getH()
            if err.nnz and np.max(np.abs(err.data))>2e-10: raise AssertionError('Non-Hermitian momentum block')
            if k==0 or 2*k==L: H=H.real
            ans.append(H)
        self._blocks[k]=(ids,ans)
        return ids,ans
    def solve(self,coeff,k,n=3):
        ids,hs=self.blocks(k)
        H=sum(float(c)*h for c,h in zip(coeff,hs))
        dim=H.shape[0]
        if dim==0: return np.empty(0),np.empty((0,0)),0
        if dim<=max(24,n+1):
            e,v=linalg.eigh(H.toarray(),check_finite=False); e=e[:n];v=v[:,:n]
        else:
            rng=np.random.default_rng(20261021+k+100*self.L)
            init=rng.normal(size=dim)
            if np.iscomplexobj(H.data):init=init+1j*rng.normal(size=dim)
            e,v=eigsh(H,k=min(n,dim-1),which='SA',tol=2e-10,maxiter=12000,v0=init,ncv=min(dim,max(30,3*n+3)))
            ix=np.argsort(e); e=e[ix];v=v[:,ix]
        residual=float(np.max(np.linalg.norm(H@v-v*e,axis=0)))
        return e,v,residual
    def expand(self,k,v):
        ids=np.flatnonzero((k*self.rperiod)%self.L==0); rv=np.zeros(len(self.reps),dtype=complex);rv[ids]=v
        psi=rv[self.rid]*np.exp(-2j*np.pi*k*self.shift/self.L)/np.sqrt(self.period.astype(float))
        if abs(np.vdot(psi,psi)-1)>1e-8:raise AssertionError('Expansion normalization')
        return psi
    def entropies(self,psi):
        ans=[]; d=self.d; L=self.L
        for ell in range(1,L//2+1):
            left=self.basis%d**ell;right=self.basis//d**ell
            lu=np.unique(left);ru=np.unique(right)
            q=np.zeros((len(lu),self.charges.shape[1]),dtype=int);t=lu.copy()
            for i in range(ell):q+=self.charges[t%d];t//=d
            _,group=np.unique(q,axis=0,return_inverse=True)
            li=np.searchsorted(lu,left); weights=[]
            for g in range(group.max()+1):
                sel=group[li]==g
                ls=left[sel];rs=right[sel];us=np.unique(ls);vs=np.unique(rs)
                mat=np.zeros((len(us),len(vs)),dtype=complex)
                mat[np.searchsorted(us,ls),np.searchsorted(vs,rs)]=psi[sel]
                sv=linalg.svdvals(mat,check_finite=False); weights.extend((sv*sv).tolist())
            w=np.array(weights);w=w[w>1e-15]
            if abs(w.sum()-1)>1e-8:raise AssertionError('Schmidt norm')
            ans.append(float(-np.dot(w,np.log(w))))
        return ans

def entropy_fit(S,L,cut=1):
    ell=np.arange(1,len(S)+1);mask=ell>=cut
    x=np.log((L/np.pi)*np.sin(np.pi*ell[mask]/L));y=np.array(S)[mask]
    if len(y)<2:return None
    a=np.stack([x,np.ones(len(x))],axis=1)
    fit=np.linalg.lstsq(a,y,rcond=None)[0];res=y-a@fit
    se=None
    if len(y)>2:se=float(3*np.sqrt((res@res)/(len(y)-2)*np.linalg.inv(a.T@a)[0,0]))
    return {'c':float(3*fit[0]),'regression_se':se,'rms':float(np.sqrt(np.mean(res*res))),'cut':cut}

def measure(sec,coeff,label,levels=3,entropy=True,momenta=None):
    start=time.monotonic();L=sec.L
    ks=range(L//2+1) if momenta is None else momenta
    spectra={};solutions={};residuals=[]
    for k in ks:
        e,v,res=sec.solve(coeff,k,levels)
        if len(e):spectra[str(k)]=e.tolist();solutions[k]=v;residuals.append(res)
        if not sec.cache_blocks:sec._blocks.clear()
    k0=min((int(k) for k in spectra),key=lambda k:spectra[str(k)][0]);e0=spectra[str(k0)][0]
    def fold(k):k=k%L;return min(k,L-k)
    adjacent=sorted(set(fold(k0+s) for s in (-1,1)))
    gaps=[spectra[str(k)][0]-e0 for k in adjacent if str(k) in spectra]
    gap=float(min(gaps)) if gaps else None
    v=float(L*gap/(2*np.pi)) if gap is not None else None
    out={'family':sec.family,'L':L,'label':label,'coeff':list(map(float,coeff)),'charge':sec.target.tolist(),
         'local_charges':sec.charges.tolist(),'sector_dimension':len(sec.basis),'ground_k':k0,'E0':e0,'spectra':spectra,'adjacent_gap':gap,'v_current_assumed':v,
         'max_residual':max(residuals),'momentum_ground_degeneracy':sum(abs(e[0]-e0)<1e-8 for e in spectra.values())}
    if entropy:
        S=sec.entropies(sec.expand(k0,solutions[k0][:,0]));out['entropy']=S
        out['entropy_fits']=[entropy_fit(S,L,c) for c in (1,2,3)]
    if v is not None and v>1e-7:
        out['tower_x_assumed']={k:[float((e-e0)*L/(2*np.pi*v)) for e in es] for k,es in spectra.items()}
    out['seconds']=time.monotonic()-start
    return out

def append(path,result):
    with open(path,'a') as f:f.write(json.dumps(result)+'\n')
    ent=result.get('entropy_fits',[None,None]);cf=ent[1] if len(ent)>1 else None
    print(json.dumps({k:result[k] for k in ('label','L','ground_k','E0','v_current_assumed','seconds')})+f' c_ent2={None if cf is None else cf["c"]}',flush=True)

def sphere(n):
    # Orthonormal coordinates perpendicular to the energy-shift direction.
    B=np.array([[1,-1,0,0],[1,1,-2,0],[1,1,1,-3]],float)
    B/=np.linalg.norm(B,axis=1)[:,None]
    for i in range(n):
        z=1-2*(i+.5)/n;p=i*np.pi*(3-np.sqrt(5));r=np.sqrt(1-z*z)
        yield (np.array([r*np.cos(p),r*np.sin(p),z])@B).tolist()

def cases(family,n):
    if family.startswith('uls') or family.startswith('tb'): return [(family,[1.])]
    if family=='sp4':return [(f'sp4_t{i:03d}',[np.cos(2*np.pi*i/n),np.sin(2*np.pi*i/n)]) for i in range(n)]
    if family=='su2':return [(f'su2_s{i:03d}',c) for i,c in enumerate(sphere(n))]
    if family=='orbital':
        vals=[-.75,-.25,0.,.25,.5,1.]
        return [(f'orb_x{x:g}_y{y:g}',[1.,y,x]) for i,x in enumerate(vals) for y in vals[i:]]+ [('orb_decoupled',[0,1,1]),('orb_decoupled_unequal',[0,1,2])]
    if family=='spin1':
        out=[(f'bb_t{i:03d}',[np.cos(2*np.pi*i/n),np.sin(2*np.pi*i/n),0,0]) for i in range(n)]
        for t in (np.pi/4,np.pi/3,5*np.pi/12):
            for delta,D in ((-.5,0),(.5,0),(0,-.5),(0,.5),(-.5,.5),(.5,-.5),(1,1),(-1,1)):
                out.append((f's1_t{t:.5f}_z{delta:g}_D{D:g}',[np.cos(t),np.sin(t),delta,D]))
        return out
    raise ValueError(family)

def selftest():
    checks={}
    for s in (1,1.5):
        x,ps=projectors(s);d=int(2*s+1)
        checks[f'projector_sum_{s}']=float(np.linalg.norm(sum(ps)-np.eye(d*d)))
        checks[f'projector_orthogonal_{s}']=max(float(np.linalg.norm(a@b-(a if i==j else 0))) for i,a in enumerate(ps) for j,b in enumerate(ps))
        checks[f'projector_ranks_{s}']=[int(round(np.trace(p))) for p in ps]
    for fam in ('uls3','sp4','su2','orbital','spin1'):
        L=4;sec=Sector(fam,L);rng=np.random.default_rng(4);c=rng.normal(size=len(sec.hs));h=sum(a*b for a,b in zip(c,sec.hs))
        b=sec.basis;H=np.zeros((len(b),len(b)))
        for col,code in enumerate(b):
            for i in range(L):
                j=(i+1)%L;a=code//sec.d**i%sec.d;bb=code//sec.d**j%sec.d
                for op in np.flatnonzero(np.abs(h[:,a*sec.d+bb])>1e-12):
                    cc,ee=divmod(int(op),sec.d);new=code+(cc-a)*sec.d**i+(ee-bb)*sec.d**j
                    H[np.searchsorted(b,new),col]+=h[op,a*sec.d+bb]
        ev=[];err=0
        for k in range(L):
            ids,hs=sec.blocks(k);hk=sum(a*hh for a,hh in zip(c,hs)).toarray()
            ev.extend(linalg.eigvalsh(hk).tolist())
            U=np.stack([sec.expand(k,np.eye(len(ids))[:,i]) for i in range(len(ids))],axis=1)
            err=max(err,float(np.max(np.abs(U.conj().T@H@U-hk))))
        checks[f'full_momentum_matrix_{fam}']=err
        checks[f'full_spectrum_{fam}']=float(np.max(np.abs(np.sort(ev)-linalg.eigvalsh(H))))
    assert max(v for v in checks.values() if isinstance(v,float))<1e-10,checks
    return checks

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--family',choices=['uls3','uls4','tb1','tb32','sp4','su2','orbital','spin1'])
    p.add_argument('--sizes',type=int,nargs='+',default=[8]);p.add_argument('--n',type=int,default=24)
    p.add_argument('--output',default='results.jsonl');p.add_argument('--levels',type=int,default=3)
    p.add_argument('--no-cache-momenta',action='store_true');p.add_argument('--sp4-cartan',action='store_true');p.add_argument('--charge',type=int,nargs='+');p.add_argument('--coeff',type=float,nargs='+');p.add_argument('--label',default='custom')
    p.add_argument('--case-indices',type=int,nargs='+');p.add_argument('--continue-on-error',action='store_true');p.add_argument('--no-entropy',action='store_true');p.add_argument('--selftest',action='store_true')
    a=p.parse_args()
    if a.selftest:
        checks=selftest();Path(a.output).write_text(json.dumps(checks,indent=2)+'\n');print(json.dumps(checks,indent=2));return
    if not a.family:p.error('--family required')
    cs=[(a.label,a.coeff)] if a.coeff is not None else cases(a.family,a.n)
    if a.case_indices is not None:cs=[cs[i] for i in a.case_indices]
    for L in a.sizes:
        sec=Sector(a.family,L,a.charge,a.sp4_cartan,not a.no_cache_momenta)
        print(f'Built {a.family} L={L} sector dimension {len(sec.basis)}',flush=True)
        for label,coeff in cs:
            try:
                append(a.output,measure(sec,coeff,label,a.levels,not a.no_entropy))
            except Exception as exc:
                if not a.continue_on_error: raise
                failure={'family':a.family,'L':L,'label':label,'coeff':list(map(float,coeff)),'error':repr(exc)}
                with open(a.output+'.errors','a') as f:f.write(json.dumps(failure)+'\n')
                print('FAILED '+json.dumps(failure),flush=True)
if __name__=='__main__':main()
