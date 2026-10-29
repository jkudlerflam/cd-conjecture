"""C94 finite quantum-clock and response checks. Outputs only in this folder."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
import json, math, itertools
import numpy as np
import scipy.sparse as sp
from scipy.linalg import eigh
from scipy.special import logsumexp, gammaln
from inherited_ed import bath, selector, clock_model
from inherited_second_order import source_vectors, susceptibility
HERE=Path(__file__).resolve().parent

def sources(lengths, mu):
    n=sum(lengths); states=list(itertools.product(*(range(b) for b in lengths)))
    ix={s:i for i,s in enumerate(states)}; d=len(states)
    hs=np.zeros((n,d,d)); off=0
    for c,b in enumerate(lengths):
        mark=off+b-1
        hs[(mark-1)%n]+=mu/2*np.eye(d);hs[mark]+=mu/2*np.eye(d)
        for q in range(b-1):
            for s,i in ix.items():
                if s[c]==q:
                    z=list(s);z[c]+=1;k=ix[tuple(z)]
                    hs[off+q,i,k]-=1;hs[off+q,k,i]-=1
        off+=b
    return hs

def coeff(lengths, pins, mu, chi):
    _,y=source_vectors(lengths,pins,mu)
    return float(np.sum(chi*(y.T@y)))

def canonical():
    rows=[]; kernels=[]; cache={}
    plans=[(2,[2,2],4,[([4],[False]),([2,2],[True,False])]),
           (3,[3,3],6,[([2,4],[False,False]),([2,2,2],[False]*3),([6],[False]),([3,3],[True,False])]),
           (3,[3,3,3],4,[([2,4,3],[False]*3),([6,3],[False]*2)]),
           (15,[15,15],1,[([14,16],[False]*2),([30],[False])])]
    for b,good,v,defects in plans:
        n=sum(good);mu=selector(b);chi,res=susceptibility(n,v)
        dist=np.minimum(np.arange(n),n-np.arange(n))
        kernels.append(dict(N=n,V=v,linear_residual=res,row=chi[0].tolist(),
            norm0=float(np.abs(chi[0]).sum()),norm1=float(np.dot(1+dist,np.abs(chi[0]))),
            eigen_min=float(np.linalg.eigvalsh(chi)[0])))
        bg=coeff(good,[False]*len(good),mu,chi)
        for j in [.04,.02,.01,.005,.0025]:
            ref=clock_model(good,[False]*len(good),v,j,mu,cache)
            for lengths,pins in defects:
                out=clock_model(lengths,pins,v,j,mu,cache)
                w=sum(l for l in lengths if l!=b)+sum(p for l,p in zip(lengths,pins) if l==b)
                db=coeff(lengths,pins,mu,chi)-bg
                dr=out['dressing_residual']-ref['dressing_residual']
                rows.append(dict(b=b,N=n,V=v,J=j,lengths=lengths,pins=pins,W=w,
                    delta_B=db,delta_R=dr,epsilon_observed=max(0.,-dr)/(j*w),
                    absolute_normalized=abs(dr)/(j*w),quadratic_error=(dr+j*j*db)/(j*j*w),
                    energy_cost=out['energy']-ref['energy'],dimension=max(ref['dimension'],out['dimension']),
                    residual=max(ref['eigenpair_residual'],out['eigenpair_residual'])))
            print('canonical',b,n,v,j,flush=True)
    return rows,kernels

def bath_spectral(n,v):
    h,a,psi,e,alpha,res=bath(n,v)
    en,un=eigh(h.toarray())
    mats=np.array([un.T@(a[:,i,None]*un) for i in range(n)])
    return en,mats

def thermal_free(spectra, temp, tilt):
    vals=np.concatenate([ev+tilt*v for v,ev in enumerate(spectra)])
    return -temp*logsumexp(-vals/temp)

def response_zero(n, sectors, vmax, temp, tilt):
    # At u=0, D is the integrated connected bath kernel contracted
    # with the normalized material Hilbert-Schmidt Gram matrix.
    spec=[bath_spectral(n,v) for v in range(vmax+1)]
    f=thermal_free([p[0] for p in spec],temp,tilt)
    a0=0.; k=np.zeros((n,n))
    for v,(en,aa) in enumerate(spec):
        shifted=en+tilt*v
        w=np.exp(-(shifted-f)/temp)
        a0+=np.dot(w,np.diagonal(aa[0]))
        de=en[None,:]-en[:,None]
        dw=w[:,None]-w[None,:]
        near=np.abs(de)<1e-11
        weight=np.zeros_like(de)
        np.divide(dw,de,out=weight,where=~near)
        weight[near]=np.broadcast_to(w[:,None]/temp,de.shape)[near]
        k+=np.einsum('ab,iab,jab->ij',weight,aa,aa)
    k-=a0*a0/temp
    return [float(np.sum(k*np.einsum('iab,jba->ij',hs,hs)/hs.shape[1])) for hs in sectors],a0

def thermal():
    n=4;b=4;mu=selector(b);lengthsets=[[4],[2,2]];vmax=9
    hs=[sources(ls,mu) for ls in lengthsets]
    js=[.04,.02,.01,.005,.0025];s=.65
    allspec={j:[[],[]] for j in js};free=[];contacts=[]
    eigres=0.;maxdim=0
    for v in range(vmax+1):
        h,a,_,_,_,_=bath(n,v); free.append(eigh(h.toarray(),eigvals_only=True));contacts.append(a)
        for sector,hh in enumerate(hs):
            dim=hh.shape[1];base=sp.kron(h,sp.eye(dim),format='csr')
            interaction=sum((sp.kron(sp.diags(a[:,i]),sp.csr_matrix(hh[i]),format='csr') for i in range(n)))
            for j in js:
                mat=(base+j*interaction).toarray()
                ev=eigh(mat,eigvals_only=True)
                allspec[j][sector].append(ev);maxdim=max(maxdim,len(ev))
        print('thermal blocks',v,flush=True)
    rows=[]
    for j in js:
        for ratio in [.1,1.]:
            temp=j*ratio;fg=thermal_free(free,temp,s)
            # Free A expectation uses exact free eigenvectors in each vacancy block.
            a0=0.
            for v in range(vmax+1):
                ev,aa=bath_spectral(n,v)
                a0+=np.dot(np.exp(-(ev+s*v-fg)/temp),np.diagonal(aa[0]))
            rr=[];tails=[]
            first=vmax+1
            logfree=gammaln(first+n)-gammaln(first+1)-gammaln(n)-(s*first-2*n-fg)/temp
            free_tail=math.exp(logfree)/(1-(first+n)/(first+1)*math.exp(-s/temp))
            for sec,hh in enumerate(hs):
                fm=thermal_free(allspec[j][sec],temp,s)
                fmat=-temp*logsumexp(-j*a0*np.linalg.eigvalsh(hh.sum(axis=0))/temp)
                rr.append(fm-fg-fmat)
                # Tail bound uses ||V|| <= sum_i ||h_i||, E_gap >= -2N.
                norm=sum(np.linalg.norm(x,2) for x in hh)
                first=vmax+1
                logfirst=math.log(hh.shape[1])+gammaln(first+n)-gammaln(first+1)-gammaln(n)
                logfirst-=(s*first-2*n-j*norm-fm)/temp
                ratio_tail=(first+n)/(first+1)*math.exp(-s/temp)
                tails.append(float(math.exp(logfirst)/(1-ratio_tail)))
            dr=rr[1]-rr[0]
            rows.append(dict(J=j,T=temp,T_over_J=ratio,tilt=s,Vmax=vmax,N=n,W=n,
                delta_R=dr,epsilon_observed=max(0.,-dr)/(j*n),absolute_normalized=abs(dr)/(j*n),
                delta_R_over_J2W=dr/(j*j*n),relative_partition_tail_bound=max(tails),free_relative_partition_tail_bound=free_tail,
                delta_R_cutoff_error_bound=temp*sum(tails)+j*sum(np.linalg.norm(h.sum(axis=0),2) for h in hs)*free_tail))
    # A free vacancy crossing in C92's tilt interval.
    def e(v): return -2*math.sin(math.pi*n/(n+v))/math.sin(math.pi/(n+v))
    crossing=e(4)-e(5)
    _,_,_,_,a4,_=bath(4,4);_,_,_,_,a5,_=bath(4,5)
    moments=[float(np.trace(x.sum(axis=0)@x.sum(axis=0))/x.shape[1]) for x in hs]
    slope=(a4-a5)**2/4*(moments[1]-moments[0])
    cross=[]
    for temp in [.01,.003,.001,.0003,.0001]:
        ds,a0=response_zero(n,hs,7,temp,crossing)
        cross.append(dict(T=temp,tilt=crossing,D_difference=ds[1]-ds[0],T_D_difference=temp*(ds[1]-ds[0]),predicted_limit=slope))
    return rows,cross,dict(max_dimension=maxdim,crossing_slope=slope,crossing_tilt=crossing,mu=mu)

def main():
    rows,kernels=canonical()
    (HERE/'canonical_results.json').write_text(json.dumps(dict(rows=rows,kernels=kernels),indent=2)+'\n')
    therm,cross,meta=thermal()
    (HERE/'thermal_results.json').write_text(json.dumps(dict(rows=therm,crossing=cross,metadata=meta),indent=2)+'\n')
    assert max(r['residual'] for r in rows)<2e-9
    assert min(k['eigen_min'] for k in kernels)>-1e-9
    assert all(r['energy_cost']>0 for r in rows)
    assert abs(cross[-1]['T_D_difference']/meta['crossing_slope']-1)<.1
    print('ALL CHECKS PASS',flush=True)
if __name__=='__main__':main()
