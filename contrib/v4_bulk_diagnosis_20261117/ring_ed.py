#!/usr/bin/env python3
"""Literal periodic V=4 spectra and translation-restored static correlations."""
import argparse,time,resource,sys
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh
from threadpoolctl import threadpool_limits
from phase_common import HERE,Z,write_json


def lattice(L,P):
    states=np.array([s for s in range(2**L) if s.bit_count()%2==P],dtype=np.int64)
    inv=np.full(2**L,-1,dtype=np.int64);inv[states]=np.arange(len(states))
    bits=np.array([(states>>(L-1-j))&1 for j in range(L)])
    h=Z.pauli_bond(Jx=.5,Jy=-.5,Jz=1.,D=1/16);data=[];rows=[];cols=[]
    diag=-1.75*np.sum(1-2*bits,axis=0).astype(complex);terms=[]
    for j in range(L):
        k=(j+1)%L;s=2*bits[j]+bits[k];d=np.diag(h)[s];a=h[3-s,s]
        dest=inv[states^((1<<(L-1-j))|(1<<(L-1-k)))];assert min(dest)>=0
        diag+=d;data.append(a);rows.append(dest);cols.append(np.arange(len(states)))
        terms.append((d,a,dest))
    data.append(diag);rows.append(np.arange(len(states)));cols.append(np.arange(len(states)))
    H=coo_matrix((np.concatenate(data),(np.concatenate(rows),np.concatenate(cols))),shape=(len(states),len(states))).tocsr()
    assert np.max(abs((H-H.getH()).data),initial=0)<1e-13
    return H,states,inv,bits,terms


def restored(raw,means):
    L=len(means);loc=raw-np.outer(means,means)
    rr=np.array([np.mean([raw[j,(j+r)%L] for j in range(L)]) for r in range(L)])
    cc=rr-np.mean(means)**2
    qs=2*np.pi*np.arange(L//2+1)/L
    sc=np.real(np.exp(1j*np.outer(qs,np.arange(L)))@cc)
    # Explicit decomposition before restoration, which need not have flat means.
    phase=np.exp(1j*np.outer(qs,np.arange(L)))
    sl=np.real(np.einsum('qi,ij,qj->q',phase.conj(),loc,phase,optimize=True))/L
    sd=abs(phase@(means-means.mean()))**2/L
    assert np.max(abs(sc-(sl+sd)))<1e-8
    k=1+np.argmax(sc[1:])
    return dict(profile=means.tolist(),raw=raw.tolist(),connected=loc.tolist(),
        translation_averaged_raw=rr.tolist(),translation_restored_connected=cc.tolist(),
        q_grid=qs.tolist(),S_restored=sc.tolist(),S_connected_before_restoring=sl.tolist(),
        S_disconnected_before_restoring=sd.tolist(),peak=dict(q=float(qs[k]),q_over_pi=float(qs[k]/np.pi),
        S=float(sc[k]),S_over_L=float(sc[k]/L)))


def run(L,P):
    start=time.time();H,states,inv,bits,terms=lattice(L,P)
    rng=np.random.default_rng(19+P);v0=rng.normal(size=len(states))
    vals,vec=eigsh(H,k=4,which='SA',tol=1e-11,v0=v0,maxiter=20000,ncv=32)
    order=np.argsort(vals);vals=vals[order];vec=vec[:,order]
    residual=[float(np.linalg.norm(H@vec[:,i]-vals[i]*vec[:,i])) for i in range(4)]
    assert max(residual)<1e-8
    v=vec[:,0];prob=abs(v)**2;dens=bits@prob;raw_n=(bits*prob)@bits.T
    out=dict(L=L,P=P,boundary='periodic',energies=vals.tolist(),gaps_within_parity=(vals-vals[0]).tolist(),
        residuals=residual,density=restored(raw_n,dens))
    for name in ('bond','kinetic'):
        vv=[]
        for d,a,dest in terms:
            w=np.zeros_like(v);w[dest]=a*v
            if name=='bond':w+=d*v
            vv.append(w)
        vv=np.array(vv);means=np.real(vv.conj()@v);raw=np.real(vv.conj()@vv.T)
        out[name]=restored(raw,means)
    pv=[]
    for j,(_,_,dest) in enumerate(terms):
        w=np.zeros_like(v);w[dest]=2*(bits[j]!=bits[(j+1)%L])*v;pv.append(w)
    pv=np.array(pv);pm=np.real(pv.conj()@v);pr=np.real(pv.conj()@pv.T)
    out['valley_bond']=restored(pr,pm)
    out['valley_bond']['uniform_order_weight']=float(out['valley_bond']['S_restored'][0]/L+np.mean(pm)**2)
    av=[]
    for j in range(L):
        k=(j+2)%L;sel=(bits[j]==1)&(bits[k]==1)
        dest=inv[states[sel]^((1<<(L-1-j))|(1<<(L-1-k)))]
        w=np.zeros_like(v);w[dest]=v[sel];av.append(w)
    av=np.array(av);am=np.array([np.vdot(v,w) for w in av]);araw=av.conj()@av.T
    aconn=araw-np.outer(am.conj(),am)
    out['pair']=dict(mean_re=am.real.tolist(),mean_im=am.imag.tolist(),raw_re=araw.real.tolist(),raw_im=araw.imag.tolist(),
        connected_re=aconn.real.tolist(),connected_im=aconn.imag.tolist(),
        distance_average_connected=[[float(np.mean([aconn[j,(j+r)%L] for j in range(L)]).real),
            float(np.mean([aconn[j,(j+r)%L] for j in range(L)]).imag)] for r in range(L)])
    translated=((states<<1)&(2**L-1))|(states>>(L-1));tv=np.zeros_like(v);tv[inv[translated]]=v
    ov=np.vdot(v,tv);out['translation_overlap']=[float(ov.real),float(ov.imag)]
    overlaps=[];momenta=[]
    for j in range(4):
        t=np.zeros_like(v);t[inv[translated]]=vec[:,j];o=np.vdot(vec[:,j],t)
        assert np.linalg.norm(H@t-vals[j]*t)<1e-8
        overlaps.append([float(o.real),float(o.imag)])
        momenta.append(float(np.angle(o)) if abs(abs(o)-1)<1e-6 else None)
    out['translation_overlap_by_level']=overlaps;out['momentum_if_translation_eigenstate']=momenta
    out['ground_degeneracy_within_four']=int(np.sum(vals-vals[0]<1e-9))
    out['runtime_seconds']=time.time()-start
    out['peak_RSS_GB']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024)/1e9
    out['scope']='One lowest eigenvector per parity; static results explicitly translation-averaged. Small-ring gaps do not determine a bulk velocity or c.'
    write_json(HERE/'periodic'/f'L{L}_P{P}.json',out)
    print('done',L,P,vals.tolist(),out['density']['peak'],'seconds',out['runtime_seconds'],flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--L',type=int,nargs='+',default=[8,12,16]);args=ap.parse_args()
    with threadpool_limits(limits=1):
        for L in args.L:
            for P in (0,1):run(L,P)
