"""C13. Python standard library, numpy, scipy only. All outputs local.
Physical projected mobile-cell graph, not the microscopic wrapper.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import itertools as it, json, time, argparse
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh
from scipy.linalg import svdvals

def masks(M,N):
    return [sum(1<<j for j in c) for c in it.combinations(range(M),N)]
def bits(x,M): return np.array([(x>>j)&1 for j in range(M)],dtype=int)
def patterns(b,M,r):
    L=b*M+r; out=set()
    for eta in masks(M,r):
        d=b+bits(eta,M); xs=np.r_[0,np.cumsum(d)[:-1]]
        for x in range(L): out.add(tuple(sorted((xs+x)%L)))
    return sorted(out)
def sea_amp(mask,M):
    pos=np.flatnonzero(bits(mask,M)); N=len(pos)
    if N==0: return 1.
    ks=2*np.pi*(np.arange(N)-(N-1)/2)/M
    return abs(np.linalg.det(np.exp(1j*np.outer(pos,ks))/np.sqrt(M)))
def basis(b,M,r,Ns):
    pats=patterns(b,M,r); states=[]
    for xs in pats:
        for zs in it.product(*[masks(M,N) for N in Ns]): states.append((xs,zs))
    return states

def matrix(b,M,r,Ns,generic):
    states=basis(b,M,r,Ns); lookup={s:i for i,s in enumerate(states)}
    L=b*M+r; n=len(Ns); hops=[1.,.73,.51][:n+1]
    rows=[];cols=[];vals=[]; phase=np.zeros(len(states)); psi=np.zeros(len(states))
    for col,(xs,zs) in enumerate(states):
        eta=(np.roll(xs,-1)-xs)%L-b
        occ=[bits(z,M) for z in zs]
        phase[col]=.37*sum(np.dot(eta,occ[a])/(a+1) for a in range(n))
        emask=sum(int(e)<<j for j,e in enumerate(eta))
        psi[col]=np.sqrt(M/L)*sea_amp(emask,M)*np.prod([sea_amp(z,M) for z in zs])
        diag=0.
        if generic:
            for a in range(n):
                diag+=.07/(a+1)*np.dot(eta-.37,occ[a]-[.27,.16][a])
                diag+=.04*np.dot(occ[a],np.roll(occ[a],-1))
            diag+=.05*np.dot(eta,np.roll(eta,-1))
            if n==2: diag+=.035*np.dot(occ[0],occ[1])
        rows.append(col);cols.append(col);vals.append(diag)
        for j in range(M):
            jm=(j-1)%M
            # Move marker j right or left, preserving cell order and logical labels.
            for direction in [-1,1]:
                allowed=(eta[jm]==0 and eta[j]==1) if direction==1 else (eta[jm]==1 and eta[j]==0)
                if not allowed: continue
                raw=list(xs);raw[j]=(raw[j]+direction)%L
                order=np.argsort(raw); nx=tuple(raw[k] for k in order)
                nz=tuple(sum(int(occ[a][k])<<l for l,k in enumerate(order)) for a in range(n))
                amp=-hops[0]*(1+.025*sum((occ[a][jm]+occ[a][j]-2*[.27,.16][a])/(a+1) for a in range(n))) if generic else -hops[0]
                rows.append(lookup[(nx,nz)]);cols.append(col);vals.append(amp)
            k=(j+1)%M
            for a in range(n):
                if occ[a][j]==occ[a][k]: continue
                nz=list(zs);nz[a]^=(1<<j)|(1<<k)
                amp=-hops[a+1]*(1+.03*(eta[j]+eta[k]-2*.37)) if generic else -hops[a+1]
                rows.append(lookup[(xs,tuple(nz))]);cols.append(col);vals.append(amp)
    H=coo_matrix((vals,(rows,cols)),shape=(len(states),len(states))).tocsr()
    return states,H,psi,phase,hops

def entropy(states,vec,L,ell,n):
    # Explicit physical local alphabet: empty=0, marker with payload=1+payload bits.
    # Block by conserved marker and payload counts in A.
    groups={}; w=n+1
    for (xs,zs),amp in zip(states,vec):
        A=B=0; qa=[0]*(n+1)
        for j,x in enumerate(xs):
            val=1+sum(((z>>j)&1)<<a for a,z in enumerate(zs))
            if x<ell:
                A|=val<<(w*x);qa[0]+=1
                for a,z in enumerate(zs): qa[a+1]+=(z>>j)&1
            else: B|=val<<(w*(x-ell))
        groups.setdefault(tuple(qa),[]).append((A,B,amp))
    S=0.; tr2=0.;maxblock=0
    for entries in groups.values():
        aa={a:i for i,a in enumerate(sorted({v[0] for v in entries}))}
        bb={b:i for i,b in enumerate(sorted({v[1] for v in entries}))}
        C=np.zeros((len(aa),len(bb)),dtype=complex)
        for a,b,z in entries:C[aa[a],bb[b]]=z
        p=svdvals(C)**2;p=p[p>1e-15]
        S-=np.dot(p,np.log(p));tr2+=np.dot(p,p);maxblock=max(maxblock,min(C.shape))
    return float(S),float(-np.log(tr2)),maxblock

def graph_tests():
    out=[]
    for n in [1,2]:
        for M,r in [(6,2),(7,3),(8,3)]:
            Ns=[2,1][:n]; b=2;L=b*M+r
            for generic in [False,True]:
                t=time.time();states,H,psi,phase,hops=matrix(b,M,r,Ns,generic)
                es,vs=eigsh(H,k=3,which='SA',tol=2e-11,v0=np.random.default_rng(13013).normal(size=len(psi)))
                order=np.argsort(es);es=es[order];v=vs[:,order[0]]
                exact=-sum(2*h*np.sin(np.pi*N/M)/np.sin(np.pi/M) for h,N in zip(hops,[r]+Ns))
                ell=list(range(max(2,L//4),L//2+1)); ent=[]; entU=[]
                for l in ell:
                    ent.append(entropy(states,v,L,l,n))
                    if not generic: entU.append(entropy(states,v*np.exp(1j*phase),L,l,n))
                x=np.log((L/np.pi)*np.sin(np.pi*np.array(ell)/L))
                cS=3*np.polyfit(x,np.array(ent)[:,0],1)[0]
                cR=4*np.polyfit(x,np.array(ent)[:,1],1)[0]
                row=dict(n=n,b=b,M=M,r=r,Ns=Ns,L=L,generic=generic,dim=len(states),energies=es.tolist(),gap=float(es[1]-es[0]),exact_free_energy=exact,free_residual=float(np.linalg.norm(H@psi-exact*psi)) if not generic else None,ell=ell,entropies=ent,c_entropy=float(cS),c_renyi2_diagnostic=float(cR),elapsed=time.time()-t)
                if entU: row.update(unitary_entropies=entU,c_unitary_entropy=float(3*np.polyfit(x,np.array(entU)[:,0],1)[0]))
                out.append(row);print(json.dumps({k:row[k] for k in ['n','M','L','generic','dim','gap','c_entropy','elapsed']}),flush=True)
                Path('graph_results.json').write_text(json.dumps(out,indent=2))
    return out

from pathlib import Path
def correlation_entropy(M,N,ell):
    ks=2*np.pi*(np.arange(N)-(N-1)/2)/M
    v=np.exp(1j*np.outer(np.arange(ell),ks))/np.sqrt(M)
    p=np.clip(np.linalg.eigvalsh(v@v.conj().T),1e-15,1-1e-15)
    return float(-np.sum(p*np.log(p)+(1-p)*np.log(1-p)))
def controls():
    out={}; b=4; fills=np.array([np.sqrt(2)/4,np.sqrt(3)/7,np.sqrt(5)/11]); ts=np.array([1,.73,.51])
    for n in [0,1,2]:
        rows=[]
        for M in [40,64,96,128,192]:
            Ns=np.rint(M*fills[:n+1]).astype(int);r=Ns[0];L=b*M+r; a=L/M
            vs=2*ts[:n+1]*np.sin(np.pi*Ns/M)*a
            E=-sum(2*ts[:n+1]*np.sin(np.pi*Ns/M)/np.sin(np.pi/M))
            bulk=-sum(2*ts[:n+1]*np.sin(np.pi*Ns/M))*M/np.pi
            ccas=(n+1)*(-6*L*(E-bulk)/(np.pi*sum(vs)))
            ells=np.arange(M//4,M//2+1)
            ent=np.array([sum(correlation_entropy(M,int(N),int(l)) for N in Ns) for l in ells])
            x=np.log(M/np.pi*np.sin(np.pi*ells/M));cs=3*np.polyfit(x,ent,1)[0]
            rows.append(dict(M=M,L=int(L),counts=Ns.tolist(),velocities=vs.tolist(),c_casimir=float(ccas),c_cell_entropy=float(cs)))
        out[str(n)]=rows
    # Minimize ALL charge sectors for every ring, including primes. Exact formula.
    unique=[]
    for n in [1,2]:
        fs=fills[:n+1];t=ts[:n+1];mu=2*t[1:]*np.cos(np.pi*fs[1:]);nu=fs[0]
        # grand payload energy per cell minimized over payload fillings
        epay=sum(-2*t[1:]*np.sin(np.pi*fs[1:])/np.pi+mu*fs[1:])
        kappa=2*t[0]/np.pi*(np.sin(np.pi*nu)-np.pi*(b+nu)*np.cos(np.pi*nu))-epay
        for L in range(40,401):
            candidates=[]
            for M in range((L+b)//(b+1),L//b+1):
                if M<2:continue
                r=L-b*M
                if not 0<=r<=M: continue
                ee=kappa*M-2*t[0]*np.sin(np.pi*r/M)/np.sin(np.pi/M)
                species=[]
                for g,u in zip(t[1:],mu):
                    vals=[(-2*g*np.sin(np.pi*N/M)/np.sin(np.pi/M)+u*N,N) for N in range(M+1)]
                    species.append(sorted(vals)[:2])
                for choice in it.product(*species):
                    candidates.append((ee+sum(z[0] for z in choice),M,r,*[z[1] for z in choice]))
            candidates.sort(); best=candidates[0];second=candidates[1]
            unique.append(dict(n=n,L=L,sector=list(best[1:]),sector_gap=second[0]-best[0],interior=0<best[2]<best[1]))
    out['all_ring_charge_minima']=unique
    Path('controls.json').write_text(json.dumps(out,indent=2))
    print('CONTROL SUMMARY',json.dumps({str(n):out[str(n)][-1] for n in [0,1,2]}),flush=True)
    print('minimum charge-sector gap',min(z['sector_gap'] for z in unique),flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--controls-only',action='store_true');a=p.parse_args()
    controls()
    if not a.controls_only:graph_tests()
