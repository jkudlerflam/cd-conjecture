"""Small ED for H=XX-YY+u ZZ+b Z, after alternate-site X rotation.
Rotated H=XX+YY-u ZZ+(-1)^j b Z. Pauli normalization, periodic even L.
Only numpy, scipy, and the Python standard library are used.
"""
import argparse, json, time
from pathlib import Path
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh


def sector(L, N, u, b):
    allstates = np.arange(1 << L, dtype=np.int64)
    states = allstates[np.bitwise_count(allstates) == N]
    z = 1-2*((states[:, None] >> np.arange(L)) & 1)
    diag = -u*np.sum(z*np.roll(z,-1,axis=1),axis=1)+b*(z@((-1.)**np.arange(L)))
    rows=[np.arange(len(states))]; cols=[rows[0]]; data=[diag]
    for j in range(L):
        k=(j+1)%L
        ids=np.flatnonzero(z[:,j]!=z[:,k])
        dest=states[ids]^((1<<j)|(1<<k))
        rows.append(np.searchsorted(states,dest)); cols.append(ids)
        data.append(np.full(len(ids),2.))
    H=coo_matrix((np.concatenate(data),(np.concatenate(rows),np.concatenate(cols))),shape=(len(states),len(states))).tocsr()
    return H, states


def low(H,k):
    if H.shape[0]<=max(64,k+2):
        e,v=np.linalg.eigh(H.toarray());e,v=e[:k],v[:,:k]
    else:
        rng=np.random.default_rng(731)
        e,v=eigsh(H,k=min(k,H.shape[0]-2),which='SA',tol=2e-11,v0=rng.normal(size=H.shape[0]))
        ix=np.argsort(e);e,v=e[ix],v[:,ix]
    res=np.max(np.linalg.norm(H@v-v*e,axis=0))
    return e,v,float(res)


def entropy_profile(psi,states,L):
    full=np.zeros(1<<L);full[states]=psi
    profiles=[]
    for start in [0,1]:
        # Cyclic relabeling is only a permutation of tensor factors.
        ids=np.arange(1<<L,dtype=np.int64)
        if start: ids=((ids<<1)&((1<<L)-1))|(ids>>(L-1))
        vec=full[ids]
        ss=[]
        for ell in range(1,L//2+1):
            mat=vec.reshape((1<<(L-ell),1<<ell))
            vals=np.linalg.eigvalsh(mat.T@mat)
            vals=vals[vals>1e-15]
            ss.append(float(-np.sum(vals*np.log(vals))))
        profiles.append(ss)
    ss=np.mean(profiles,axis=0)
    ell=np.arange(1,L//2+1); use=ell>=max(2,L//4)
    x=np.log((L/np.pi)*np.sin(np.pi*ell[use]/L))
    fit=np.polyfit(x,ss[use],1)
    return ss.tolist(),float(3*fit[0]),float(np.max(np.abs(np.polyval(fit,x)-ss[use])))


def momenta(e,v,states,L):
    shift=((states<<2)&((1<<L)-1))|(states>>(L-2))
    ids=np.searchsorted(states,shift)
    out=[];start=0
    while start<len(e):
        end=start+1
        while end<len(e) and abs(e[end]-e[start])<2e-7: end+=1
        vv=v[:,start:end]
        # T versus T inverse changes only the momentum sign.
        T=vv.T@vv[ids,:]
        zz=np.linalg.eigvals(T)
        for z in zz:
            out.append({'energy':float(np.mean(e[start:end])), 'phase_T2':float(np.angle(z)), 'mod_T2':float(abs(z))})
        start=end
    return out


def run(L,u,b,k=12):
    t=time.time();H,s=sector(L,L//2,u,b);e,v,r=low(H,k)
    ent,cS,efit=entropy_profile(v[:,0],s,L)
    modes=momenta(e,v,s,L)
    p0=modes[0]['phase_T2']; target=4*np.pi/L
    phonons=[m['energy']-e[0] for m in modes if m['mod_T2']>.999 and abs(abs(np.angle(np.exp(1j*(m['phase_T2']-p0))))-target)<1e-5]
    speed=L*min(phonons)/(2*np.pi) if phonons else None
    charges={}
    for q in [1,2]:
        HH,_=sector(L,L//2+q,u,b);ee,_,rr=low(HH,1)
        charges[str(q)]=float(ee[0]-e[0]);r=max(r,rr)
    same=[m['energy']-e[0] for m in modes if m['energy']-e[0]>1e-8 and m['mod_T2']>.999 and abs(np.angle(np.exp(1j*(m['phase_T2']-p0))))<1e-5]
    return dict(L=L,u=u,b=b,E0=float(e[0]),c_chord=cS,entropy_fit_max_error=efit,S=ent,
                gaps=(e-e[0]).tolist(),modes=modes,v_phonon=speed,charge_gaps=charges,
                neutral_same_T2=min(same) if same else None,residual=r,seconds=time.time()-t)


def check():
    X=np.array([[0.,1.],[1.,0.]]);Y=np.array([[0.,-1j],[1j,0.]]);Z=np.diag([1.,-1.]);I=np.eye(2)
    def op(terms,L):
        out=np.array([[1.]])
        for j in reversed(range(L)):out=np.kron(out,terms.get(j,I))
        return out
    L=6;u=.83;b=1.2
    H=np.zeros((1<<L,1<<L),complex)
    Hrot=np.zeros_like(H)
    for j in range(L):
        k=(j+1)%L
        H+=op({j:X,k:X},L)-op({j:Y,k:Y},L)+u*op({j:Z,k:Z},L)+b*op({j:Z},L)
        Hrot+=op({j:X,k:X},L)+op({j:Y,k:Y},L)-u*op({j:Z,k:Z},L)+(-1)**j*b*op({j:Z},L)
    es=[];matrix_error=0.
    for n in range(L+1):
        HH,s=sector(L,n,u,b);matrix_error=max(matrix_error,float(np.max(np.abs(HH.toarray()-Hrot[np.ix_(s,s)]))));es.extend(np.linalg.eigvalsh(HH.toarray()))
    spectral_error=float(np.max(np.abs(np.sort(es)-np.linalg.eigvalsh(H))))
    # Entire spectrum and exact upper-wing ground energy and gap.
    rr=np.sqrt(1+b*b/4);uu=rr+.2;es=[]
    for n in range(L+1):
        HH,_=sector(L,n,uu,b);es.extend(np.linalg.eigvalsh(HH.toarray()))
    es=np.sort(es);wing_error=max(abs(es[0]+uu*L),abs(es[2]-es[0]-.8))
    # XX ring ground energy for Pauli XX+YY.
    HH,_=sector(8,4,0,0);ee,_,_=low(HH,1);xx_error=abs(ee[0]+4/np.sin(np.pi/8))
    report=dict(matrix_error=matrix_error,rotation_spectrum_error=spectral_error,wing_error=float(wing_error),xx_energy_error=float(xx_error))
    assert max(report.values())<1e-10,report
    return report

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--out',default='scan.jsonl');ap.add_argument('--L',nargs='+',type=int,default=[8,10,12,14]);ap.add_argument('--points',nargs='+',default=['0,0']);ap.add_argument('--k',type=int,default=12);args=ap.parse_args()
    print(json.dumps({'checks':check()}),flush=True)
    with Path(args.out).open('a') as f:
        for p in args.points:
            u,b=map(float,p.split(','))
            for L in args.L:
                row=run(L,u,b,args.k);f.write(json.dumps(row)+'\n');f.flush()
                print(json.dumps({key:row[key] for key in ['L','u','b','E0','c_chord','v_phonon','charge_gaps','neutral_same_T2','residual','seconds']}),flush=True)
