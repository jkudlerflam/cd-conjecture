#!/usr/bin/env python3
"""Closed-book nearest-neighbor qubit-chain ED. Pauli, not spin, couplings.

All generated files are relative to the working directory. See report.md for
assumptions and limitations. Set OPENBLAS_NUM_THREADS=1 for predictable runtime.
"""
import os
for _key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_key] = '1'
import argparse, json, time, math
from pathlib import Path
from functools import lru_cache
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh
from scipy.optimize import minimize
from numba import njit

PAULI = np.array([[[0,1],[1,0]], [[0,-1j],[1j,0]], [[1,0],[0,-1]]], complex)

def density(J, field, field_split=.5):
    """Symmetric field density, sum_j h_j is the requested Hamiltonian."""
    h = np.zeros((4,4), complex)
    for a in range(3):
        h += field[a]*(field_split*np.kron(PAULI[a],np.eye(2))+(1-field_split)*np.kron(np.eye(2),PAULI[a]))
        for b in range(3):
            h += J[a,b]*np.kron(PAULI[a],PAULI[b])
    return h

@njit(cache=True)
def _orbits(L):
    N=1<<L
    oid=np.full(N,-1,np.int32); shift=np.zeros(N,np.int8)
    reps=np.empty(N,np.int64); periods=np.empty(N,np.int8)
    nr=0; mask=N-1
    for s in range(N):
        if oid[s]>=0: continue
        t=s; p=0
        while True:
            oid[t]=nr; shift[t]=p; p+=1
            t=((t<<1)&mask)|(t>>(L-1))
            if t==s: break
        reps[nr]=s; periods[nr]=p; nr+=1
    return reps[:nr].copy(),periods[:nr].copy(),oid,shift

@lru_cache(maxsize=8)
def orbits(L): return _orbits(L)

@njit(cache=True)
def _popcounts(reps):
    out=np.zeros(len(reps),np.int8)
    for i in range(len(reps)):
        s=reps[i]
        while s:
            out[i]+=s&1; s>>=1
    return out

def basis(L,m,nup=None):
    reps,p,oid,shift=orbits(L)
    allowed=(p.astype(np.int64)*m)%L==0
    if nup is not None: allowed &= _popcounts(reps)==nup
    active=np.flatnonzero(allowed)
    inv=np.full(len(reps),-1,np.int32); inv[active]=np.arange(len(active))
    return active,inv

@njit(cache=True)
def _sector_entries(L,h,reps,p,oid,shift,src,inv,mdst,n):
    cap=len(src)*L*4
    rows=np.empty(cap,np.int32); cols=np.empty(cap,np.int32)
    vals=np.empty(cap,np.complex128); count=0
    for col in range(len(src)):
        r=src[col]; s=reps[r]
        for j in range(L):
            jp=(j+1)%L; bj=(s>>j)&1; bp=(s>>jp)&1; a=2*bj+bp
            phase=np.exp(2j*np.pi*n*j/L)
            for b in range(4):
                v=h[b,a]
                if abs(v)<1e-14: continue
                sp=s^((bj^(b>>1))<<j)^((bp^(b&1))<<jp)
                r2=oid[sp]; row=inv[r2]
                if row<0: continue
                rows[count]=row; cols[count]=col
                vals[count]=v*phase*np.sqrt(p[r]/p[r2])*np.exp(2j*np.pi*mdst*shift[sp]/L)
                count+=1
    return rows[:count],cols[:count],vals[:count]

def sector_operator(L,h,m,n=0,nup=None):
    """h_n maps translation eigenvalue exp(2 pi i m/L) to m-n."""
    mdst=(m-n)%L
    src,_=basis(L,m,nup); dst,inv=basis(L,mdst,nup)
    reps,p,oid,shift=orbits(L)
    row,col,val=_sector_entries(L,h,reps,p,oid,shift,src,inv,mdst,n)
    H=sparse.coo_matrix((val,(row,col)),shape=(len(dst),len(src))).tocsr()
    if n==0 and (not len(H.data) or np.max(np.abs(H.data.imag))<1e-12): H=H.real
    return H

@njit(cache=True)
def _full_entries(L,h):
    N=1<<L; cap=N*L*4
    rows=np.empty(cap,np.int32); cols=np.empty(cap,np.int32)
    vals=np.empty(cap,np.complex128); count=0
    for s in range(N):
        for j in range(L):
            jp=(j+1)%L; bj=(s>>j)&1; bp=(s>>jp)&1; a=2*bj+bp
            for b in range(4):
                v=h[b,a]
                if abs(v)<1e-14: continue
                sp=s^((bj^(b>>1))<<j)^((bp^(b&1))<<jp)
                rows[count]=sp; cols[count]=s; vals[count]=v; count+=1
    return rows[:count],cols[:count],vals[:count]

def full_operator(L,h):
    row,col,val=_full_entries(L,h)
    H=sparse.coo_matrix((val,(row,col)),shape=(1<<L,1<<L)).tocsr()
    return H.real if np.max(np.abs(H.data.imag),initial=0)<1e-13 else H

def solve(H,k=4,tol=2e-10,v0=None):
    if H.shape[0]<=max(k+2,24):
        e,u=np.linalg.eigh(H.toarray()); return e[:k],u[:,:k]
    if v0 is None:
        rng=np.random.default_rng(127)
        v0=rng.normal(size=H.shape[0])
        if np.iscomplexobj(H.data): v0=v0+1j*rng.normal(size=H.shape[0])
    e,u=eigsh(H,k=k,which='SA',tol=tol,v0=v0,maxiter=20000,ncv=max(2*k+1,24))
    ix=np.argsort(e); return e[ix],u[:,ix]

def expand(L,m,u,nup=None):
    src,inv=basis(L,m,nup); reps,p,oid,shift=orbits(L)
    idx=inv[oid]; good=idx>=0
    psi=np.zeros(1<<L,complex)
    psi[good]=u[idx[good]]*np.exp(-2j*np.pi*m*shift[good]/L)/np.sqrt(p[oid[good]].astype(float))
    return psi

def translate(psi,L):
    s=np.arange(1<<L); dest=((s<<1)&((1<<L)-1))|(s>>(L-1))
    out=np.empty_like(psi); out[dest]=psi
    return out

def resolve_ground(e,u,L):
    """Pick a translation eigenstate in any ground eigenspace returned by ARPACK."""
    ix=np.flatnonzero(e-e[0]<1e-8); U=u[:,ix].astype(complex)
    U=np.linalg.qr(U)[0]
    TU=np.column_stack([translate(U[:,i],L) for i in range(U.shape[1])])
    w,v=np.linalg.eig(U.conj().T@TU)
    j=np.argmax(abs(w)); psi=U@v[:,j]; psi/=np.linalg.norm(psi)
    z=np.vdot(psi,translate(psi,L)); m=int(round(np.angle(z)*L/(2*np.pi)))%L
    return psi,m,float(1-abs(z)),len(ix)

def entropies(psi,L,nup=None):
    if np.max(np.abs(psi.imag))<1e-13: psi=psi.real
    result=[]
    for ell in range(1,L//2+1):
        A=psi.reshape(1<<(L-ell),1<<ell)
        if nup is None:
            s=np.linalg.svd(A,compute_uv=False)
        else:
            # Charge blocks avoid factoring the many structural zeros in A.
            rr=np.fromiter((i.bit_count() for i in range(A.shape[0])),np.int8)
            cc=np.fromiter((i.bit_count() for i in range(A.shape[1])),np.int8)
            parts=[]
            for charge in range(ell+1):
                sub=A[np.ix_(rr==nup-charge,cc==charge)]
                if min(sub.shape)>0: parts.append(np.linalg.svd(sub,compute_uv=False))
            s=np.concatenate(parts)
        q=s*s; q=q[q>1e-15]; result.append(float(-q@np.log(q)))
    return np.array(result)

def entropy_fit(S,L,min_cut=2):
    ell=np.arange(min_cut,L//2+1)
    x=np.log(L/np.pi*np.sin(np.pi*ell/L))/3
    y=np.asarray(S)[ell-1]
    if len(x)<2: return {'c':None}
    X=np.column_stack([x,np.ones(len(x))]); beta=np.linalg.lstsq(X,y,rcond=None)[0]
    residual=y-X@beta; rms=np.sqrt(np.mean(residual**2))
    se=float(np.sqrt(np.sum(residual**2)/max(1,len(x)-2)*np.linalg.inv(X.T@X)[0,0]))
    return dict(c=float(beta[0]),intercept=float(beta[1]),fit_se=se,rms=float(rms),min_cut=min_cut)

def quick(J,field,L,k=4):
    h=density(np.array(J),np.array(field)); H=full_operator(L,h)
    e,u=solve(H,k=k); psi,m,tdef,deg=resolve_ground(e,u,L)
    S=entropies(psi,L); fit=entropy_fit(S,L)
    residual=np.linalg.norm(H@psi-e[0]*psi)
    positive=e-e[0]; positive=positive[positive>1e-7]
    return dict(L=L,E0=float(e[0]),levels=e.tolist(),gap=float(e[1]-e[0]),
                gap_nondegenerate=float(positive[0]) if len(positive) else None,
                momentum=m,translation_defect=tdef,ground_degeneracy_returned=deg,
                residual=float(residual),entropy=S.tolist(),**fit)

def detailed(J,field,L,m=None,levels=10,verify_sectors=False,nup=None,field_split=.5):
    start=time.time(); h=density(np.array(J),np.array(field),field_split)
    all_momenta=(verify_sectors or m is None)
    if nup is not None:
        charge=np.diag([0.,1.,1.,2.])
        if np.linalg.norm(h@charge-charge@h)>1e-11:
            raise ValueError('Fixed-charge restriction requested for a non-U(1) density')
    energies=[]
    if verify_sectors or m is None:
        best=None
        for mm in range(L):
            H=sector_operator(L,h,mm,nup=nup); ee,uu=solve(H,k=1)
            energies.append(float(ee[0]))
            if best is None or ee[0]<best[0]-1e-9: best=(ee[0],mm,uu[:,0])
            if L>=20 and mm%4==0:
                print('ground-sector progress',L,mm,'best_m',best[1],flush=True)
        _,m,psi=best
    H=sector_operator(L,h,m,nup=nup); e,u=solve(H,k=4); psi=u[:,0]; E0=e[0]
    ground_residual=float(np.linalg.norm(H@psi-E0*psi)); del H
    full=expand(L,m,psi,nup); S=entropies(full,L,nup); fit=entropy_fit(S,L)
    del full
    Hp=sector_operator(L,h,(m+1)%L,nup=nup); ep,_=solve(Hp,k=2)
    del Hp
    Hm=sector_operator(L,h,(m-1)%L,nup=nup); em,_=solve(Hm,k=2)
    del Hm
    vplus=(ep[0]-E0)*L/(2*np.pi); vminus=(em[0]-E0)*L/(2*np.pi)
    stress=[]
    for n in (-2,2):
        target=(m-n)%L
        Ht=sector_operator(L,h,target,nup=nup)
        O=sector_operator(L,h,m,n=n,nup=nup); w=O@psi
        del O
        et,ut=solve(Ht,k=levels,v0=w/np.linalg.norm(w) if np.linalg.norm(w)>1e-12 else None)
        target_residual=float(np.max(np.linalg.norm(Ht@ut-ut*et,axis=0)))
        del Ht
        amps=ut.conj().T@w; weights=abs(amps)**2
        # Degenerate eigensolvers can return arbitrary nonorthogonal vectors.
        # Orthonormalize each cluster and sum its spectral weight. The normalized
        # projection of h_n|0> is the largest-overlap state in that eigenspace.
        clusters=[]; i=0
        while i<len(et):
            j=i+1
            while j<len(et) and abs(et[j]-et[i])<1e-7: j+=1
            Q=np.linalg.qr(ut[:,i:j])[0]
            weight=float(np.linalg.norm(Q.conj().T@w)**2)
            gap=float(np.mean(et[i:j])-E0)
            clusters.append(dict(energy=float(np.mean(et[i:j])),gap=gap,weight=weight,
                                 multiplicity=j-i,c_MV=8*weight/gap**2 if abs(gap)>1e-10 else None))
            i=j
        valid=[z for z in clusters if z['gap']>1e-7]
        selected=max(valid,key=lambda z:z['weight'])
        stress.append(dict(n=n,selected=selected,clusters=clusters,
                           target_eigen_residual=target_residual,
                           captured_weight=float(sum(z['weight'] for z in clusters)),
                           total_weight=float(np.vdot(w,w).real),
                           velocity=selected['gap']*L/(4*np.pi)))
    alt=entropy_fit(S,L,3)
    return dict(L=L,momentum=m,nup=nup,field_split=field_split,all_ground_momenta_checked=all_momenta,E0=float(E0),levels_same_momentum=e.tolist(),
                sector_ground_energies=energies,v_q1_plus=float(vplus),v_q1_minus=float(vminus),
                stress=stress,entropy=S.tolist(),entropy_cut3=alt,entropy_bulk=entropy_fit(S,L,max(2,L//4)),
                residual=ground_residual,seconds=time.time()-start,**fit)

CALIBRATIONS={
    'Ising':dict(J=np.diag([-1.,0,0]).tolist(),field=[0,0,-1.],expected=.5),
    'XX':dict(J=np.diag([1.,1.,0]).tolist(),field=[0,0,0],expected=1.),
    'Heisenberg':dict(J=np.eye(3).tolist(),field=[0,0,0],expected=1.),
    'XXZ_0.5':dict(J=np.diag([1.,1.,.5]).tolist(),field=[0,0,0],expected=1.),
}

def dump(path,obj):
    Path(path).parent.mkdir(parents=True,exist_ok=True)
    tmp=Path(str(path)+'.tmp'); tmp.write_text(json.dumps(obj,indent=2)+'\n'); tmp.replace(path)

def energy_fit(rows):
    L=np.array([r['L'] for r in rows],float); E=np.array([r['E0'] for r in rows])
    # E/L = e_infty + a/L^2 + b/L^4; a = -pi sum(c_a v_a)/6.
    X=np.column_stack([np.ones(len(L)),L**-2,L**-4])
    coef=np.linalg.lstsq(X,E/L,rcond=None)[0]
    v=np.mean([z['velocity'] for z in rows[-1]['stress']])
    return dict(e_infty=float(coef[0]),casimir_coefficient=float(coef[1]),
                sum_c_v=float(-6*coef[1]/np.pi),velocity_used=float(v),
                c_energy=float(-6*coef[1]/(np.pi*v)),sizes=L.astype(int).tolist())

def calibrate(sizes,output):
    result={}
    if Path(output).exists(): result=json.loads(Path(output).read_text())
    for name,par in CALIBRATIONS.items():
        rows=result.setdefault(name,dict(parameters=par,rows=[]))['rows']
        for L in sizes:
            if any(r['L']==L for r in rows): continue
            m=0 if name=='Ising' or L%4==0 else L//2
            row=detailed(par['J'],par['field'],L,m,verify_sectors=L==8,nup=None if name=='Ising' else L//2)
            rows.append(row); rows.sort(key=lambda r:r['L'])
            print(name,L,'cS',round(row['c'],6),'cMV',round(row['stress'][0]['selected']['c_MV'],6),
                  'v1',round(row['v_q1_plus'],6),'seconds',round(row['seconds'],2),flush=True)
            if len(rows)>=3: result[name]['energy_fit']=energy_fit(rows)
            dump(output,result)

def normalized(J,field):
    J=np.array(J,float); field=np.array(field,float)
    scale=np.sqrt(np.sum(J*J)+np.sum(field*field))
    return (J/scale).tolist(),(field/scale).tolist()

def parameters(seed=20261021,random_points=80):
    points=[]
    def add(family,J,field,label):
        J,field=normalized(J,field)
        points.append(dict(id=len(points),family=family,label=label,J=J,field=field))
    # XY / XYZ longitudinal-to-z field. Include Lifshitz and saturation boundaries.
    for delta in [-.9,-.5,0,.5,1,1.5]:
        for hz in [0,.5,1.,1.5,2.,2.5]:
            add('XXZ_field',np.diag([1,1,delta]),[0,0,hz],f'delta={delta},hz={hz}')
    for gamma in [.2,.6,1.]:
        for delta in [0,.4,1.]:
            for hz in [0,1,2]:
                add('XYZ_field',np.diag([1+gamma,1-gamma,delta]),[.13,.17,hz],f'g={gamma},d={delta},hz={hz}')
    for D in [.3,.7,1.2,2.]:
        for gamma in [0,.4,1.]:
            for hz in [0,.8,1.6,2.4]:
                J=np.diag([1+gamma,1-gamma,.3]); J[0,1]=D; J[1,0]=-D
                add('XYZ_DM',J,[0,0,hz],f'D={D},g={gamma},hz={hz}')
    for K in [-2.,-1.,0.,1.,2.]:
        for G in [-1.,0.,1.]:
            J=np.diag([1.,1.,1+K]); J[0,1]=G; J[1,0]=G
            add('uniform_KHG',J,[.3,.2,.1],f'K={K},Gamma={G}')
    for hx in [0,.03,.15,.5]:
        for hz in [.7,1.,1.3]:
            add('Ising_fields',np.diag([-1.,0,0]),[hx,0,hz],f'hx={hx},hz={hz}')
    rng=np.random.default_rng(seed)
    for i in range(random_points):
        # Any real J has symmetric part diagonalizable by SO(3).
        # This canonical chart has 3 diagonal + 3 DM + 3 field coordinates,
        # minus scale, so eight continuous parameters.
        diag=rng.normal(size=3); D=rng.normal(size=3)*rng.choice([.25,1.,2.])
        J=np.diag(diag)+np.array([[0,D[2],-D[1]],[-D[2],0,D[0]],[D[1],-D[0],0]])
        field=rng.normal(size=3)*rng.choice([.1,.4,1.])
        add('random_canonical',J,field,f'random={i}')
    return points

def score(rows):
    a,b=rows[-2:]; g1=a['gap_nondegenerate']; g2=b['gap_nondegenerate']
    z=-np.log(g2/g1)/np.log(b['L']/a['L']) if g1 and g2 else 10.
    c=b['c']; drift=abs(c-a['c'])
    penalty=.12*min(abs(z-1),5)+.6*drift+10*b['translation_defect']
    return dict(c_eff=c,c_drift=drift,z_gap=float(z),score=float(c-penalty))

def scan(points,sizes,output):
    done=[]
    if Path(output).exists(): done=json.loads(Path(output).read_text())
    ids={p['id'] for p in done}
    for p in points:
        if p['id'] in ids: continue
        t=time.time()
        try:
            rows=[quick(p['J'],p['field'],L) for L in sizes]
            r=dict(**p,rows=rows,**score(rows),seconds=time.time()-t)
        except Exception as exc:
            r=dict(**p,error=repr(exc),seconds=time.time()-t)
        done.append(r); dump(output,done)
        print(p['id'],p['family'],{k:round(r[k],5) for k in ('c_eff','z_gap','score') if k in r},round(r['seconds'],2),flush=True)
    return done

def optimize(scanfile,output,maxiter=24,seeds=3,sizes=(10,12),seed_ids=None,prefix='opt',velocity_penalty=False):
    data=json.loads(Path(scanfile).read_text()); candidates=sorted([r for r in data if 'score' in r],key=lambda r:r['score'],reverse=True)
    trace=[]; used=[]
    for p in candidates:
        if seed_ids is not None:
            if str(p['id']) in seed_ids: used.append(p)
        elif p['family'] not in [x['family'] for x in used]: used.append(p)
        if len(used)>=seeds: break
    for seed in used:
        J=np.array(seed['J']); f=np.array(seed['field']); S=(J+J.T)/2
        _,R=np.linalg.eigh(S); R[:,0]*=np.linalg.det(R)
        J=R.T@J@R; f=R.T@f
        x0=np.r_[np.diag(J),J[1,2],J[2,0],J[0,1],f]
        def fun(x):
            A=np.diag(x[:3])+np.array([[0,x[5],-x[4]],[-x[5],0,x[3]],[x[4],-x[3],0]])
            A,f=normalized(A,x[6:]); rows=[quick(A,f,L) for L in sizes]
            r=dict(id=prefix+str(len(trace)),family='local_optimization',seed=seed['id'],J=A,field=f,rows=rows,**score(rows))
            if velocity_penalty:
                h=density(np.array(A),np.array(f)); velocities=[]
                for row in rows:
                    L=row['L']; vv=[]
                    for sign in (-1,1):
                        B=sector_operator(L,h,(row['momentum']+sign)%L)
                        ee,_=solve(B,k=1); vv.append(float((ee[0]-row['E0'])*L/(2*np.pi)))
                    velocities.append(vv)
                v1,v2=np.mean(velocities,axis=1)
                zq=1-np.log(v2/v1)/np.log(sizes[-1]/sizes[-2]) if min(v1,v2)>1e-9 else 10.
                r['q1_velocities']=velocities; r['z_q1']=float(zq)
                r['score']-=.4*min(abs(zq-1),5)
            trace.append(r); dump(output,trace)
            print(r['id'],r['seed'],round(r['c_eff'],6),round(r['score'],6),flush=True)
            return -r['score']
        minimize(fun,x0,method='Nelder-Mead',options=dict(maxiter=maxiter,xatol=2e-3,fatol=2e-3))

def selftest():
    rng=np.random.default_rng(5); L=6; J=rng.normal(size=(3,3)); field=rng.normal(size=3); h=density(J,field)
    H=full_operator(L,h); full_e=np.linalg.eigvalsh(H.toarray()); all_e=[]
    errors=[]; fourier_errors=[]
    for m in range(L):
        A=sector_operator(L,h,m); all_e.extend(np.linalg.eigvalsh(A.toarray()))
        u=rng.normal(size=A.shape[0])+1j*rng.normal(size=A.shape[0]); u/=np.linalg.norm(u)
        psi=expand(L,m,u)
        errors.append(np.linalg.norm(H@psi-expand(L,m,A@u)))
        n=-2; w=np.zeros_like(psi)
        # Direct h_n in the full computational basis, independent of orbit formula.
        for s in range(1<<L):
            for j in range(L):
                jp=(j+1)%L; bj=(s>>j)&1; bp=(s>>jp)&1
                for b in range(4):
                    sp=s^((bj^(b>>1))<<j)^((bp^(b&1))<<jp)
                    w[sp]+=np.exp(2j*np.pi*n*j/L)*h[b,2*bj+bp]*psi[s]
        O=sector_operator(L,h,m,n)
        fourier_errors.append(np.linalg.norm(w-expand(L,(m-n)%L,O@u)))
    result=dict(spectrum_error=float(np.max(abs(np.sort(all_e)-full_e))),
                projection_error=float(max(errors)),fourier_error=float(max(fourier_errors)),
                hermiticity_error=float(sparse.linalg.norm(H-H.getH())))
    assert max(result.values())<1e-10,result
    dump('results/selftest.json',result); print(result)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    sub=ap.add_subparsers(dest='command',required=True)
    sub.add_parser('selftest')
    p=sub.add_parser('calibrate'); p.add_argument('--sizes',nargs='+',type=int,default=[8,10,12,14,16,18,20]); p.add_argument('--output',default='results/calibration.json')
    p=sub.add_parser('points'); p.add_argument('--random',type=int,default=80); p.add_argument('--output',default='results/points.json')
    p=sub.add_parser('scan'); p.add_argument('--points',default='results/points.json'); p.add_argument('--sizes',nargs='+',type=int,default=[10,12]); p.add_argument('--output',default='results/scan.json'); p.add_argument('--index',type=int)
    p=sub.add_parser('optimize'); p.add_argument('--input',default='results/scan.json'); p.add_argument('--output',default='results/optimization.json'); p.add_argument('--maxiter',type=int,default=24); p.add_argument('--seeds',type=int,default=3); p.add_argument('--sizes',nargs=2,type=int,default=[10,12]); p.add_argument('--seed-ids',nargs='+'); p.add_argument('--prefix',default='opt'); p.add_argument('--velocity-penalty',action='store_true')
    p=sub.add_parser('refine'); p.add_argument('--input',default='results/scan.json'); p.add_argument('--ids',nargs='+',required=True); p.add_argument('--sizes',nargs='+',type=int,default=[14,16,18,20]); p.add_argument('--output',default='results/refined.json'); p.add_argument('--levels',type=int,default=12); p.add_argument('--momentum',type=int); p.add_argument('--af-momentum',action='store_true'); p.add_argument('--half-filling',action='store_true'); p.add_argument('--field-split',type=float,default=.5)
    args=ap.parse_args()
    # A PID record makes long runs interruptible without enumerating processes.
    if hasattr(args,'output'):
        Path(args.output).parent.mkdir(parents=True,exist_ok=True)
        Path(args.output+'.pid').write_text(str(os.getpid())+'\n')
    if args.command=='selftest': selftest()
    elif args.command=='calibrate': calibrate(args.sizes,args.output)
    elif args.command=='points': dump(args.output,parameters(random_points=args.random))
    elif args.command=='scan':
        points=json.loads(Path(args.points).read_text())
        if args.index is not None: points=[points[args.index]]
        scan(points,args.sizes,args.output)
    elif args.command=='optimize': optimize(args.input,args.output,args.maxiter,args.seeds,args.sizes,args.seed_ids,args.prefix,args.velocity_penalty)
    elif args.command=='refine':
        data=json.loads(Path(args.input).read_text()); out={}
        if Path(args.output).exists(): out=json.loads(Path(args.output).read_text())
        for p in data:
            if str(p['id']) not in args.ids: continue
            entry=out.setdefault(str(p['id']),dict(parameters={k:p[k] for k in ('id','family','J','field')},rows=[]))
            for L in args.sizes:
                if any(r['L']==L for r in entry['rows']): continue
                m=(0 if L%4==0 else L//2) if args.af_momentum else args.momentum
                row=detailed(p['J'],p['field'],L,m,levels=args.levels,nup=L//2 if args.half_filling else None,field_split=args.field_split)
                entry['rows'].append(row); entry['rows'].sort(key=lambda r:r['L'])
                if len(entry['rows'])>=3: entry['energy_fit']=energy_fit(entry['rows'])
                dump(args.output,out)
                print(p['id'],L,'m',row['momentum'],'cS',row['c'],'cMV',row['stress'][0]['selected']['c_MV'],'secs',row['seconds'],flush=True)

if __name__=='__main__': main()
