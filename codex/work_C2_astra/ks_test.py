#!/usr/bin/env python3
"""Closed-book Koo-Saleur test. numpy/scipy only; periodic spin chains.
Run: OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 ks_test.py
Outputs results.json and progress.txt. To reproduce the full reported sweep, use
  python3 ks_test.py --maxN 22 --output results_extended.json
The default maximum is 20. Ising/OF densities are site-centered; the literal
Ising density is included as an additional control. Spin-chain densities are
the displayed forward-bond summands. Inverse-current moments include all
non-ground states and are obtained by a projected positive resolvent. Momentum eigenvalue is exp(2pi*i*k/N).
All symmetry sectors are searched (both Ising parities; Sz=0 for SU(2)).
SU(2) guarantees every energy multiplet of even N has an Sz=0 member.
"""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[key]='1'
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh, cg, LinearOperator
import time, json, argparse
START=time.time()

def parity(a):
    a=np.asarray(a,dtype=np.uint32).copy()
    a ^= a>>16; a ^= a>>8; a ^= a>>4
    return (np.uint32(0x6996)>>(a&15))&1

def add(d,key,c):
    d[key]=d.get(key,0)+c

def translate(d,j,N):
    mask=(1<<N)-1; j%=N
    def rot(x): return ((x<<j)|(x>>(N-j)))&mask if j else x
    return {(rot(x),rot(z)):c for (x,z),c in d.items()}

def product(a,b):
    d={}
    for (x,z),c in a.items():
        for (u,v),e in b.items():
            add(d,(x^u,z^v),c*e*(-1 if (z&u).bit_count()%2 else 1))
    return d

def clean(d): return {k:v for k,v in d.items() if abs(v)>1e-13}

def density(N,model,coupling=0.,literal=False):
    d={}
    def p(s,c):
        x=z=0
        for j,op in s:
            bit=1<<(j%N)
            if op in 'XY': x ^=bit
            if op in 'ZY': z ^=bit
            if op=='Y': c*=1j
        add(d,(x,z),c)
    if model=='ising':
        p([(0,'Z')],-1.)
        if literal: p([(0,'X'),(1,'X')],-1.)
        else:
            p([(-1,'X'),(0,'X')],-.5); p([(0,'X'),(1,'X')],-.5)
        if coupling:
            p([(-1,'X'),(0,'X'),(1,'Z')],coupling)
            p([(-1,'Z'),(0,'X'),(1,'X')],coupling)
    else:
        for r,c in [(1,1.),(2,coupling)]:
            for op in 'XYZ': p([(0,op),(r,op)],c*.25)
    return clean(d)

def total(d,N,mode=0):
    out={}
    for j in range(N):
        f=np.exp(-2j*np.pi*mode*j/N) if mode else 1.
        for key,c in translate(d,j,N).items(): add(out,key,c*f)
    return clean(out)

def current(d,N):
    seed={}
    for r in (1,2):
        dr=translate(d,r,N)
        for key,c in product(d,dr).items(): add(seed,key,1j*r*c)
        for key,c in product(dr,d).items(): add(seed,key,-1j*r*c)
    return total(clean(seed),N)

class Orbits:
    def __init__(self,N,sector):
        self.N=N; mask=(1<<N)-1
        s=np.arange(1<<N,dtype=np.uint32)
        if sector in ('even','odd'): s=s[parity(s)==(sector=='odd')]
        else:
            counts=np.zeros(len(s),np.uint8)
            for j in range(N): counts+=(s>>j)&1
            s=s[counts==N//2]
        r=s.copy(); cur=s.copy(); sh=np.zeros(len(s),np.int16)
        for t in range(1,N):
            cur=((cur<<1)|(cur>>(N-1)))&mask
            take=cur<r; r[take]=cur[take]; sh[take]=N-t
        reps,inv,L=np.unique(r,return_inverse=True,return_counts=True)
        self.reps=reps; self.L=L
        self.oid=np.full(1<<N,-1,np.int32); self.oid[s]=inv
        self.shift=np.zeros(1<<N,np.int16); self.shift[s]=sh
        self.bases={}
    def basis(self,k):
        k%=self.N
        if k not in self.bases:
            ids=np.flatnonzero((k*self.L)%self.N==0)
            idx=np.full(len(self.reps),-1,np.int32); idx[ids]=np.arange(len(ids))
            self.bases[k]=(ids,idx)
        return self.bases[k]
    def translation_error(self,k,vec):
        ids,imap=self.basis(k)
        states=np.flatnonzero(self.oid>=0).astype(np.uint32)
        oi=self.oid[states]; keep=imap[oi]>=0; states=states[keep]; oi=oi[keep]
        amplitudes=vec[imap[oi]]*np.exp(-2j*np.pi*k*self.shift[states]/self.N)/np.sqrt(self.L[oi])
        prev=(states>>1)|((states&1)<<(self.N-1))
        po=self.oid[prev]
        shifted=vec[imap[po]]*np.exp(-2j*np.pi*k*self.shift[prev]/self.N)/np.sqrt(self.L[po])
        return float(np.linalg.norm(shifted-np.exp(2j*np.pi*k/self.N)*amplitudes))
    def matrix(self,terms,ks,kd=None):
        if kd is None: kd=ks
        si,_=self.basis(ks); di,dmap=self.basis(kd)
        s=self.reps[si]; cols0=np.arange(len(si),dtype=np.int32)
        rows=[]; cols=[]; vals=[]
        for (x,z),c in terms.items():
            dest=s^x; oi=self.oid[dest]
            valid=oi>=0
            row=np.full(len(s),-1,np.int32); row[valid]=dmap[oi[valid]]
            valid &= row>=0
            a=c*(1.-2.*parity(s&z))*np.exp(2j*np.pi*kd*self.shift[dest]/self.N)
            a[valid]*=np.sqrt(self.L[si[valid]]/self.L[oi[valid]])
            rows.append(row[valid]); cols.append(cols0[valid]); vals.append(a[valid])
        M=coo_matrix((np.concatenate(vals),(np.concatenate(rows),np.concatenate(cols))),shape=(len(di),len(si))).tocsr()
        M.sum_duplicates(); M.eliminate_zeros()
        if np.max(np.abs(M.data.imag),initial=0)<1e-12: M=M.real
        return M

def lowest(H,n=1,v0=None):
    if H.shape[0]<=max(n+1,15):
        e,v=np.linalg.eigh(H.toarray()); return e[:n],v[:,:n]
    e,v=eigsh(H,k=n,which='SA',tol=2e-10,maxiter=5000,v0=v0,ncv=max(25,2*n+5))
    order=np.argsort(e); return e[order],v[:,order]

def run(N,model,coupling,label,check_all=True):
    tick=time.time(); h=density(N,model,coupling); ht=total(h,N)
    candidates=[]; blocks={}; bases={}
    # k and -k have equal spectra because H is real.
    for sector in (['even','odd'] if model=='ising' else ['sz0']):
        orb=Orbits(N,sector); bases[sector]=orb
        for k in (range(N//2+1) if check_all else sorted(set([0,N//2]))):
            H=orb.matrix(ht,k); e,v=lowest(H,2 if k in (0,N//2) else 1)
            candidates.append((float(e[0]),sector,k))
            if k in (0,N//2): blocks[sector,k]=(H,e,v)
    candidates.sort(); E,sector,k=candidates[0]
    orb=bases[sector]
    if (sector,k) in blocks: H,e,v=blocks[sector,k]
    else: H=orb.matrix(ht,k); e,v=lowest(H,2)
    gs=v[:,0]; E=float(e[0]); gapcheck=min(float(e[1])-E,candidates[1][0]-E)
    if k not in (0,N//2): gapcheck=0. # distinct conjugate ground momentum
    if gapcheck<1e-8:
        return dict(N=N,label=label,degenerate=True,gapcheck=gapcheck)
    kd=(k+2)%N; Ht=orb.matrix(ht,kd); et,vt=lowest(Ht,10)
    out=dict(N=N,label=label,E0=E,k=k,sector=sector,gapcheck=gapcheck,unique_all_sectors=check_all,
             ground_residual=float(np.linalg.norm(H@gs-E*gs)),translation_residual=orb.translation_error(k,gs))
    for kind,d in [('centered',h)]+([('literal',density(N,model,coupling,literal=True))] if model=='ising' and coupling==0 else []):
        hv=orb.matrix(total(d,N,2),k,kd)@gs
        m0=float(np.vdot(hv,hv).real); m1=float(np.vdot(hv,Ht@hv-E*hv).real)
        ovs=np.abs(vt.conj().T@hv)**2; it=int(np.argmax(ovs)); gap=float(et[it]-E)
        J=orb.matrix(current(d,N),k); jv=J@gs; jmean=np.vdot(gs,jv)
        var=float((np.vdot(jv,jv)-abs(jmean)**2).real/N)
        spectral={}
        if var>1e-15:
            rhs=jv-gs*jmean
            op=LinearOperator(H.shape,matvec=lambda a:H@a-E*a+gs*np.vdot(gs,a),dtype=complex)
            sol,info=cg(op,rhs,rtol=1e-10,atol=1e-13,maxiter=2000)
            A2=float(np.vdot(sol,sol).real/N); A1=float(np.vdot(rhs,sol).real/N)
            spectral=dict(Aminus2=A2,Aminus1=A1,c2_asymptotic_prefactor=128*np.pi**2*A2**3/A1**2,
                          resolvent_info=int(info),resolvent_residual=float(np.linalg.norm(op@sol-rhs)))
        out[kind]=dict(m0=m0,m1=m1,c2=8*m0**3/m1**2,m0N2=m0*N*N,m1N3=m1*N**3,
                       JvarN=var,cMV=8*float(ovs[it])/gap**2,Tgap=gap,Tindex=it,low_gaps=(et-E).tolist(),low_weights=ovs.tolist(),
                       low_weight_fraction=float(ovs.sum()/m0),Jmean=float(jmean.real),**spectral)
    out['seconds']=time.time()-tick
    return out

def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--maxN',type=int,default=20); parser.add_argument('--quick',action='store_true'); parser.add_argument('--minN',type=int,default=8); parser.add_argument('--output',default='results.json')
    args=parser.parse_args(); results=[]
    for N in range(args.minN,args.maxN+1,2):
        for model,c,label in [('ising',0.,'TFIM'),('spin',0.,'Heisenberg'),('ising',.2,'OF0.2'),('ising',.35,'OF0.35'),('spin',.2,'J1J2')]:
            out=run(N,model,c,label,not args.quick); results.append(out)
            with open(args.output,'w') as f: json.dump(results,f,indent=2)
            line=f'{label:10} N={N:2} '+json.dumps(out)
            print(line,flush=True)
            with open('progress.txt','a') as f: f.write(line+'\n')
    print('TOTAL SECONDS',time.time()-START,flush=True)
if __name__=='__main__': main()
