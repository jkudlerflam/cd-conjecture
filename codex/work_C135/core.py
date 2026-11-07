"""C135 small computations. Numerical certificates are floating point, not intervals."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
import json, time
from pathlib import Path
from functools import lru_cache
import numpy as np
import scipy.linalg as la
import scipy.sparse as sp
import scipy.sparse.linalg as sla
import cvxpy as cp

I=np.eye(2); P=np.array([I,[[0,1],[1,0]],[[0,-1j],[1j,0]],[[1,0],[0,-1]]],complex)
S=P[1:]/2
B=np.array([np.kron(s,s) for s in S]+[np.kron(S[a],S[b])-np.kron(S[b],S[a]) for a,b in [(1,2),(2,0),(0,1)]]+[(np.kron(s,I)+np.kron(I,s))/2 for s in S])
G=np.array([np.kron(s,I)-np.kron(I,s) for s in S])
def xxz(d,b=0,dm=0): return B[0]+B[1]+d*B[2]+dm*B[5]+b*B[8]
def tfim(g): return -4*B[0]-2*g*B[8]
def gauge(h):
    o=cp.Variable(3); a=cp.Variable(); m=cp.Variable()
    q=h+a*np.eye(4)+sum(o[j]*G[j] for j in range(3))
    p=cp.Problem(cp.Minimize(m),[m*np.eye(4)-q>>0,m*np.eye(4)+q>>0])
    try:
        p.solve(solver='CLARABEL',tol_gap_abs=1e-10,tol_feas=1e-10,tol_gap_rel=1e-10)
    except cp.error.SolverError:
        p.solve(solver='SCS',eps=1e-9,max_iters=20000)
    if p.status!='optimal':
        p.solve(solver='SCS',eps=1e-9,max_iters=20000)
    rep=h+sum(o.value[j]*G[j] for j in range(3))
    ev=la.eigvalsh(rep); upper=(ev[-1]-ev[0])/2
    dual=herm(p.constraints[0].dual_value-p.constraints[1].dual_value)
    dual-=np.trace(dual)*np.eye(4)/4
    for g in G:dual-=np.trace(dual@g).real/np.trace(g@g).real*g
    dual/=max(1.,float(abs(la.eigvalsh(dual)).sum()))
    lower=float(np.trace(h@dual).real)
    return dict(value=float(m.value),lower=lower,upper=float(upper),gap=float(upper-lower),status=p.status,o=o.value.tolist())
def herm(a): return (a+a.conj().T)/2
def left(r):
    d=r.shape[0]//2
    return r[:d,:d]+r[d:,d:]
def right(r): return r[::2,::2]+r[1::2,1::2]
def marginal(r,n,k=2):
    for j in range(n-k): r=right(r)
    return r

@lru_cache(None)
def pauli_orbits(n):
    groups={}
    for idx in range(4**n):
        digs=[]; v=idx
        for _ in range(n): digs.append(v%4); v//=4
        digs=digs[::-1]; non=[j for j,x in enumerate(digs) if x]
        key=tuple(digs[non[0]:non[-1]+1]) if non else ()
        groups.setdefault(key,[]).append(idx)
    return tuple(np.array(v) for v in groups.values())
def pauli_transform(r,n):
    t=r.reshape([2]*(2*n)).transpose(sum(([j,n+j] for j in range(n)),[])).reshape([4]*n)
    # Basis rows are P[a].T flattened, so this gives Tr(P r)/2^n.
    u=P.transpose(0,2,1).reshape(4,4)/2
    for axis in range(n): t=np.moveaxis(np.tensordot(u,t,axes=(1,axis)),0,axis)
    return t
def pauli_inverse(t,n):
    u=P.reshape(4,4).T
    for axis in range(n): t=np.moveaxis(np.tensordot(u,t,axes=(1,axis)),0,axis)
    return t.reshape([2]*(2*n)).transpose(list(range(0,2*n,2))+list(range(1,2*n,2))).reshape(2**n,2**n)
def repair(r,n):
    c=pauli_transform(herm(r),n).real.flatten()
    for inds in pauli_orbits(n): c[inds]=np.mean(c[inds])
    c[0]=1/2**n
    q=herm(pauli_inverse(c.reshape([4]*n),n))
    neg=max(0.,-float(la.eigvalsh(q,subset_by_index=[0,0])[0]))+2e-13
    q=(q+neg*np.eye(2**n))/(1+neg*2**n)
    return q,float(neg*2**n),float(la.norm(left(q)-right(q)))

class LTI:
    def __init__(self,n,real=True,sym='none',distance=False):
        self.n=n; D=2**n; self.real=real
        if sym=='u1': keys=np.array([i.bit_count() for i in range(D)])
        elif sym=='parity': keys=np.array([i.bit_count()%2 for i in range(D)])
        else: keys=np.zeros(D,int)
        terms=[]; self.blocks=[]; cons=[]
        for key in sorted(set(keys)):
            inds=np.flatnonzero(keys==key); k=len(inds)
            v=cp.Variable((k,k),symmetric=True) if real else cp.Variable((k,k),hermitian=True)
            emb=sp.csc_matrix((np.ones(k),(inds,np.arange(k))),shape=(D,k))
            terms.append(emb@v@emb.T); cons.append(v>>0); self.blocks.append((inds,v))
        r=sum(terms); self.r=r
        self.norm=cp.trace(r)==1; cons.append(self.norm)
        # Tr_n - Tr_1, so its multiplier is X in the requested dual sign convention.
        self.eq=r[::2,::2]+r[1::2,1::2]-r[:D//2,:D//2]-r[D//2:,D//2:]
        self.ec=self.eq==0; cons.append(self.ec)
        self.hp=cp.Parameter((4,4),symmetric=True) if real else cp.Parameter((4,4),hermitian=True)
        obj=cp.trace(cp.kron(self.hp,np.eye(D//4))@r)
        self.problem=cp.Problem(cp.Minimize(obj if real else cp.real(obj)),cons)
        self.constraints=cons
    def solve(self,h,solver=None,eps=2e-7,iters=15000):
        start=time.time(); self.hp.value=h.real if self.real else h
        solver=solver or ('CLARABEL' if self.n<=5 else 'SCS')
        opts=dict(tol_gap_abs=1e-9,tol_gap_rel=1e-9,tol_feas=1e-9,max_iter=150) if solver=='CLARABEL' else dict(eps=eps,max_iters=iters,acceleration_lookback=10)
        try:
            self.problem.solve(solver=solver,**opts)
        except cp.error.SolverError:
            if solver!='CLARABEL':raise
            solver='SCS';self.problem.solve(solver=solver,eps=eps,max_iters=iters)
        r=herm(np.array(self.r.value)); x=herm(np.array(self.ec.dual_value)); D=len(r)
        H=np.kron(h,np.eye(D//4)); K=herm(H+np.kron(x,I)-np.kron(I,x))
        # This eigenvalue is a feasible dual objective without trusting solver feasibility.
        lower=float(la.eigvalsh(K,subset_by_index=[0,0])[0])-2e-12*max(1,la.norm(K,2))
        q,mix,eqfix=repair(r,self.n)
        upper=float(np.trace(H@q).real)
        stat=self.problem.solver_stats
        rec=dict(n=self.n,value=float(self.problem.value),lower=lower,upper=upper,status=self.problem.status,solver=solver,seconds=time.time()-start,iterations=stat.num_iters,trace_residual=float(abs(np.trace(r)-1)),lti_residual=float(la.norm(left(r)-right(r))),min_rho=float(la.eigvalsh(r,subset_by_index=[0,0])[0]),repair_mix=mix,repaired_lti_residual=eqfix,dual_gap=upper-lower,slackness=float(la.norm((K-self.problem.value*np.eye(D))@r)))
        if solver=='SCS':
            info=stat.extra_stats.get('info',{})
            rec['scs']={k:info[k] for k in ['res_pri','res_dual','gap','status'] if k in info}
        return rec,q,x

def ring_matrix(h,L):
    D=2**L; states=np.arange(D,dtype=np.int64); rows=[];cols=[];vals=[]
    for j in range(L):
        k=(j+1)%L; a=L-1-j; b=L-1-k
        old=2*((states>>a)&1)+((states>>b)&1); base=states&~((1<<a)|(1<<b))
        for out in range(4):
            v=h[out,old]; use=np.flatnonzero(abs(v)>1e-14)
            rows.append(base[use]|((out>>1)<<a)|((out&1)<<b)); cols.append(states[use]); vals.append(v[use])
    H=sp.coo_matrix((np.concatenate(vals),(np.concatenate(rows),np.concatenate(cols))),shape=(D,D)).tocsr()
    return H.real if np.max(abs(H.data.imag),initial=0)<1e-14 else H
def ring(h,L,kmax=6):
    H=ring_matrix(h,L); D=2**L
    if L<=7: ev,v=la.eigh(H.toarray())
    else:
        rng=np.random.default_rng(135+L); v0=rng.normal(size=D)
        if np.iscomplexobj(H.data): v0=v0+1j*rng.normal(size=D)
        ev,v=sla.eigsh(H,k=min(kmax,D-2),which='SA',tol=2e-11,v0=v0)
        ix=np.argsort(ev);ev=ev[ix];v=v[:,ix]
    sel=ev<ev[0]+2e-8; vs=v[:,sel]
    residual=float(max(la.norm(H@vs[:,j]-ev[0]*vs[:,j]) for j in range(vs.shape[1])))
    return dict(L=L,E=float(ev[0]),e=float(ev[0]/L),residual=residual,multiplicity_sampled=int(sel.sum()),groundspace_complete=bool(L<=7 or not sel[-1])),vs
def ti_marginal(v,L,l):
    ans=np.zeros((2**l,2**l),complex)
    for j in range(v.shape[1]):
        t=v[:,j].reshape([2]*L)
        for shift in range(L):
            q=t.transpose(list(range(shift,L))+list(range(shift))).reshape(2**l,-1)
            ans+=q@q.conj().T/(L*v.shape[1])
    return herm(ans)
def power(r,p):
    w,v=la.eigh(herm(r)); f=np.zeros_like(w); use=w>max(w.max()*1e-11,1e-13);f[use]=w[use]**p
    return (v*f)@v.conj().T
def petz(r):
    n=round(np.log2(len(r))); b=right(r); s=power(r,.5); ib=power(b,-.5)
    z=np.kron(I,ib)@r@np.kron(I,ib)
    z=np.kron(I,s)@np.kron(z,I)@np.kron(I,s)
    out=marginal(z,n+1)
    return herm(out),dict(trace_residual=float(abs(np.trace(out)-1)),lti_input_residual=float(la.norm(left(r)-right(r))),min_output=float(la.eigvalsh(out)[0]))
def write_json(path,x): Path(path).write_text(json.dumps(x,indent=2)+'\n')
def pack(a): return dict(real=a.real.tolist(),imag=a.imag.tolist())
def unpack(a): return np.array(a['real'])+1j*np.array(a['imag'])
