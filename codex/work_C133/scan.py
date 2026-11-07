"""C133 small full-ring ED. numpy/scipy only. All output stays beside this file."""
import argparse, json, time
from pathlib import Path
import numpy as np
from scipy.sparse import coo_matrix, csr_matrix
from scipy.sparse.linalg import eigsh, LinearOperator

I=np.eye(2); X=np.array([[0.,1.],[1.,0.]]); Y=np.array([[0.,-1j],[1j,0.]]); Z=np.diag([1.,-1.])

def spin(d):
    s=(d-1)/2; m=np.arange(s,-s-1,-1)
    plus=np.diag(np.sqrt(s*(s+1)-m[1:]*(m[1:]+1)),1)
    return [(plus+plus.T)/2,(plus-plus.T)/(2j),np.diag(m)]

def bonds():
    out=[]
    for delta in [-.9,0.,.5,.9,1.,1.05]:
        h=(np.kron(X,X)+np.kron(Y,Y)+delta*np.kron(Z,Z))/4
        out.append(('xxz',delta,2,h,[4,6,8,10,12,16]))
    xa=np.kron(X,I); xb=np.kron(I,X); za=np.kron(Z,I); zb=np.kron(I,Z); i4=np.eye(4)
    for ep in [-.5,0.,.5,1.]:
        field=xa+xb+ep*xa@xb
        h=-np.kron(za,za)-np.kron(zb,zb)-ep*np.kron(za@zb,za@zb)-.5*(np.kron(field,i4)+np.kron(i4,field))
        out.append(('ashkin_teller',ep,4,h,[2,3,4,5,6,8]))
    ss=spin(2)
    aa=sum(np.kron(np.kron(s,I),np.kron(s,I)) for s in ss)
    bb=sum(np.kron(np.kron(I,s),np.kron(I,s)) for s in ss)
    for g in [0.,1.,2.,4.]:
        out.append(('spin_orbital',g,4,aa+bb+g*aa@bb,[4,8]))
    s3=spin(3); sdot=sum(np.kron(s,s) for s in s3)
    for theta in [-np.pi/4,np.pi/4]:
        ls=[4,6,8,10] if theta<0 else [3,6,9]
        out.append(('spin1_blbq',float(theta),3,np.cos(theta)*sdot+np.sin(theta)*sdot@sdot,ls))
    for q in range(2,7):
        clock=np.diag(np.exp(2j*np.pi*np.arange(q)/q)); shift=np.roll(np.eye(q),1,axis=0)
        h=np.zeros((q*q,q*q),complex)
        for k in range(1,q):
            zk=np.linalg.matrix_power(clock,k); xk=np.linalg.matrix_power(shift,k)
            h-=(np.kron(zk.conj().T,zk)+.5*(np.kron(xk,np.eye(q))+np.kron(np.eye(q),xk)))/np.sin(np.pi*k/q)
        ls={2:[4,6,8,12,16],3:[3,4,6,8,10],4:[2,3,4,6,8],5:[2,3,4,6],6:[2,3,4,6]}[q]
        out.append(('fz_parafermion',q,q,h,ls))
    return out

def ring(h,d,L):
    h=np.real_if_close(h,tol=1000)
    ids=np.arange(d**L,dtype=np.int64); powers=d**np.arange(L-1,-1,-1,dtype=np.int64)
    digits=[(ids//p)%d for p in powers]
    diag=np.zeros(len(ids),dtype=h.dtype); rows=[]; cols=[]; vals=[]
    transitions=[[(r,h[r,c]) for r in np.flatnonzero(np.abs(h[:,c])>1e-13)] for c in range(d*d)]
    for j in range(L):
        k=(j+1)%L; pair=digits[j]*d+digits[k]
        diag+=h.diagonal()[pair]
        for c, ts in enumerate(transitions):
            inds=ids[pair==c]
            a,b=divmod(c,d)
            for r,v in ts:
                if r==c: continue
                ap,bp=divmod(r,d)
                rows.append(inds+(ap-a)*powers[j]+(bp-b)*powers[k]);cols.append(inds); vals.append(np.full(len(inds),v,dtype=h.dtype))
    rows.append(ids);cols.append(ids);vals.append(diag)
    return coo_matrix((np.concatenate(vals),(np.concatenate(rows),np.concatenate(cols))),shape=(len(ids),len(ids))).tocsr()

def ground(h,d,L,matrixfree=False):
    n=d**L
    if n <= 500_000 and not matrixfree:
        H=ring(h,d,L)
    else:
        hr=np.real_if_close(h,tol=1000)
        axes=[([j,(j+1)%L]+[k for k in range(L) if k not in (j,(j+1)%L)]) for j in range(L)]
        inverses=[np.argsort(a) for a in axes]
        def mv(v):
            t=v.reshape((d,)*L); out=np.zeros_like(t)
            for a,inv in zip(axes,inverses):
                z=(hr@t.transpose(a).reshape(d*d,-1)).reshape((d,)*L)
                out+=z.transpose(inv)
            return out.ravel()
        H=LinearOperator((n,n),matvec=mv,dtype=hr.dtype)
    if n<=128 and not matrixfree:
        e,V=np.linalg.eigh(H.toarray()); e=e[0];v=V[:,0]
    else:
        x=np.arange(n,dtype=float);v0=np.cos(x*np.sqrt(2))+np.sin((x+1)*np.sqrt(3))
        es,vs=eigsh(H,k=1,which='SA',tol=2e-11,v0=v0,maxiter=12000,ncv=30);e=float(es[0]);v=vs[:,0]
    residual=float(np.linalg.norm(H@v-e*v))
    return float(e),residual,n

def run(which=None,large=False):
    dest=Path(('ed_large_'+(which or 'all')+'.jsonl') if large else 'ed_results.jsonl')
    done=set()
    if dest.exists():
        for line in dest.read_text().splitlines():
            a=json.loads(line);done.add((a['family'],a['parameter'],a['L']))
    for family,p,d,h,ls in bonds():
        if which and family!=which:continue
        if large:
            ls={'xxz':[18,20,24], 'ashkin_teller':[10,12], 'spin_orbital':[12], 'spin1_blbq':([12,16] if p<0 else [12,18]),'fz_parafermion':{2:[20,24],3:[12,16],4:[10,12],5:[8,10],6:[8,10]}[d] if family=='fz_parafermion' else []}.get(family,ls)
        # Reflection averaging proves that the ungauged half-spread is optimal.
        sw=np.arange(d*d).reshape(d,d).T.ravel()
        assert np.max(np.abs(h-h.conj().T))<1e-11
        assert np.max(np.abs(h-h[sw][:,sw]))<1e-11
        ev=np.linalg.eigvalsh(h); M=float((ev[-1]-ev[0])/2)
        for L in ls:
            if (family,p,L) in done:continue
            if large and d**L>20_000_000:
                print(f'Skipping {family} L={L}: dimension exceeds configured 20000000 limit',flush=True)
                continue
            t=time.time();e,res,n=ground(h,d,L)
            a=dict(family=family,parameter=p,d=d,L=L,E=e,M=M,residual=res,dimension=n,seconds=time.time()-t,tag='numerical')
            with dest.open('a') as f:f.write(json.dumps(a)+'\n')
            print(f'{family} {p:g} L={L} E={e:.12g} M={M:.8g} residual={res:.2g} seconds={a["seconds"]:.1f}',flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--family');ap.add_argument('--large',action='store_true');args=ap.parse_args();run(args.family,args.large)
