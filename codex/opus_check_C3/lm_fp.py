import numpy as np, scipy.linalg as la
rng=np.random.default_rng(1)
h14=-0.25*np.log(0.25)-0.75*np.log(0.75)
I2=np.eye(2)
def ptr(r,keep,dims):
    n=len(dims); r=r.reshape(dims+dims)
    for i in reversed(range(n)):
        if i not in keep: r=np.trace(r,axis1=i,axis2=i+r.ndim//2)
    k=int(np.prod([dims[i] for i in keep])); return r.reshape(k,k)
def logm(r):
    w,v=la.eigh((r+r.conj().T)/2); return (v*np.log(np.clip(w,1e-300,None)))@v.conj().T
def step(s):
    t=ptr(s,[0],[2,2]); ls=logm(s); lt=logm(t)
    X=np.kron(ls,I2)+np.kron(I2,ls)-0.25*(np.kron(lt,np.eye(4))+np.kron(np.eye(4),lt))-0.75*np.kron(np.kron(I2,lt),I2)
    w,v=la.eigh((4/3)*X); m=w.max(); E=(v*np.exp(w-m))@v.conj().T; Z=np.trace(E).real
    lam=-(3/4)*(np.log(Z)+m); r3=E/Z
    a=ptr(r3,[0,1],[2,2,2]); b=ptr(r3,[1,2],[2,2,2])
    return r3,a,b,lam
def solve(s,iters=4000,damp=0.5):
    for k in range(iters):
        r3,a,b,lam=step(s); new=(a+b)/2
        new=(new+new.conj().T)/2
        res=max(np.abs(a-s).max(),np.abs(b-s).max())
        if res<1e-12: break
        s=(1-damp)*s+damp*new
    return s,lam,res,k
def rand_sigma(kind):
    if kind=='haar':
        G=rng.normal(size=(4,4))+1j*rng.normal(size=(4,4)); s=G@G.conj().T
    elif kind=='real':
        G=rng.normal(size=(4,4)); s=G@G.T
    elif kind=='diag':
        s=np.diag(rng.random(4))
    s=s/np.trace(s)
    s=(s+ptr(s,[0,1],[2,2]).reshape(2,2,2,2).transpose(1,0,3,2).reshape(4,4))/2  # symmetrize under swap for equal marginals
    return s/np.trace(s)
out=[]
for kind in ['haar','real','diag']:
    for t in range(60):
        s0=rand_sigma(kind); s,lam,res,k=solve(s0)
        ev=np.sort(la.eigvalsh(s))
        out.append((kind,res,lam,3*lam/h14,ev))
conv=[o for o in out if o[1]<1e-9]
print("converged",len(conv),"of",len(out))
cs=np.array([o[3] for o in conv]); print("c_Delta range",cs.min(),cs.max())
for o in conv[:5]+conv[-5:]: print(o[0],round(o[3],6),np.round(o[4],6))
uniq=sorted(set(np.round(cs,5))); print("distinct c_Delta:",uniq[:20])
