"""C126 ED and permutation trial checks. Python, NumPy and SciPy only.
Large runs require --cluster within a SLURM allocation. All outputs stay here.
"""
from pathlib import Path
import argparse, json, os, time
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh, LinearOperator
from scipy.optimize import minimize
ROOT=Path(__file__).resolve().parent

def basis(d,n):
    target=n//d
    codes=np.zeros(1,dtype=np.int64)
    counts=np.zeros((1,d-1),dtype=np.int8)
    for pos in range(n):
        cc=(codes[:,None]*d+np.arange(d)).reshape(-1)
        nn=np.repeat(counts,d,axis=0)
        for a in range(d-1):nn[a::d,a]+=1
        last=pos+1-nn.sum(axis=1)
        rem=n-pos-1
        good=((nn<=target)&(nn+rem>=target)).all(axis=1)&(last<=target)&(last+rem>=target)
        codes,counts=cc[good],nn[good]
    return codes

def ground(d,n,delta=1.,cluster=False):
    codes=basis(d,n);dim=len(codes);powers=d**np.arange(n,dtype=np.int64)
    inv=np.full(d**n,-1,dtype=np.int32);inv[codes]=np.arange(dim,dtype=np.int32)
    perms=[];zz=np.zeros(dim)
    for x in range(n):
        y=(x+1)%n;a=codes//powers[x]%d;b=codes//powers[y]%d
        perms.append(inv[codes+(b-a)*(powers[x]-powers[y])])
        if d==2:zz+=np.where(a==b,1.,-1.)
    del inv
    factor=.5 if d==2 else 1.
    diag=-n/4+(delta-1)*zz/4 if d==2 else np.zeros(dim)
    def mv(v):
        out=diag*v
        for p in perms:out+=factor*v[p]
        return out
    if cluster and dim>500000:
        H=LinearOperator((dim,dim),matvec=mv,dtype=np.float64)
    else:
        rows=np.tile(np.arange(dim,dtype=np.int32),n+1)
        cols=np.concatenate(perms+[np.arange(dim,dtype=np.int32)])
        data=np.concatenate([np.full(n*dim,factor),diag])
        H=coo_matrix((data,(rows,cols)),shape=(dim,dim)).tocsr()
    v0=np.sin(np.arange(dim)*1.131+0.37)+np.cos(np.arange(dim)*.731)
    ev,vec=eigsh(H,k=1,which='SA',tol=2e-12,v0=v0,maxiter=10000)
    v=vec[:,0];e=float(ev[0]);res=float(np.linalg.norm(H@v-e*v))
    return e,codes,v,res,dim

def bond(d,delta):
    if d==3:return np.eye(9).reshape(3,3,3,3).transpose(1,0,2,3).reshape(9,9)
    sx=np.array([[0,1],[1,0]])/2;sy=np.array([[0,-1j],[1j,0]])/2;sz=np.diag([.5,-.5])
    return (np.kron(sx,sx)+np.kron(sy,sy)+delta*np.kron(sz,sz)).real

def marginal(v,d,n,r):
    mat=v.reshape(d**r,d**(n-r))
    return mat@mat.conj().T

def trial_parameters(v,d,n,h):
    rho=marginal(v,d,n,2);sigma=marginal(v,d,n,1)
    e=float(np.trace(rho@h).real);p=float(np.trace(np.kron(sigma,sigma)@h).real)
    # Basis order a1,b1,a2,b2. Swap a1,a2 and b1,b2.
    ids=np.arange(d**4).reshape((d,)*4)
    sx=np.eye(d**4)[ids.transpose(2,1,0,3).reshape(-1)]
    sy=np.eye(d**4)[ids.transpose(0,3,2,1).reshape(-1)]
    ident=np.eye(d**4);ops=[ident,sx,sy,sx@sy]
    H=np.kron(h,np.eye(d*d))+np.kron(np.eye(d*d),h)
    R=np.kron(rho,rho)
    K=np.array([[np.trace(R@a.T@H@b) for b in ops] for a in ops])
    def direct(t,u):
        coeff=np.array([np.cos(t)*np.cos(u),1j*np.sin(t)*np.cos(u),1j*np.cos(t)*np.sin(u),-np.sin(t)*np.sin(u)])
        return float(np.vdot(coeff,K@coeff).real)
    a=p-e;b=direct(np.pi/4,np.pi/4)-(e+p)
    err=0.
    for t,u in [(.17,.21),(-.35,.91),(1.2,-.73),(.53,1.93)]:
        err=max(err,abs(direct(t,u)-(e+p+(e-p)*np.cos(2*t)*np.cos(2*u)+b*np.sin(2*t)*np.sin(2*u))))
    return a,b,e,p,err

def optimize_angles(n,a,b):
    # phi=2 theta. Twisted seam has antiperiodic vector (cos phi,sin phi).
    signs=np.ones(n);signs[-1]=-1
    def fg(phi):
        z=np.cos(phi);w=np.sin(phi);zn=np.roll(z,-1);wn=np.roll(w,-1)
        f=n*a+np.sum(signs*(-a*z*zn+b*w*wn))
        g=signs*(a*w*zn+b*z*wn)+np.roll(signs*(a*z*wn+b*w*zn),1)
        return f,g
    starts=[np.zeros(n),np.pi*np.arange(n)/n,np.pi*(-1.)**np.arange(n)/2]
    for center in [0.,.5,1.5]:
        for width in [.5,1.,2.]:
            z=np.tanh((np.arange(n)-(n-1)/2-center)/width)
            starts.append((-1.)**np.arange(n)*np.arccos(z))
    for j in range(5):starts.append(np.sin(np.arange(n)*(.61+j*.21)+j)*2)
    results=[minimize(fg,s,jac=True,method='BFGS',options={'gtol':2e-10,'maxiter':1500}) for s in starts]
    best=min(results,key=lambda r:r.fun)
    uniform=min(fg(np.pi*(np.arange(n)+shift)/n)[0] for shift in [0.,.5])
    return float(best.fun),float(uniform),best.x.tolist(),float(np.max(np.abs(fg(best.x)[1])))

def coherent(v,d,n,a):
    pur=[1.]
    for r in range(1,n):
        rr=min(r,n-r);rho=marginal(v,d,n,rr);pur.append(float(np.trace(rho@rho).real))
    return 2*a/sum(pur),pur

def run(cluster=False):
    start=time.time();records=[];energies={}
    models=[('Heisenberg',2,1.),('XXZ_0.5',2,.5),('XXZ_2',2,2.),('SU3',3,0.)]
    for name,d,delta in models:
        lengths=([2,4,6,8,10] if d==2 else ([3,6,9] if cluster else [3,6]))
        needed=sorted(set(lengths+[2*l for l in lengths]));states={}
        for n in needed:
            now=time.time();e,c,v,res,dim=ground(d,n,delta,cluster)
            energies[f'{name}_{n}']={'energy':e,'residual':res,'dimension':dim,'seconds':time.time()-now}
            if n in lengths:
                full=np.zeros(d**n);full[c]=v;states[n]=full
            print(name,n,e,'residual',res,'dim',dim,flush=True)
        h=bond(d,delta)
        for n in lengths:
            a,b,e,p,err=trial_parameters(states[n],d,n,h)
            opt,uniform,angles,grad=optimize_angles(n,a,b)
            coh,pur=coherent(states[n],d,n,a)
            defect=energies[f'{name}_{2*n}']['energy']-2*energies[f'{name}_{n}']['energy']
            M=1. if d==3 else (.5 if delta<=1 else (delta+1)/4)
            rec=dict(model=name,L=n,E_L=energies[f'{name}_{n}']['energy'],E_2L=energies[f'{name}_{2*n}']['energy'],D=defect,M=M,a=a,b=b,
                bare=2*a,fractional=opt,uniform=uniform,coherent=coh,
                ratio_fractional=opt/defect,ratio_uniform=uniform/defect,ratio_coherent=coh/defect,
                angles=angles,gradient=grad,formula_error=err,purities=pur)
            assert err<1e-10, (name,n,'fractional contraction',err)
            assert min(2*a,opt,uniform,coh)>=defect-1e-9, (name,n,'variational bound')
            if name=='Heisenberg' and n>=4:
                assert opt>=3/32-1e-10
                assert abs(b-4*a*a/3)<1e-10
            records.append(rec)
    assert max(v['residual'] for v in energies.values())<1e-8
    out=dict(status='numerical',seconds=time.time()-start,cluster=cluster,energies=energies,records=records)
    fname='checks_cluster.json' if cluster else 'checks_local.json'
    (ROOT/fname).write_text(json.dumps(out,indent=2)+'\n')
    lines=['# C126 numerical tables','', '[Numerical] Energies use balanced-charge ED. Fractional means the best multistart variational profile found, not a certified global minimum. Coherent uses the proved C114 optimum in its seam span.','', '| Model | L | D(L) | Bare cost | Uniform cost | Optimized fractional / D | Coherent / D |','| :--- | ---: | ---: | ---: | ---: | ---: | ---: |']
    for r in records:lines.append(f"| {r['model']} | {r['L']} | {r['D']:.10g} | {r['bare']:.8g} | {r['uniform']:.8g} | {r['ratio_fractional']:.7g} | {r['ratio_coherent']:.7g} |")
    lines+=['',f"[Numerical] Maximum eigenpair residual {max(v['residual'] for v in energies.values()):.3g}. Maximum contraction formula error {max(r['formula_error'] for r in records):.3g}. Elapsed seconds {out['seconds']:.2f}."]
    (ROOT/('TABLES_CLUSTER.md' if cluster else 'TABLES.md')).write_text('\n'.join(lines)+'\n')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--cluster',action='store_true');args=p.parse_args()
    if args.cluster and not os.environ.get('SLURM_JOB_ID'):raise SystemExit('Large computation requires a SLURM allocation.')
    run(args.cluster)
