"""Closed-book C3 computations. Natural logarithms; Pauli normalization Tr(P_i P_j)=4 delta_ij.
Run with OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 c3.py <mode>.
Only small matrix calculations and explicitly bounded finite-chain ED are used.
"""
import os, sys, json, time
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import least_squares, root, brentq
I=np.eye(2,dtype=complex)
P=np.array([[[0,1],[1,0]],[[0,-1j],[1j,0]],[[1,0],[0,-1]]],complex)
def kron(*args):
    a=np.array([[1]],complex)
    for b in args:a=np.kron(a,b)
    return a
B1=np.array([kron(p,I)+kron(I,p) for p in P])
BT=np.array([kron(p,q) for p in P for q in P])
B=np.concatenate([B1,BT])
PA=np.array([kron(p,I) for p in P]+[kron(I,p) for p in P]+list(BT))
D3=np.array([kron(b,I)+kron(I,b) for b in B])
ON=np.array([kron(p,I,I)/3+kron(I,p,I)+kron(I,I,p)/3 for p in P])
h=lambda x:-x*np.log(x)-(1-x)*np.log1p(-x)
def logm(a):
    e,u=eigh(a,check_finite=False)
    if e[0]<=0:raise ValueError('not positive')
    return (u*np.log(e))@u.conj().T

def gibbs(a):
    e,u=eigh(a,check_finite=False); z=np.logaddexp.reduce(e)
    return (u*np.exp(e-z))@u.conj().T,float(z)

def ptr(a,keep,n):
    keep=list(keep); rest=[j for j in range(n) if j not in keep]
    ar=a.reshape([2]*(2*n)).transpose(keep+rest+[j+n for j in keep]+[j+n for j in rest])
    d=2**len(keep); r=2**len(rest)
    return np.einsum('abcb->ac',ar.reshape(d,r,d,r))

def coeff(s):
    v=np.real(np.einsum('aij,ji->a',PA,s))
    return np.r_[(v[:3]+v[3:6])/2,v[6:]]

def sigma_x(x):return (np.eye(4)+np.einsum('a,aij->ij',x,B))/4

def sigma_g(y):return gibbs(np.einsum('a,aij->ij',y,B))[0]

def y_from_s(s):
    l=logm(s)
    return np.r_[np.real(np.einsum('aij,ji->a',B1,l))/8,np.real(np.einsum('aij,ji->a',BT,l))/4]

def map3(s):
    t=(ptr(s,[0],2)+ptr(s,[1],2))/2
    ls=logm(s);lt=logm(t)
    a=4/3*(kron(ls,I)+kron(I,ls))-(kron(lt,I,I)+kron(I,I,lt))/3-kron(I,lt,I)
    r,z=gibbs(a)
    return r,-.75*z

def residual_x(x,penalty=True):
    s=sigma_x(x);e=np.linalg.eigvalsh(s)
    if e[0]<=1e-14:
        return np.full(30,1e1+1e2*max(0,-e[0]))
    r,la=map3(s)
    return np.r_[np.real(np.einsum('aij,ji->a',PA,ptr(r,[0,1],3)-s)),np.real(np.einsum('aij,ji->a',PA,ptr(r,[1,2],3)-s))]

# Positive Gibbs parametrization with all 15 traceless Hermitian generators.
# Equal marginal constraints are included explicitly, not presumed for this parametrization.
def sigma_full(y):return gibbs(np.einsum('a,aij->ij',y,PA))[0]
def residual_g(y):
    s=sigma_full(y)
    try:
        r,la=map3(s)
    except ValueError:return np.full(33,1e2)
    return np.r_[np.real(np.einsum('aij,ji->a',PA,ptr(r,[0,1],3)-s)),np.real(np.einsum('aij,ji->a',PA,ptr(r,[1,2],3)-s)),np.real(np.einsum('aij,ji->a',P,ptr(s,[0],2)-ptr(s,[1],2)))]

def full_y(s):return np.real(np.einsum('aij,ji->a',PA,logm(s)))/4

def entropy(s):
    e=np.linalg.eigvalsh(s);e=e[e>0];return float(-e@np.log(e))

def describe(s):
    r,la=map3(s); x=coeff(s);v=x[:3];T=x[3:].reshape(3,3)
    return dict(x=x.tolist(),vnorm=float(np.linalg.norm(v)),T_eig=np.linalg.eigvalsh((T+T.T)/2).tolist(),T_singular=np.linalg.svd(T,compute_uv=False).tolist(),T_asym=float(np.linalg.norm(T-T.T)),lambda_=la,c=3*la/h(.25),min_eig=float(np.linalg.eigvalsh(s)[0]),residual=float(np.linalg.norm(residual_x(x))),S1=entropy(ptr(s,[0],2)),S2=entropy(s),S3=entropy(r))

def werner():
    def f(k):
        t=(np.exp(3*k)-1)/(3*np.exp(3*k)+1)
        w=np.exp(np.array([2*k,0,-4*k])-max(2*k,0,-4*k))
        return (2/3)*(w[0]-w[2])/(2*w[0]+w[1]+w[2])-t
    out=[]
    for lo,hi in [(-5,-.001),(.001,5)]:
        k=brentq(f,lo,hi,xtol=1e-14);t=(np.exp(3*k)-1)/(3*np.exp(3*k)+1)
        s=sigma_x(np.r_[np.zeros(3),(t*np.eye(3)).ravel()]);out.append(s)
        print('Werner',k,t,json.dumps(describe(s)),flush=True)
    return out

def scale(s):
    r3,la=map3(s);t=ptr(s,[0],2);L1=logm(t);L2=logm(s);L3=logm(r3)
    # Triple (1,1,2): AB=12, BC=234, A=1, C=34, B=2.
    a=1.5*(kron(L2,I,I)+kron(I,L3))-.5*(kron(L1,I,I,I)+kron(I,I,L2))-kron(I,L1,I,I)
    r4,z=gibbs(a);lam4=-2*z/3
    # Direct K_Delta operator for each 4-site partition.
    L4=logm(r4)
    ks=[]
    for sizes in [(1,1,2),(2,1,1),(1,2,1)]:
        aa,b,c=sizes;eta=aa*c/((aa+b)*(b+c))
        def emb(start,length):
            loc={1:L1,2:L2,3:L3,4:L4}[length]
            return kron(np.eye(2**start),loc,np.eye(2**(4-start-length)))
        K=-emb(0,aa+b)-emb(aa,b+c)+eta*(emb(0,aa)+emb(aa+b,c))+(1-eta)*(emb(aa,b)+L4)
        target=la*h(eta)/h(.25)
        ks.append(dict(sizes=sizes,eta=eta,target=target,expect=float(np.trace(r4@K).real),operator_fro=float(np.linalg.norm(K-target*np.eye(16))),operator_spread=float(np.ptp(np.linalg.eigvalsh(K)))))
    return dict(lambda112=lam4,marg123=float(np.linalg.norm(ptr(r4,[0,1,2],4)-r3)),marg234=float(np.linalg.norm(ptr(r4,[1,2,3],4)-r3)),triples=ks)

def save_records(name,records):
    with open(name,'w') as f:json.dump(records,f,indent=2)


def random_start(rng,kind):
    if kind=='generic':
        y=rng.normal(size=15)*rng.uniform(.1,1.8)
        return sigma_full(y)
    if kind=='near_pure':
        z=rng.normal(size=4)+1j*rng.normal(size=4);z/=np.linalg.norm(z)
        eps=10**rng.uniform(-4,-1)
        return (1-eps)*np.outer(z,z.conj())+eps*np.eye(4)/4
    if kind=='bell':
        w=rng.dirichlet(np.full(4,.4));
        bs=[np.array([1,0,0,1])/np.sqrt(2),np.array([1,0,0,-1])/np.sqrt(2),np.array([0,1,1,0])/np.sqrt(2),np.array([0,1,-1,0])/np.sqrt(2)]
        s=sum(v*np.outer(b,b) for v,b in zip(w,bs));return .98*s+.02*np.eye(4)/4
    if kind=='werner':
        t=rng.uniform(-.95,.32);s=sigma_x(np.r_[np.zeros(3),(t*np.eye(3)).ravel()])
        perturb=sigma_full(rng.normal(size=15));return .98*s+.02*perturb
    if kind=='ising':
        # Critical H=-sum XX-sum Z: mz=2/pi, xx=2/pi, yy=-2/(3pi), zz=16/(3pi^2).
        x=np.r_[0,0,2/np.pi,np.diag([2/np.pi,-2/(3*np.pi),16/(3*np.pi**2)]).ravel()]
        return .96*sigma_x(x)+.04*sigma_full(rng.normal(size=15))
    if kind=='heisenberg':
        t=(1-4*np.log(2))/3
        return .97*sigma_x(np.r_[np.zeros(3),(t*np.eye(3)).ravel()])+.03*sigma_full(rng.normal(size=15))

def roots(n=240,seed=8291):
    rng=np.random.default_rng(seed);records=[];start=time.time()
    kinds=['generic','near_pure','bell','werner','ising','heisenberg']
    for j in range(n):
        kind=kinds[j%len(kinds)];s=random_start(rng,kind)
        sol=least_squares(residual_g,full_y(s),xtol=2e-11,ftol=2e-11,gtol=2e-11,max_nfev=350,diff_step=1e-5)
        s=sigma_full(sol.x); rr=float(np.linalg.norm(residual_g(sol.x)))
        rec=dict(index=j,kind=kind,nfev=sol.nfev,success=bool(sol.success),residual_full=rr,y=sol.x.tolist())
        if np.linalg.eigvalsh(s)[0]>1e-11:rec.update(describe(s))
        else:rec.update(min_eig=float(np.linalg.eigvalsh(s)[0]))
        records.append(rec)
        if j%10==0 or j==n-1:
            save_records('c3_roots.json',records)
            print(j+1,'elapsed',round(time.time()-start,1),'res',rr,'min',rec['min_eig'],'c',rec.get('c'),'T',rec.get('T_singular'),'v',rec.get('vnorm'),flush=True)
    save_records('c3_roots.json',records)


def pfaff(a):
    if a.shape[0]==0:return 1.
    return sum((-1)**(j+1)*a[0,j]*pfaff(np.delete(np.delete(a,[0,j],0),[0,j],1)) for j in range(1,a.shape[0]))

def ising_rho(n):
    from itertools import combinations
    ma=[]
    for j in range(n):
        ma.extend([kron(*([P[2]]*j+[p]+[I]*(n-j-1))) for p in P[:2]])
    G=np.zeros((2*n,2*n))
    for j in range(n):
        for k in range(n):
            G[2*j,2*k+1]=2/(np.pi*(2*(j-k)-1));G[2*k+1,2*j]=-G[2*j,2*k+1]
    rho=np.zeros((2**n,2**n),complex)
    for size in range(0,2*n+1,2):
        for subset in combinations(range(2*n),size):
            Bm=np.eye(2**n,dtype=complex)
            for j in subset:Bm=Bm@ma[j]
            phase=(1j)**(size*(size-1)//2);Bm*=phase
            ex=phase*pfaff(-1j*G[np.ix_(subset,subset)])
            rho+=ex*Bm/(2**n)
    return (rho+rho.conj().T)/2

def heisenberg_rho():
    from scipy.special import zeta
    g1=(1-4*np.log(2))/3;g2=1/3-16*np.log(2)/3+3*zeta(3,1)
    r=np.eye(8,dtype=complex)
    for p in P:r+=g1*(kron(p,p,I)+kron(I,p,p))+g2*kron(p,I,p)
    return r/8,g1,g2

def diagnose_rho3(r):
    s=ptr(r,[0,1],3);t=ptr(r,[0],3);L1=logm(t);L2=logm(s);L3=logm(r)
    K=-kron(L2,I)-kron(I,L2)+.25*(kron(L1,I,I)+kron(I,I,L1))+.75*(kron(I,L1,I)+L3)
    lam=float(np.trace(r@K).real);R=K-lam*np.eye(8)
    rr,lfp=map3(s)
    return dict(S1=entropy(t),S2=entropy(s),S3=entropy(r),lambda_=lam,c=3*lam/h(.25),variance=float(np.trace(r@R@R).real),weighted_residual=float(np.sqrt(np.trace(r@R@R).real)),operator_fro=float(np.linalg.norm(R)),operator_norm=float(np.linalg.norm(R,2)),K_eigenvalues=np.linalg.eigvalsh(K).tolist(),rho_eigenvalues=np.linalg.eigvalsh(r).tolist(),fixed_point_marginal_residual=float(np.linalg.norm(ptr(rr,[0,1],3)-s)),rho3_map_distance=float(np.linalg.norm(rr-r)))


def reduced_roots(n=1600,seed=29183):
    rng=np.random.default_rng(seed);records=[];start=time.time()
    def fun(j):
        y=np.r_[np.zeros(6),j];return residual_g(y)[:30]
    for k in range(n):
        j=rng.normal(size=9)*10**rng.uniform(-1,0.4)
        if k%4==0:j=np.diag(rng.normal(size=3)*rng.uniform(.1,2)).ravel()
        sol=least_squares(fun,j,xtol=1e-12,ftol=1e-12,gtol=1e-12,max_nfev=250,diff_step=2e-5)
        s=sigma_full(np.r_[np.zeros(6),sol.x]);rec=describe(s);rec.update(index=k,nfev=sol.nfev,j=sol.x.tolist())
        records.append(rec)
        if k%100==0 or k==n-1:
            save_records('c3_reduced_roots.json',records)
            print(k+1,'seconds',time.time()-start,'res',rec['residual'],'c',rec['c'],'sing',rec['T_singular'],flush=True)
    save_records('c3_reduced_roots.json',records)


def square_fun(y):
    s=sigma_full(y)
    try:r,_=map3(s)
    except ValueError:return np.full(15,1e2)
    return np.real(np.einsum('aij,ji->a',PA,(ptr(r,[0,1],3)+ptr(r,[1,2],3))/2-s))

def square_roots(n=500,seed=801):
    rng=np.random.default_rng(seed);records=[];start=time.time();kinds=['generic','near_pure','bell','werner','ising','heisenberg']
    for j in range(n):
        kind=kinds[j%6];s=random_start(rng,kind)
        sol=root(square_fun,full_y(s),method='hybr',options={'xtol':1e-10,'maxfev':1500})
        s=sigma_full(sol.x);rec=dict(index=j,kind=kind,nfev=int(sol.nfev),root_success=bool(sol.success),square_residual=float(np.linalg.norm(square_fun(sol.x))),min_eig=float(np.linalg.eigvalsh(s)[0]),y=sol.x.tolist())
        try:rec.update(describe(s));rec['residual_full']=float(np.linalg.norm(residual_g(sol.x)))
        except (ValueError,np.linalg.LinAlgError):pass
        records.append(rec)
        if j%50==0 or j==n-1:
            save_records('c3_square_roots.json',records);print(j+1,'elapsed',time.time()-start,'sq',rec['square_residual'],'full',rec.get('residual_full'),'c',rec.get('c'),flush=True)
    save_records('c3_square_roots.json',records)

def iterate(n=120,seed=5928):
    rng=np.random.default_rng(seed);records=[];start=time.time();kinds=['generic','near_pure','bell','werner','ising','heisenberg']
    for j in range(n):
        kind=kinds[j%6];s=random_start(rng,kind);alpha=[.15,.4,.7][(j//6)%3]
        for it in range(5000):
            try:r,_=map3(s)
            except ValueError:break
            sf=(ptr(r,[0,1],3)+ptr(r,[1,2],3))/2
            change=np.linalg.norm(sf-s);s=(1-alpha)*s+alpha*sf
            if change<1e-11 or np.linalg.eigvalsh(s)[0]<1e-12:break
        rec=dict(index=j,kind=kind,iterations=it+1,alpha=alpha,change=float(change),min_eig=float(np.linalg.eigvalsh(s)[0]))
        try:rec.update(describe(s));rec['residual_full']=float(np.linalg.norm(residual_g(full_y(s))))
        except (ValueError,np.linalg.LinAlgError):pass
        records.append(rec)
        if j%20==0 or j==n-1:
            save_records('c3_iterations.json',records);print(j+1,'elapsed',time.time()-start,'iters',it+1,'res',rec.get('residual_full'),'mineig',rec['min_eig'],'c',rec.get('c'),flush=True)
    save_records('c3_iterations.json',records)


def j_sigma(J):return sigma_full(np.r_[np.zeros(6),np.asarray(J).ravel()])

def canonical_exotic(v):
    a,b,d=v;return np.array([[0,a,0],[b,0,0],[0,0,d]])

def exotic_fun(v):
    J=canonical_exotic(v);s=j_sigma(J);A=np.einsum('a,aij->ij',J.ravel(),BT);r,_=gibbs(4/3*(kron(A,I)+kron(I,A)))
    diff=coeff(ptr(r,[0,1],3)-s)[3:].reshape(3,3)
    return diff[[0,1,2],[1,0,2]]

def numerical_jac(fun,x,eps=1e-5):
    return np.column_stack([(fun(x+eps*np.eye(len(x))[k])-fun(x-eps*np.eye(len(x))[k]))/(2*eps) for k in range(len(x))])

def canonical():
    ss={'tracial':np.eye(4)/4}
    for name,k in [('Werner_AF',-.6178646150823423),('Werner_F',.9740798988555522)]:ss[name]=j_sigma(3*k/4*np.eye(3))
    j=brentq(lambda j:np.tanh(j)-np.tanh(4*np.sqrt(2)*j/3)/np.sqrt(2),.1,2,xtol=1e-14)
    ss['XX_rotation']=j_sigma(np.diag([j,j,0]));ss['XX_reflection']=j_sigma(np.diag([j,-j,0]))
    J=np.zeros((3,3));J[0,1]=j;ss['nilpotent']=j_sigma(J)
    exotic=root(exotic_fun,[.691,.250,-.375],tol=1e-12).x
    ss['exotic_plus']=j_sigma(canonical_exotic(exotic));ss['exotic_minus']=j_sigma(canonical_exotic(exotic*np.array([1,-1,-1])))
    records={}
    for name,s in ss.items():
        rec=describe(s);rec['J']=np.array(full_y(s)[6:]).reshape(3,3).tolist();rec['scale']=scale(s)
        y=full_y(s);sv=np.linalg.svd(numerical_jac(residual_g,y),compute_uv=False)
        rec['jacobian_singular_values']=sv.tolist();rec['nullity_1e-6']=int(np.sum(sv<1e-6))
        rec['rho3_diagnosis']=diagnose_rho3(map3(s)[0]);records[name]=rec
        print(name,'lambda',rec['lambda_'],'c',rec['c'],'J',rec['J'],'T',np.array(rec['x'][3:]).reshape(3,3),'null',rec['nullity_1e-6'],'scale',rec['scale'],flush=True)
    save_records('c3_canonical.json',records)
    return ss


def map4(s,r3=None,la=None):
    if la is None:r3,la=map3(s)
    t=(ptr(s,[0],2)+ptr(s,[1],2))/2;L1=logm(t);L2=logm(s)
    Q=1.5*(kron(L2,I,I)+kron(I,I,L2))+2*kron(I,L2,I)-.5*(kron(L1,I,I,I)+kron(I,I,I,L1))-1.5*(kron(I,L1,I,I)+kron(I,I,L1,I))
    r4,z=gibbs(Q);return r4,-2/3*(z+2*la)

P3=np.array([kron(a,b,c) for a in [I,*P] for b in [I,*P] for c in [I,*P]][1:])
def combined_fun(y):
    s=sigma_full(y)
    try:
        r3,la=map3(s);r4,l4=map4(s,r3,la)
        f=residual_g(y)
        return np.r_[f,np.real(np.einsum('aij,ji->a',P3,ptr(r4,[0,1,2],4)-r3)),np.real(np.einsum('aij,ji->a',P3,ptr(r4,[1,2,3],4)-r3)),l4-la*h(1/3)/h(.25)]
    except (ValueError,np.linalg.LinAlgError):return np.full(160,100.)

def combined(n=180,seed=1708):
    rng=np.random.default_rng(seed);records=[];start=time.time();kinds=['generic','near_pure','bell','werner','ising','heisenberg']
    for j in range(n):
        kind=kinds[j%6];s=random_start(rng,kind)
        sol=least_squares(combined_fun,full_y(s),xtol=1e-10,ftol=1e-10,gtol=1e-10,max_nfev=250,diff_step=1e-5)
        s=sigma_full(sol.x);rec=dict(index=j,kind=kind,nfev=sol.nfev,combined_residual=float(np.linalg.norm(combined_fun(sol.x))),min_eig=float(np.linalg.eigvalsh(s)[0]),y=sol.x.tolist())
        try:rec.update(describe(s))
        except (ValueError,np.linalg.LinAlgError):pass
        records.append(rec)
        if j%15==0 or j==n-1:
            save_records('c3_combined.json',records);print(j+1,'elapsed',time.time()-start,'res',rec['combined_residual'],'min',rec['min_eig'],'c',rec.get('c'),flush=True)
    save_records('c3_combined.json',records)


def heisenberg_ed(lengths=(10,12,14,16,18)):
    from itertools import combinations
    from scipy.sparse import coo_matrix
    from scipy.sparse.linalg import eigsh
    records=[]
    for L in lengths:
        begin=time.time();states=np.array([sum(1<<j for j in cc) for cc in combinations(range(L),L//2)],dtype=np.int64);index={int(s):i for i,s in enumerate(states)}
        n=len(states);rows=[];cols=[];data=[];diag=np.zeros(n)
        for j in range(L):
            k=(j+1)%L;different=((states>>j)^(states>>k))&1;diag+=1-2*different
            ids=np.flatnonzero(different);flips=states[ids]^((1<<j)|(1<<k))
            rows.extend(ids.tolist());cols.extend(index[int(s)] for s in flips);data.extend([2.]*len(ids))
        rows.extend(range(n));cols.extend(range(n));data.extend(diag.tolist())
        H=coo_matrix((data,(rows,cols)),shape=(n,n)).tocsr()
        vals,vecs=eigsh(H,k=1,which='SA',tol=2e-12,v0=np.random.default_rng(L).normal(size=n))
        psi=np.zeros(2**L);psi[states]=vecs[:,0];a=psi.reshape(-1,8);rho=a.T@a
        rec=diagnose_rho3(rho);rec.update(L=L,dimension=n,energy_per_bond=float(vals[0]/L),seconds=time.time()-begin,g1=float(np.trace(rho@kron(P[2],P[2],I)).real),g2=float(np.trace(rho@kron(P[2],I,P[2])).real))
        records.append(rec);print(json.dumps(rec),flush=True)
    save_records('c3_heisenberg_ed.json',records)


def symmetric_iterations(n=80,seed=5932):
    rng=np.random.default_rng(seed);out=[]
    for idx in range(n):
        kind=['Werner_AF','Werner_F','XX','nilpotent'][idx%4];alpha=.35
        t=rng.uniform(-.9,-.01) if kind=='Werner_AF' else rng.uniform(.01,.32) if kind=='Werner_F' else rng.uniform(.01,.9)
        for it in range(3000):
            if kind.startswith('Werner'):
                s=sigma_x(np.r_[np.zeros(3),(t*np.eye(3)).ravel()]);r,_=map3(s);tf=coeff(ptr(r,[0,1],3))[3]
            elif kind=='XX':
                s=j_sigma(np.diag([np.arctanh(t),np.arctanh(t),0]));r,_=map3(s);tf=coeff(ptr(r,[0,1],3))[3]
            else:
                s=sigma_x(np.r_[np.zeros(3),np.array([[0,t,0],[0,0,0],[0,0,0]]).ravel()]);r,_=map3(s);tf=coeff(ptr(r,[0,1],3))[4]
            if abs(tf-t)<2e-13:break
            t=(1-alpha)*t+alpha*tf
        out.append(dict(kind=kind,iterations=it+1,parameter=float(t),residual=describe(s)['residual'],c=describe(s)['c']))
    save_records('c3_iterations_symmetry.json',out);print('symmetric iteration',len(out),'max residual',max(r['residual'] for r in out),'iterations max',max(r['iterations'] for r in out),flush=True)

def checks(seed=301):
    from scipy.spatial.transform import Rotation
    rng=np.random.default_rng(seed);ss=canonical();errs=[];deltas=[]
    for name,s in ss.items():
        J=full_y(s)[6:].reshape(3,3);base=describe(s)
        for _ in range(12):
            R=Rotation.random(random_state=rng).as_matrix();Q=R@J@R.T
            if name.startswith('Werner'):Q=Q@Rotation.random(random_state=rng).as_matrix()
            if name=='XX_rotation':
                phi=rng.uniform(0,2*np.pi);V=np.array([[np.cos(phi),-np.sin(phi),0],[np.sin(phi),np.cos(phi),0],[0,0,1]])
                Q=R@(J@V)@R.T
            z=j_sigma(Q);d=describe(z);errs.append(d['residual']);deltas.append(abs(d['lambda_']-base['lambda_']))
    # Check partial-trace H^3 formula used in analytic reduction.
    J=rng.normal(size=(3,3));A=np.einsum('a,aij->ij',J.ravel(),BT);H=kron(A,I)+kron(I,A);H3=H@H@H
    diff=(ptr(H3,[0,1],3)-ptr(H3,[1,2],3))/2
    measured=np.real(np.einsum('aij,ji->a',BT,diff))/4
    predicted=2*(J@J@J.T-J.T@J@J)
    cubic_error=float(np.linalg.norm(measured.reshape(3,3)-predicted))
    # Exact Ising two-site correlations, checked from the Wick construction.
    ri=ising_rho(3);x=coeff(ptr(ri,[0,1],3));ex=np.r_[0,0,2/np.pi,np.diag([2/np.pi,-2/(3*np.pi),16/(3*np.pi**2)]).ravel()]
    ising_error=float(np.linalg.norm(x-ex))
    # Pure infinite-line Slater construction: its local correlation matrix matches XX rho3.
    j=.7277678875729927;r=np.tanh(j);C=.5*np.eye(3)+r/2*(np.diag(np.ones(2),1)+np.diag(np.ones(2),-1))
    ec,uc=eigh(C);hc=(uc*np.log(ec/(1-ec)))@uc.T
    annih=[kron(*([P[2]]*k+[(P[0]+1j*P[1])/2]+[I]*(2-k))) for k in range(3)]
    Hc=sum(hc[k,l]*annih[k].conj().T@annih[l] for k in range(3) for l in range(3))
    rg,_=gibbs(Hc);pure_xx_error=float(np.linalg.norm(rg-map3(ss['XX_rotation'])[0]))
    psi=np.zeros(64);psi[[9*k for k in range(8)]]=1/np.sqrt(8)
    rb=np.outer(psi,psi);bell_marginal_error=float(np.linalg.norm(ptr(rb,[0,1,2],6)-np.eye(8)/8))
    shifted=np.array([psi[((k<<1)&63)|(k>>5)] for k in range(64)])
    bell_translation_error=float(np.linalg.norm(shifted-psi))
    report=dict(twisted_states=96,max_fixed_point_residual=max(errs),max_lambda_change=max(deltas),H3_formula_error=cubic_error,ising_two_site_formula_error=ising_error,pure_XX_local_marginal_error=pure_xx_error,pure_Bell_ring_marginal_error=bell_marginal_error,pure_Bell_ring_translation_error=bell_translation_error)
    assert pure_xx_error<1e-12 and bell_marginal_error<1e-12 and bell_translation_error<1e-12
    assert max(errs)<1e-11 and max(deltas)<1e-11 and cubic_error<1e-11 and ising_error<1e-12
    save_records('c3_checks.json',report);print(json.dumps(report),flush=True)


def canonical_search(n=1200,seed=29302):
    rng=np.random.default_rng(seed);records=[];start=time.time()
    def diagonal_fun(v):
        J=np.diag(v);s=j_sigma(J);A=np.einsum('a,aij->ij',J.ravel(),BT);r,_=gibbs(4/3*(kron(A,I)+kron(I,A)));return np.diag(coeff(ptr(r,[0,1],3)-s)[3:].reshape(3,3))
    for i in range(n):
        kind='diagonal' if i%2==0 else 'offdiagonal';fun=diagonal_fun if i%2==0 else exotic_fun
        x=rng.normal(size=3)*10**rng.uniform(-1,.5)
        sol=least_squares(fun,x,xtol=1e-12,gtol=1e-12,ftol=1e-12,max_nfev=250,diff_step=1e-5)
        J=np.diag(sol.x) if i%2==0 else canonical_exotic(sol.x);s=j_sigma(J)
        try:rec=describe(s)
        except (ValueError,np.linalg.LinAlgError):rec=dict(min_eig=float(np.linalg.eigvalsh(s)[0]),residual=100.)
        rec.update(index=i,kind=kind,parameters=sol.x.tolist());records.append(rec)
        if i%100==0 or i==n-1:
            save_records('c3_canonical_search.json',records);print(i+1,'elapsed',time.time()-start,'res',rec['residual'],'c',rec.get('c'),flush=True)
    save_records('c3_canonical_search.json',records)


def polish_square():
    records=json.load(open('c3_square_roots.json'));polished=[]
    for rec in records:
        if rec.get('residual_full',100)<1e-7 and rec['min_eig']>1e-7:
            sol=least_squares(residual_g,np.array(rec['y']),xtol=1e-14,ftol=1e-14,gtol=1e-14,max_nfev=100,diff_step=1e-5)
            z=describe(sigma_full(sol.x));z.update(index=rec['index'],kind=rec['kind'],residual_full=float(np.linalg.norm(residual_g(sol.x))))
            polished.append(z)
    save_records('c3_square_polished.json',polished)
    print('polished',len(polished),'maxres',max(r['residual_full'] for r in polished),'maxv',max(r['vnorm'] for r in polished),flush=True)


def critical_report():
    rh,g1,g2=heisenberg_rho()
    out={'Ising':diagnose_rho3(ising_rho(3)),'Heisenberg':diagnose_rho3(rh),'Heisenberg_g1':g1,'Heisenberg_g2':g2}
    save_records('c3_critical.json',out);print(json.dumps(out,indent=2),flush=True)

def summarize():
    from collections import Counter
    out={}
    for label,filename in [('full','c3_roots.json'),('unpolarized','c3_reduced_roots.json'),('square','c3_square_roots.json'),('canonical','c3_canonical_search.json'),('iteration','c3_iterations.json')]:
        rs=json.load(open(filename));good=[r for r in rs if r.get('residual_full',r.get('residual',100))<1e-7 and r['min_eig']>1e-7]
        out[label]=dict(total=len(rs),accepted=len(good),boundary=sum(r['min_eig']<=1e-7 for r in rs),counts=dict(Counter(f"{r['c']:.6f}" for r in good)),maximum_residual=max((r.get('residual_full',r.get('residual',100)) for r in good),default=None),maximum_magnetization=max((r.get('vnorm',0) for r in good),default=None))
    save_records('c3_summary.json',out);print(json.dumps(out,indent=2),flush=True)


def main():
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode',choices=['werner','roots','reduced','square','iterate','symiterate','canonical','canonical_search','combined','critical','ed','checks','polish','summary'])
    parser.add_argument('n',nargs='?',type=int)
    args=parser.parse_args()
    if args.mode=='werner':
        for s in [np.eye(4)/4]+werner():print(json.dumps(scale(s)))
    elif args.mode=='roots':roots(args.n or 2400)
    elif args.mode=='reduced':reduced_roots(args.n or 2000)
    elif args.mode=='square':square_roots(args.n or 1000)
    elif args.mode=='iterate':iterate(args.n or 120)
    elif args.mode=='canonical_search':canonical_search(args.n or 1200)
    elif args.mode=='combined':combined(args.n or 180)
    else:
        {'symiterate':symmetric_iterations,'canonical':canonical,'critical':critical_report,'ed':heisenberg_ed,'checks':checks,'polish':polish_square,'summary':summarize}[args.mode]()

if __name__=='__main__':main()
