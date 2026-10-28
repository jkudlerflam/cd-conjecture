"""Small C15 ED checks. NumPy/SciPy only; writes only next to this file.
Run with OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 checks.py.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
import json, time
from pathlib import Path
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh
from scipy.linalg import eigvalsh
ROOT = Path(__file__).resolve().parent

def basis(L, parity=None, N=None):
    x = np.arange(1 << L, dtype=np.int64)
    count = np.zeros(len(x), dtype=np.int8)
    for j in range(L): count += ((x >> j) & 1).astype(np.int8)
    if N is not None: return x[count == N]
    return x if parity is None else x[count % 2 == parity]

def ham(L, Jx, Jy, Jz, D, g, parity=None, N=None):
    """Periodic SPIN chain, no fermion boundary sign ambiguity."""
    b = basis(L, parity, N)
    rows, cols, vals = [], [], []
    diag = np.zeros(len(b)); ii = np.arange(len(b))
    for j in range(L):
        k = (j+1) % L
        zj = 1-2*((b >> j) & 1); zk = 1-2*((b >> k) & 1)
        diag += Jz*zj*zk - g*zj
        amp = Jx-Jy*zj*zk + 1j*D*(zk-zj)
        mask = np.abs(amp) > 1e-15
        target = b[mask] ^ ((1 << j) | (1 << k))
        rows.append(np.searchsorted(b, target)); cols.append(ii[mask]); vals.append(amp[mask])
    rows.append(ii); cols.append(ii); vals.append(diag)
    H = coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                   shape=(len(b), len(b))).tocsr()
    if D == 0: H = H.real
    return H, b

def low(H, k=3):
    if H.shape[0] <= 10:
        e, u = np.linalg.eigh(H.toarray()); return e[:k], u[:,:k]
    # Deterministic generic starting vector.
    v = np.cos(np.arange(H.shape[0]) * np.sqrt(2.0))
    e,u = eigsh(H, k=k, which='SA', tol=2e-10, v0=v, maxiter=12000)
    order = np.argsort(e); return e[order],u[:,order]

def entropy(psi,L,l):
    s = np.linalg.svd(psi.reshape(1 << (L-l), 1 << l), compute_uv=False)
    p=s*s; p=p[p>1e-15]
    return float(-p@np.log(p))

def full_case(name,L,params):
    es=[]; gs=[]
    for p in (0,1):
        H,b=ham(L,**params,parity=p); e,u=low(H)
        es.append(e)
        gs.append((b,u[:,0],float(np.linalg.norm(H@u[:,0]-e[0]*u[:,0]))))
    p=int(es[1][0]<es[0][0]); b,u,res=gs[p]
    psi=np.zeros(1 << L,dtype=complex); psi[b]=u
    cuts=np.arange(max(2,L//4),L//2+1)
    ss=np.array([entropy(psi,L,int(l)) for l in cuts])
    x=np.log(L/np.pi*np.sin(np.pi*cuts/L))
    c=float(3*np.polyfit(x,ss,1)[0])
    return dict(name=name,L=L,params=params,E0=float(es[p][0]),parity=p,
                levels=[a.tolist() for a in es],residual=res,c_chord=c,S_half=float(ss[-1]))

def few_particle(L,V,N):
    # H=-hopping+V nn; subtract centered-density constant in spin form.
    H,b=ham(L,-.5,-.5,V/4,0,0,N=N)
    e,u=low(H,k=1)
    E=float(e[0]-V*L/4+V*N)
    return E,float(np.linalg.norm(H@u[:,0]-e[0]*u[:,0]))

def constrained(L,N):
    # Periodic hard-exclusion gas, same matrix elements as t=1 spin hopping.
    b=basis(L,N=N)
    shift=((b << 1)&((1 << L)-1)) | (b >> (L-1))
    b=b[(b & shift)==0]; index={int(v):i for i,v in enumerate(b)}
    rows=[]; cols=[]
    for j in range(L):
        k=(j+1)%L
        mask=((b >> j)&1)!=((b >> k)&1)
        for c,target in zip(np.where(mask)[0],b[mask]^((1 << j)|(1 << k))):
            r=index.get(int(target))
            if r is not None: rows.append(r);cols.append(int(c))
    H=coo_matrix((-np.ones(len(rows)),(rows,cols)),shape=(len(b),len(b))).tocsr()
    e,u=low(H,k=3); prob=u[:,0]**2
    q=2*np.pi/L
    nq=sum(((b >> j)&1)*np.exp(1j*q*j) for j in range(L))
    K=float(2*np.pi/q*(prob@np.abs(nq)**2)/L)
    # c_E uses the thermodynamic energy and analytic velocity of squeezed free fermions.
    rho=N/L; kf=np.pi*rho/(1-rho)
    einf=-2*(1-rho)/np.pi*np.sin(kf); vel=2*np.sin(kf)/(1-rho)
    cE=float(-6*L/(np.pi*vel)*(e[0]-L*einf))
    psi=np.zeros(1 << L);psi[b]=u[:,0]
    cuts=np.arange(max(2,L//4),L//2+1)
    ss=[entropy(psi,L,int(l)) for l in cuts]
    cS=float(3*np.polyfit(np.log(L/np.pi*np.sin(np.pi*cuts/L)),ss,1)[0])
    return dict(L=L,N=N,dim=len(b),E0=float(e[0]),K_structure=K,K_exact=(1-rho)**2,
                c_energy=cE,c_chord=cS,S_half=float(ss[-1]),
                residual=float(np.linalg.norm(H@u[:,0]-e[0]*u[:,0])))

def main():
    out={'cases':[],'few_particle':[],'constrained':[]}
    # Independent dense Kronecker validation of signs and both parity sectors.
    X=np.array([[0,1],[1,0]]);Y=np.array([[0,-1j],[1j,0]]);Z=np.diag([1,-1]);I=np.eye(2)
    L=4; pars=dict(Jx=-.65,Jy=-.35,Jz=.17,D=.23,g=.41)
    Hd=np.zeros((1 << L,1 << L),complex)
    def op(ops):
        a=np.array([[1.]])
        for j in reversed(range(L)): a=np.kron(a,ops.get(j,I))
        return a
    for j in range(L):
        k=(j+1)%L
        for A,B,w in [(X,X,pars['Jx']),(Y,Y,pars['Jy']),(Z,Z,pars['Jz']),
                      (X,Y,pars['D']),(Y,X,-pars['D'])]: Hd+=w*op({j:A,k:B})
        Hd-=pars['g']*op({j:Z})
    H,b=ham(L,**pars)
    out['dense_check']=float(np.max(np.abs(H.toarray()-Hd)))
    cases={
      'ising':dict(Jx=-1,Jy=0,Jz=0,D=0,g=1),
      'xx':dict(Jx=-.5,Jy=-.5,Jz=0,D=0,g=0),
      'free_collision':dict(Jx=-.75,Jy=-.25,Jz=0,D=.25,g=1),
      'collision_interacting_probe':dict(Jx=-.75,Jy=-.25,Jz=.1,D=.25,g=1.2),
      'attractive_pair_probe':dict(Jx=-.55,Jy=-.45,Jz=-.55,D=.1,g=0),
    }
    for name,params in cases.items():
        for L in (10,12,14,16):
            start=time.time();rec=full_case(name,L,params);out['cases'].append(rec)
            print(name,L,'cS',round(rec['c_chord'],6),'S',round(rec['S_half'],6),
                  'sec',round(time.time()-start,2),flush=True)
            (ROOT/'checks.json').write_text(json.dumps(out,indent=2))
    for V in (-1.8,-2.,-2.2,-4.):
        vals=[]
        for N in (1,2,3,4):
            E,res=few_particle(24,V,N);vals.append(E)
            out['few_particle'].append(dict(L=24,V=V,N=N,E=E,E_per_particle=E/N,residual=res))
        print('binding',V,vals,flush=True)
    for L,N in ((12,4),(15,5),(18,6)):
        rec=constrained(L,N);out['constrained'].append(rec);print('constrained',rec,flush=True)
    out['free_slow_branch']=[]
    for L in (32,64,128,256):
        q=2*np.pi/L
        # Rationalized small eigenvalue, avoids cancellation.
        a=2*(1-np.cos(q));b=np.sin(q)
        gap=a*a/(np.sqrt(a*a+b*b)+b)
        out['free_slow_branch'].append(dict(L=L,gap=float(gap),Lgap=float(L*gap),L3gap=float(L**3*gap)))
    (ROOT/'checks.json').write_text(json.dumps(out,indent=2))
    print('dense check',out['dense_check'],flush=True)
if __name__=='__main__': main()
