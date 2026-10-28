"""Full microscopic three-state vacancy/XX control; no code-space projection.

Run with numpy/scipy. All generated files stay beside this script.
"""
from pathlib import Path
from itertools import combinations
import argparse, json, time
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh
from scipy.linalg import svdvals

HERE=Path(__file__).resolve().parent

def basis(L, up, down):
    powers=3**np.arange(L,dtype=np.int64)
    out=[]
    for occ in combinations(range(L),up+down):
        baseline=int(sum(powers[list(occ)]))
        for ups in combinations(occ,up):
            out.append(baseline+int(sum(powers[list(ups)])))
    return np.sort(np.asarray(out,dtype=np.int64))

def build(L,up,down):
    codes=basis(L,up,down); dim=len(codes)
    powers=3**np.arange(L,dtype=np.int64)
    digits=(codes[:,None]//powers[None,:])%3
    rows=[[],[]]; cols=[[],[]]
    for x in range(L):
        y=(x+1)%L; a=digits[:,x]; b=digits[:,y]
        for kind,mask in enumerate(((a==0)^(b==0),(a+b==3)&(a!=0)&(b!=0))):
            ii=np.flatnonzero(mask)
            target=codes[ii]+(b[ii]-a[ii])*(powers[x]-powers[y])
            jj=np.searchsorted(codes,target)
            assert np.all(codes[jj]==target)
            rows[kind].append(ii); cols[kind].append(jj)
    mats=[]
    for r,c in zip(rows,cols):
        r=np.concatenate(r); c=np.concatenate(c)
        mat=coo_matrix((np.ones(len(r)),(r,c)),shape=(dim,dim)).tocsr()
        assert (mat-mat.T).nnz==0
        mats.append(mat)
    return codes,mats

def ground(H,v0=None,k=1):
    # Deterministic non-symmetry-adapted start avoids missing ground momentum.
    if v0 is None: v0=np.random.default_rng(1739).normal(size=H.shape[0])
    val,vec=eigsh(H,k=k,which='SA',v0=v0,tol=2e-11,ncv=max(24,2*k+1),maxiter=15000)
    ix=np.argsort(val); val=val[ix]; vec=vec[:,ix]
    res=[float(np.linalg.norm(H@vec[:,i]-val[i]*vec[:,i])) for i in range(k)]
    return val,vec,res

def entropy(codes,psi,L,ell,up,down):
    a=codes%(3**ell); b=codes//(3**ell)
    d=(a[:,None]//(3**np.arange(ell)))%3
    au=np.sum(d==2,axis=1); ad=np.sum(d==1,axis=1)
    eigen=[]
    for u in range(up+1):
        for v in range(down+1):
            ix=np.flatnonzero((au==u)&(ad==v))
            if not len(ix): continue
            _,ai=np.unique(a[ix],return_inverse=True)
            _,bi=np.unique(b[ix],return_inverse=True)
            mat=np.zeros((int(ai.max())+1,int(bi.max())+1))
            mat[ai,bi]=psi[ix]
            ss=svdvals(mat,check_finite=False)
            eigen.extend(ss[ss>1e-15]**2)
    p=np.array(eigen)
    assert abs(p.sum()-1)<1e-9
    return float(-np.dot(p,np.log(p)))

def independent_small_check(L,up,down,codes,mats):
    if L>8:return None
    # Independently assemble by literal tensor-product two-site 9x9 bond.
    h=np.zeros((9,9))
    for spin in [1,2]: h[spin,3*spin]=h[3*spin,spin]=-1
    h[5,7]=h[7,5]=-.173
    index={int(v):i for i,v in enumerate(codes)}
    dense=np.zeros((len(codes),len(codes)))
    for i,code in enumerate(codes):
        ds=[(int(code)//3**x)%3 for x in range(L)]
        for x in range(L):
            y=(x+1)%L; src=3*ds[x]+ds[y]
            for dest in np.flatnonzero(h[:,src]):
                a,b=divmod(int(dest),3)
                target=int(code)+(a-ds[x])*3**x+(b-ds[y])*3**y
                dense[index[target],i]+=h[dest,src]
    delta=dense-(-mats[0]-.173*mats[1]).toarray()
    return float(np.max(np.abs(delta)))

def run(L,js):
    started=time.time(); up=down=L//4
    assert L%4==0
    codes,mats=build(L,up,down)
    out={'L':L,'N':up+down,'Nup':up,'dimension':len(codes),'hopping_nnz':mats[0].nnz,
         'exchange_nnz':mats[1].nnz,'literal_tensor_check':independent_small_check(L,up,down,codes,mats),'runs':[]}
    psi=None
    for J in js:
        H=-mats[0]-J*mats[1]
        ev,vs,res=ground(H,psi,k=2 if L<=12 else 1); psi=vs[:,0]
        ents={str(ell):entropy(codes,psi,L,ell,up,down) for ell in range(1,L//2+1)}
        fits={}
        for cut in [1,2,3]:
            ells=np.arange(cut,L//2+1)
            if len(ells)<3:continue
            x=np.log(L/np.pi*np.sin(np.pi*ells/L)); y=np.array([ents[str(e)] for e in ells])
            slope,intercept=np.polyfit(x,y,1)
            fits[str(cut)]={'c':float(3*slope),'rms':float(np.sqrt(np.mean((y-slope*x-intercept)**2)))}
        # Uniform density and XX bond expectations are obtained directly by HF.
        exch=float(-psi@(mats[1]@psi)); kin=float(-psi@(mats[0]@psi))
        # Exact finite free charge sea and squeezed XX ground energy at first order.
        echarge=-2/np.sin(np.pi/L)
        alphaL=.5-2/(L*L*np.sin(np.pi/L)**2)
        espin=-2/np.sin(np.pi/(L//2))
        leading=echarge+J*alphaL*espin
        out['runs'].append({'J':J,'energies':ev.tolist(),'residuals':res,'kinetic_expectation':kin,
            'exchange_expectation':exch,'entropy':ents,'entropy_fits':fits,
            'first_order_energy':leading,'first_order_alpha':alphaL,'free_charge_energy':echarge,
            'squeezed_xx_energy':espin,'energy_remainder_over_J2':float((ev[0]-leading)/J**2),
            'contact_exchange_ratio':float(exch/espin),'elapsed_seconds':time.time()-started})
        (HERE/f'result_L{L}.json').write_text(json.dumps(out,indent=2)+'\n')
        print(json.dumps({'L':L,'J':J,'dim':len(codes),'E0':ev[0],'c_fits':fits,'remainder/J2':(ev[0]-leading)/J**2,'seconds':time.time()-started}),flush=True)
    return out

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--L',type=int,default=12);ap.add_argument('--J',type=float,nargs='+',default=[.2,.1,.05]);a=ap.parse_args()
    run(a.L,a.J)
