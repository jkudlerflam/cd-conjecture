"""C25 exact tagged-vacancy / C18 Laplacian composition. numpy/scipy only.
Run with python3 -B composition.py. All outputs stay beside this script.
Read-only import of C18's finite rule table, bytecode disabled.
"""
import sys
sys.dont_write_bytecode = True
import importlib.util
import itertools
import json
import time
from pathlib import Path
import numpy as np
from scipy.linalg import eigh, cho_factor, cho_solve

ROOT=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('c18clock',ROOT/'../work_C18/clock_checks.py')
clock=importlib.util.module_from_spec(spec); spec.loader.exec_module(clock)
C=30; K=4

def clock_data():
    words=[('PL','I','Z','PR')]; prev=None
    while True:
        cand=[x for x in clock.neighbors(words[-1]) if x!=prev]
        if not cand: break
        assert len(cand)==1
        prev=words[-1]; words.append(cand[0])
    assert len(words)==C
    parts=np.zeros((3,C,C)); parity=np.zeros(C,dtype=int); swaps=[]
    for t,(a,b) in enumerate(zip(words[:-1],words[1:])):
        matches=[]
        for src,dst in clock.RULES:
            width=len(src)
            for pos in range(K-width+1):
                if a[pos:pos+width] not in (src,dst): continue
                other=dst if a[pos:pos+width]==src else src
                if a[:pos]+other+a[pos+width:]==b:
                    matches.append((pos,width))
        assert len(matches)==1
        pos,width=matches[0]
        coeff=np.zeros(3)
        if width==2: coeff[pos]=1
        else:
            assert pos in (1,2)
            coeff[pos-1:pos+1]=.5
        lap=np.zeros((C,C)); lap[t,t]=lap[t+1,t+1]=1
        lap[t,t+1]=lap[t+1,t]=-1
        parts+=coeff[:,None,None]*lap
        swap=a[1].startswith('D') and b[1].startswith('A')
        parity[t+1]=parity[t]^swap
        if swap: swaps.append(t)
    assert swaps==[4,10,18,24]
    assert np.linalg.norm(parts.sum(0)@np.ones(C))==0
    return words,parts,parity

WORDS,PARTS,PARITY=clock_data()
PHI=np.ones(C)/np.sqrt(C)
# i=0 gives both work bits zero; i=1 starts with its occupied bit at left.
POT={}
for side in ('L','R'):
    POT[side,0]=np.full(C,1.25)
    bit=(1-PARITY) if side=='L' else PARITY
    POT[side,1]=1.25-.5*bit


def charge(p,w,U=1.):
    N=K+p; configs=list(itertools.combinations(range(N),p))
    index={s:i for i,s in enumerate(configs)}; d=len(configs)
    H=np.zeros((d,d)); D=np.zeros((d,3)); occ={x:np.zeros(d) for x in ('L','R')}
    for a,s in enumerate(configs):
        vac=set(s); bits=[int(j in vac) for j in range(N)]
        H[a,a]=U*sum(bits[j]!=bits[j+1] for j in range(N-1))
        occupied=[j for j in range(N) if j not in vac]
        D[a]=[occupied[j+1]==occupied[j]+1 for j in range(3)]
        occ['L'][a]=0 not in vac; occ['R'][a]=N-1 not in vac
        for j in range(N-1):
            if bits[j]!=bits[j+1]:
                new=tuple(sorted(vac.symmetric_difference({j,j+1})))
                H[index[new],a]=-w
    # Diagonalize reflection sectors separately: the very small even/odd
    # splitting otherwise contaminates the ground vector at small w.
    reflect=[index[tuple(sorted(N-1-j for j in x))] for x in configs]
    even=[]; odd=[]; seen=set()
    for j,jr in enumerate(reflect):
        if j in seen: continue
        seen.update((j,jr)); v=np.zeros(d)
        if j==jr: v[j]=1; even.append(v)
        else:
            v[j]=v[jr]=1/np.sqrt(2); even.append(v)
            v=v.copy(); v[jr]*=-1; odd.append(v)
    all_e=[]; all_v=[]
    for basis in (even,odd):
        B=np.array(basis).T; ee,VV=eigh(B.T@H@B)
        all_e.extend(ee); all_v.extend((B@VV).T)
    order=np.argsort(all_e); ev=np.array(all_e)[order]
    V=np.array(all_v).T[:,order]; chi=V[:,0]
    if chi.sum()<0: chi=-chi; V[:,0]=chi
    assert min(chi)>0
    # Symmetric matrix and nondegenerate ground state imply reflection equality.
    r=float(np.dot(chi**2,occ['L']))
    assert abs(r-np.dot(chi**2,occ['R']))<2e-9, (p,w,r,np.dot(chi**2,occ['R']),ev[:3],np.linalg.norm(H-H.T))
    return dict(p=p,N=N,d=d,configs=configs,H=H-ev[0]*np.eye(d),chi=chi,
                gap=ev[1]-ev[0],charge_ev=ev-ev[0],charge_vec=V,D=D,occ=occ,r=r,
                cut_mass=float(np.dot(chi**2,np.any(D==0,axis=1))),
                avgD=chi**2@D)


def core_matrix(ch,side,i,s):
    # Exact conserved Bell block, initialized gauged input (i,0).
    H=np.kron(ch['H'],np.eye(C))
    for j in range(3): H+=s*np.kron(np.diag(ch['D'][:,j]),PARTS[j])
    edge=0 if side=='L' else 2
    H+=s*np.diag((ch['D'][:,edge,None]*POT[side,i]).ravel())
    return H


def compensator(r1,r2,phi,kind='positive'):
    plus=lambda theta: np.array([1,np.exp(1j*theta)])/np.sqrt(2)
    u1=plus(-4*phi); u2=plus(-5*phi)
    rho1=np.outer(u1,u1.conj()); rho2=np.outer(u2,u2.conj())
    x=np.sin(phi/2)**2; m1,m2=1/r1,1/r2
    if kind=='positive': M=(m1*(np.eye(2)-rho2)+m2*(np.eye(2)-rho1))/x
    else: M=(m1+m2)/2*np.eye(2)+(m1-m2)/(2*x)*(rho1-rho2)
    for u,m in ((u1,m1),(u2,m2)):
        assert abs(np.vdot(u,M@u)-m)<1e-10
    return M


def helix_measure(N,M,phi):
    # Left endpoint, B track pinned at right. Ground phase = -(N-1)phi.
    theta=-(N-1)*phi
    plus=np.array([1,np.exp(1j*theta)])/np.sqrt(2)
    minus=np.array([1,-np.exp(1j*theta)])/np.sqrt(2)
    m=np.vdot(plus,M@plus); a=np.vdot(minus,M@plus)
    H=np.zeros((N,N))
    for j in range(N-1):
        H[j,j]+=.5; H[j+1,j+1]+=.5; H[j,j+1]-=.5; H[j+1,j]-=.5
    H[-1,-1]+=1
    e,V=eigh(H)
    weight=abs(a)**2*abs(V[0,:])**2
    return np.r_[0,e],np.r_[abs(m)**2,weight],float(abs(a)**2)


def measure(ch,side,i,s,hm):
    # Stable block diagonalization. The charge splitting can be O(w^8),
    # so raw full-matrix eigh would lose the O(s) eigenvalues.
    d=ch['d']; T=ch['charge_vec']; dim=d*C
    W=np.zeros((dim,dim))
    for j in range(3):
        W+=np.kron(T.T@(ch['D'][:,j,None]*T),PARTS[j])
    edge=0 if side=='L' else 2
    W+=np.kron(T.T@(ch['D'][:,edge,None]*T),np.diag(POT[side,i]))
    A=W[:C,:C]; B=W[C:,:C]
    D=W[C:,C:]+np.diag(np.repeat(ch['charge_ev'][1:]/s,C))
    cf=cho_factor(D,check_finite=False)
    X=-cho_solve(cf,B,check_finite=False)
    for it in range(200):
        Xnew=cho_solve(cf,-B+X@(A+B.T@X),check_finite=False)
        if np.linalg.norm(Xnew-X)<1e-14:
            X=Xnew; break
        X=Xnew
    else: raise RuntimeError('Riccati iteration failed')
    def invsqrt(a):
        e,v=eigh(a,check_finite=False)
        return (v/np.sqrt(e))@v.T
    N=invsqrt(np.eye(C)+X.T@X)
    # The high normalization differs from identity only on rank <= C.
    uu,ss,vh=np.linalg.svd(X,full_matrices=False)
    Nh=np.eye(dim-C)+(uu*(1/np.sqrt(1+ss**2)-1))@uu.T
    L=N@(A+B.T@X+X.T@B+X.T@D@X)@N
    HH=Nh@(D-B@X.T-X@B.T+X@A@X.T)@Nh
    el,vl=eigh(L,check_finite=False); eh,vh=eigh(HH,check_finite=False)
    residual=np.linalg.norm(B+D@X-X@(A+B.T@X))
    assert residual<2e-11, residual
    src=np.kron(T.T@(ch['occ'][side]*ch['chi']),PHI)
    src_low=N@(src[:C]+X.T@src[C:])
    src_high=Nh@(src[C:]-X@src[:C])
    e=s*np.r_[el,eh]
    weights=np.r_[abs(vl.T@src_low)**2,abs(vh.T@src_high)**2]
    assert min(e)>0
    assert abs(weights.sum()-np.dot(src,src))<1e-10
    he,ha,_=hm
    return (e[:,None]+he).ravel(),(weights[:,None]*ha).ravel()


def contraction(mu,nu):
    e,a=mu; f,b=nu
    val=0.
    for lo in range(0,len(e),96):
        val+=np.dot(a[lo:lo+96], (1/(e[lo:lo+96,None]+f[None,:]))@b)
    return -3*val


def coefficients(diag):
    return (np.array([[1,1,1,1],[1,1,-1,-1],[1,-1,1,-1],[1,-1,-1,1]])@diag/4).tolist()


def reference(ch=None):
    ms={}
    for side in ('L','R'):
        for i in (0,1):
            weights=np.ones(3) if ch is None else ch['avgD']
            H=np.einsum('j,jab->ab',weights,PARTS)
            H+=np.diag(weights[0 if side=='L' else 2]*POT[side,i])
            e,V=eigh(H); ms[side,i]=(e,abs(V.T@PHI)**2)
    return ms


def diag_from(ms,ns):
    return np.array([contraction(ms['R',i],ns['L',j]) for i in (0,1) for j in (0,1)])


def main():
    start=time.monotonic()
    bare=reference(); baseline=diag_from(bare,bare)
    out=dict(status='numerical',k=K,C=C,physical_lengths=[5,6],
             baseline_diagonal=baseline.tolist(),baseline_coefficients=coefficients(baseline),runs=[])
    for w in (.16,.10,.06):
        charges={p:charge(p,w) for p in (1,2)}
        gap=min(x['gap'] for x in charges.values())
        avg={p:reference(ch) for p,ch in charges.items()}
        avdiag={str(p)+str(q):diag_from(avg[p],avg[q]) for p in (1,2) for q in (1,2)}
        for phi,kind in ((np.pi/2,'positive'),(.2,'positive'),(.2,'small')):
            M=compensator(charges[1]['r'],charges[2]['r'],phi,kind)
            hms={p:helix_measure(ch['N'],M,phi) for p,ch in charges.items()}
            for ratio in (.1,.01,.001):
                s=ratio*gap
                ms={p:{(side,i):measure(ch,side,i,s,hms[p])
                       for side in ('L','R') for i in (0,1)} for p,ch in charges.items()}
                ds={str(p)+str(q):s*diag_from(ms[p],ms[q]) for p in (1,2) for q in (1,2)}
                # Compare logical (traceless) norm to unwrapped logical norm.
                center=lambda x:x-np.mean(x)
                norm=max(abs(center(baseline)))
                mismatch=max(max(abs(center(ds[key]-ds['11']))) for key in ds)/norm
                total=max(max(abs(ds[key]-ds['11'])) for key in ds)
                limit_error=max(max(abs(center(ds[key]-avdiag[key]))) for key in ds)/norm
                target_error=max(max(abs(center(ds[key]-baseline))) for key in ds)/norm
                rec=dict(w_over_U=w,phi=phi,compensator=kind,source_over_charge_gap=ratio,
                         s=s,charge_gap=gap,M_norm=float(np.linalg.norm(M,2)),
                         source_dimensions={str(p):len(ms[p]['L',0][0]) for p in (1,2)},
                         r=[charges[p]['r'] for p in (1,2)],
                         cut_mass=[charges[p]['cut_mass'] for p in (1,2)],
                         avgD=[charges[p]['avgD'].tolist() for p in (1,2)],
                         scaled_diagonals={key:val.tolist() for key,val in ds.items()},
                         scaled_coefficients={key:coefficients(val) for key,val in ds.items()},
                         logical_length_mismatch=mismatch,absolute_length_mismatch=total,
                         logical_limit_error=limit_error,logical_target_error=target_error,
                         averaged_diagonals={key:val.tolist() for key,val in avdiag.items()})
                out['runs'].append(rec)
                print(w,round(phi,4),kind,ratio,'mismatch',mismatch,'limit_error',limit_error,flush=True)
                (ROOT/'composition_results.json').write_text(json.dumps(out,indent=2)+'\n')
    out['seconds']=time.monotonic()-start
    (ROOT/'composition_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print('seconds',out['seconds'],flush=True)

if __name__=='__main__': main()
