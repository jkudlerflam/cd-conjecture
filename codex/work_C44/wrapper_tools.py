# C44 snapshot of ../work_C32/mobile.py. The C32 core constructor is
# retained for provenance but is not used. C44 supplies C39 core matrices.
"""C32 positive folded bus and ordered-coordinate mobile wrapper.

Only numpy/scipy. All output is produced by check_mobile.py in this folder.
The literal core alphabet has 28 states. Source calculations retain all
legal core and vacancy excitations; neutral filter factors stay in their
exact ground states. See RESULT.md for the scope of the mobile ED.
"""
import itertools
import numpy as np
from scipy.linalg import eigh, cho_factor, cho_solve

SP = ('v', '0', '+', '-')
EX = ('0', '+', '-')
CH = {'v': 0, '0': 0, '+': 1, '-': -1}
LET = ([('L', None, a) for a in (0, 1)]
       + [('R', None, a) for a in (0, 1)]
       + [('d', s, a) for s in SP for a in (0, 1)]
       + [('x', s, a) for s in EX for a in (0, 1)]
       + [(end, None, None) for end in ('A', 'B')]
       + [(end, s, None) for end in ('A', 'B') for s in SP])
assert len(LET) == 28
FMAP = {'0': {'v': '0', '0': 'v'},
        '+': {'v': '+', '-': 'v'}, '-': {'v': '-', '+': 'v'}}
TERMS = [('0', '0', 1.), ('+', '-', 1j), ('-', '+', -1j)]


def core(n, g=.12, mass=.6, express=3., loop=.003, phase=.4):
    """Each hopping is its positive two-state edge matrix, with original
    off-diagonal amplitude. Onsite mass plus conversion is positive.
    The local first/last penalties selecting A/B are zero in this sector.
    """
    k = n+2
    states = []
    positions = [('A', 0)] + [('d', j+1) for j in range(n)]
    positions += [('x', j+1) for j in range(n)] + [('B', k-1)]
    for kind, pos in positions:
        for s in (EX if kind == 'x' else SP):
            for bits in itertools.product((0, 1), repeat=n):
                word = [('A', s if kind == 'A' else None, None)]
                for j, a in enumerate(bits, 1):
                    if kind in ('d', 'x') and j == pos:
                        word.append((kind, s, a))
                    else:
                        word.append(('L' if j < pos else 'R', None, a))
                word.append(('B', s if kind == 'B' else None, None))
                states.append(tuple(word))
    idx = {w:i for i,w in enumerate(states)}
    d = len(states)
    parts = np.zeros((k-1, d, d), complex)
    ends = np.zeros((d,d), complex)
    rules = []
    def edge(a,b,c): rules.append((a,b,c))
    for a,b in itertools.product((0,1), repeat=2):
        for s in SP:
            z = (-1)**a if s == '0' else 1
            edge((('d',s,a),('R',None,b)),
                 (('L',None,a),('d',s,b)), -z)
        for s in EX:
            z = (-1)**a if s == '0' else 1
            edge((('x',s,a),('R',None,b)),
                 (('L',None,a),('x',s,b)),
                 -z*np.exp(1j*phase*(s=='0')))
    for a in (0,1):
        for s in SP:
            edge((('A',s,None),('R',None,a)),
                 (('A',None,None),('d',s,a)), -1.)
        edge((('d','v',a),('B',None,None)),
             (('L',None,a),('B','v',None)), -1.)
        for s in EX:
            z = (-1)**a if s == '0' else 1
            edge((('x',s,a),('B',None,None)),
                 (('L',None,a),('B',s,None)),
                 -z*np.exp(1j*phase*(s=='0')))
            if loop:
                edge((('d',s,a),('B',None,None)),
                     (('L',None,a),('B',s,None)),
                     -loop*z*np.exp(.9j*(s=='0')))
            for sA in (None,)+SP:
                edge((('A',sA,None),('d',s,a)),
                     (('A',sA,None),('x',s,a)), -1.)
    for i,w in enumerate(states):
        for p in range(k-1):
            for aa,bb,c in rules:
                if w[p:p+2] == aa:
                    w2=w[:p]+bb+w[p+2:]
                    if w2 not in idx: continue
                    j=idx[w2]
                    parts[p,i,i]+=abs(c); parts[p,j,j]+=abs(c)
                    parts[p,j,i]+=c; parts[p,i,j]+=np.conj(c)
        # onsite terms distributed as half per incident bond, including
        # the positive endpoint halves in a_pay and b_pay.
        for p,(kind,s,a) in enumerate(w):
            if s not in EX: continue
            val=mass+express*(kind in ('A','B','x'))
            entries=[(i,val)]
            if kind=='d' and s=='0':
                s2,a2=('+',0) if a else ('-',1)
                ww=w[:p]+(('d',s2,a2),)+w[p+1:]
                entries.append((idx[ww],g))
            for j,v in entries:
                onsite=np.zeros((d,d),complex)
                onsite[j,i]=v
                if j!=i: onsite[i,j]=np.conj(v)
                if p>0: parts[p-1]+=onsite/2
                else: ends+=onsite/2
                if p<k-1: parts[p]+=onsite/2
                else: ends+=onsite/2
    H=parts.sum(axis=0)+ends
    Q=np.array([sum((a or 0)+CH.get(s,0) for _,s,a in w) for w in states])
    ex=np.array([any(s in EX for _,s,_ in w) for w in states])
    G=np.zeros((d,2**n),complex)
    for i,w in enumerate(states):
        if not ex[i]:
            bits=[a for _,_,a in w if a is not None]
            col=sum(a*2**(n-1-j) for j,a in enumerate(bits))
            G[i,col]=1/np.sqrt(n+2)
    F={}
    for end,p in (('L',0),('R',k-1)):
        F[end]={}
        for f,mp in FMAP.items():
            ff=np.zeros((d,d),complex)
            for i,w in enumerate(states):
                kind,s,a=w[p]
                if s in mp:
                    ww=w[:p]+((kind,mp[s],a),)+w[p+1:]
                    ff[idx[ww],i]=1
            F[end][f]=ff
    return dict(n=n,k=k,d=d,states=states,H=H,parts=parts,ends=ends,
                Q=Q,ex=ex,G=G,F=F)


def charge(k,p,w,U=1.):
    N=k+p; cfg=list(itertools.combinations(range(N),p)); ix={a:i for i,a in enumerate(cfg)}
    d=len(cfg); H=np.zeros((d,d)); D=np.zeros((d,k-1)); occ=np.zeros((2,d))
    for i,s in enumerate(cfg):
        vac=set(s); bits=[j in vac for j in range(N)]
        H[i,i]=U*sum(bits[j]!=bits[j+1] for j in range(N-1))
        pos=[j for j in range(N) if j not in vac]
        D[i]=[pos[j+1]==pos[j]+1 for j in range(k-1)]
        occ[:,i]=[0 not in vac,N-1 not in vac]
        for j in range(N-1):
            if bits[j]!=bits[j+1]:
                H[ix[tuple(sorted(vac.symmetric_difference({j,j+1})))],i]=-w
    # Reflection resolution prevents loss of the tiny symmetric splitting.
    seen=set(); bases=[[],[]]
    for i,s in enumerate(cfg):
        if i in seen: continue
        j=ix[tuple(sorted(N-1-x for x in s))]; seen.update((i,j))
        v=np.zeros(d)
        if i==j: v[i]=1; bases[0].append(v)
        else:
            v[i]=v[j]=1/np.sqrt(2); bases[0].append(v)
            v=v.copy(); v[j]*=-1; bases[1].append(v)
    ee=[]; vv=[]
    for bb in bases:
        if not bb: continue
        B=np.array(bb).T; e,v=eigh(B.T@H@B)
        ee.extend(e); vv.extend((B@v).T)
    order=np.argsort(ee); ev=np.array(ee)[order]; V=np.array(vv).T[:,order]
    if V[:,0].sum()<0: V[:,0]*=-1
    chi=V[:,0]; r=float(chi**2@occ[0])
    return dict(k=k,p=p,N=N,d=d,cfg=cfg,H=H,ev=ev-ev[0],V=V,chi=chi,
                D=D,occ=occ,r=r,gap=ev[1]-ev[0],avgD=chi**2@D)


def measure(c,ch=None,s=1.):
    """Exact source eigensystem of H_wrapped/s, stabilized by a low graph.
    Returns all source amplitudes, not just their scalar squared weights.
    """
    ex=c['ex']; de=int(ex.sum()); G=c['G']
    if ch is None:
        e,V=eigh(c['H'][np.ix_(ex,ex)])
        return e,{side:{f:V.conj().T@(A@G)[ex] for f,A in c['F'][side].items()}
                  for side in ('L','R')},0.
    T=ch['V']; d=ch['d']; dim=d*de
    W=np.kron(np.eye(d),c['ends'][np.ix_(ex,ex)])
    for j,B in enumerate(c['parts']):
        W+=np.kron(T.T@(ch['D'][:,j,None]*T),B[np.ix_(ex,ex)])
    A=W[:de,:de]; B=W[de:,:de]
    D=W[de:,de:]+np.diag(np.repeat(ch['ev'][1:]/s,de))
    cf=cho_factor(D,check_finite=False)
    X=-cho_solve(cf,B,check_finite=False)
    for it in range(300):
        Xn=cho_solve(cf,-B+X@(A+B.conj().T@X),check_finite=False)
        if np.linalg.norm(Xn-X)<1e-13: X=Xn; break
        X=Xn
    else: raise RuntimeError('graph iteration')
    e,v=eigh(np.eye(de)+X.conj().T@X)
    N=(v/np.sqrt(e))@v.conj().T
    u,sv,_=np.linalg.svd(X,full_matrices=False)
    Nh=np.eye(dim-de)+(u*(1/np.sqrt(1+sv**2)-1))@u.conj().T
    L=N@(A+B.conj().T@X+X.conj().T@B+X.conj().T@D@X)@N
    HH=Nh@(D-B@X.conj().T-X@B.conj().T+X@A@X.conj().T)@Nh
    el,vl=eigh(L,check_finite=False); eh,vh=eigh(HH,check_finite=False)
    res=np.linalg.norm(B+D@X-X@(A+B.conj().T@X))
    out={}
    for side,no in (('L',0),('R',1)):
        out[side]={}
        for f,F in c['F'][side].items():
            src=np.kron((T.T@(ch['occ'][no]*ch['chi']))[:,None],(F@G)[ex])
            low=N@(src[:de]+X.conj().T@src[de:])
            high=Nh@(src[de:]-X@src[:de])
            out[side][f]=np.vstack((vl.conj().T@low,vh.conj().T@high))
    return np.r_[el,eh],out,float(res)


def compensator(r1,r2,b,phi):
    plus=lambda x:np.array([1,np.exp(1j*x)])/np.sqrt(2)
    u,v=plus(-(b-1)*phi),plus(-b*phi)
    chi=np.sin(phi/2)**2; m1,m2=1/r1,1/r2
    return (m1+m2)/2*np.eye(2)+(m1-m2)/(2*chi)*(np.outer(u,u.conj())-np.outer(v,v.conj()))


def helix(N,M,phi,s):
    u=np.array([1,np.exp(-1j*(N-1)*phi)])/np.sqrt(2)
    v=u*np.array([1,-1]); m=np.vdot(u,M@u); a=np.vdot(v,M@u)
    H=np.zeros((N,N))
    for j in range(N-1):
        H[j,j]+=.5; H[j+1,j+1]+=.5; H[j,j+1]-=.5; H[j+1,j]-=.5
    H[-1,-1]+=1
    e,V=eigh(H)
    return np.r_[0,e/s],np.r_[abs(m)**2,abs(a)**2*V[0]**2]


def bond(me1,me2,hm1=None,hm2=None):
    e,aa,_=me1; f,bb,_=me2; A=aa['R']; B=bb['L']
    C=np.zeros((len(e),len(f)))
    if hm1 is None: hm1=(np.array([0.]),np.array([1.]))
    if hm2 is None: hm2=(np.array([0.]),np.array([1.]))
    for h,a in zip(*hm1):
        for hh,b in zip(*hm2): C+=a*b/(e[:,None]+f[None,:]+h+hh)
    N=A['0'].shape[1]; H=np.zeros((N*N,N*N),complex)
    for ar,al,v in TERMS:
        for br,bl,w in TERMS:
            X=np.einsum('ei,ej->eij',A[ar].conj(),A[br]).reshape(len(e),-1)
            Y=np.einsum('ei,ej->eij',B[al].conj(),B[bl]).reshape(len(f),-1)
            H-=np.conj(v)*w*(X.T@C@Y).reshape(N,N,N,N).transpose(0,2,1,3).reshape(N*N,N*N)
    return H


def transport(ch1,ch2):
    ix={x:i for i,x in enumerate(ch2['cfg'])}; k=ch1['k']
    return sum(ch1['chi'][i]*ch2['chi'][ix[(j,k+1)]] for i,(j,) in enumerate(ch1['cfg']))**2


def flat_model(c,ch):
    out=dict(c)
    out['H']=c['ends']+np.einsum('j,jab->ab',ch['avgD'],c['parts'])
    return out
