"""15 letters: the vacuum walker is the L/R domain wall.
Finite penalties, including every off-grammar transition, are kept.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import splu
import bond18 as b18
import bond17 as b17

def merge(l):
    l=b17.merge(l)
    return ('R','',l[2]) if l[:2]==('d','v') else l
LET=list(dict.fromkeys(merge(l) for l in b18.LET));IX={l:i for i,l in enumerate(LET)};d=len(LET);assert d==15
def allowed(l,r):
    if l[0] in ['e','m']:return True
    if r[0]=='e':return True
    if r[0]=='m':return l[0]=='L'
    if l[0]=='L':return r[0] in ['L','R','d','x']
    if l[0] in ['R','d','x']:return r[0]=='R'
    raise ValueError((l,r))
def parts(p):
    hk,_,oo,vv=b18.parts(p)
    h=np.zeros((d*d,d*d),complex);o=np.zeros((d,d),complex);v=np.zeros_like(h)
    for src,dst in [(hk,h),(vv,v)]:
        seen={}
        for j,i in zip(*np.nonzero(src)):
            jj=IX[merge(b18.LET[j//18])]*d+IX[merge(b18.LET[j%18])]
            ii=IX[merge(b18.LET[i//18])]*d+IX[merge(b18.LET[i%18])]
            if (jj,ii) in seen:assert abs(seen[jj,ii]-src[j,i])<1e-14
            else:seen[jj,ii]=src[j,i]
        for (j,i),x in seen.items():dst[j,i]=x
    for j,i in zip(*np.nonzero(oo)):
        jj=IX[merge(b18.LET[j])];ii=IX[merge(b18.LET[i])]
        if o[jj,ii]!=0:assert o[jj,ii]==oo[j,i]
        o[jj,ii]=oo[j,i]
    pen=np.diag([p['K']*(not allowed(l,r)) for l in LET for r in LET])
    return h,pen,o,v
def bond(p,lam=1.):
    h,pen,o,v=parts(p);return h+pen+.5*(np.kron(o,np.eye(d))+np.kron(np.eye(d),o))+lam*v

def cluster(p):
    """Literal two-cell invariant closure. Full L/R vacuum code is used."""
    n=p['n'];L=2*(n+1);hk,pen,o,v=parts(p);hk+=pen
    Ws=b18.words(n);base=[tuple(IX[merge(l)] for l in x+y) for x in Ws for y in Ws]
    assert len(set(base))==len(base)
    states=list(base);ix={w:i for i,w in enumerate(states)};transitions=[];vtrans=[]
    fixed=IX['e','',None]
    hc={i:[(j,hk[j,i]) for j in np.flatnonzero(hk[:,i])] for i in range(d*d)}
    vc={i:[(j,v[j,i]) for j in np.flatnonzero(v[:,i])] for i in range(d*d)}
    oc={i:[(j,o[j,i]) for j in np.flatnonzero(o[:,i])] for i in range(d)}
    cursor=0
    while cursor<len(states):
        w=states[cursor]
        def emit(z,val,vs=False):
            z=tuple(z)
            if z not in ix:ix[z]=len(states);states.append(z)
            (vtrans if vs else transitions).append((ix[z],cursor,val))
        for pos in range(L):
            for row,val in oc[w[pos]]:
                z=list(w);z[pos]=row;emit(z,val)
        for pos in range(-1,L-1):
            a=fixed if pos==-1 else w[pos];bb=w[pos+1]
            for row,val in hc[a*d+bb]:
                z=list(w)
                if pos==-1:assert row//d==fixed
                else:z[pos]=row//d
                z[pos+1]=row%d;emit(z,val)
        pos=n
        for row,val in vc[w[pos]*d+w[pos+1]]:
            z=list(w);z[pos]=row//d;z[pos+1]=row%d;emit(z,val,True)
        cursor+=1
        if len(states)>100000:raise RuntimeError('local computation cap')
    def matrix(ts):
        rows,cols,vals=zip(*ts);return sp.coo_matrix((vals,(rows,cols)),shape=(len(states),len(states))).tocsc()
    H=matrix(transitions);V=matrix(vtrans)
    ex=np.array([sum(LET[l][1] in ['A','B'] for l in w) for w in states])
    # Each logical vacuum block has the SAME scalar orientation Hamiltonian.
    bits=np.array([sum((LET[w[cell*(n+1)+j]][2] or 0)<<(cell*n+j) for cell in range(2) for j in range(n)) for w in states])
    vacs=[np.flatnonzero((ex==0)&(bits==a)) for a in range(4**n)]
    # All no-excitation states have unoccupied markers. Their order may
    # differ between bits, so map orientations explicitly.
    def orient(w):return tuple(LET[l][0] for l in w)
    vacs=[np.array(sorted(s,key=lambda i:orient(states[i]))) for s in vacs]
    Hvac=H[vacs[0]][:,vacs[0]].toarray();ev,U=np.linalg.eigh(Hvac);e=ev[0]
    G=np.zeros((len(states),4**n),complex)
    for a,s in enumerate(vacs):
        assert np.max(abs(H[s][:,s].toarray()-Hvac))<1e-10
        G[s,a]=U[:,0]
    assert np.max(abs(H@G-e*G))<1e-9
    sel=np.flatnonzero(ex==2);A=H[sel][:,sel]-e*sp.eye(len(sel),format='csc');S=(V@G)[sel]
    H2=-S.conj().T@splu(A.tocsc()).solve(S)
    bad=np.array([sum(not allowed(LET[a],LET[b]) for a,b in zip((fixed,)+w,w)) for w in states])
    return dict(HH=H,V=V,G=G,H2=H2,E0=e/2,bad=bad,states=states,legal_dim=len(base),source_dim=len(sel),vac_gap=ev[1]-ev[0])
