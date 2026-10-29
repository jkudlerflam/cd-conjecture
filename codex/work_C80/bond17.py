# C80 local copy of C75 source, used as a reference or scalar certificate.
"""Merge empty and vacuum marker. Illegal port excursions are retained."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import itertools
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import splu,eigsh
import bond18 as b18

def merge(l):return ('e','',None) if l==('m','v',None) else l
LET=list(dict.fromkeys(merge(l) for l in b18.LET));IX={l:i for i,l in enumerate(LET)};d=len(LET);assert d==17
def allowed(l,r):
    # The merged marker accepts either incoming arrow and has outgoing L.
    if r[0]=='e':return True
    return b18.AR[l[0]][1]==b18.AR[r[0]][0]
def parts(p):
    hk,_,oo,vv=b18.parts(p)
    h=np.zeros((d*d,d*d),complex);o=np.zeros((d,d),complex);v=np.zeros_like(h)
    # Duplicate identical rules, notably switches controlled by either old
    # empty or vacuum marker, become a SINGLE rule after identification.
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
    """Full invariant closure of two legal cells, with a frozen left marker.
    Build from literal local rules, including all states reached by V and H0.
    Only this proven invariant component is diagonalized.
    """
    n=p['n'];L=2*(n+1);hk,pen,o,v=parts(p);hk+=pen
    Ws=b18.words(n);base=[tuple(IX[merge(l)] for l in x+y) for x in Ws for y in Ws]
    states=list(dict.fromkeys(base));ix={w:i for i,w in enumerate(states)}
    transitions=[];vtrans=[]
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
        if len(states)>150000:raise RuntimeError('local computation cap')
    def matrix(ts):
        rows,cols,vals=zip(*ts);return sp.coo_matrix((vals,(rows,cols)),shape=(len(states),len(states))).tocsc()
    H=matrix(transitions);V=matrix(vtrans)
    h18,g18,_,_,_=b18.cell(p);nn=len(Ws)
    G0=np.kron(g18,g18);perm=[a*2**n+b for b in range(2**n) for a in range(2**n)];G0=G0[:,perm]
    G=np.zeros((len(states),4**n),complex)
    for i,w in enumerate(base):G[ix[w]]=G0[i]
    E0=p['musel']-2*p['t']*np.cos(np.pi/(n+2))
    # V creates two excited species. H0 preserves their number, so this
    # subspace contains the ENTIRE second-order resolvent source space.
    ex=[sum(LET[l][1] in ['A','B'] for l in w) for w in states]
    sel=np.flatnonzero(np.array(ex)==2)
    A=H[sel][:,sel]-2*E0*sp.eye(len(sel),format='csc');S=(V@G)[sel]
    H2=-S.conj().T@splu(A.tocsc()).solve(S)
    bad=np.array([sum(not allowed(LET[a],LET[b]) for a,b in zip((fixed,)+w,w)) for w in states])
    assert np.max(abs(H@G-2*E0*G))<1e-10
    return dict(HH=H,V=V,G=G,H2=H2,E0=E0,bad=bad,states=states,legal_dim=len(base),source_dim=len(sel))
