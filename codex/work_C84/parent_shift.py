# C84 local source copied from C80. See PROOF.md for changes.
"""Bitless-token bonds with 14, 13 or 12 letters. C80, all moves literal."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import itertools
from functools import lru_cache
import numpy as np
import scipy.sparse as sp
import bond18 as parent

BASE=([('L',a) for a in (0,1)]+[('R',a) for a in (0,1)]
      +[('d',s) for s in 'vAB']+[('x',s) for s in 'AB']
      +[('e',a) for a in (0,1)]+[('m',s) for s in 'vAB'])

def merge(l,size):
    if size<=13 and l==('m','v'):return ('e',0)
    if size<=12 and l==('d','v'):return ('R',0)
    return l

def alphabet(size):return list(dict.fromkeys(merge(l,size) for l in BASE))
def charge(l):
    if l[0] in ['L','R','e']:return 2*l[1]
    return {'v':0,'A':-1,'B':1}[l[1]]
def excited(l):return l[0] in ['d','x','m'] and l[1] in ['A','B']
@lru_cache(None)
def allowed(l,r,size):
    # The union of the original ordered-pair grammar under identification.
    ar={'L':('L','L'),'R':('R','R'),'d':('L','R'),'x':('L','R'),'e':('R','L'),'m':('L','L')}
    return any(ar[a[0]][1]==ar[b[0]][0] for a in BASE for b in BASE
               if merge(a,size)==l and merge(b,size)==r)

def parts(p,size=12):
    ls=alphabet(size);ix={l:i for i,l in enumerate(ls)};d=len(ls)
    h=np.zeros((d*d,d*d),complex);o=np.zeros((d,d),complex);v=np.zeros_like(h)
    seen={}
    def add(out,inn,val,port=False):
        j=ix[merge(out[0],size)]*d+ix[merge(out[1],size)]
        i=ix[merge(inn[0],size)]*d+ix[merge(inn[1],size)]
        key=(port,j,i)
        if key in seen:
            assert abs(seen[key]-val)<1e-12
            return
        seen[key]=val
        mat=v if port else h
        mat[j,i]+=val;mat[i,j]+=np.conj(val)
    for a in (0,1):
        for s in 'vAB':
            z=(-1)**a if s=='A' else 1
            add((('L',a),('d',s)),(('d',s),('R',a)),-z)
            if s=='v':add((('L',a),('m','v')),(('d','v'),('e',a)),-1)
        for s in 'AB':
            z=(-1)**a if s=='A' else 1
            phi=p['flux'] if s=='A' else 0
            psi=p['n']*p['flux']+np.pi/2 if s=='A' else np.pi/2
            add((('L',a),('x',s)),(('x',s),('R',a)),-np.exp(1j*phi)*z)
            add((('L',a),('m',s)),(('x',s),('e',a)),-np.exp(1j*phi)*z)
            add((('L',a),('m',s)),(('d',s),('e',a)),-p['td']*np.exp(1j*psi)*z)
    for m in [l for l in BASE if l[0] in ['e','m']]:
        for s in 'AB':add((m,('x',s)),(m,('d',s)),-1)
    for k in ['R','e']:
        add((('d','B'),(k,0)),(('d','A'),(k,1)),p['g'])
    for s,t,amp in [('A','B',np.exp(1j*p['theta'])),('B','A',np.exp(-1j*p['theta']))]:
        add((('m',s),('d',t)),(('m','v'),('d','v')),amp,True)
    for i,l in enumerate(ls):
        o[i,i]=p.get('musel',0)*(l[0] in ['e','m'])
        if excited(l):o[i,i]+=p['Delta']+p['D']*(l[0] in ['x','m'])
        # Field on the physical bit letters. There are exactly n bits in legal words.
        if l[0] in ['L','R','e']:o[i,i]+=p.get('field',0)*l[1]
    pen=np.diag([p['K']*(not allowed(a,b,size)) for a in ls for b in ls])
    return h,pen,o,v

def bond(p,size=12,lam=1):
    h,pen,o,v=parts(p,size);d=len(o)
    return h+pen+lam*v+(np.kron(o,np.eye(d))+np.kron(np.eye(d),o))/2

def words(n,size=12):
    out=[]
    for w in parent.words(n):
        bits=[l[2] for l in w if l[2] is not None]
        pos=next(i for i,l in enumerate(w) if l[0] in ['d','x','m'])
        kind,s,_=w[pos]
        if pos==n:
            nw=[('L',a) for a in bits]+[('m',s)]
        else:
            nw=[('L',a) for a in bits[:pos]]+[(kind,s)]+[('R',a) for a in bits[pos:-1]]+[('e',bits[-1])]
        out.append(tuple(merge(l,size) for l in nw))
    assert len(set(out))==len(out)
    return out

def cluster(p,size=12,closure=True):
    n=p['n'];ls=alphabet(size);ix={l:i for i,l in enumerate(ls)};d=len(ls);L=2*(n+1)
    hk,pen,o,v=parts(p,size);hk+=pen
    ws=words(n,size);states=[tuple(ix[l] for l in x+y) for x in ws for y in ws]
    idx={w:i for i,w in enumerate(states)};legal_dim=len(states)
    fixed=ix[merge(('m','v'),size)]
    hc={i:[(j,hk[j,i]) for j in np.flatnonzero(hk[:,i])] for i in range(d*d)}
    vc={i:[(j,v[j,i]) for j in np.flatnonzero(v[:,i])] for i in range(d*d)}
    oc={i:[(j,o[j,i]) for j in np.flatnonzero(o[:,i])] for i in range(d)}
    ts=[];vs=[];cursor=0
    while cursor<len(states):
        w=states[cursor]
        def emit(z,val,port=False):
            z=tuple(z)
            if z not in idx:
                if not closure:return
                idx[z]=len(states);states.append(z)
            (vs if port else ts).append((idx[z],cursor,val))
        for pos in range(L):
            for row,val in oc[w[pos]]:
                z=list(w);z[pos]=row;emit(z,val)
        for pos in range(-1,L-1):
            a=fixed if pos==-1 else w[pos];b=w[pos+1]
            for row,val in hc[a*d+b]:
                z=list(w)
                if pos==-1:assert row//d==fixed
                else:z[pos]=row//d
                z[pos+1]=row%d;emit(z,val)
        for row,val in vc[w[n]*d+w[n+1]]:
            z=list(w);z[n]=row//d;z[n+1]=row%d;emit(z,val,True)
        cursor+=1
        if len(states)>100000:raise RuntimeError('local dimension cap')
    def matrix(items):
        r,c,x=zip(*items)
        return sp.coo_matrix((x,(r,c)),shape=(len(states),len(states))).tocsc()
    H=matrix(ts);V=matrix(vs)
    ex=np.array([sum(excited(ls[l]) for l in w) for w in states])
    bad=np.array([sum(not allowed(ls[a],ls[b],size) for a,b in zip((fixed,)+w,w)) for w in states])
    return dict(HH=H,V=V,exc=ex,bad=bad,states=states,legal_dim=legal_dim)
