"""Eleven letters. p['n'] counts physical data sites, n-1 logical bits."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import itertools
import numpy as np
import scipy.sparse as sp
import parent_shift as old

def alphabet():return [l for l in old.alphabet(12) if l!=('e',1)]
def charge(l):return old.charge(l)
def excited(l):return old.excited(l)
def marker(l):return l[0] in ['e','m']
def allowed(l,r):
 return old.allowed(l,r,12) and not (l==('L',1) and marker(r))
def parts(p):
 ls=alphabet();orig=old.alphabet(12);idx=[orig.index(l) for l in ls]
 pairs=[12*i+j for i in idx for j in idx]
 h,_,o,v=old.parts(p,12)
 h=h[np.ix_(pairs,pairs)];v=v[np.ix_(pairs,pairs)];o=o[np.ix_(idx,idx)]
 pen=np.diag([p['K']*(not allowed(l,r)) for l in ls for r in ls])
 return h,pen,o,v

def bond(p,lam=1):
 h,pen,o,v=parts(p);d=len(o)
 return h+pen+lam*v+(np.kron(o,np.eye(d))+np.kron(np.eye(d),o))/2

def words(p):
 if p==0:return [(l,) for l in alphabet() if marker(l)]
 out=[]
 for bits in itertools.product((0,1),repeat=p-1):
  for j in range(p):
   for token in [('R',0),('d','A'),('d','B'),('x','A'),('x','B')]:
    out.append(tuple([('L',a) for a in bits[:j]]+[token]+[('R',a) for a in bits[j:]]+[('e',0)]))
  for m in [('e',0),('m','A'),('m','B')]:
   out.append(tuple([('L',a) for a in bits]+[('L',0),m]))
 assert len(set(out))==len(out)==(5*p+3)*2**(p-1)
 return out

def cluster(p,closure=True):
    n=p['n'];ls=alphabet();ix={l:i for i,l in enumerate(ls)};d=len(ls);L=2*(n+1)
    hk,pen,o,v=parts(p);hk+=pen
    ws=words(n);states=[tuple(ix[l] for l in x+y) for x in ws for y in ws]
    idx={w:i for i,w in enumerate(states)};legal_dim=len(states)
    fixed=ix['e',0]
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
    bad=np.array([sum(not allowed(ls[a],ls[b]) for a,b in zip((fixed,)+w,w)) for w in states])
    return dict(HH=H,V=V,exc=ex,bad=bad,states=states,legal_dim=legal_dim)
