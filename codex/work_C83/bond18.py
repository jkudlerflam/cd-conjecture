# C83 local copy from C80. Original attribution follows.
# C80 local copy of C75 source, used as a reference or scalar certificate.
"""Literal 18-letter bond and independent finite-cell assembly."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import itertools
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import splu,eigsh

LET=([('L','',a) for a in (0,1)]+[('R','',a) for a in (0,1)]
     +[('d',s,a) for s in 'vAB' for a in (0,1)]
     +[('x',s,a) for s in 'AB' for a in (0,1)]
     +[('e','',None)]+[('m',s,None) for s in 'vAB'])
IX={l:i for i,l in enumerate(LET)};d=len(LET);assert d==18
AR={'L':('L','L'),'R':('R','R'),'d':('L','R'),'x':('L','R'),'e':('R','L'),'m':('L','L')}
MARK=[l for l in LET if l[0] in ['e','m']]
def charge(l):return 2*(l[2] or 0)+{'':0,'v':0,'A':-1,'B':1}[l[1]]
def params(n=2,**kw):
    p=dict(n=n,t=1.,D=10.,Delta=4-2*np.cos(np.pi/(n+2)),g=.03,td=1e-5,
           flux=np.pi/(2*n),theta=-np.pi/4,K=20.,musel=0.,field=0.)
    p.update(kw);return p
def parts(p):
    h=np.zeros((d*d,d*d),complex);o=np.zeros((d,d),complex)
    def add(out,inn,v):
        j=IX[out[0]]*d+IX[out[1]];i=IX[inn[0]]*d+IX[inn[1]]
        h[j,i]+=v;h[i,j]+=np.conj(v)
    for a,b in itertools.product((0,1),repeat=2):
        for s in 'vAB':
            z=(-1)**a if s=='A' else 1
            add((('L','',a),('d',s,b)),(('d',s,a),('R','',b)),-p['t']*z)
        for s in 'AB':
            z=(-1)**a if s=='A' else 1
            ph=p['flux'] if s=='A' else 0.
            add((('L','',a),('x',s,b)),(('x',s,a),('R','',b)),-np.exp(1j*ph)*z)
    for a in (0,1):
        add((('L','',a),('m','v',None)),(('d','v',a),('e','',None)),-p['t'])
        for s in 'AB':
            z=(-1)**a if s=='A' else 1;ph=p['flux'] if s=='A' else 0.
            add((('L','',a),('m',s,None)),(('x',s,a),('e','',None)),-np.exp(1j*ph)*z)
            lp=p['n']*p['flux']+np.pi/2 if s=='A' else np.pi/2
            add((('L','',a),('m',s,None)),(('d',s,a),('e','',None)),-p['td']*np.exp(1j*lp)*z)
            for m in MARK:add((m,('x',s,a)),(m,('d',s,a)),-1.)
    for i,l in enumerate(LET):
        k,s,a=l;o[i,i]=p['field']*(a or 0)+p['musel']*(k in ['e','m'])
        if s in ['A','B']:o[i,i]+=p['Delta']+p['D']*(k in ['x','m'])
    i=IX['d','A',1];j=IX['d','B',0];o[i,j]=p['g'];o[j,i]=p['g']
    pen=np.diag([p['K']*(AR[l[0]][1]!=AR[r[0]][0]) for l in LET for r in LET])
    v=np.zeros_like(h)
    for a in (0,1):
        inn=IX['m','v',None]*d+IX['d','v',a]
        for s,t,amp in [('A','B',np.exp(1j*p['theta'])),('B','A',np.exp(-1j*p['theta']))]:
            out=IX['m',s,None]*d+IX['d',t,a];v[out,inn]=amp;v[inn,out]=np.conj(amp)
    return h,pen,o,v
def bond(p,lam=1.):
    h,pen,o,v=parts(p);return h+pen+.5*(np.kron(o,np.eye(d))+np.kron(np.eye(d),o))+lam*v
def words(n):
    out=[]
    for bits in itertools.product((0,1),repeat=n):
        for j in range(n):
            for k,s in [('d','v'),('d','A'),('d','B'),('x','A'),('x','B')]:
                out.append(tuple([('L','',bits[i]) for i in range(j)]+[(k,s,bits[j])]+[('R','',bits[i]) for i in range(j+1,n)]+[('e','',None)]))
        for s in 'vAB':out.append(tuple([('L','',a) for a in bits]+[('m',s,None)]))
    return out
def cell(p):
    n=p['n'];w=words(n);ix={x:i for i,x in enumerate(w)};N=len(w)
    h,pen,o,_=parts(p);h+=pen;H=np.zeros((N,N),complex)
    for i,x in enumerate(w):
        for pos,l in enumerate(x):
            for j in np.flatnonzero(o[:,IX[l]]):
                y=list(x);y[pos]=LET[j];H[ix[tuple(y)],i]+=o[j,IX[l]]
        full=(('e','',None),)+x
        for pos in range(n+1):
            col=IX[full[pos]]*d+IX[full[pos+1]]
            for row in np.flatnonzero(h[:,col]):
                y=list(full);y[pos]=LET[row//d];y[pos+1]=LET[row%d]
                assert y[0]==full[0]
                H[ix[tuple(y[1:])],i]+=h[row,col]
    # Local creation maps at the two boundary ports.
    fs={}
    for end in ['left','right']:
        for s in 'AB':
            F=np.zeros_like(H)
            for i,x in enumerate(w):
                pos=0 if end=='left' else n;k,sp_,a=x[pos]
                if (end=='left' and k=='d' and sp_=='v') or (end=='right' and k=='m' and sp_=='v'):
                    y=list(x);y[pos]=(k,s,a);F[ix[tuple(y)],i]=1.
            fs[end,s]=sp.csr_matrix(F)
    # Little-endian logical bit order, unlike the physical tuple enumeration.
    G=np.zeros((N,2**n),complex)
    vv=np.sqrt(2/(n+2))*np.sin(np.pi*np.arange(1,n+2)/(n+2))
    for i,x in enumerate(w):
        if not any(l[1]=='v' for l in x):continue
        bits=[l[2] for l in x if l[2] is not None];b=sum(a<<j for j,a in enumerate(bits))
        j=next(j for j,l in enumerate(x) if l[1]=='v');G[i,b]=vv[j]
    exc=np.array([any(l[1] in ['A','B'] for l in x) for x in w])
    return H,G,fs,exc,w
def pair(p):
    H,G,fs,ex,w=cell(p);N=len(H)
    V=np.exp(1j*p['theta'])*sp.kron(fs['right','A'],fs['left','B'])+np.exp(-1j*p['theta'])*sp.kron(fs['right','B'],fs['left','A'])
    V=(V+V.conj().T).tocsr();HH=sp.kronsum(sp.csr_matrix(H),sp.csr_matrix(H),format='csc')
    # kron cell A then B maps bit index (a,b) to a+2^n b.
    GG=np.kron(G,G);perm=[a*(2**p['n'])+b for b in range(2**p['n']) for a in range(2**p['n'])];GG=GG[:,perm]
    E0=p['musel']-2*p['t']*np.cos(np.pi/(p['n']+2))
    sel=np.flatnonzero(np.kron(ex,ex));S=(V@GG)[sel]
    A=HH[sel][:,sel]-2*E0*sp.eye(len(sel),format='csc')
    H2=-S.conj().T@splu(A.tocsc()).solve(S)
    p0=dict(p,g=0.);H0=cell(p0)[0];C=H-H0
    AA=sp.kronsum(sp.csr_matrix(H0),sp.csr_matrix(H0),format='csc')[sel][:,sel]-2*E0*sp.eye(len(sel),format='csc')
    CC=sp.kronsum(sp.csr_matrix(C),sp.csr_matrix(C),format='csc')[sel][:,sel]/p['g']
    lu=splu(AA.tocsc());X=lu.solve(S);Y=CC@X
    H22=-Y.conj().T@lu.solve(Y)
    H20=-S.conj().T@X
    return dict(HH=HH,V=V,G=GG,H2=H2,H22=H22,H20=H20,E0=E0,cell=H,code=G,exc=ex)
