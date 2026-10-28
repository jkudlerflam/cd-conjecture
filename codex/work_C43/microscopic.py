"""C43 literal core lift, adapted from C39/literal.py's published rules.

No neighboring imports. All matrices are microscopic site-basis matrices.
"""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
import itertools as it
import numpy as np
import scipy.sparse as sp


def core(n, w=.8, v=.65, D=.2, g=.4, q=.8, a=4., K=10.,
         lam=.03, eta=.6, readout=.03):
    species = ['v'] + [(r,s) for r in (1,2) for s in (0,1,-1)]
    ex = species[1:]
    letters = [(z,None,b) for z in ('L','R') for b in (0,1)]
    letters += [('d',s,b) for s in species for b in (0,1)]
    letters += [('x',s,b) for s in ex for b in (0,1)]
    letters += [(z,s,None) for z in ('A','B') for s in [None]+species]
    ix = {x:i for i,x in enumerate(letters)}
    d = len(letters)
    assert d == 46
    Q = np.array([s in ex for z,s,b in letters],dtype=float)
    high = np.array([s in ex and z in ('A','B','x') for z,s,b in letters])
    charges = np.array([(b or 0)+(s[1] if s in ex else 0) for z,s,b in letters])
    site = sp.diags((a-3)*Q+D*high,format='lil',dtype=complex)
    for r in (1,2):
        for x,y in [(('d',(r,0),1),('d',(r,1),0)),
                    (('d',(r,0),0),('d',(r,-1),1))]:
            site[ix[x],ix[y]] = site[ix[y],ix[x]] = g
    site = site.tocsr()
    bond = sp.diags(np.repeat(Q,d)+2*np.tile(Q,d),format='lil',dtype=complex)
    def edge(x,y,z,vac=False):
        i,j = (ix[p[0]]*d+ix[p[1]] for p in (x,y))
        if vac:
            bond[i,i] += q*q
            bond[j,j] += 1
            z = -q
        bond[j,i] += z
        bond[i,j] += np.conj(z)
    for b,c in it.product((0,1),repeat=2):
        edge((('d','v',b),('R',None,c)),(('L',None,b),('d','v',c)),-q,True)
        for r,s in ex:
            edge((('d',(r,s),b),('R',None,c)),(('L',None,b),('d',(r,s),c)),
                 -w*(v if r==2 else 1)*((-1)**b if s==0 else 1))
            edge((('x',(r,s),b),('R',None,c)),(('L',None,b),('x',(r,s),c)),
                 -np.exp(1j*np.pi/(2*n)*(s==0))*((-1)**b if s==0 else 1))
    for b in (0,1):
        edge((('A','v',None),('R',None,b)),(('A',None,None),('d','v',b)),-q,True)
        edge((('d','v',b),('B',None,None)),(('L',None,b),('B','v',None)),-q,True)
        for r,s in ex:
            edge((('A',(r,s),None),('R',None,b)),(('A',None,None),('d',(r,s),b)),-1)
            edge((('x',(r,s),b),('B',None,None)),(('L',None,b),('B',(r,s),None)),
                 -np.exp(1j*np.pi/(2*n)*(s==0))*((-1)**b if s==0 else 1))
            for am in [None]+species:
                edge((('A',am,None),('d',(r,s),b)),(('A',am,None),('x',(r,s),b)),-1)
    def incoming(x):
        z,s,b=x
        if z=='B': return 'R' if s is None else 'L'
        return 'R' if z=='R' else 'L'
    def outgoing(x):
        z,s,b=x
        if z=='A': return 'L' if s is None else 'R'
        return 'L' if z=='L' else 'R'
    bad = np.array([r[0]=='A' or l[0]=='B' or outgoing(l)!=incoming(r)
                    for l in letters for r in letters])
    bond = bond.tocsr()+sp.diags(K*bad)
    eye = sp.eye(d,format='csr')
    bond += (sp.kron(site,eye)+sp.kron(eye,site))/2
    R = sp.diags([int(z=='d' and s=='v' and b==1) for z,s,b in letters],dtype=float)
    bond += readout*lam**2*(sp.kron(R,eye)+sp.kron(eye,R))/2
    def F(end,r,s):
        out=sp.lil_matrix((d,d),dtype=complex)
        def add(a,b):out[ix[(end,a,None)],ix[(end,b,None)]]+=1
        if s==0:
            add((r,0),'v');add('v',(r,0))
        elif s==1:
            add((r,1),'v');add('v',(r,-1))
        else:return F(end,r,1).getH()
        return out.tocsr()
    for r,amp,phase in [(1,1,np.pi/2),(2,np.sqrt(eta),np.pi)]:
        bond += lam*amp*(sp.kron(F('B',r,0),F('A',r,0))
             +np.exp(1j*phase)*sp.kron(F('B',r,1),F('A',r,-1))
             +np.exp(-1j*phase)*sp.kron(F('B',r,-1),F('A',r,1)))
    bond = bond.tocsr();bond.eliminate_zeros()
    return letters,charges,bond


def lift(h,t=1.,J=.07,mu=.17):
    d=int(round(np.sqrt(h.shape[0])));D=d+1
    h=h.tocoo()
    rows=list((h.row//d+1)*D+h.row%d+1)
    cols=list((h.col//d+1)*D+h.col%d+1)
    vals=list(J*h.data)
    for a in range(1,D):
        rows.extend([a,a*D]);cols.extend([a*D,a]);vals.extend([-t,-t])
    for a,b in it.product(range(D),repeat=2):
        rows.append(a*D+b);cols.append(a*D+b)
        vals.append(-mu*((a!=0)+(b!=0))/2)
    out=sp.coo_matrix((vals,(rows,cols)),shape=(D*D,D*D)).tocsr()
    out.eliminate_zeros()
    return out


def ring_full3(h,d=47):
    """All d^3 states, including every occupied-number and illegal sector."""
    h=h.tocoo();rows=[];cols=[];vals=[]
    for a,b in ((0,1),(1,2),(2,0)):
        other=3-a-b;powers=[d*d,d,1]
        for c in range(d):
            rows.append((h.row//d)*powers[a]+(h.row%d)*powers[b]+c*powers[other])
            cols.append((h.col//d)*powers[a]+(h.col%d)*powers[b]+c*powers[other])
            vals.append(h.data)
    out=sp.coo_matrix((np.concatenate(vals),(np.concatenate(rows),np.concatenate(cols))),
                       shape=(d**3,d**3)).tocsr()
    out.eliminate_zeros()
    return out


def legal_words(n,letters):
    """Complete initialized cell, then every cyclic translate, no code projection."""
    ix={z:i+1 for i,z in enumerate(letters)}
    species=['v']+[(r,s) for r in (1,2) for s in (0,1,-1)]
    words=set();N=n+2
    for s in species:
        for p in range(N if s=='v' else 2*n+2):
            if s=='v':kind='A' if p==0 else ('B' if p==N-1 else 'd');site=p
            elif p==0:kind='A';site=0
            elif p==2*n+1:kind='B';site=N-1
            elif p<=n:kind='d';site=p
            else:kind='x';site=p-n
            for bits in range(2**n):
                word=[('A',s if kind=='A' else None,None)]
                for j in range(n):
                    bit=(bits>>j)&1
                    word.append((kind,s,bit) if site==j+1 else
                                ('L' if j+1<site else 'R',None,bit))
                word.append(('B',s if kind=='B' else None,None))
                w=tuple(ix[x] for x in word)
                for r in range(N):words.add(w[r:]+w[:r])
    return sorted(words)


def restricted_physical(h,words):
    """Apply literal bond to each physical word and reject any leakage."""
    d=int(round(np.sqrt(h.shape[0])));h=h.tocsc();ix={w:i for i,w in enumerate(words)}
    rows=[];cols=[];vals=[];L=len(words[0])
    for col,w in enumerate(words):
        for x in range(L):
            y=(x+1)%L;c=w[x]*d+w[y]
            for k in range(h.indptr[c],h.indptr[c+1]):
                z=h.indices[k];v=list(w);v[x]=z//d;v[y]=z%d;v=tuple(v)
                assert v in ix,('leakage',w,v)
                rows.append(ix[v]);cols.append(col);vals.append(h.data[k])
    out=sp.coo_matrix((vals,(rows,cols)),shape=(len(words),len(words))).tocsr()
    out.eliminate_zeros();return out
