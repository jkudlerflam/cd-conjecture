"""C85 literal eleven-material, three-color vacancy bond.
All float parameters are for resolved diagnostics. The certified hierarchy
is specified symbolically in PROOF.md and finite_budget.py.
"""
import numpy as np
MAT=[('L',0),('L',1),('R',0),('R',1),('d','A'),('d','B'),
     ('x','A'),('x','B'),('e',0),('m','A'),('m','B')]
LET=MAT+[('V',0),('V',1),('V',2)]
IX={a:i for i,a in enumerate(LET)}
def outgoing(a):
    if a[0]=='V':return a[1]
    if a==('L',1):return 2
    return int(a[0] in ('R','d','x'))
def incoming(a):
    if a[0]=='V':return (a[1],)
    if a==('R',0):return (0,1,2)
    if a==('R',1):return (1,)
    if a==('e',0):return (0,1)
    if a[0]=='m':return (0,)
    return (0,2)
def allowed(a,b): return outgoing(a) in incoming(b)
def charge(a):
    if a[0]=='V': return 0
    if a[0] in ('L','R','e'): return 2*a[1]
    return -1 if a[1]=='A' else 1
def marker(a): return a[0] in ('e','m')
def params(n=3,**kw):
    """n counts physical data sites, with n-1 active bits."""
    p=dict(n=n,D=4.,Delta=4-2*np.cos(np.pi/(n+2)),g=.07,
           td=.003,K=100.,musel=0.,field=0.,flux=np.pi/(2*n),theta=-np.pi/4)
    p.update(kw)
    return p
def material_bond(p,lam):
    d=11;h=np.zeros((d*d,d*d),complex)
    def add(out,inn,z):
        i=IX[out[0]]*d+IX[out[1]];j=IX[inn[0]]*d+IX[inn[1]]
        h[i,j]+=z;h[j,i]+=np.conj(z)
    for a in (0,1):
        la,ra,ea=('L',a),('R',a),('e',a)
        add((la,('R',0)),(('R',0),ra),-1.)
        if a==0:add((la,('e',0)),(('R',0),ea),-1.)
        for s in ('A','B'):
            z=(-1)**a if s=='A' else 1
            phi=p['flux'] if s=='A' else 0.
            psi=p['n']*p['flux']+np.pi/2 if s=='A' else np.pi/2
            add((la,('d',s)),(('d',s),ra),-z)
            add((la,('x',s)),(('x',s),ra),-z*np.exp(1j*phi))
            if a==0:add((la,('m',s)),(('x',s),ea),-z*np.exp(1j*phi))
            if a==0:add((la,('m',s)),(('d',s),ea),-p['td']*z*np.exp(1j*psi))
    for m in [a for a in MAT if marker(a)]:
        for s in ('A','B'): add((m,('x',s)),(m,('d',s)),-1.)
    for k in ('R',): add((('d','B'),(k,0)),(('d','A'),(k,1)),p['g'])
    add((('m','A'),('d','B')),(('e',0),('R',0)),lam*np.exp(1j*p['theta']))
    add((('m','B'),('d','A')),(('e',0),('R',0)),lam*np.exp(-1j*p['theta']))
    o=[]
    for a in MAT:
        v=p['musel']*marker(a)
        if a[0] in ('d','x','m'): v+=p['Delta']+p['D']*(a[0] in ('x','m'))
        if a[0] in ('L','R','e'): v+=p['field']*a[1]
        o.append(v)
    h+=np.diag([(o[i]+o[j])/2+p['K']*(not allowed(a,b))
                for i,a in enumerate(MAT) for j,b in enumerate(MAT)])
    return h
def bond(p,lam=.1,J=1.,t=1.,mu=0.,colored=True):
    ls=LET if colored else MAT+[('V',0)];d=len(ls)
    h=np.zeros((d*d,d*d),complex)
    occ=[d*i+j for i in range(11) for j in range(11)]
    h[np.ix_(occ,occ)]=J*material_bond(p,lam)
    for i,a in enumerate(MAT):
        for q in incoming(a) if colored else (0,):
            src=(11+q)*d+i;dst=i*d+11+(outgoing(a) if colored else 0)
            h[dst,src]=h[src,dst]=-t
    for i,a in enumerate(ls):
        for j,b in enumerate(ls):
            if colored and (i>=11 or j>=11):
                h[d*i+j,d*i+j]+=J*p['K']*(not allowed(a,b))
            h[d*i+j,d*i+j]-=mu*((i<11)+(j<11))/2
    return h
