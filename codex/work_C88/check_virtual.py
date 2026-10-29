"""All second-order grammar terms on every legal four-site ring.
This constructs illegal intermediates directly from the literal bond.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import itertools,json
import numpy as np
import scipy.sparse as sp
import bond13 as b

def check(L=4):
    ls=b.LET;d=len(ls);p=b.params(n=2,K=0.,musel=.93,field=.021)
    lam=.17
    T=b.bond(p,lam=lam,J=0.,t=1.)
    W=b.bond(p,lam=lam,J=1.,t=0.)
    matrices=[T,W]
    tr=[{i:[(int(j),h[j,i]) for j in np.flatnonzero(h[:,i])] for i in range(d*d)} for h in matrices]
    def penalty(w):return sum(not b.allowed(ls[a],ls[w[(i+1)%L]]) for i,a in enumerate(w))
    legal=[w for w in itertools.product(range(d),repeat=L) if penalty(w)==0 and any(a<11 for a in w)]
    high={};records=[[],[]];pen=[]
    for col,w in enumerate(legal):
        for typ in range(2):
            for i,a in enumerate(w):
                j=(i+1)%L
                for q,v in tr[typ][a*d+w[j]]:
                    z=list(w);z[i],z[j]=divmod(q,d);z=tuple(z)
                    nz=penalty(z)
                    if not nz:continue
                    if z not in high:high[z]=len(high);pen.append(nz)
                    records[typ].append((high[z],col,v))
    def matrix(rec):
        r,c,v=zip(*rec)
        return sp.coo_matrix((v,(r,c)),shape=(len(high),len(legal))).tocsc()
    A,B=map(matrix,records);inv=sp.diags(1/np.array(pen,dtype=float))
    tt=A.getH()@inv@A;ww=B.getH()@inv@B;cross=A.getH()@inv@B
    ct=[];cv=[];cp=[]
    for w in legal:
        ct.append(sum(len(b.incoming(ls[a]))-1 for i,a in enumerate(w) if a<11 and w[(i+1)%L]>=11))
        vac=port=0
        for i,a in enumerate(w):
            pred=ls[w[(i-1)%L]];nxt=ls[w[(i+1)%L]];a=ls[a]
            vac+=int(b.outgoing(pred)==1 and a==('R',0) and nxt[0] in ('R','e'))
            port+=int(b.outgoing(pred)==1 and a==('e',0) and nxt==('R',0))
        cv.append(vac);cp.append(port)
    maxabs=lambda a:float(np.max(abs(a.data),initial=0))
    errt=maxabs(tt-sp.diags(ct,dtype=float));errw=maxabs(ww-sp.diags(np.array(cv)+2*lam**2*np.array(cp)));errx=maxabs(cross)
    assert max(errt,errw,errx)<1e-12
    assert set(pen)=={1}
    out=dict(L=L,legal_dimension=len(legal),illegal_intermediates=len(high),
      tt_diagonal_formula_error=errt,ww_diagonal_formula_error=errw,mixed_tW_error=errx,
      illegal_penalty_values=sorted(set(pen)),
      ambiguity={str(a):len(b.incoming(a))-1 for a in b.MAT})
    print(json.dumps(out),flush=True)
    return out
if __name__=='__main__':
    json.dump([check(3),check(4),check(5)],open('virtual_checks.json','w'),indent=2)
