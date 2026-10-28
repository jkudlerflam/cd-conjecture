"""Full periodic-ring check of the nine-letter clock selector plus #.
All clock words and all connected components, at L=3,4,5.
"""
from itertools import product
import json
import numpy as np
from scipy.linalg import eigh
import binary_clock_checks as bc

ABC=bc.ALPHABET+('#',)
PAIRS=bc.PAIRS|{('#',c) for c in bc.LEFT}|{(c,'#') for c in bc.RIGHT}
RULES=bc.RULES+[(('#','R'+b),('#','C'+b),'restart') for b in bc.BITS]+[(('I','#'),('R0','#'),'activate')]


def defects(w):return tuple((w[j],w[(j+1)%len(w)]) not in PAIRS for j in range(len(w)))


def forward(w):
    L=len(w);out=[]
    for j in range(L):
        k=(j+1)%L
        for a,b,label in RULES:
            if (w[j],w[k])==a:
                z=list(w);z[j],z[k]=b;out.append(tuple(z))
        if w[j]=='C0':
            z=list(w);z[j]='R1';out.append(tuple(z))
    return out


def check(L):
    n=L-1
    def a(m):return 4*np.sin(np.pi/(2*(m+4*(2.**m-1)+1)))**2
    lower=(n+2)*a(n)-(n+1)*a(n+1)
    upper=(n+1)*a(n-1)-n*a(n)
    delta=(lower+upper)/2;mu=2-delta;K=8.
    words=list(product(ABC,repeat=L));idx={w:i for i,w in enumerate(words)};adj=[[] for w in words]
    for i,w in enumerate(words):
        dv=defects(w)
        for v in forward(w):
            assert defects(v)==dv
            assert tuple(s=='#' for s in v)==tuple(s=='#' for s in w)
            j=idx[v];adj[i].append(j);adj[j].append(i)
    seen=set();minE=np.inf;deg=0;minbad=np.inf;maxdim=0;components=0
    expected=a(n)-delta
    for i,w in enumerate(words):
        if i in seen:continue
        todo=[i];seen.add(i);comp=[]
        while todo:
            j=todo.pop();comp.append(j)
            for k in adj[j]:
                if k not in seen:seen.add(k);todo.append(k)
        components+=1;maxdim=max(maxdim,len(comp));ix={v:j for j,v in enumerate(comp)}
        A=np.zeros((len(comp),len(comp)))
        for j in comp:
            for k in adj[j]:A[ix[j],ix[k]]-=1
        v=sum(defects(w));M=w.count('#');shift=K*v+mu*M
        eig=eigh(A,eigvals_only=True)+shift
        E=float(eig[0])
        if v:
            minbad=min(minbad,E)
            assert E >= expected+(K-4)*v-1e-11
        if E<minE-1e-10:minE=E;deg=int(np.count_nonzero(abs(eig-E)<1e-10))
        elif abs(E-minE)<1e-10:deg+=int(np.count_nonzero(abs(eig-minE)<1e-10))
    assert abs(minE-expected)<1e-11 and deg==L
    return dict(L=L,selected_interior_length=n,full_dimension=len(words),components=components,largest_component=maxdim,mu=mu,ground_energy=minE,ground_multiplicity=deg,predicted_energy=expected,lowest_illegal_energy=minbad)


if __name__=='__main__':
    out=[check(L) for L in (3,4,5)]
    with open('clock_ring_checks.json','w') as f:json.dump(out,f,indent=2)
    print(json.dumps(out,indent=2))
