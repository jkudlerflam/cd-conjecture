"""Full three-cell contraction at n=1 and exact rational Sidon certificates."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import itertools,json
from fractions import Fraction
from pathlib import Path
import numpy as np
import sidon_cell as sc
from checks_sidon import spectral,norm,tl,exact_graph,pair_engine

def exact_margins(n):
    d=sc.sidon_data(n);p=[Fraction(6*(a+d['B'])+1,6*d['T']) for a in d['A']]
    circ=lambda x:min(x%1,(-x)%1)
    d2=min([circ(2*x) for x in p]+[circ(p[i]+s*p[j]) for i,j in itertools.combinations(range(2*n),2) for s in [-1,1]])
    ds=[circ(sum((s*p[i] for s,i in zip(ss,inds)),Fraction())) for inds in itertools.combinations(range(2*n),4) for ss in itertools.product([-1,1],repeat=4)]
    d4=min(ds) if ds else None
    assert d2>0 and (d4 is None or d4>0)
    return dict(n=n,delta2_over_2pi=str(d2),delta4_over_2pi=str(d4) if d4 is not None else None)

c=sc.reflection_cell(1,True);N=c['N']
spec=[spectral(c,s) for s in ['R',None,'L']]
es=[s[0] for s in spec];fs=[s[1] for s in spec];dims=[len(e) for e in es];dim=int(np.prod(dims));low=N**3
idx=np.array([(i*dims[1]+j)*dims[2]+k for i in range(N) for j in range(N) for k in range(N)])
g=np.zeros((dim,low),complex);g[idx,np.arange(low)]=1
energy=(es[0][:,None,None]+es[1][None,:,None]+es[2][None,None,:]).ravel()
inv=np.zeros(dim);np.divide(1,energy,out=inv,where=energy>0)
def axis(a,x,j):
    v=np.tensordot(a,x,axes=(1,j));return np.moveaxis(v,0,j)
def V(x):
    x=x.reshape(*dims,-1);out=np.zeros_like(x)
    for j in [0,1]:
        for r,l in zip(fs[j]['R'],fs[j+1]['L']):out+=axis(r,axis(l,x,j+1),j)
    return out.reshape(dim,-1)
b=V(g);w=inv[:,None]*b;y=V(w);A=b.conj().T@w;T=w.conj().T@w
k2=-A;k3=w.conj().T@y;k4=-y.conj().T@(inv[:,None]*y)+(A@T+T@A)/2
stored=np.load('three_cell_coefficients_n1.npz')['K4']
e=dict(V=V,g=g,idx=idx,inv=inv,b=b,energy=energy,K2=k2,K3=k3,K4=k4)
pair=pair_engine(c)
result=dict(exact_Sidon=[exact_margins(n) for n in [1,2]],three_cell_dimension=dim,
    direct_vs_streamed_K4=norm(k4-stored),
    K2_additivity_error=norm(k2-np.kron(pair['K2'],np.eye(N))-np.kron(np.eye(N),pair['K2'])),
    K3_additivity_error=norm(k3-np.kron(pair['K3'],np.eye(N))-np.kron(np.eye(N),pair['K3'])),
    finite_coupling=[exact_graph(e,lam) for lam in [.04,.02]])
assert result['direct_vs_streamed_K4']<1e-12
assert result['K2_additivity_error']<1e-12 and result['K3_additivity_error']<1e-12
Path('independent_results.json').write_text(json.dumps(result,indent=2))
print(json.dumps(result,indent=2))
