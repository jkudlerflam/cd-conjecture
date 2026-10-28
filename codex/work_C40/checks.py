"""Exact finite checks for C40. No microscopic phase simulation."""
from fractions import Fraction as F
from itertools import combinations, product
from math import comb
import json
from pathlib import Path

out = {"tag": "proved by exact integer/rational enumeration"}
rows=[]
for n in range(1,7):
    N=2*n
    counts={1:0,2:0}
    # All neutral nonzero integer vectors with squared norm <=4 have
    # two or four nonzero entries, all +/-1. A +/-2 cannot be neutral.
    vectors=[]
    for k in (2,4):
        for support in combinations(range(N), k):
            for signs in product((-1,1), repeat=k):
                if sum(signs): continue
                z=[0]*N
                for j,s in zip(support,signs): z[j]=s
                vectors.append(tuple(z))
                counts[k//2]+=1
    assert counts[1]==N*(N-1)
    assert counts[2]==6*comb(N,4)
    payload=sum(counts.values())//2
    mobile=payload+1+N*(N-1)
    assert payload==comb(N,2)+3*comb(N,4)
    assert mobile==1+3*comb(N,2)+3*comb(N,4)
    rows.append(dict(n=n, branches=N, oriented_vectors=len(vectors),
                     payload_conjugate_pairs=payload,
                     mobile_conjugate_pairs=mobile))
out['vertex_counts']=rows
# Weak relative to the fast velocity need not preserve positivity.
vslow=F(1,10**12)
u=F(1,10**5)
bad_det=vslow-u*u
assert bad_det<0
u_good=F(1,10**7)
good_det=vslow-u_good*u_good
assert good_det>0
out['hierarchy_hessian']={'v_fast':'1','v_slow':str(vslow),
    'bad_off_diagonal':str(u),'bad_determinant':str(bad_det),
    'good_off_diagonal':str(u_good),'good_determinant':str(good_det)}
# At eta=1/4, kappa^2=3/5 and (3*kappa)^2>4.
eta=F(1,4)
kappa2=(1-eta)/(1+eta)
assert 9*kappa2>4
assert (1-F(5,13))/(1+F(5,13))==F(4,9)
out['stiffness']={'eta':str(eta),'kappa_squared':str(kappa2),
                  'dimension_three_squared':str(9*kappa2)}
# Exact majorant solution x(t)=x0/(1-A*x0*t).
A=F(7); ell=F(1000); x0=F(1,100000)
xend=x0/(1-A*x0*ell)
assert A*x0*ell<F(1,2) and xend<2*x0
out['rg_majorant']={k:str(v) for k,v in dict(A=A,ell=ell,x0=x0,xend=xend).items()}
Path('checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
