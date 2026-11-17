"""Generic real-q sixth constraint, exact symbolic interval calculation."""
import json
import time
from pathlib import Path
import numpy as np
import sympy as sp

import exact_local as ex
from model import embed


def main():
    start=time.perf_counter();ex.canon=lambda x:sp.cancel(x)
    q=sp.symbols('q',real=True);L=9;dim=1<<L
    v0=np.array([0,0,0,1]);v1=np.array([0,-1,1,0]);v2=np.array([1,0,0,0]);o=np.outer
    H=[ex.exact(embed(L,o(v0,v0),False)),ex.exact(embed(L,o(v0,v1)+o(v1,v0),False)),
       ex.exact(embed(L,o(v1,v1),False))+q*ex.exact(embed(L,o(v0,v2)+o(v2,v0),False)),
       q*ex.exact(embed(L,o(v1,v2)+o(v2,v1),False)),q*q*ex.exact(embed(L,o(v2,v2),False)),
       sp.zeros(dim),sp.zeros(dim)]
    words=list(range(dim))
    for n in (0,2,4):
        keep=[i for i in range(len(words)) if sp.cancel(H[0][i,i])==0]
        H=ex.schur(H,keep);words=[words[i] for i in keep]
        if n<4:H=H[2:]
    cs={5:(q*q-10*q+5)**2/3,
        6:(3*q-1)**2*(q*q-6*q+25)/12,
        7:(q**3-21*q*q+35*q-7)**2/(4*(5*q*q-18*q+21)),
        8:(q*q-6*q+1)**2/8}
    expected=[]
    for state in words:
        n=[(state>>j)&1 for j in range(L)]
        expected.append(sum(c*sum(n[j]*n[j+d] for j in range(L-d)) for d,c in cs.items()))
    diff=H[2]-sp.diag(*expected)
    assert all(sp.cancel(v)==0 for v in diff.todok().values())
    out={'L':L,'chart_excludes_q':[1,3],
         'coefficients':{str(d):str(sp.factor(c)) for d,c in cs.items()},
         'full_symbolic_matrix_identity':True,
         'resonances':{str(m):{str(d):str(sp.simplify(c.subs(q,qq))) for d,c in cs.items()}
                       for m,qq in ((5,5-2*sp.sqrt(5)),(6,sp.Rational(1,3)))},
         'runtime_seconds':time.perf_counter()-start}
    Path(__file__).with_name('sixth_constraint_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':main()
