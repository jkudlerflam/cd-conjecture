"""Displace q by alpha*r inside the exact parent; retain all companion terms."""
import json
import time
from pathlib import Path
import numpy as np
import sympy as sp

import exact_local as ex
from model import embed,hierarchy,qvalue,maxabs,occupations


def exact_case(m):
    q=5-2*sp.sqrt(5) if m==5 else sp.Rational(1,3)
    a=sp.symbols('alpha',real=True);L=m+1;dim=1<<L
    bare=[np.array(v) for v in ([0,0,0,1],[0,-1,1,0],[1,0,0,0],[1,0,0,0])]
    pref=[1,1,q,a];H=[]
    for n in range(7):
        h=sp.zeros(dim)
        for i in range(4):
            if 0<=n-i<4:
                h+=pref[i]*pref[n-i]*ex.exact(embed(L,np.outer(bare[i],bare[n-i]),False))
        H.append(h)
    H += [sp.zeros(dim),sp.zeros(dim)];words=list(range(dim))
    for n in (0,2,4,6):
        H=[h.applyfunc(ex.canon) for h in H]
        keep=[i for i in range(len(words)) if H[0][i,i]==0]
        assert all(i==j or v==0 for (i,j),v in H[0].todok().items())
        H=ex.schur(H,keep);words=[words[i] for i in keep]
        assert all(v==0 for v in H[1].todok().values())
        if n<6:H=H[2:]
    t=20992-sp.Rational(28160,3)*sp.sqrt(5) if m==5 else sp.Rational(832,2187)
    c=(1600-704*sp.sqrt(5))/3 if m==5 else sp.Rational(208,81)
    v=sp.Rational(80,3) if m==5 else sp.Rational(52,3)
    expected=sp.zeros(len(words));ids=[words.index(s) for s in (0,1,1<<m,1+(1<<m))]
    block=sp.Matrix([[0,0,0,0],[0,t,-t,c*a],[0,-t,t,-c*a],[0,c*a,-c*a,v*a*a]])
    for i,ii in enumerate(ids):
        for j,jj in enumerate(ids):expected[ii,jj]=block[i,j]
    assert all(ex.canon(x)==0 for x in (H[2]-expected).todok().values())
    assert ex.canon(c*c-t*v)==0
    return {'m':m,'L':L,'periodic':False,'symbolic_alpha':True,'full_matrix_identity':True,
            't':str(t),'c':str(c),'v':str(v),'positive_rank_one_identity_c_squared_equals_tv':True}


def larger_case(m,L,periodic,alpha):
    words,K,_=hierarchy(L,qvalue(m),periodic=periodic,q_slope=alpha)
    base,B,_=hierarchy(L,qvalue(m),periodic=periodic)
    assert np.array_equal(words,base)
    c=(1600-704*np.sqrt(5))/3 if m==5 else 208/81
    v=80/3 if m==5 else 52/3
    ns=occupations(L)[words];positions={int(s):i for i,s in enumerate(words)}
    expected=np.zeros(K[2].shape,complex)
    for col,state in enumerate(words):
        for j in range(L if periodic else L-m):
            k=(j+m)%L
            expected[col,col]+=v*alpha*alpha*ns[col,j]*ns[col,k]
            for site,weight in ((k,ns[col,j]),(j,-ns[col,k])):
                dest=int(state)^(1<<site)
                if dest in positions:expected[positions[dest],col]+=c*alpha*weight
    residual=maxabs(K[2]-B[2]-expected)
    assert residual<1e-9 and maxabs(K[1])<1e-9
    return {'m':m,'L':L,'periodic':periodic,'alpha':alpha,'dimension':len(words),
            'whole_matrix_displacement_error':residual,'seventh_order_norm':maxabs(K[1]),
            'number_changing_norm':maxabs(expected*(ns.sum(axis=1)[:,None]!=ns.sum(axis=1)[None,:]))}


def main():
    start=time.perf_counter();exact=[];larger=[]
    for m in (5,6):
        exact.append(exact_case(m));print(exact[-1],flush=True)
    for args in ((5,13,False,-.4),(5,15,True,.7),(6,16,False,.7),(6,12,True,-.4)):
        larger.append(larger_case(*args));print(larger[-1],flush=True)
    out={'exact_interval_results':exact,'larger_full_sector_checks':larger,
         'runtime_seconds':time.perf_counter()-start}
    Path(__file__).with_name('resonance_displacement_results.json').write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':main()
