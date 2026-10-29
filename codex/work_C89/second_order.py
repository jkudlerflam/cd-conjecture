"""Finite-volume charge resolvent, with operator-valued clock sources."""
from pathlib import Path
import json
import itertools
import math
import numpy as np
import scipy.sparse.linalg as sla
from check_defects import bath,selector

HERE=Path(__file__).resolve().parent

def source_vectors(lengths,pins,mu):
    dims=[1 if p else b for b,p in zip(lengths,pins)]
    states=list(itertools.product(*(range(d) for d in dims)))
    idx={s:i for i,s in enumerate(states)}
    psi=np.ones(1)
    for d in dims:
        w=np.sin(np.pi*np.arange(1,d+1)/(d+1)); w/=np.linalg.norm(w)
        psi=np.kron(psi,w)
    N=sum(lengths); Y=np.zeros((len(psi),N)); start=0
    for c,(b,p) in enumerate(zip(lengths,pins)):
        marker=start+b-1
        Y[:,(marker-1)%N]+=mu*psi/2
        Y[:,marker]+=mu*psi/2
        if not p:
            for q in range(b-1):
                for s,i in idx.items():
                    if s[c]==q:
                        z=list(s);z[c]+=1;k=idx[tuple(z)]
                        Y[i,start+q]-=psi[k];Y[k,start+q]-=psi[i]
        start+=b
    return psi,Y

def susceptibility(N,V):
    H,A,psi,E0,alpha,res=bath(N,V)
    def mv(v): return H@v-E0*v+psi*np.dot(psi,v)
    op=sla.LinearOperator(H.shape,matvec=mv,dtype=float)
    source=(A[:,0]-alpha)*psi
    v,info=sla.cg(op,source,rtol=1e-12,atol=1e-14,maxiter=10000)
    assert info==0
    row=(A-alpha).T@(psi*v)
    chi=np.array([[row[(j-i)%N] for j in range(N)] for i in range(N)])
    assert np.max(abs(chi-chi.T)) < 1e-10
    return chi,float(np.linalg.norm(mv(v)-source))

def run():
    ed=json.loads((HERE/'ed_results.json').read_text())
    rows=[]; cache={}; coefficients={}
    for c in ed['clock']:
        vals=[]
        for which in ['reference','defect']:
            r=c[which];key=(tuple(r['lengths']),tuple(r['pins']),r['V'])
            if key not in coefficients:
                bk=(r['N'],r['V'])
                if bk not in cache: cache[bk]=susceptibility(*bk)
                chi,res=cache[bk]
                psi,Y=source_vectors(r['lengths'],r['pins'],c['mu_sel'])
                B=float(np.sum(chi*(Y.T@Y)))
                assert B > -1e-10
                coefficients[key]=(B,res)
            vals.append(coefficients[key][0])
        J=c['reference']['J']
        deltaB=vals[1]-vals[0]
        rows.append(dict(b=c['b'],N=c['reference']['N'],V=c['reference']['V'],
                    lengths=c['defect']['lengths'],pins=c['defect']['pins'],J=J,
                    B_reference=vals[0],B_defect=vals[1],delta_B=deltaB,
                    predicted_residual_over_J2=-deltaB,
                    measured_residual_over_J2=c['residual_difference']/J**2,
                    second_order_cost=c['leading_cost']-J**2*deltaB,
                    exact_cost=c['cost']))
    out=dict(tag='numerical_finite_volume_quantum_source_susceptibility',
             largest_linear_solve_residual=max(r[1] for r in cache.values()),rows=rows)
    (HERE/'second_order_results.json').write_text(json.dumps(out,indent=2)+'\n')
    for r in rows:
        if r['J']==min(z['J'] for z in rows if z['N']==r['N']):
            print(r)
    print('max linear residual',out['largest_linear_solve_residual'])

if __name__=='__main__': run()
