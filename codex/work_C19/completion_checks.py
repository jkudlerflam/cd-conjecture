"""Analytic reflection checks and higher-order Gaussian diagnostic."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import json
import numpy as np
import scipy.linalg as la
from uniform_cell import *
from checks import algebra,ed

def gaussian_residual(h,n,M):
    gs,_,P=modes(n);N=2**n
    globalg=[kron_all([P]*j+[g]+[np.eye(N)]*(M-j-1)) for j in range(M) for g in gs]
    d=N**M
    fit=np.trace(h)/d*np.eye(d,dtype=complex)
    for a in range(len(globalg)):
        for b in range(a):
            op=1j*globalg[a]@globalg[b]
            fit+=np.trace(op@h)/d*op
    return float(la.norm(h-fit,2))

def third_order(c):
    H,G=c['H'],c['G'];D,N=c['D'],c['N']
    ev,U=la.eigh(H)
    L=[U.conj().T@f@U for f in c['FL']]
    R=[U.conj().T@f@U for f in c['FR']]
    gg=U.conj().T@G
    b=sum((np.kron(r@gg,l@gg) for r,l in zip(R,L)),np.zeros((D*D,N*N),complex))
    den=(ev[:,None]+ev[None,:]).ravel()
    inv=np.zeros_like(den);inv[den>.5]=1/den[den>.5]
    w=(inv[:,None]*b).reshape(D,D,N*N)
    vw=sum((np.einsum('ab,cd,bdk->ack',r,l,w,optimize=True) for r,l in zip(R,L)),np.zeros_like(w))
    out=w.reshape(D*D,N*N).conj().T@vw.reshape(D*D,N*N)
    return dict(norm=float(la.norm(out,2)),gaussian_residual=gaussian_residual(out,c['n'],2),
                hermitian_error=float(la.norm(out-out.conj().T)))

def main():
    out=dict(reflection=[],third_order=[],reflection_ed=[])
    for n in [1,2]:
        c=reflection_cell(n);a=algebra(c);a['reflection_error']=c['reflection_error']
        out['reflection'].append(a);print(json.dumps(a),flush=True)
        for name,cc in [('reflection',c),('compact',make_cell(n))]:
            a=third_order(cc);a.update(n=n,completion=name);out['third_order'].append(a);print(json.dumps(a),flush=True)
        if n==1:out['reflection_ed']+=ed(c,2,[.04,.02])
    with open('completion_results.json','w') as f:json.dump(out,f,indent=2)

if __name__ == "__main__":main()
