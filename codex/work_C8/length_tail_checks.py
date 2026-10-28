"""Finite scalar-cell control of the factorial grand-canonical length tail."""
import json
from pathlib import Path
import numpy as np
from scipy.linalg import eigh,norm
from local_dressing_checks import embed


def main():
    b=3;lengths=np.arange(2,12);d=len(lengths);M=3;t=.02;gamma=1.
    energy=np.array([0. if ell in (b,b+1) else gamma*(1+abs(ell-b)) for ell in lengths])
    shift=np.diag(np.ones(d-1),-1)
    hopping=-t*(np.kron(shift,shift.T)+np.kron(shift.T,shift))
    v=2*float(norm(hopping,2))
    H=sum(embed(np.diag(energy),j,1,M,d) for j in range(M)).astype(complex)
    for j in range(M):
        for left,right in ((shift,shift.T),(shift.T,shift)):
            term=np.array([[1.]])
            for k in range(M):term=np.kron(term,left if k==j else right if k==(j+1)%M else np.eye(d))
            H-=t*term
    e,U=eigh(H,subset_by_index=(0,0));psi=U[:,0]
    weights=(abs(psi)**2).reshape((d,)*M)
    local=weights.sum(axis=(1,2))
    rows=[];prev=float(np.sqrt(local[lengths>=b+1].sum()))
    for R in (b+2,b+3,b+4):
        actual=float(np.sqrt(local[lengths>=R].sum()))
        step=v/(gamma*(R-b)-2*v)*prev
        from math import factorial
        bound=(v/gamma)**(R-b-1)/factorial(R-b-1)
        assert actual<=step+1e-13 and actual<=bound+1e-13
        rows.append(dict(cutoff=R,tail_amplitude=actual,recursion_bound=step,factorial_bound=bound))
        prev=actual
    out=dict(status='finite check of unconstrained length-tail theorem',cells=M,
             site_dimension=d,dimension=d**M,b=b,hopping=t,local_norm=v,
             ground=float(e[0]),rows=rows)
    Path('length_tail_checks.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':main()
