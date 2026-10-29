"""C88 second-order local penalty normal form on a literal small ring."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import json
import numpy as np
from scipy.linalg import expm
from check_lift import closure
from bond13 import params

def main():
    H,pen,bad,states=closure([('L',0),('e',0),('V',0)],params(n=1,K=0.),.1,True)
    B=H.toarray();diff=pen[:,None]-pen[None,:];mask=diff!=0
    B0,off=B*(~mask),B*mask;out=[]
    comm=lambda a,b:a@b-b@a
    for K in (40.,80.,160.,320.):
        s1=np.zeros_like(B);s1[mask]=B[mask]/(K*diff[mask])
        f2=comm(s1,B0)+comm(s1,off)/2
        s2=np.zeros_like(B);s2[mask]=f2[mask]/(K*diff[mask])
        S=s1+s2;assert np.max(abs(S+S.conj().T))<1e-14
        U=expm(S);target=K*np.diag(pen)+B0+f2*(~mask)
        err=float(np.linalg.norm(U@(K*np.diag(pen)+B)@U.conj().T-target,2))
        good=np.flatnonzero(pen==0);other=np.flatnonzero(pen>0)
        expected=B[np.ix_(good,good)]-(B[np.ix_(good,other)]/(K*pen[other]))@B[np.ix_(other,good)]
        signerr=float(np.max(abs(target[np.ix_(good,good)]-expected),initial=0))
        assert signerr<1e-12
        rec=dict(K=K,dimension=len(pen),remainder_norm=err,remainder_times_K_squared=err*K*K,
                 second_order_compression_error=signerr,
                 unitarity_error=float(np.linalg.norm(U@U.conj().T-np.eye(len(pen)),2)))
        out.append(rec);print(json.dumps(rec),flush=True)
    assert out[-1]['remainder_norm']<out[-2]['remainder_norm']/3
    with open('sw_checks.json','w') as f:json.dump(out,f,indent=2)
if __name__=='__main__':main()
