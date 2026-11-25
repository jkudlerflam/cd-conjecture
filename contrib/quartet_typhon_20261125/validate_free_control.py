"""Dense independent check of the finite XX control's Wick expressions."""
import argparse
import numpy as np
from common import write_json


def check():
    L=8;Id=np.eye(2);sm=np.array([[0.,0.],[1.,0.]]);sp=sm.T
    def op(terms):
        x=np.ones((1,1))
        for i in range(L):
            y=Id
            for k,a in terms:
                if i==k:y=y@a
            x=np.kron(x,y)
        return x
    H=np.zeros((2**L,2**L))
    for i in range(L-1):H-=op([(i,sm),(i+1,sp)])+op([(i,sp),(i+1,sm)])
    E,V=np.linalg.eigh(H);v=V[:,0]
    sites=np.arange(1,L+1);modes=np.arange(1,L//2+1)
    U=np.sqrt(2/(L+1))*np.sin(np.pi*np.outer(sites,modes)/(L+1));C=U@U.T
    err_pair=err_transverse=0.
    for i in range(L-3):
        for j in range(i+2,L-1):
            dense=v@op([(i,sm),(i+1,sm),(j,sp),(j+1,sp)])@v
            wick=C[i,j]*C[i+1,j+1]-C[i,j+1]*C[i+1,j]
            err_pair=max(err_pair,abs(dense-wick))
    for i in range(L-1):
        for j in range(i+1,L):
            a=np.arange(i,j);b=np.arange(i+1,j+1)
            G=2*C[np.ix_(a,b)]-(a[:,None]==b[None,:])
            dense=v@op([(i,sm),(j,sp)])@v
            err_transverse=max(err_transverse,abs(dense-.5*np.linalg.det(G)))
    assert max(err_pair,err_transverse)<1e-12
    return {'passed':True,'L':L,'ground_energy':float(E[0]),'adjacent_pair_wick_error':float(err_pair),
            'transverse_determinant_error':float(err_transverse),'two_chain_rule':'independent copies factor, so square each one-chain correlation'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);a=p.parse_args()
    r=check();write_json(a.output,r);print(r)
