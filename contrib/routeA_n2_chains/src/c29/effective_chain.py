"""C29 n=2 effective open chain, in units J=lambda^2 alpha_2^4 g^2.
leading: normal free two-flavor reference.
full: exact order-lambda^2 bond at fixed g, including non-Gaussian terms.
Neither includes order-lambda^3 corrections. NumPy/SciPy only.
"""
import json,argparse
import numpy as np
from scipy import sparse as sp
from scipy.linalg import eigh
from model import data_ops,target_data,combined_bond


def tensor(xs):
    ans=sp.csr_matrix([[1.]])
    for x in xs:ans=sp.kron(ans,x,format='csr')
    return ans


def logical_bond(g=.2,level='leading'):
    T,mu,W=target_data(2,g)
    if level=='full':
        K=combined_bond(2,g)/g**2
        K=K-np.trace(K)*np.eye(16)/16
        return (K+K.conj().T)/2,W
    cs,_=data_ops(4);K=np.zeros((16,16),complex)
    for i in range(2):
        for j in range(2):
            A=T[i,j]*(cs[i].conj().T@cs[j+2]);K+=A+A.conj().T
    return K,W


def one_body(cells):
    T,mu,_=target_data(2);H=np.kron(np.eye(cells),np.diag(mu)).astype(complex)
    for x in range(cells-1):
        H[2*x:2*x+2,2*x+2:2*x+4]=T
        H[2*x+2:2*x+4,2*x:2*x+2]=T.conj().T
    return H


def chain(cells=4,g=.2,level='leading'):
    K,W=logical_bond(g,level);D=4**cells
    H=sum((tensor([sp.eye(4**j),K,sp.eye(4**(cells-j-2))]) for j in range(cells-1)),sp.csr_matrix((D,D)))
    H+=sum((tensor([sp.eye(4**j),W,sp.eye(4**(cells-j-1))]) for j in range(cells)),sp.csr_matrix((D,D)))
    return H

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--cells',type=int,default=4);p.add_argument('--level',choices=['leading','full'],default='leading');p.add_argument('--g',type=float,default=.2);a=p.parse_args()
    if not 1<=a.cells<=4:raise ValueError('This local dense diagnostic is capped at four cells; import chain for sparse work.')
    H=chain(a.cells,a.g,a.level);vals=eigh(H.toarray(),eigvals_only=True)
    out=dict(cells=a.cells,logical_dimension=H.shape[0],level=a.level,g=a.g,hermiticity=float(sp.linalg.norm(H-H.conj().T)),lowest_eigenvalues=vals[:16].tolist())
    if a.level=='leading':
        es=eigh(one_body(a.cells),eigvals_only=True);sums=np.array([0.])
        for e in es:sums=np.r_[sums,sums+e]
        out['free_spectrum_residual']=float(np.max(np.abs(vals-np.sort(sums))))
    print(json.dumps(out,indent=2))
    with open(f'chain_{a.level}_M{a.cells}.json','w') as f:json.dump(out,f,indent=2)
