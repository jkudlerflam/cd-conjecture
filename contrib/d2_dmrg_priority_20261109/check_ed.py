"""Independent sparse Pauli assembly checks of the full-chain DMRG output."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ.setdefault(k,'2')
import argparse,json,time
from pathlib import Path
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh
X=np.array([[0,1],[1,0]],complex);Y=np.array([[0,-1j],[1j,0]],complex);Z=np.diag([1.,-1.]);I=np.eye(2)

def build(L,cp,periodic=False):
    h=cp['Jx']*np.kron(X,X)+cp['Jy']*np.kron(Y,Y)+cp['Jz']*np.kron(Z,Z)
    h+=cp['D']*(np.kron(X,Y)-np.kron(Y,X))+cp.get('Kxz',0)*(np.kron(X,Z)-np.kron(Z,X))
    H=sp.csr_matrix((2**L,2**L),dtype=complex)
    for j in range(L-1):H+=sp.kron(sp.kron(sp.eye(2**j,format='csr'),sp.csr_matrix(h)),sp.eye(2**(L-j-2),format='csr'),format='csr')
    if periodic:
        wrap=lambda a,b:sp.kron(sp.kron(sp.csr_matrix(a),sp.eye(2**(L-2),format="csr")),sp.csr_matrix(b),format="csr")
        H+=cp["Jx"]*wrap(X,X)+cp["Jy"]*wrap(Y,Y)+cp["Jz"]*wrap(Z,Z)
        H+=cp["D"]*(wrap(Y,X)-wrap(X,Y))+cp.get("Kxz",0)*(wrap(Z,X)-wrap(X,Z))
    hz=cp['hz']*Z
    for j in range(L):H+=sp.kron(sp.kron(sp.eye(2**j,format='csr'),sp.csr_matrix(hz)),sp.eye(2**(L-j-1),format='csr'),format='csr')
    return H

def main():
    a=argparse.ArgumentParser();a.add_argument('input');a.add_argument('--output',required=True);a.add_argument('--periodic',action='store_true');a.add_argument('--L',type=int);args=a.parse_args()
    d=json.loads(Path(args.input).read_text());L=args.L or d['L'];H=build(L,d['couplings'],args.periodic);t=time.time()
    sector=d.get('parity')
    ids=np.array([s for s in range(1<<L) if s.bit_count()%2==(1 if sector=='odd' else 0)]) if sector else np.arange(1<<L)
    Hsector=H[ids][:,ids]
    e,vs=eigsh(Hsector,k=2,which='SA',tol=1e-12,maxiter=100000,v0=np.random.default_rng(41).normal(size=Hsector.shape[0]));ii=np.argsort(e);e=e[ii];v=np.zeros(1<<L,complex);v[ids]=vs[:,ii[0]]
    mag=[]
    for op in (X,Y,Z):
        m=[]
        for j in range(L):
            psi=v.reshape(2**j,2,2**(L-j-1));rho=np.einsum('asb,atb->st',psi,psi.conj());m.append(float(np.trace(rho@op).real))
        mag.append(m)
    S=[]
    for j in range(1,L):
        sv=np.linalg.svd(v.reshape(2**j,-1),compute_uv=False);p=sv*sv;p=p[p>1e-15];S.append(float(-p@np.log(p)))
    res=dict(L=L,parity=sector,couplings=d['couplings'],E=float(e[0]),gap=float(e[1]-e[0]),
             dmrg_energy_error=(float(d['E']-e[0]) if not args.periodic else None),max_entropy_error=(float(np.max(np.abs(np.array(S)-d['S']))) if not args.periodic else None),
             max_X_error=(float(np.max(np.abs(np.array(mag[0])-d['X']))) if not args.periodic else None),max_Y_error=(float(np.max(np.abs(np.array(mag[1])-d['Y']))) if not args.periodic else None),
             max_Z_error=(float(np.max(np.abs(np.array(mag[2])-d['Z']))) if not args.periodic else None),X=mag[0],Y=mag[1],Z=mag[2],S=S,seconds=time.time()-t)
    res['periodic']=args.periodic
    res['Xmean']=float(np.mean(mag[0]));res['Ymean']=float(np.mean(mag[1]));res['N']=float(np.sum((1-np.array(mag[2]))/2))
    if args.periodic:
        for k in ('dmrg_energy_error','max_entropy_error','max_X_error','max_Y_error','max_Z_error'):res[k]=None
    Path(args.output).write_text(json.dumps(res,indent=2)+'\n');print({k:v for k,v in res.items() if k not in ['X','Y','Z','S','couplings']})
if __name__=='__main__':main()
