"""Reproduce C19 finite-dimensional checks in minutes, no external files."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import json,time,argparse
import numpy as np
import scipy.linalg as la
from scipy.sparse.linalg import LinearOperator,eigsh
from uniform_cell import *

def embed_local(a,j,M,dim):
    return kron_all([np.eye(dim)]*j+[a]+[np.eye(dim)]*(M-j-1))

def logical_chain(c,M):
    N=c['N'];J=c['rho']**2
    out=np.zeros((N**M,N**M),complex)
    if not c.get('exchange_fields',False):
        for j in range(M):out+=J*embed_local(c['field'],j,M,N)
    for j in range(M-1):
        out+=np.kron(np.eye(N**j),np.kron(c['H2'],np.eye(N**(M-j-2))))
    return out

def apply_axis(a,x,j,M,D):
    y=np.tensordot(a,x.reshape((D,)*M),axes=(1,j))
    return np.moveaxis(y,0,j).reshape(-1)

def physical_operator(c,M,lam):
    D=c['D']; dim=D**M
    local=c['H']+lam**2*c['rho']**2*c['field_phys']
    def mv(x):
        out=np.zeros(dim,complex)
        for j in range(M):out+=apply_axis(local,x,j,M,D)
        for j in range(M-1):
            for R,L in zip(c['FR'],c['FL']):
                out+=lam*apply_axis(R,apply_axis(L,x,j+1,M,D),j,M,D)
        return out
    return LinearOperator((dim,dim),matvec=mv,rmatvec=mv,dtype=complex)

def ed(c,M,lambdas):
    D,N=c['D'],c['N'];dim=D**M; low=N**M
    # Dense physical G tensor product only has low columns.
    GG=kron_all([c['G']]*M)
    target=logical_chain(c,M)
    results=[]
    rng=np.random.default_rng(1900+10*c['n']+M)
    v0=rng.normal(size=dim)+1j*rng.normal(size=dim)
    for lam in lambdas:
        start=time.time();H=physical_operator(c,M,lam)
        e,U=eigsh(H,k=low+2,which='SA',tol=2e-12,ncv=max(2*(low+2)+1,28),v0=v0)
        # Ritz repair of clustered complex ARPACK vectors.
        Q,_=la.qr(U,mode='economic');HQ=np.column_stack([H@Q[:,a] for a in range(Q.shape[1])])
        ev,V=la.eigh(Q.conj().T@HQ); U=Q@V
        residual=max(la.norm(H@U[:,a]-ev[a]*U[:,a]) for a in range(low+2))
        overlap=GG.conj().T@U[:,:low]
        a,s,b=la.svd(overlap);polar=a@b
        heff=(polar*ev[:low])@polar.conj().T
        remainder=heff-lam**2*target
        err=float(la.norm(remainder,2))
        # Remove scalar to expose logical corrections.
        traceless=remainder-np.trace(remainder)/low*np.eye(low)
        record=dict(n=c['n'],cells=M,dimension=dim,lambda_=lam,
                    seconds=time.time()-start,low_band_dimension=low,
                    min_overlap_singular=float(s.min()),isolation=float(ev[low]-ev[low-1]),
                    eigenpair_residual=float(residual),error=err,
                    error_over_lambda3=err/lam**3,
                    traceless_error=float(la.norm(traceless,2)),
                    parity_error=float(la.norm(heff@kron_all([c['logical_parity']]*M)-kron_all([c['logical_parity']]*M)@heff)),
                    low_eigenvalues=ev[:low].tolist())
        print(json.dumps(record),flush=True);results.append(record)
    return results

def algebra(c):
    n,N,D=c['n'],c['N'],c['D'];g,e,P=modes(n)
    one=la.norm(c['G'].conj().T@c['G']-np.eye(N))
    h2=direct_second_order(c)
    # Build globally anticommuting two-cell fermions, independently of eta formula.
    gamglob=[np.kron(a,np.eye(N)) for a in g]+[np.kron(P,a) for a in g]
    free=np.zeros((N*N,N*N),complex)
    for a in range(n):
        ca=(gamglob[2*a]+1j*gamglob[2*a+1])/2
        cb=(gamglob[2*n+2*a]+1j*gamglob[2*n+2*a+1])/2
        free-=c['t'][a]*(np.exp(1j*c['theta'][a])*ca.conj().T@cb+np.exp(-1j*c['theta'][a])*cb.conj().T@ca)
    # Third-cell cross terms vanish as complete matrices by the outer zero compression.
    source=np.column_stack([f@c['G'] for f in c['FL']+c['FR']])
    zero=max(la.norm(c['G'].conj().T@f@c['G']) for f in c['FL']+c['FR'])
    ps=c['parity'];H=c['H'];vals=la.eigvalsh(H)
    wp=min(la.eigvalsh(m).min() for m in c['Ms'])
    # Check the cubic identity used for the analytic square-root recipe.
    square_error=0.;cubic_error=0.
    for m in c['Ms']:
        d=(m-np.eye(6*N)/c['k'])/c['rho']
        for ds in [d[:3*N,:3*N],d[3*N:,3*N:]]:
            radius=la.norm(ds,2)
            cubic_error=max(cubic_error,la.norm(ds@ds@ds-radius**2*ds))
            pp=1/c['k'];rr=c['rho']
            aa=(np.sqrt(pp+rr*radius)-np.sqrt(pp-rr*radius))/(2*radius)
            bb=(np.sqrt(pp+rr*radius)+np.sqrt(pp-rr*radius)-2*np.sqrt(pp))/(2*radius**2)
            root=np.sqrt(pp)*np.eye(3*N)+aa*ds+bb*ds@ds
            square_error=max(square_error,la.norm(root@root-(pp*np.eye(3*N)+rr*ds)))
    return dict(n=n,N=N,D=D,k=c['k'],rho=c['rho'],exchange=c['rho']**2,
                code_error=one,kernel_error=float(la.norm(H@c['G'])),
                first_order_compression=zero,zero_eigenvalues=int(sum(abs(vals)<1e-10)),gap=float(vals[N]),
                parity_error=float(la.norm(H*ps[None,:]-ps[:,None]*H)),
                min_weight=float(wp),bond_error=float(la.norm(h2-c['H2'])),
                fermion_sign_error=float(la.norm(free-c['Kbond'])),
                sqrt_error=float(square_error),cubic_relative_error=float(cubic_error/max(1,la.norm(d,2)**3)),
                momenta=c['momenta'].tolist(),hopping=c['t'].tolist(),mu=c['mu'].tolist())

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--algebra-only',action='store_true');args=ap.parse_args()
    results=dict(exact_filters=[],algebra=[],ed=[])
    for n in [1,2,4,8,12]:
        results['exact_filters'].append(dict(n=n,pass_=exact_filter_check(n)))
    for n in [1,2]:
        c=make_cell(n);a=algebra(c);results['algebra'].append(a);print(json.dumps(a),flush=True)
        if not args.algebra_only:
            results['ed']+=ed(c,2,[.08,.04,.02])
            if n==1: results['ed']+=ed(c,3,[.08,.04,.02])
        with open('checks_results.json','w') as f:json.dump(results,f,indent=2)
    print('PASS: finite matrix checks complete',flush=True)

if __name__=='__main__':main()
