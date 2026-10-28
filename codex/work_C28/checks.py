"""Finite checks independent of the Bloch fitting and ED comparison."""
import json
import numpy as np
from scipy import sparse as sp
from scipy.sparse.linalg import eigsh
from scipy.linalg import eigh
from bus_base import h2,data_ops
from ff_bus import literal,cell,legal_embedding,PHASES,field
from explore import coefficients,normal_order

def field_rank(n):
    # All diagonal repeated two-site terms and arbitrary diagonal endpoint terms.
    A=np.zeros((n,17**2+34))
    N=2**n
    for bits in range(N):
        data=[(bits>>(n-1-i))&1 for i in range(n)]
        zs=1-2*np.array(data)
        for j in range(n+1):
            word=[a if i<j else 4+a if i==j else 2+a for i,a in enumerate(data)]
            word.append(12 if j<n else 13)
            for a,b in zip(word[:-1],word[1:]):A[:,17*a+b]+=zs/(N*(n+1))
            A[:,289+word[0]]+=zs/(N*(n+1))
            A[:,306+word[-1]]+=zs/(N*(n+1))
    return int(np.linalg.matrix_rank(A,tol=1e-10))

def main():
    out={};h,a,b=literal()
    out['local_minima']=[float(np.linalg.eigvalsh(x)[0]) for x in (h,a,b)]
    literal_rows=[]
    for n in (1,2):
        H=sum((sp.kron(sp.eye(17**j),sp.kron(sp.csr_matrix(h),sp.eye(17**(n-j-1))),format='csr') for j in range(n)))
        H+=sp.kron(sp.csr_matrix(a),sp.eye(17**n))+sp.kron(sp.eye(17**n),sp.csr_matrix(b))
        E=legal_embedding(n);C,G,_,_=cell(n)
        if n==1:
            ev=eigh(H.toarray(),eigvals_only=True,subset_by_index=(0,2**n+1))
        else:
            ev=eigsh(H,k=2**n+2,which='SA',tol=1e-12,return_eigenvectors=False,
                     v0=np.random.default_rng(28).normal(size=H.shape[0]),ncv=60)
        ev=np.sort(ev)
        assert np.max(abs(ev[:2**n]))<1e-10
        assert ev[2**n]>.9
        literal_rows.append(dict(n=n,dimension=17**(n+1),restriction=float(np.linalg.norm((E.T@H@E).toarray()-C)),
                                 invariance=float(sp.linalg.norm(H@E-E@sp.csr_matrix(C))),low_eigenvalues=ev.tolist()))
    out['literal']=literal_rows
    out['diagonal_field_ranks']={str(n):field_rank(n) for n in range(1,9)}
    comparisons=[]
    for n in (1,2):
        K=coefficients(n,phases=PHASES,ratio=1.)
        errors=[]
        for g in (.2,.1,.05):
            exact,_=h2(n,g=g,phases=PHASES)
            err=exact-K[0]-g*g*K[2]-g**4*K[4]
            errors.append(float(np.linalg.norm(err,2)))
        ratio=coefficients(n,phases=PHASES,ratio=.45,order=2)[2]
        hfull=normal_order(K[2],2*n)[1];hasym=normal_order(ratio,2*n)[1]
        comparisons.append(dict(n=n,g6_errors=errors,halving_exponents=np.log2(np.array(errors[:-1])/errors[1:]).tolist(),
                                conversion_hopping_error=float(np.linalg.norm(hasym[:n,n:]-.45*hfull[:n,n:])),
                                asymmetric_onsite_norm=float(np.linalg.norm(hasym[:n,:n])+np.linalg.norm(hasym[n:,n:]))))
    out['resolvent_expansion']=comparisons
    with open('checks_results.json','w') as f:json.dump(out,f,indent=2)
    print(json.dumps(out,indent=2))

if __name__=='__main__':main()
