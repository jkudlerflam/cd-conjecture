"""Checks for boundary_rank_obstruction.md and qubit_ff_witness.md.
No continuum or central-charge claim is made by this script.
"""
import json
import numpy as np
from scipy.linalg import null_space
from scipy.sparse import csr_matrix, eye, kron
from scipy.sparse.linalg import eigsh

RNG = np.random.default_rng(821104)


def random_complex(shape):
    return RNG.normal(size=shape) + 1j*RNG.normal(size=shape)


def projector(columns):
    Q, _ = np.linalg.qr(np.column_stack(columns))
    return Q@Q.conj().T


def chain(h, q, length, sparse=False):
    if sparse:
        H = csr_matrix((q**length, q**length), dtype=complex)
        for j in range(length-1):
            H += kron(kron(eye(q**j), csr_matrix(h)), eye(q**(length-j-2)))
    else:
        H = np.zeros((q**length, q**length), complex)
        for j in range(length-1):
            H += np.kron(np.kron(np.eye(q**j), h), np.eye(q**(length-j-2)))
    return H


def low_rank_case(q, r, N):
    W, _ = np.linalg.qr(random_complex((q*q,r)))
    h = W@W.conj().T
    H = chain(h,q,N)
    es, U = np.linalg.eigh(H)
    V = U[:,es<1e-9]
    extrema_error=0.; product_error=0.
    for _ in range(6):
        F=random_complex((q,q)); F=(F+F.conj().T)/2
        fs, fv=np.linalg.eigh(F)
        A=V.conj().T@np.kron(F,np.eye(q**(N-1)))@V
        ae=np.linalg.eigvalsh(A)
        extrema_error=max(extrema_error,abs(ae[0]-fs[0]),abs(ae[-1]-fs[-1]))
        for a in [fv[:,0],fv[:,-1]]:
            letters=[a]
            for _j in range(N-1):
                constraints=np.einsum('abr,a->rb',W.conj().reshape(q,q,r),letters[-1])
                b=null_space(constraints)[:,0]
                letters.append(b)
            psi=letters[0]
            for b in letters[1:]: psi=np.kron(psi,b)
            product_error=max(product_error,np.linalg.norm(H@psi))
    assert max(extrema_error, product_error)<1e-10
    return dict(q=q,rank=r,N=N,ground_dimension=V.shape[1],
                extrema_error=float(extrema_error),product_error=float(product_error))


def real_image_rank(q,s,E):
    # An explicit rank-s vector, with distinct nonzero Schmidt weights.
    psi=np.zeros((q,E),complex)
    vals=np.arange(1,s+1,dtype=float); vals/=np.linalg.norm(vals)
    for j in range(s): psi[j,j]=vals[j]
    herm=[]
    for i in range(q):
        M=np.zeros((q,q),complex); M[i,i]=1; herm.append(M)
        for j in range(i):
            M=np.zeros((q,q),complex); M[i,j]=M[j,i]=1; herm.append(M)
            M=np.zeros((q,q),complex); M[i,j]=1j; M[j,i]=-1j; herm.append(M)
    images=np.stack([(F@psi).reshape(-1) for F in herm],axis=1)
    rank=int(np.linalg.matrix_rank(np.vstack([images.real,images.imag]),tol=1e-10))
    assert rank==2*q*s-s*s
    return dict(q=q,s=s,image_rank=rank,predicted=2*q*s-s*s)


def path_witness(q):
    eyeq=np.eye(q); v=np.ones(q)/np.sqrt(q)
    columns=[np.kron(eyeq[j],eyeq[j+1] if j<q-1 else v) for j in range(q)]
    h=np.eye(q*q)-projector(columns)
    rows=[]
    for N in [q+1,q+2]:
        H=chain(h,q,N,sparse=True)
        e=np.sort(eigsh(H,k=3,which='SA',return_eigenvectors=False,tol=1e-11))
        rows.append(dict(N=N,lowest=e.real.tolist()))
    assert abs(rows[0]['lowest'][0])<1e-9 and rows[0]['lowest'][1]>1e-7
    assert rows[1]['lowest'][0]>1e-7
    return dict(q=q,spectra=rows)


def qubit_family_cases():
    e=np.eye(2); plus=np.ones(2)/np.sqrt(2)
    examples={
        'secant_no_edge':[np.kron(e[0],plus),np.kron(e[1],(e[0]+1j*e[1])/np.sqrt(2))],
        'secant_path':[np.kron(e[0],e[1]),np.kron(e[1],plus)],
        'secant_cycle':[np.kron(e[0],e[1]),np.kron(e[1],e[0])],
        'tangent_matching':[np.kron(e[0],e[0]),np.kron(e[0],e[1])+np.kron(e[1],e[0])],
        'tangent_mismatch':[np.kron(e[0],e[1]),np.kron(e[0],e[0])+np.kron(e[1],e[1])],
        'common_factor':[np.kron(e[0],e[0]),np.kron(e[0],e[1])],
    }
    expected={'secant_no_edge':[0,0,0], 'secant_path':[1,0,0],
              'secant_cycle':[2,2,2], 'tangent_matching':[2,2,2],
              'tangent_mismatch':[0,0,0], 'common_factor':[2,2,2]}
    out=[]
    for name,columns in examples.items():
        # A common invertible map tests a nonorthogonal physical basis too.
        for trial in range(3):
            A=np.eye(2) if trial==0 else np.eye(2)+.2*random_complex((2,2))
            h=np.eye(4)-projector([np.kron(A,A)@v for v in columns])
            dims=[]
            for N in [3,4,5]:
                es=np.linalg.eigvalsh(chain(h,2,N))
                dims.append(int(np.sum(es<1e-9)))
            assert dims==expected[name],(name,dims)
            out.append(dict(family=name,trial=trial,dimensions_N3_N4_N5=dims))
    return out


def main():
    out=dict(low_rank=[low_rank_case(q,q-1,4) for q in [2,3,4]],
             real_image=[real_image_rank(q,s,q+1) for q in [2,3,4,5] for s in range(1,q+1)],
             qubit_families=qubit_family_cases(),
             sharp_witness=[path_witness(q) for q in [2,3,4]])
    with open('boundary_ff_checks.json','w') as f: json.dump(out,f,indent=2)
    print(json.dumps(out,indent=2))

if __name__=='__main__': main()
