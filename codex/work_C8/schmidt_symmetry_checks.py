"""Check symmetry representations on the connected Schmidt spaces."""
import json
import numpy as np
from scipy.linalg import null_space


def hermitian_basis(d):
    out=[]
    for j in range(d):
        for i in range(j):
            X=np.zeros((d,d),complex); X[i,j]=X[j,i]=1/np.sqrt(2); out.append(X)
            Y=np.zeros((d,d),complex); Y[i,j]=-1j/np.sqrt(2); Y[j,i]=1j/np.sqrt(2); out.append(Y)
    for j in range(1,d):
        Z=np.zeros((d,d)); Z[:j,:j]+=np.eye(j)/np.sqrt(j*(j+1)); Z[j,j]=-j/np.sqrt(j*(j+1)); out.append(Z)
    return out


def analyze(name,b,v,expected_sym,expected_moving,expected_mult):
    d=len(b); basis=hermitian_basis(d); I=np.eye(d)
    coeff=np.array([[np.trace(np.kron(A,B)@v).real for B in basis] for A in basis])
    U,s,Vt=np.linalg.svd(coeff); r=int(np.sum(s>1e-9)); s=s[:r]
    left=[sum(U[j,a]*basis[j] for j in range(len(basis))) for a in range(r)]
    eq=[]
    for T in basis:
        Q=np.kron(T,I)+np.kron(I,T)
        x=np.concatenate([(1j*(Q@v-v@Q)).reshape(-1),(1j*(T@b-b@T)).reshape(-1)])
        eq.append(np.concatenate([x.real,x.imag]))
    symmetry_coeff=null_space(np.stack(eq,axis=1),rcond=1e-10)
    generators=[sum(symmetry_coeff[j,a]*basis[j] for j in range(len(basis))) for a in range(symmetry_coeff.shape[1])]
    actions=[]; leak=0.
    for T in generators:
        action=np.array([[np.trace(A@(1j*(T@B-B@T))).real for B in left] for A in left])
        for a,B in enumerate(left):
            leak=max(leak,np.linalg.norm(1j*(T@B-B@T)-sum(action[i,a]*left[i] for i in range(r))))
        assert np.linalg.norm(action+action.T)<1e-10
        assert np.linalg.norm(action@np.diag(s*s)-np.diag(s*s)@action)<1e-10
        actions.append(action.reshape(-1))
    moving=int(np.linalg.matrix_rank(np.stack(actions,axis=1),tol=1e-9)) if actions else 0
    mult=[]
    for x in s:
        if mult and abs(x-mult[-1][0])<1e-8: mult[-1][1]+=1
        else: mult.append([float(x),1])
    multiplicities=[x[1] for x in mult]
    bound=sum(n*(n-1)//2 for n in multiplicities)
    assert len(generators)==expected_sym and moving==expected_moving,(name,len(generators),moving)
    assert multiplicities==expected_mult
    assert moving<=bound and leak<1e-10
    return dict(name=name,d=d,r=r,singular_values=s.tolist(),multiplicities=multiplicities,
                symmetry_dimension=len(generators),moving_image_dimension=moving,
                local_kernel_dimension=len(generators)-moving,dimension_bound=bound,
                torus_bound=sum(n//2 for n in multiplicities),invariance_error=float(leak))


def main():
    X=np.array([[0,1],[1,0]],complex); Y=np.array([[0,-1j],[1j,0]]); Z=np.diag([1.,-1.]); I=np.eye(2)
    XX=np.kron(X,X); YY=np.kron(Y,Y); ZZ=np.kron(Z,Z)
    zero=np.zeros((2,2)); out=[]
    out.append(analyze('classical_Ising',zero,ZZ,1,0,[1]))
    out.append(analyze('critical_transverse_Ising',X,ZZ,0,0,[1]))
    out.append(analyze('XX',zero,XX+YY,1,1,[2]))
    out.append(analyze('XXZ',zero,XX+YY+.3*ZZ,1,1,[2,1]))
    out.append(analyze('XYZ',zero,XX+.7*YY+.3*ZZ,0,0,[1,1,1]))
    out.append(analyze('Heisenberg',zero,XX+YY+ZZ,3,3,[3]))
    v4=sum(np.kron(np.kron(F,I),np.kron(F,I))+.37*np.kron(np.kron(I,F),np.kron(I,F)) for F in [X,Y])
    out.append(analyze('two_XX_unequal',np.zeros((4,4)),v4,2,2,[2,2]))
    with open('schmidt_symmetry_checks.json','w') as f: json.dump(out,f,indent=2)
    print(json.dumps(out,indent=2))

if __name__=='__main__': main()
