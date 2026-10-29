"""Independent checks of the new analytic identities, finite matrices only."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import json
import numpy as np
from scipy.linalg import eigh, norm
from scipy.integrate import quad_vec
from explore import coefficient, path


def filters(H, source):
    e,U=eigh(H)
    # A[p,i] = <source|p><p|i>.
    A=(source.conj()@U)[:,None]*U.conj().T
    def f(t):
        E=np.empty((len(e),len(e)))
        for p in range(len(e)):
            for q in range(len(e)):
                if abs(e[p]-e[q])<1e-12:
                    E[p,q]=t*np.exp(-e[p]*t)
                else:
                    E[p,q]=(np.exp(-e[q]*t)-np.exp(-e[p]*t))/(e[p]-e[q])
        return np.einsum('pi,qi,pq->i',A,A.conj(),E).real
    return f


def main():
    out=[]
    for s in [4,5,6]:
        H=path(s,mass=6)
        for i in range(s-2):
            H[i,i+2]=(1+1j*(2*i-s+3))/8
            H[i+2,i]=H[i,i+2].conjugate()
        J=np.eye(s)[::-1];H=(H+J@H@J)/2
        left=np.eye(s,dtype=complex)[0];left[1]=(1+1j)/2
        right=J@left
        F,M=coefficient(H,left,right)
        f=filters(H,right);g=filters(H,left)
        integral,err=quad_vec(lambda t:np.outer(f(t),g(t))/2,0,np.inf,epsabs=1e-20,epsrel=1e-11)
        gram_res=float(norm(F-integral)/norm(F))
        assert gram_res<1e-9
        sym=float(norm(J@F@J-F.T)/norm(F))
        onsite=float(norm(J@M@J-M)/max(norm(M),1e-100))
        assert max(sym,onsite)<1e-9
        pencil=[]
        for k in [.17,.71,1.39,2.47]:
            h=lambda q:M-2*F*np.exp(1j*q)-2*F.T*np.exp(-1j*q)
            pencil.append(float(norm(h(-k)-J@h(k)@J)/norm(h(k))))
        assert max(pencil)<1e-9
        out.append(dict(status='numerical_consistency_check',vertices=s,
            gram_integral_relative_residual=gram_res,quadrature_error_estimate=float(err),
            hopping_reflection_relative_residual=sym,onsite_reflection_relative_residual=onsite,
            bloch_reflection_relative_residual=max(pencil)))
    # A common source gives Hermitian C and symmetric Re C, without reflection.
    H=path(5,end_mass=3);H[0,2]=1j/3;H[2,0]=-1j/3
    source=np.array([1,1j/2,0,0,0],complex)
    F,M=coefficient(H,source,source)
    common=float(norm(F-F.T)/norm(F));assert common<1e-12
    report=dict(status='numerical_checks_of_proved_identities',reflection=out,
        equal_source_symmetry_relative_residual=common,
        note='Quadrature and eigenvalues are consistency checks, not acceptance certificates.')
    json.dump(report,open('consistency_checks.json','w'),indent=2)
    print(json.dumps(report))


if __name__=='__main__':
    main()
