"""Dense Pauli/Kronecker checks, independent of the sparse ED assembly."""
import json
import time
import numpy as np
from ed import Chain, HERE
from threadpoolctl import threadpool_limits

I = np.eye(2); X = np.array([[0,1],[1,0]])
Y = np.array([[0,-1j],[1j,0]]); Z = np.diag([1,-1])


def op(L, terms):
    out = np.array([[1.]])
    for j in reversed(range(L)):
        out = np.kron(out, terms.get(j,I))
    return out


def literal(L,h,periodic):
    H = sum(-h*op(L,{j:Z}) for j in range(L)).astype(complex)
    for j in range(L if periodic else L-1):
        k=(j+1)%L
        H += 9/16*op(L,{j:X,k:Y})+7/16*op(L,{j:Y,k:X})+op(L,{j:Z,k:Z})
    return H


def run():
    start=time.monotonic(); checks=[]
    for L in (4,6,8):
        for periodic in (False,True):
            h=1.7466; H=literal(L,h,periodic)
            P=op(L,{j:Z for j in range(L)})
            A=op(L,{j:Z for j in range(1,L,2)})
            C=op(L,{j:X for j in range(L)})
            F=op(L,{j:X if j%2==0 else Y for j in range(L)})
            errors=dict(parity=float(np.max(abs(P@H@P-H))),
                        antiunitary_A=float(np.max(abs(A@H.conj()@A-H))),
                        h_reflection=float(np.max(abs(C@H.conj()@C-literal(L,-h,periodic)))))
            errors['unitary_field_reversal']=float(np.max(abs(F@H@F-literal(L,-h,periodic))))
            # Independent Jordan-Wigner matrices, with both open-end degrees.
            cs=[op(L,{**{k:Z for k in range(j)},j:(X+1j*Y)/2}) for j in range(L)]
            ns=[c.conj().T@c for c in cs]
            HF=np.eye(1<<L,dtype=complex)*((L if periodic else L-1)-h*L)
            for j,n in enumerate(ns):
                degree=2 if periodic or 0<j<L-1 else 1
                HF+=(2*h-2*degree)*n
            for j in range(L if periodic else L-1):
                k=(j+1)%L; cj,ck=cs[j],cs[k]
                pair=cj.conj().T@ck.conj().T
                transverse=1j/8*(ck.conj().T@cj-cj.conj().T@ck)+1j*(pair-pair.conj().T)
                if periodic and j==L-1: transverse=-P@transverse
                HF+=transverse+4*ns[j]@ns[k]
            errors['jordan_wigner_full_field_and_parity_seam']=float(np.max(abs(HF-H)))
            refl=np.zeros_like(H)
            for s in range(1<<L):
                dest=sum(((s>>j)&1)<<(L-1-j) for j in range(L))
                refl[dest,s]=1
            R=op(L,{j:np.diag([np.exp(-1j*np.pi/4),np.exp(1j*np.pi/4)]) for j in range(L)})
            B=refl@R
            errors['antiunitary_inversion']=float(np.max(abs(B@H.conj()@B.conj().T-H)))
            for parity in (0,1):
                c=Chain(L,parity,periodic)
                sub=H[np.ix_(c.states,c.states)]
                rotated=c.phase.conj()[:,None]*sub*c.phase[None,:]
                errors[f'real_frame_P{parity}']=float(np.max(abs(rotated-c.H(h).toarray())))
                rng=np.random.default_rng(L+parity)
                v=rng.normal(size=c.n)+1j*rng.normal(size=c.n); v/=np.linalg.norm(v)
                for name,mat in [('kappa',op(L,{0:X,2:Y})-op(L,{0:Y,2:X})),
                                 ('B',op(L,{0:X,2:X})-op(L,{0:Y,2:Y}))]:
                    errors[f'{name}_action_P{parity}']=float(np.max(abs(c.local(name,0,v)-mat[np.ix_(c.states,c.states)]@v)))
                if periodic:
                    T=np.zeros_like(sub); T[c.trans_dest,np.arange(c.n)]=1
                    errors[f'translation_P{parity}']=float(np.max(abs(T@sub-sub@T)))
                    As=A[np.ix_(c.states,c.states)]
                    errors[f'projective_A_T_P{parity}']=float(np.max(abs(T@As-(-1)**parity*As@T)))
            kap=op(L,{0:X,2:Y})-op(L,{0:Y,2:X})
            errors['kappa_A_odd']=float(np.max(abs(A@kap.conj()@A+kap)))
            errors['kappa_C_even']=float(np.max(abs(C@kap.conj()@C-kap)))
            errors['kappa_F_odd']=float(np.max(abs(F@kap@F+kap)))
            assert max(errors.values())<2e-12, errors
            checks.append(dict(L=L,periodic=periodic,errors=errors))
    out=dict(checks=checks,max_error=max(max(x['errors'].values()) for x in checks),
             runtime_seconds=time.monotonic()-start)
    (HERE/'results'/'exact_checks.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':
    with threadpool_limits(limits=1):
        run()
