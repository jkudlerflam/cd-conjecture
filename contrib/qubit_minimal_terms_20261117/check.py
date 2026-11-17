"""Symbolic bond identities and independent complete finite-spin checks.

No particle-number restriction, entanglement fit or inferred CFT count.
"""
import json
import time
from itertools import product
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.spatial.transform import Rotation

I=np.eye(2,dtype=complex)
X=np.array([[0,1],[1,0]],complex)
Y=np.array([[0,-1j],[1j,0]],complex)
Z=np.diag([1.,-1.]).astype(complex)


def kron_all(items):
    out=np.ones((1,1),complex)
    for item in items:out=np.kron(out,item)
    return out


def site_ops(L):
    return [[kron_all([p if k==j else I for k in range(L)]) for p in (X,Y,Z)] for j in range(L)]


def spin_hamiltonian(ops,a,b,c,h,periodic,closing_sign=1):
    L=len(ops);out=np.zeros_like(ops[0][0])
    for j in range(L if periodic else L-1):
        k=(j+1)%L;sgn=closing_sign if k==0 else 1
        x,y,z=ops[j];u,v,w=ops[k]
        out+=sgn*(a*(x@u-y@v)+b*(x@v-y@u))+c*z@w
    fields=np.full(L,h) if np.isscalar(h) else h
    for j in range(L):out+=fields[j]*ops[j][2]
    return out


def xxz(ops,jxy,jz,fields,periodic,closing_sign=1):
    out=np.zeros_like(ops[0][0]);L=len(ops)
    for j in range(L if periodic else L-1):
        k=(j+1)%L;sgn=closing_sign if k==0 else 1
        x,y,z=ops[j];u,v,w=ops[k]
        out+=sgn*jxy*(x@u+y@v)+jz*z@w
    for j in range(L):out+=fields[j]*ops[j][2]
    return out


def norm(A):return float(np.max(np.abs(A)))


def symbolic_checks():
    x=sp.Matrix([[0,1],[1,0]]);y=sp.Matrix([[0,-sp.I],[sp.I,0]]);z=sp.diag(1,-1);ii=sp.eye(2)
    k=sp.kronecker_product;a,b,c,h,l=sp.symbols('a b c h lambda',real=True)
    hb=a*(k(x,x)-k(y,y))+b*(k(x,y)-k(y,x))+c*k(z,z)
    P=(x+y)/sp.sqrt(2);Q=(x-y)/sp.sqrt(2);W=[ii,P,z,Q]
    target=b*(k(x,x)-k(y,y))+a*(k(x,y)-k(y,x))-c*k(z,z)
    def zero(A):return all(sp.simplify(v)==0 for v in A)
    for j in range(4):
        u=k(W[j],W[(j+1)%4]);assert zero(u*hb*u.H-target)
        assert zero(W[j]*z*W[j].H-(-1)**j*z)
    u=sp.diag(sp.exp(-sp.I*sp.pi/8),sp.exp(sp.I*sp.pi/8));U=k(u,u)
    assert zero(U*hb*U.H-((a+b)*k(x,y)+(a-b)*k(y,x)+c*k(z,z)))
    zplus=k(z,ii)+k(ii,z);zminus=k(z,ii)-k(ii,z)
    assert zero((hb*zplus-zplus*hb).subs(a,0))
    assert zero((hb*zminus-zminus*hb).subs(b,0))
    J=sp.Matrix([[a,b,0],[-b,-a,0],[0,0,c]]);S=(J+J.T)/2;A=(J-J.T)/2
    assert sp.simplify((sp.trace(S*S)-sp.trace(S)**2)/2-a*a)==0
    assert sp.simplify(sp.trace(A.T*A)/2-b*b)==0
    assert sp.factor((l*sp.eye(3)-J*J.T).det()-(l-(a+b)**2)*(l-(a-b)**2)*(l-c*c))==0
    return {'four_symbolic_duality_bonds':True,'common_rotation':True,
            'uniform_and_staggered_charge_commutators':True,'normal_invariants_and_singular_values':True}


def majoranas(ops):
    prefix=np.eye(ops[0][0].shape[0],dtype=complex);aa=[];bb=[]
    for x,y,z in ops:
        aa.append(prefix@x);bb.append(prefix@y);prefix=prefix@z
    return aa,bb


def matrix_checks():
    rows=[]
    for L,periodic in ((4,False),(6,False),(4,True),(6,True),(8,True)):
        ops=site_ops(L);dim=1<<L;a,b,c,h=.7,1.1,-.4,.3
        seam=1 if not periodic or L%4==0 else -1
        G=kron_all([np.diag(np.exp(-.5j*(np.pi/4-j*np.pi/2)*np.array([1.,-1.]))) for j in range(L)])
        source=spin_hamiltonian(ops,0,b,c,h,periodic)
        spiral=norm(G@source@G.conj().T-xxz(ops,b,c,np.full(L,h),periodic,seam))
        V=kron_all([I if j%2==0 else X for j in range(L)])
        source=spin_hamiltonian(ops,a,0,c,h,periodic)
        staggered=norm(V@source@V.conj().T-xxz(ops,a,-c,h*(-1.)**np.arange(L),periodic))
        ww=[I,(X+Y)/np.sqrt(2),Z,(X-Y)/np.sqrt(2)]
        W=kron_all([ww[j%4] for j in range(L)])
        source=spin_hamiltonian(ops,a,b,c,h,periodic)
        target=spin_hamiltonian(ops,b,a,-c,h*(-1.)**np.arange(L),periodic,seam)
        dual=norm(W@source@W.conj().T-target)
        source0=spin_hamiltonian(ops,a,b,c,0,periodic)
        aa,bb=majoranas(ops);As=[(u+v)/np.sqrt(2) for u,v in zip(aa,bb)];Bs=[(u-v)/np.sqrt(2) for u,v in zip(aa,bb)]
        PA=(1j)**(L//2)*np.eye(dim,dtype=complex);PB=PA.copy()
        for A0,B0 in zip(As,Bs):PA=PA@A0;PB=PB@B0
        parity_error=max(norm(source0@PA-PA@source0),norm(source0@PB-PB@source0))
        data={'L':L,'periodic':periodic,'closing_xy_sign':seam,
              'spiral_identity_error':spiral,'staggered_identity_error':staggered,
              'duality_with_field_and_seam_error':dual,'separate_majorana_parity_error':parity_error}
        if not periodic:
            AT=np.zeros((dim,dim),complex)
            for j in range(L-1):
                ea=1j*As[j]@As[j+1];eb=1j*Bs[j]@Bs[j+1]
                AT+=-(a+b)*ea-(b-a)*eb-c*ea@eb
            for A0,B0 in zip(As,Bs):AT+=1j*h*A0@B0
            data['full_majorana_matrix_error']=norm(AT-source)
        if periodic and L%4==2:
            wrong=spin_hamiltonian(ops,b,a,-c,h*(-1.)**np.arange(L),True)
            data['error_if_twist_is_omitted']=norm(W@source@W.conj().T-wrong)
            assert data['error_if_twist_is_omitted']>1
        assert max(spiral,staggered,dual,parity_error)<1e-11
        if not periodic:assert data['full_majorana_matrix_error']<1e-11
        rows.append(data)
    return rows


def free_spectra():
    rows=[];L=6;ops=site_ops(L)
    for a,b,h in ((1.,.5,.3),(.5,1.,.5),(1.,1.,0.),(0.,1.,2.),(1.,0.,0.)):
        M=np.zeros((2*L,2*L))
        def add(i,j,value):M[i,j]+=value;M[j,i]-=value
        for j in range(L):add(2*j,2*j+1,-2*h)
        for j in range(L-1):
            for p,q,v in ((2*j+1,2*j+2,-2*a),(2*j,2*j+3,-2*a),
                          (2*j,2*j+2,-2*b),(2*j+1,2*j+3,-2*b)):add(p,q,v)
        energies=np.linalg.eigvalsh(1j*M)[L:]
        expected=np.sort([np.dot(energies,np.array(bits)-.5) for bits in product((0,1),repeat=L)])
        actual=np.linalg.eigvalsh(spin_hamiltonian(ops,a,b,0.,h,False))
        error=norm(actual-expected);assert error<1e-10
        rows.append({'L':L,'a':a,'b':b,'h':h,'complete_open_spectrum_error':error})
    return rows


def charge_nullities():
    basis=[sp.Matrix([1,0,0]),sp.Matrix([0,1,0]),sp.Matrix([0,0,1])]
    def cross(v):return sp.Matrix([[0,-v[2],v[1]],[v[2],0,-v[0]],[-v[1],v[0],0]])
    rows=[];L=4
    for a,b,c,h,expected in ((0,2,3,1,1),(2,0,3,1,1),(1,2,4,0,0),(1,2,3,0,0),(1,2,1,0,0),(1,1,2,0,0),(1,1,0,0,0)):
        J=sp.Matrix([[a,b,0],[-b,-a,0],[0,0,c]]);B=sp.Matrix([0,0,h]);columns=[]
        for site in range(L):
            for v in basis:
                col=[]
                for j in range(L):
                    block=sp.zeros(3)
                    if j==site:block+=cross(v)*J
                    if (j+1)%L==site:block-=J*cross(v)
                    col+=list(block)
                for j in range(L):col+=list(v.cross(B) if j==site else sp.zeros(3,1))
                columns.append(col)
        mat=sp.Matrix.hstack(*[sp.Matrix(col) for col in columns]);nullity=3*L-mat.rank()
        assert nullity==expected
        rows.append({'L':L,'a':a,'b':b,'c':c,'h':h,'exact_additive_charge_dimension':nullity})
    return rows


def invariants_and_dimensions():
    a,b,c,h=.5,1/160,10.,-19.9760472506
    J=np.array([[a,b,0],[-b,-a,0],[0,0,c]])
    B=np.array([0.,0.,h]);errors=[]
    for R in Rotation.random(12,random_state=20261117).as_matrix():
        T=R@J@R.T;S=(T+T.T)/2;A=(T-T.T)/2;F=R@B
        recovered=np.array([np.trace(S),(np.trace(S@S)-np.trace(S)**2)/2,np.sum(A*A)/2,F@F])
        errors.append(norm(recovered-[c,a*a,b*b,h*h]))
    assert max(errors)<1e-10
    ks=[.55,.75,1.,1.5,1.95];operators=[]
    for K in ks:
        lattice=[(p*p*K+q*q/K,p,q) for p in range(-5,6) for q in range(-5,6)
                 if (p or q) and (p+q)%2==0]
        minimum=min(t[0] for t in lattice);expected=min(K+1/K,4*K,4/K)
        assert abs(minimum-expected)<1e-12
        operators.append({'K':K,'minimum_allowed_vertex_dimension':minimum})
    return {'maximum_rotation_invariant_error':max(errors),'reference_ratio_a_over_b':a/b,
            'reference_ratio_c_over_b':c/b,'original_literal_spread':abs(h)/b,
            'rotated_literal_spread':abs(h)/abs(a-b),'operator_lattice_checks':operators}


def main():
    start=time.perf_counter()
    result={'symbolic':symbolic_checks(),'full_spin_matrix_checks':matrix_checks(),
            'free_many_body_spectra':free_spectra(),'exact_charge_nullities':charge_nullities(),
            'invariants_and_operator_dimensions':invariants_and_dimensions(),
            'numpy_version':np.__version__,'sympy_version':sp.__version__}
    result['runtime_seconds']=time.perf_counter()-start
    Path(__file__).with_name('checks.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
