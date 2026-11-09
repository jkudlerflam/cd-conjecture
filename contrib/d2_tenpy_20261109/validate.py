"""Independent matrix construction versus TeNPy's MPO and ground-state solver."""
import json
from pathlib import Path
import numpy as np
from scipy.linalg import eigh
from tenpy.algorithms.exact_diag import ExactDiag
from run import PauliChain, couplings, run

X=np.array([[0,1],[1,0]],complex);Y=np.array([[0,-1j],[1j,0]]);Z=np.diag([1.,-1.]);I=np.eye(2)


def literal_matrix(L,p):
    H=np.zeros((2**L,2**L),complex)
    def term(ops):
        a=np.ones((1,1),complex)
        for j in range(L):a=np.kron(a,ops.get(j,I))
        return a
    d=p.get('distance',1)
    for j in range(L-d):
        for a,b,v in [(X,X,p['Jx']),(Y,Y,p['Jy']),(Z,Z,p['Jz']),
                      (X,Y,p['Dz']),(Y,X,-p['Dz']),(X,Z,p['Dy']),(Z,X,-p['Dy'])]:
            H+=v*term({j:a,j+d:b})
    for j in range(L):H+=p['hz']*term({j:Z})
    return H


if __name__=='__main__':
    results=[]
    for name in ['c183','dyadic','xx','two_xx']:
        L=8;p=couplings(name);model=PauliChain(dict(L=L,bc_MPS='finite',bc_x='open',couplings=p))
        ed=ExactDiag(model);ed.build_full_H_from_mpo()
        A=literal_matrix(L,p);shift=model.energy_offset
        Hmpo=ed.full_H.to_ndarray()
        matrix_err=float(np.max(np.abs(Hmpo-(A-shift*np.eye(2**L)))))
        en,vec=eigh(A,subset_by_index=[0,1]);v=vec[:,0];S=[]
        for cut in range(1,L):
            s=np.linalg.svd(v.reshape(2**cut,-1),compute_uv=False);p2=s*s;p2=p2[p2>1e-25]
            S.append(float(-p2@np.log(p2)))
        r,_,_=run(dict(name=name,L=L,chi=32,sweeps=12,seed=1))
        row=dict(name=name,L=L,matrix_error=matrix_err,energy_error=abs(r['E_unshifted']-en[0]),
                 entropy_error=float(np.max(np.abs(np.array(S)-r['entropy']))),gap=float(en[1]-en[0]))
        assert row['matrix_error']<1e-12 and row['energy_error']<1e-8 and row['entropy_error']<1e-6,row
        results.append(row)
    for parity,bit in [('even',0),('odd',1)]:
        name='c181153';L=8;p=couplings(name);A=literal_matrix(L,p)
        inds=[j for j in range(2**L) if j.bit_count()%2==bit]
        en,vs=eigh(A[np.ix_(inds,inds)],subset_by_index=[0,0])
        v=np.zeros(2**L,complex);v[inds]=vs[:,0];S=[]
        for cut in range(1,L):
            sv=np.linalg.svd(v.reshape(2**cut,-1),compute_uv=False);p2=sv*sv;p2=p2[p2>1e-25];S.append(float(-p2@np.log(p2)))
        r,_,_=run(dict(name=name,L=L,chi=32,sweeps=12,seed=1,parity=parity))
        row=dict(name=name,L=L,parity=parity,energy_error=abs(r['E_unshifted']-en[0]),entropy_error=float(np.max(np.abs(np.array(S)-r['entropy']))))
        assert row['energy_error']<1e-8 and row['entropy_error']<1e-6,row
        results.append(row)
    Path(__file__).with_name('validation.json').write_text(json.dumps(results,indent=2)+'\n')
    print(json.dumps(results,indent=2))
