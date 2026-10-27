"""An explicit sparse two-site Hamiltonian for the C10 short program.

Local basis: selector (#,L,X,R), computation (M,A1..A4,R1..R4,H0..H2),
then data (0,1). This implements the partial construction in RESULT.md.
No universal compiler is claimed. All matrices below are nearest-neighbor
and the same bond is used at every position.
"""
import json
import numpy as np
from scipy import sparse as sp
from scipy.linalg import block_diag, eigvalsh

SEL = ['#', 'L', 'X', 'R']
COMP = ['M','A1','A2','A3','A4','R1','R2','R3','R4','H0','H1','H2']
ALLOWED = [('M','A1'),('R4','M'),('A1','A2'),('A2','A3'),('A3','A4'),
           ('R1','R2'),('R2','R3'),('R3','R4'),
           ('A1','H0'),('H0','R3'),('A2','H1'),('H1','R4'),('A1','H2'),('H2','R3')]
SELECTOR_ALLOWED = [('#','L'),('#','X'),('L','L'),('L','X'),('X','R'),('X','#'),('R','R'),('R','#')]


def unit(dim, row, col=None):
    if col is None:
        col = row
    return sp.csr_matrix(([1.], ([row], [col])), shape=(dim, dim), dtype=complex)


def tensor(*xs):
    out = sp.csr_matrix([[1.]], dtype=complex)
    for x in xs:
        out = sp.kron(out, sp.csr_matrix(x), format='csr')
    return out


def make_bond(J=100., kappa=20., K=5., t=1., mu=0.9, Delta=1., g=0.2, marker_data_penalty=20., exact_counterterm=False):
    """Return (full h, selector h, computation h, metadata), all CSR matrices.

    mu is the original selector chemical potential. The bare marker
    coefficient is mu-rho. By default rho=C=3|g| is an explicit upper bound
    on the valid-cell ground energy e0 before the counterterm. This preserves
    the selector ordering by the density comparison in RESULT.md, without
    requiring a computed eigenvalue as a coupling. The optional exact
    counterterm instead uses rho=e0, which makes the valid cell energy zero.
    """
    isel, icomp, idata = np.eye(4), np.eye(12), np.eye(2)
    d = 96
    ident = sp.eye(d, format='csr', dtype=complex)
    psel = {s: unit(4, SEL.index(s)) for s in SEL}
    pc = {a: unit(12, COMP.index(a)) for a in COMP}
    chi = sp.eye(4, format='csr')-psel['#']
    chiM = pc['M']
    def s_op(a,b=None):
        if b is None: b=a
        return tensor(unit(4,SEL.index(a),SEL.index(b)),icomp,idata)
    def c_op(a,b=None):
        if b is None: b=a
        return tensor(isel,unit(12,COMP.index(a),COMP.index(b)),idata)
    def track_data(rowpair,colpair,D):
        # Product selector projectors chi on both endpoints prohibit all
        # propagation and time energies across an actual selector marker.
        out = sp.csr_matrix((d*d,d*d),dtype=complex)
        for u in range(4):
            for v in range(4):
                if D[u,v] == 0: continue
                a,b = divmod(u,2)
                c,e = divmod(v,2)
                left = tensor(chi,unit(12,COMP.index(rowpair[0]),COMP.index(colpair[0])),unit(2,a,c))
                right = tensor(chi,unit(12,COMP.index(rowpair[1]),COMP.index(colpair[1])),unit(2,b,e))
                out += D[u,v]*sp.kron(left,right,format='csr')
        return out
    def onsite_to_bond(A):
        return (sp.kron(A,ident,format='csr')+sp.kron(ident,A,format='csr'))/2

    hsel = sp.csr_matrix((d*d,d*d),dtype=complex)
    for a in SEL:
        for b in SEL:
            if (a,b) not in SELECTOR_ALLOWED:
                hsel += K*sp.kron(s_op(a),s_op(b),format='csr')
    hop = sp.kron(s_op('L','X'),s_op('X','R'),format='csr')
    hsel += -t*(hop+hop.getH()) + onsite_to_bond(mu*s_op('#'))

    hcomp = sp.csr_matrix((d*d,d*d),dtype=complex)
    for a in COMP:
        for b in COMP:
            if (a,b) not in ALLOWED:
                hcomp += kappa*sp.kron(c_op(a),c_op(b),format='csr')
    mismatch = tensor(psel['#'],np.eye(12)-chiM,idata) + tensor(chi,chiM,idata)
    onsite = kappa*mismatch
    onsite += marker_data_penalty*tensor(psel['#'],icomp,unit(2,1))

    X=np.array([[0,1],[1,0]],complex)
    Y=np.array([[0,-1j],[1j,0]],complex)
    omega=np.array([1,0,0,1],complex)/np.sqrt(2)
    PO=np.outer(omega,omega.conj())
    SW=np.array([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]],complex)
    eye4=np.eye(4)
    rules=[(('H0','R3'),('A2','H1'),SW), (('A2','H1'),('H2','R3'),eye4)]
    for a,b,W in rules:
        cross=track_data(b,a,W)
        hcomp += (J/2)*(track_data(a,a,eye4)+track_data(b,b,eye4)-cross-cross.getH())
    energy_terms=[(('A1','H0'),3*Delta*(eye4-PO)),
                  (('H0','R3'),1.5*g*(np.kron(X,X)-np.kron(Y,Y))),
                  (('A1','H2'),1.5*g*(np.kron(X,X)+np.kron(Y,Y)))]
    for pair,D in energy_terms:
        hcomp += track_data(pair,pair,D)
    C=3*abs(g)
    onsite += C*tensor(chi,pc['H0']+pc['H1']+pc['H2'],idata)

    # e0: the legal cell ground energy including the uniform per-head shift.
    P=np.kron(PO,np.eye(4))
    Ws=np.kron(np.eye(2),np.kron(SW,np.eye(2)))
    U=block_diag(np.eye(16),Ws,Ws)
    lap=np.array([[1,-1,0],[-1,2,-1],[0,-1,1]],float)/2
    data0=3*Delta*(np.eye(16)-P)+1.5*g*np.kron(np.eye(2),np.kron(np.kron(X,X)-np.kron(Y,Y),np.eye(2)))
    data2=1.5*g*np.kron(np.kron(X,X)+np.kron(Y,Y),np.eye(4))
    legal=J*U@np.kron(lap,np.eye(16))@U.conj().T+block_diag(data0,np.zeros((16,16)),data2)
    e0=float(eigvalsh(legal)[0]+C)
    rho=e0 if exact_counterterm else C
    onsite += -rho*s_op('#')
    hcomp += onsite_to_bond(onsite)
    h=(hsel+hcomp).tocsr()
    h.eliminate_zeros()
    meta=dict(d0=d,k_c=12,q=2,b=5,J=J,kappa=kappa,K=K,t=t,mu=mu,Delta=Delta,g=g,
              marker_data_penalty=marker_data_penalty,C=C,e0=e0,rho=rho,
              legal_cell_energy_after_counterterm=e0-rho,
              bare_marker_chemical_potential=mu-rho,
              illegal_cell_lower_bound=kappa/2-rho,
              conservative_illegal_cell_lower_bound=kappa/2-C)
    return h,hsel,hcomp,meta


if __name__ == '__main__':
    h,hs,hc,meta=make_bond()
    herm=h-h.getH()
    comm=hs@hc-hc@hs
    meta.update(shape=list(h.shape),nnz=int(h.nnz),
                hermiticity_max_error=float(np.max(np.abs(herm.data))) if herm.nnz else 0.,
                same_bond_selector_commutator_max_error=float(np.max(np.abs(comm.data))) if comm.nnz else 0.)
    assert meta['hermiticity_max_error']<1e-12
    assert meta['same_bond_selector_commutator_max_error']<1e-12
    with open('bond_verification.json','w') as f:
        json.dump(meta,f,indent=2)
    print(json.dumps(meta,indent=2))
