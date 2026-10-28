"""Strict all-sectors port-idle version of the explicit C10 short-program bond.

Local order is selector s (4), clock letter a (12), workspace z (2), port p (2).
The computation alphabet consists of (a,z), hence k_c=24 and d0=192.
Propagation acts on z registers and leaves every p register idle, including
in illegal sectors. The final marker counterterm is the explicit rho=3|g|.
"""
import json
import numpy as np
from scipy import sparse as sp
from scipy.linalg import block_diag, eigvalsh
from build_bond_c10 import SEL, COMP, ALLOWED, SELECTOR_ALLOWED, unit, tensor


def make_strict_bond(J=100.,kappa=20.,K=5.,t=1.,mu=.9,Delta=1.,g=.2,nu=20.):
    ds=4
    da=12
    d=192
    I2=np.eye(2,dtype=complex)
    I4=np.eye(4,dtype=complex)
    I16=np.eye(16,dtype=complex)
    ident=sp.eye(d,format='csr',dtype=complex)
    ps={s:unit(ds,SEL.index(s)) for s in SEL}
    pa={a:unit(da,COMP.index(a)) for a in COMP}
    chi=sp.eye(ds,format='csr')-ps['#']
    def ss(a,b=None):
        return tensor(unit(ds,SEL.index(a),SEL.index(a if b is None else b)),np.eye(da),I4)
    def aa(a,b=None):
        return tensor(np.eye(ds),unit(da,COMP.index(a),COMP.index(a if b is None else b)),I4)
    def onsite(A):
        return (sp.kron(A,ident,format='csr')+sp.kron(ident,A,format='csr'))/2
    def term(row,col,D):
        out=sp.csr_matrix((d*d,d*d),dtype=complex)
        rows,cols=np.nonzero(D)
        for u,v in zip(rows,cols):
            a,b=divmod(int(u),4)
            c,e=divmod(int(v),4)
            left=tensor(chi,unit(da,COMP.index(row[0]),COMP.index(col[0])),unit(4,a,c))
            right=tensor(chi,unit(da,COMP.index(row[1]),COMP.index(col[1])),unit(4,b,e))
            out += D[u,v]*sp.kron(left,right,format='csr')
        return out

    hs=sp.csr_matrix((d*d,d*d),dtype=complex)
    for a in SEL:
        for b in SEL:
            if (a,b) not in SELECTOR_ALLOWED:
                hs += K*sp.kron(ss(a),ss(b),format='csr')
    hop=sp.kron(ss('L','X'),ss('X','R'),format='csr')
    hs += -t*(hop+hop.getH())+onsite(mu*ss('#'))
    hc=sp.csr_matrix((d*d,d*d),dtype=complex)
    for a in COMP:
        for b in COMP:
            if (a,b) not in ALLOWED:
                hc += kappa*sp.kron(aa(a),aa(b),format='csr')
    sync=tensor(ps['#'],np.eye(da)-pa['M'],I4)+tensor(chi,pa['M'],I4)
    local=kappa*sync
    Q=np.diag([0.,1.])
    # Both registers on actual markers are pinned. These are energies.
    local += nu*tensor(ps['#'],np.eye(da),np.kron(Q,I2)+np.kron(I2,Q))
    position={f'A{j}':j for j in range(1,5)}
    position.update({f'R{j}':j for j in range(1,5)})
    position.update(H0=2,H1=3,H2=2)
    for a,j in position.items():
        unused=np.kron(Q,I2) if j in (1,4) else np.kron(I2,Q)
        local += nu*tensor(chi,pa[a],unused)
    C=3*abs(g)
    local += C*tensor(chi,pa['H0']+pa['H1']+pa['H2'],I4)-C*ss('#')
    hc += onsite(local)

    # Four local qubits in two-site order (z_left,p_left,z_right,p_right).
    SWz=np.zeros((16,16),complex)
    for z1 in range(2):
        for p1 in range(2):
            for z2 in range(2):
                for p2 in range(2):
                    col=((z1*2+p1)*2+z2)*2+p2
                    row=((z2*2+p1)*2+z1)*2+p2
                    SWz[row,col]=1
    hp=sp.csr_matrix((d*d,d*d),dtype=complex)
    for a,b,W in [(('H0','R3'),('A2','H1'),SWz),(('A2','H1'),('H2','R3'),I16)]:
        cross=term(b,a,W)
        hp += (J/2)*(term(a,a,I16)+term(b,b,I16)-cross-cross.getH())
    hc += hp
    X=np.array([[0,1],[1,0]],complex)
    Y=np.array([[0,-1j],[1j,0]],complex)
    omega=np.array([1,0,0,1],complex)/np.sqrt(2)
    PO=np.outer(omega,omega.conj())
    def kk(a,b,c,e):
        return np.kron(np.kron(a,b),np.kron(c,e))
    D_Bell=3*Delta*np.kron(I2,np.kron(I4-PO,I2))
    D_ws=1.5*g*(kk(X,I2,X,I2)-kk(Y,I2,Y,I2))
    D_port=1.5*g*(kk(I2,X,X,I2)+kk(I2,Y,Y,I2))
    for pair,D in [(('A1','H0'),D_Bell),(('H0','R3'),D_ws),(('A1','H2'),D_port)]:
        hc += term(pair,pair,D)
    h=(hs+hc).tocsr()
    h.eliminate_zeros()
    meta=dict(k_c=24,q=2,d0=d,b=5,J=J,kappa=kappa,nu=nu,K=K,t=t,mu=mu,
              Delta=Delta,g=g,C=C,rho=C,bare_marker_chemical_potential=mu-C,
              illegal_shared_boundary_bound=kappa/2-C,unused_register_gap=nu)
    return h,hs,hc,hp,meta


def maximum(A):
    return float(np.max(np.abs(A.data))) if A.nnz else 0.


if __name__=='__main__':
    h,hs,hc,hp,meta=make_strict_bond()
    d=meta['d0']
    I2=np.eye(2,dtype=complex)
    X=np.array([[0,1],[1,0]],complex)
    Y=np.array([[0,-1j],[1j,0]],complex)
    portX=tensor(np.eye(4),np.eye(12),I2,X)
    portY=tensor(np.eye(4),np.eye(12),I2,Y)
    ident=sp.eye(d,format='csr',dtype=complex)
    port_errors=[]
    for F in (portX,portY):
        for Fbond in (sp.kron(F,ident,format='csr'),sp.kron(ident,F,format='csr')):
            port_errors.append(maximum(hp@Fbond-Fbond@hp))
    # Freeze unused registers and independently project the actual bond.
    clocks=[('A1','H0','R3','R4'),('A1','A2','H1','R4'),('A1','H2','R3','R4')]
    basis=[]
    import itertools
    for word in clocks:
        for p,s,l,r in itertools.product(range(2),repeat=4):
            zvals=(0,s,l,0)
            pvals=(p,0,0,r)
            physical=[(((1*12+COMP.index(a))*2+z)*2+v) for a,z,v in zip(word,zvals,pvals)]
            basis.append(tuple([0]+physical+[0]))
    actual=np.zeros((48,48),complex)
    for ai,a in enumerate(basis):
        for bi,b in enumerate(basis):
            for bond in range(5):
                if any(a[s]!=b[s] for s in range(6) if s not in (bond,bond+1)):
                    continue
                actual[ai,bi]+=hc[a[bond]*d+a[bond+1],b[bond]*d+b[bond+1]]
    SW=np.array([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]],complex)
    W=np.kron(I2,np.kron(SW,I2))
    Ud=block_diag(np.eye(16),W,W)
    lap=np.array([[1,-1,0],[-1,2,-1],[0,-1,1]],float)/2
    omega=np.array([1,0,0,1],complex)/np.sqrt(2)
    PO=np.outer(omega,omega.conj())
    D0=3*(np.eye(16)-np.kron(PO,np.eye(4)))+.3*np.kron(I2,np.kron(np.kron(X,X)-np.kron(Y,Y),I2))
    D2=.3*np.kron(np.kron(X,X)+np.kron(Y,Y),np.eye(4))
    expected=100*Ud@np.kron(lap,np.eye(16))@Ud.conj().T+block_diag(D0,np.zeros((16,16)),D2)
    projection_error=float(np.linalg.norm(actual-expected,2))
    meta.update(shape=list(h.shape),nnz=int(h.nnz),
                hermiticity_max_error=maximum(h-h.getH()),
                same_bond_selector_commutator_max_error=maximum(hs@hc-hc@hs),
                all_sector_propagation_port_commutator_errors=port_errors,
                legal_projection_error=projection_error,
                legal_ground_energy=float(eigvalsh(actual)[0]),
                scope='Actual strict two-site bond and its 48-dimensional initialized legal-cell compression. Full illegal quantum Hilbert space of this enlarged variant was not diagonalized; the positivity proof supplies its bounds.')
    assert meta['hermiticity_max_error']<1e-12
    assert meta['same_bond_selector_commutator_max_error']<1e-12
    assert max(port_errors)<1e-12
    assert projection_error<1e-11
    with open('strict_bond_verification.json','w') as f:
        json.dump(meta,f,indent=2)
    print(json.dumps(meta,indent=2))
