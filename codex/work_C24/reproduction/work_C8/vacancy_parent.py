"""Local vacancy extension of a frustration-free open-chain code.

The demo uses the symmetric three-qubit ferromagnetic code. It does not
provide a high-channel Clifford code or a critical microscopic chain.
"""
import itertools
import json
import math
import numpy as np
from scipy.sparse import csr_matrix, kron


def charge_hamiltonian(k, p, U, w):
    N = k+p
    states = list(itertools.combinations(range(N), p))
    index = {s:i for i, s in enumerate(states)}
    H = np.zeros((len(states), len(states)))
    for i, holes in enumerate(states):
        occ = [int(j in holes) for j in range(N)]
        H[i,i] = U*sum(occ[j] != occ[j+1] for j in range(N-1))
        for j in range(N-1):
            if occ[j] != occ[j+1]:
                new = list(holes)
                old = j if occ[j] else j+1
                new[new.index(old)] = j+1 if occ[j] else j
                H[index[tuple(sorted(new))], i] -= w
    e, v = np.linalg.eigh(H)
    chi = v[:,0]
    if chi.sum() < 0:
        chi = -chi
    # The tiny even/odd tunneling splitting makes an unconstrained eigenvector
    # sensitive to roundoff. Project to the rigorously known even sector.
    reflect = [index[tuple(sorted(N-1-j for j in s))] for s in states]
    chi = chi + chi[reflect]
    chi /= np.linalg.norm(chi)
    assert np.linalg.norm(H@chi-e[0]*chi) < 1e-12
    assert chi.min() > 0
    return states, H, chi, e


def symmetric_code(k):
    V = np.zeros((2**k, k+1))
    for i in range(2**k):
        n = i.bit_count()
        V[i,n] = 1/np.sqrt(math.comb(k,n))
    return V


def lifted_code(k, holes, chi, V):
    N = k+len(holes[0])
    rows, cols, vals = [], [], []
    physical_states = []
    for a, hs in enumerate(holes):
        for z in range(2**k):
            bits = np.unravel_index(z, (2,)*k)
            it = iter(bits)
            letters = tuple(0 if j in hs else next(it)+1 for j in range(N))
            physical_states.append(letters)
            row = np.ravel_multi_index(letters, (3,)*N)
            for d in range(V.shape[1]):
                if V[z,d] != 0:
                    rows.append(row); cols.append(d); vals.append(chi[a]*V[z,d])
    W = csr_matrix((vals,(rows,cols)), shape=(3**N,V.shape[1]))
    return W, physical_states


def physical_hamiltonian(states, U, w, g):
    index = {s:i for i, s in enumerate(states)}
    H = np.zeros((len(states), len(states)))
    for a, s in enumerate(states):
        for j in range(len(s)-1):
            x,y = s[j:j+2]
            if (x == 0) != (y == 0):
                H[a,a] += U
                target = list(s); target[j],target[j+1] = y,x
                H[index[tuple(target)],a] -= w
            elif x and y and x != y:
                # g times the singlet projector (I-swap)/2.
                H[a,a] += g/2
                target = list(s); target[j],target[j+1] = y,x
                H[index[tuple(target)],a] -= g/2
    return H


def compression(W, F, length, endpoint):
    v = W.toarray().reshape((3,)*length+(W.shape[1],))
    fv = np.tensordot(F,v,axes=(1,endpoint))
    fv = np.moveaxis(fv,0,endpoint)
    return W.conj().T@fv.reshape(W.shape)


def one_case(w):
    k=3; U=1.; g=.7
    V = symmetric_code(k); D=V.shape[1]
    X=np.array([[0,1],[1,0.]])
    Y=np.array([[0,-1j],[1j,0]])
    Z=np.diag([1.,-1.])
    paulis=[X,Y,Z]
    As=[V.T@np.kron(f,np.eye(2**(k-1)))@V for f in paulis]
    W=[]; rows=[]; chis=[]; all_holes=[]
    for p in [1,2]:
        holes,Hc,chi,ec=charge_hamiltonian(k,p,U,w)
        code,states=lifted_code(k,holes,chi,V)
        Hp=physical_hamiltonian(states,U,w,g)
        ev=np.linalg.eigvalsh(Hp)
        physical_ids=[np.ravel_multi_index(s,(3,)*(k+p)) for s in states]
        frame=code.toarray()[physical_ids,:]
        residual=np.linalg.norm(Hp@frame-ec[0]*frame)
        ground_split=ev[D-1]-ev[0]
        assert residual < 1e-10 and abs(ev[0]-ec[0]) < 1e-10 and ground_split < 1e-10
        occ=1-sum(chi[i]**2 for i,s in enumerate(holes) if 0 in s)
        comp_err=0.
        for f,A in zip(paulis,As):
            F=np.zeros((3,3),complex); F[1:,1:]=f
            for end in [0,k+p-1]:
                comp_err=max(comp_err,np.linalg.norm(compression(code,F,k+p,end)-occ*A))
        assert comp_err < 1e-10
        rows.append(dict(holes=p,length=k+p,charge_dimension=len(holes),
                         physical_sector_dimension=len(states),
                         charge_gap=float(ec[1]-ec[0]),
                         full_gap=float(ev[D]-ev[D-1]),
                         exact_code_residual=float(residual),ground_split=float(ground_split),
                         endpoint_occupation=float(occ),port_error=float(comp_err)))
        W.append(code); chis.append(chi); all_holes.append(holes)
    before=kron(W[0],W[1],format='csr'); after=kron(W[1],W[0],format='csr')
    actual=(after.conj().T@before).toarray()
    lookup={s:i for i,s in enumerate(all_holes[1])}
    overlap=sum(chis[0][i]*chis[1][lookup[(s[0],k+1)]]
                for i,s in enumerate(all_holes[0]))
    beta=float(overlap**2)
    err=np.linalg.norm(actual-beta*np.eye(D**2))
    assert err < 1e-10
    return dict(k=k,U=U,w=w,g=g,core_dimension=D,cell_checks=rows,
                boundary_transport_amplitude=beta,transport_error=float(err))


def main():
    out=[one_case(w) for w in [.15,.3,.5]]
    with open('vacancy_parent_checks.json','w') as f:
        json.dump(out,f,indent=2)
    print(json.dumps(out,indent=2))


if __name__ == '__main__':
    main()
