"""Exact-parity version of the C9 spectral gadget; NumPy only.

All physical cells below are arbitrary finite matrices, not microscopic
nearest-neighbour realizations.  Full-resolvent checks use two and three
cells; no central charge is estimated.
"""
import json
import time
from pathlib import Path

import numpy as np


I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], complex)
Y = np.array([[0, -1j], [1j, 0]], complex)
Z = np.diag([1., -1.]).astype(complex)


def norm(a):
    return float(np.linalg.norm(a))


def null(a, tol=1e-11):
    _, s, vh = np.linalg.svd(a, full_matrices=True)
    rank = np.count_nonzero(s > tol)
    return vh[rank:].conj().T


def kron_all(xs):
    out = np.ones((1, 1), complex)
    for x in xs:
        out = np.kron(out, x)
    return out


def spectral_weights(As, Bs, parities, coupling, P):
    N = len(P)
    m = len(As)
    k = m + 2
    energies = 4.**np.arange(k)
    p = np.ones(k)/k
    C = 1/(energies[:, None] + energies[None, :])
    B = null(np.array([np.ones(k), C @ p]))
    vals, vecs = np.linalg.eigh(B.T @ C @ B)
    Zc = B @ ((vecs/np.sqrt(vals)) @ vecs.T)
    Te = np.diag([1., 0., 0.])
    To = np.array([[0., 1., 0.], [1., 0., 0.], [0., 0., 0.]])/np.sqrt(2)
    Ts = [Te, To]
    if any(abs(coupling[i, j]) > 1e-13 and parities[i] != parities[j]
           for i in range(m) for j in range(m)):
        raise ValueError('Target is not globally parity even')
    R = [sum(Zc[e, i]*np.kron(Ts[parities[i]], As[i])
             for i in range(m)) for e in range(k)]
    L = [sum((Zc @ coupling)[e, j]*np.kron(Ts[parities[j]], Bs[j])
             for j in range(m)) for e in range(k)]
    largest = max(max(abs(np.linalg.eigvalsh(a))) for a in R + L)
    rho = min(p)/(2*largest)
    weights = []
    for e in range(k):
        ml = p[e]*np.eye(3*N) - rho*L[e]
        mr = p[e]*np.eye(3*N) + rho*R[e]
        w = np.zeros((6*N, 6*N), complex)
        w[:3*N, :3*N] = ml
        w[3*N:, 3*N:] = mr
        weights.append(w)
    target = sum(coupling[i, j]*np.kron(As[i], Bs[j])
                 for i in range(m) for j in range(m))
    scalar = -3*(p @ C @ p)
    expected = scalar*np.eye(N*N) + rho*rho*target
    direct = np.zeros((N*N, N*N), complex)
    for e in range(k):
        for f in range(k):
            for a in range(3):
                for b in range(3):
                    mr = weights[e][(3+a)*N:(4+a)*N, (3+b)*N:(4+b)*N]
                    ml = weights[f][a*N:(a+1)*N, b*N:(b+1)*N]
                    direct -= C[e, f]*np.kron(mr, ml)
    sp = np.tile(np.repeat([1., -1., -1.], N)*np.tile(np.diag(P).real, 3), 2)
    checks = dict(
        N=N, m=m, k=k, rho=float(rho), rho_squared=float(rho*rho),
        scalar=float(scalar), min_weight_eigenvalue=float(min(
            min(np.linalg.eigvalsh(w)) for w in weights)),
        weight_sum_error=norm(sum(weights)-np.eye(6*N)),
        source_parity_commutator=max(norm(sp[:, None]*w-w*sp[None, :]) for w in weights),
        C_metric_error=norm(Zc.T @ C @ Zc-np.eye(m)),
        linear_cancellation_error=norm(p @ C @ Zc),
        target_parity_commutator=norm(np.kron(P,P)@target-target@np.kron(P,P)),
        weight_formula_error=norm(direct-expected),
        weight_formula_error_relative_to_exchange=norm(direct-expected)/(rho*rho),
    )
    return weights, energies, sp, expected, checks


def physical_cell(weights, energies, sp, P):
    N = len(P)
    k = len(energies)
    source_dim = 6*N
    D = 4*((max(16*N, (6*k+1)*N)+3)//4)
    aux = D//4
    G = np.zeros((D, N), complex)
    aux_parity = np.ones(aux)
    for a in range(2):
        for b in range(2):
            for i in range(N):
                idx = (2*a+b)*N+i
                G[(2*a+b)*aux+idx, i] = .5
                aux_parity[idx] = (-1)**(a+b)*P[i, i].real
    up = np.kron(np.kron([1., -1.], [1., -1.]), aux_parity)
    left = [np.kron(np.kron(f, I2), np.eye(aux)) for f in (Z, X, Y)]
    right = [np.kron(np.kron(I2, f), np.eye(aux)) for f in (Z, X, Y)]
    S = np.hstack([f @ G for f in left+right])
    roots = []
    for w in weights:
        vals, vecs = np.linalg.eigh(w)
        roots.append((vecs*np.sqrt(np.maximum(vals,0))) @ vecs.conj().T)
    W = np.vstack(roots)
    dp = np.tile(sp, k)
    J = np.zeros((D, k*source_dim), complex)
    for sign in [1., -1.]:
        phys = np.flatnonzero(up == sign)
        dil = np.flatnonzero(dp == sign)
        src = np.flatnonzero(sp == sign)
        code = np.flatnonzero(np.diag(P).real == sign)
        ws = W[np.ix_(dil, src)]
        ss = S[np.ix_(phys, src)]
        gs = G[np.ix_(phys, code)]
        wb = np.hstack([ws, null(ws.conj().T)])
        sb = np.hstack([ss, null(np.hstack([gs, ss]).conj().T)])
        J[np.ix_(phys, dil)] = sb[:, :len(dil)] @ wb.conj().T
    projectors = [J[:, e*source_dim:(e+1)*source_dim] @
                  J[:, e*source_dim:(e+1)*source_dim].conj().T
                  for e in range(k)]
    remainder = np.eye(D)-G @ G.conj().T-J @ J.conj().T
    H = sum(e*p for e,p in zip(energies, projectors)) + 4*energies[-1]*remainder
    H = (H+H.conj().T)/2
    eig, U = np.linalg.eigh(H)
    checks = dict(
        D=D, parity_plus=int(sum(up>0)), parity_minus=int(sum(up<0)),
        G_isometry_error=norm(G.conj().T @ G-np.eye(N)),
        syndrome_isometry_error=norm(S.conj().T @ S-np.eye(source_dim)),
        code_intertwining_error=norm(up[:, None]*G-G@P),
        dilation_intertwining_error=norm(up[:, None]*J-J*dp[None,:]),
        realization_error=max(norm(S.conj().T @ pe @ S-w) for pe,w in zip(projectors,weights)),
        cell_parity_commutator=norm(up[:,None]*H-H*up[None,:]),
        cell_ground_error=norm(H@G), cell_gap=float(eig[N]),
        first_order_port_norm=max(norm(G.conj().T @ f @ G) for f in left+right),
    )
    return H,G,left,right,eig,U,checks


def apply_axis(x, matrix, axis):
    return np.moveaxis(np.tensordot(matrix, x, axes=(1, axis)), 0, axis)


def cluster_check(cell, expected, count):
    H,G,left,right,energies,U,_ = cell
    D,N = G.shape
    GG = kron_all([G]*count)
    def perturb(x):
        out = np.zeros_like(x)
        for bond in range(count-1):
            for a,b in zip(right,left):
                out += apply_axis(apply_axis(x,a,bond),b,bond+1)
        return out
    vgs = np.stack([perturb(GG[:,i].reshape((D,)*count)).ravel()
                    for i in range(N**count)], axis=1)
    den = np.zeros((D,)*count)
    for axis in range(count):
        shape = [1]*count
        shape[axis] = D
        den += energies.reshape(shape)
    inv = np.zeros_like(den)
    inv[den>.25] = 1/den[den>.25]
    rgs = []
    for i in range(N**count):
        x = vgs[:,i].reshape((D,)*count)
        for axis in range(count):
            x = apply_axis(x,U.conj().T,axis)
        x *= inv
        for axis in range(count):
            x = apply_axis(x,U,axis)
        rgs.append(x.ravel())
    first = GG.conj().T @ vgs
    second = -vgs.conj().T @ np.stack(rgs,axis=1)
    predicted = np.zeros_like(second)
    for bond in range(count-1):
        predicted += np.kron(np.eye(N**bond),
                             np.kron(expected,np.eye(N**(count-bond-2))))
    return dict(cells=count, physical_dimension=D**count,
                first_order_norm=norm(first),
                full_resolvent_second_order_error=norm(second-predicted))


def clifford_square_root_check(n):
    gammas=[]
    for a in range(n):
        gammas.extend([kron_all([Z]*a+[pauli]+[I2]*(n-a-1)) for pauli in [X,Y]])
    parity=kron_all([Z]*n)
    coeff=np.sin(np.arange(1,2*n+1)*.71)+.17
    vector=sum(c*g for c,g in zip(coeff,gammas))
    length=np.linalg.norm(coeff)
    p=.23
    rho=.17*np.sqrt(2)/length
    eta=length/np.sqrt(2)
    plus,minus=np.sqrt(p+rho*eta),np.sqrt(p-rho*eta)
    aa=(plus+minus)/2
    bb=(plus-minus)/(2*eta)
    ts=np.array([[0.,1.,0.],[1.,0.,0.],[0.,0.,0.]])/np.sqrt(2)
    ee=np.diag([1.,1.,0.])
    weight=p*np.eye(3*2**n)+rho*np.kron(ts,vector)
    root=np.kron(np.sqrt(p)*(np.eye(3)-ee)+aa*ee,np.eye(2**n))+bb*np.kron(ts,vector)
    eig,u=np.linalg.eigh(weight)
    generic=(u*np.sqrt(eig))@u.conj().T
    rotation=np.eye(2**n,dtype=complex)
    reduced=coeff.copy()
    for j in range(2*n-1,0,-1):
        angle=np.arctan2(reduced[j],reduced[j-1])
        step=np.cos(angle/2)*np.eye(2**n)+np.sin(angle/2)*(gammas[j-1]@gammas[j])
        rotation=step@rotation
        reduced[j-1]=np.hypot(reduced[j-1],reduced[j])
        reduced[j]=0
    report=dict(n=n,rotation_count=2*n-1,
                square_root_error=norm(root@root-weight),
                square_root_vs_eigendecomposition=norm(root-generic),
                alignment_error=norm(rotation@vector@rotation.conj().T-length*gammas[0]),
                alignment_parity_commutator=norm(rotation@parity-parity@rotation))
    assert max(report[key] for key in report if key.endswith('error') or key.endswith('commutator'))<1e-11
    return report


def main():
    start = time.time()
    output = {'scope':'Abstract cells only; no microscopic realization or CFT estimate.', 'cases':[]}
    annihilation=(X+1j*Y)/2
    rx,ry=1j*X@Z,1j*Y@Z
    theta=.37
    hopping=np.exp(1j*theta)*np.kron(annihilation.conj().T@Z,annihilation)
    hopping=hopping+hopping.conj().T
    bilinear=.5*(np.sin(theta)*(np.kron(rx,X)+np.kron(ry,Y))
                 +np.cos(theta)*(np.kron(rx,Y)-np.kron(ry,X)))
    output['fermion_hopping_identity_error']=norm(hopping-bilinear)
    assert output['fermion_hopping_identity_error']<1e-13
    output['clifford_square_root_checks']=[clifford_square_root_check(n) for n in [1,2,3,4]]
    cases = [('odd_rank_one', [X], [Y], [1], np.ones((1,1)), Z, True),
             ('mixed_parity', [X,Z], [Y,Z], [1,0], np.diag([1.,.3]), Z, True)]
    for n in [1,2,3]:
        gammas=[]
        for a in range(n):
            gammas.extend([kron_all([Z]*a+[pauli]+[I2]*(n-a-1)) for pauli in [X,Y]])
        parity=kron_all([Z]*n)
        density=sum((a+.3)*kron_all([I2]*a+[Z]+[I2]*(n-a-1)) for a in range(n))
        As=[1j*g@parity for g in gammas]+[density,np.eye(2**n)]
        Bs=gammas+[density,np.eye(2**n)]
        coupling=np.zeros((2*n+2,2*n+2))
        for a in range(n):
            theta=.37+.43*a
            coupling[2*a:2*a+2,2*a:2*a+2]=np.array(
                [[np.sin(theta),np.cos(theta)],[-np.cos(theta),np.sin(theta)]])
        coupling[-2,-1]=coupling[-1,-2]=.5
        cases.append((f'free_target_n{n}',As,Bs,[1]*(2*n)+[0,0],coupling,parity,False))
    for name,As,Bs,ps,K,P,full in cases:
        weights,energies,sp,expected,report=spectral_weights(As,Bs,ps,K,P)
        report['name']=name
        if full:
            cell=physical_cell(weights,energies,sp,P)
            report['physical_cell']=cell[-1]
            report['clusters']=[cluster_check(cell,expected,2)]
            if name=='odd_rank_one':
                report['clusters'].append(cluster_check(cell,expected,3))
            assert cell[-1]['cell_gap']>.999999
            assert cell[-1]['cell_parity_commutator']<1e-10
            assert cell[-1]['realization_error']<1e-10
            assert all(c['full_resolvent_second_order_error']<1e-10 for c in report['clusters'])
        assert report['source_parity_commutator']<1e-12
        assert report['min_weight_eigenvalue']>0
        assert report['weight_sum_error']<1e-10
        assert report['weight_formula_error_relative_to_exchange']<1e-5
        output['cases'].append(report)
    output['runtime_seconds']=time.time()-start
    path=Path(__file__).with_name('checks_results.json')
    path.write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps(output,indent=2))


if __name__=='__main__':
    main()
