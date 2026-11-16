"""Reproducible algebraic and band-counting checks, not proof of c <= 2.

Only numpy is required. The projector calculation is independent of the
low-energy extraction in attack_nested.py; the final ring check compares them.
"""
import json
import time
from pathlib import Path

import numpy as np
from attack_nested import operators


def cusp_coefficients(L, V, p, g):
    """Independent bit-basis second-order elimination, retaining every P sector.

    p creates adjacent pairs; g(n_left-n_right) creates a particle.
    The particle-number-one compression is used only as a band diagnostic.
    """
    dim = 1 << L
    E = np.zeros(dim)
    perturbation = np.zeros((dim,dim),dtype=complex)
    for s in range(dim):
        ns = [(s >> j) & 1 for j in range(L)]
        E[s] = V*sum(ns[j]*ns[(j+1)%L] for j in range(L))
        for j in range(L):
            k = (j+1)%L
            if ns[j] == ns[k]:
                perturbation[s^(1 << j)^(1 << k),s] += p if ns[j] == 0 else np.conj(p)
            amp = g if ns[j] == 0 else np.conj(g)
            perturbation[s^(1 << j),s] += amp*(ns[(j-1)%L]-ns[(j+1)%L])
    P = np.flatnonzero(E == 0); Q = np.flatnonzero(E > 0)
    assert np.max(abs(perturbation-perturbation.conj().T)) == 0
    assert np.max(abs(perturbation[np.ix_(P,P)])) == 0
    A = perturbation[np.ix_(Q,P)]
    K = -A.conj().T @ (A/E[Q,None])
    pos = {s:i for i,s in enumerate(P)}
    actual = np.array([K[pos[2],pos[1]],K[pos[4],pos[1]],K[pos[8],pos[1]],
                       K[pos[1],pos[1]]-K[pos[0],pos[0]]])
    expected = np.array([abs(g)**2/V,-abs(p)**2/(2*V),0,
                         (3*abs(p)**2-2*abs(g)**2)/V])
    return float(np.max(abs(actual-expected)))


def uniform_u1_check(L):
    """Kernel of [H,sum Z] over all nine exchange and three uniform-field terms."""
    I = np.eye(2); X = np.array([[0,1],[1,0]])
    Y = np.array([[0,-1j],[1j,0]]); Z = np.diag([1,-1])
    pauli = [X,Y,Z]
    def product(ops):
        out = np.ones((1,1),complex)
        for j in range(L): out = np.kron(out,ops.get(j,I))
        return out
    basis = [sum((product({j:a,(j+1)%L:b}) for j in range(L))) for a in pauli for b in pauli]
    basis += [sum((product({j:a}) for j in range(L))) for a in pauli]
    charge = basis[-1]
    C = np.column_stack([(H@charge-charge@H).reshape(-1) for H in basis])
    real_map = np.concatenate([C.real,C.imag],axis=0)
    _,sv,vh = np.linalg.svd(real_map,full_matrices=False)
    rank = int(np.count_nonzero(sv>1e-9))
    known = np.zeros((12,4))
    known[0,0]=known[4,0]=1/np.sqrt(2)  # XX+YY
    known[1,1]=1/np.sqrt(2); known[3,1]=-1/np.sqrt(2)  # XY-YX
    known[8,2]=known[11,3]=1  # ZZ, Z
    projector_error = float(np.max(abs(vh[rank:].T@vh[rank:]-known@known.T)))
    assert rank == 8 and projector_error < 1e-12
    return {'L':L,'kernel_dimension':12-rank,'kernel_projector_error':projector_error}


def count_minima(hoppings):
    """epsilon(k)=2 Re sum_r t_r exp(i r k); polynomial roots of epsilon'."""
    R = len(hoppings)
    poly = np.zeros(2*R+1, dtype=complex)  # ascending powers of z
    for r, t in enumerate(hoppings, 1):
        poly[R+r] += 1j*r*t
        poly[R-r] -= 1j*r*np.conj(t)
    roots = np.roots(np.trim_zeros(poly, 'b')[::-1])
    angles = [float(np.angle(z)) for z in roots if abs(abs(z)-1) < 1e-7]
    minima = []
    max_residual = 0.0
    for k in angles:
        first = 2*np.real(sum(1j*r*t*np.exp(1j*r*k) for r,t in enumerate(hoppings,1)))
        second = 2*np.real(sum(-r*r*t*np.exp(1j*r*k) for r,t in enumerate(hoppings,1)))
        max_residual = max(max_residual,abs(float(first)))
        if second > 1e-7:
            minima.append(k)
    return len(minima), max_residual


def projector_bond(V, lam):
    # Stable version of (V-sqrt(V^2+8 lam^2))/4.
    w = -2*lam*lam/(V+np.sqrt(V*V+8*lam*lam))
    n = np.diag([0.,1.])
    bdag = np.array([[0.,0.],[1.,0.]])
    b = bdag.T
    X = bdag+b
    I = np.eye(2)
    mu = 2*w
    h = (V*np.kron(n,n)+lam*(np.kron(n,X)-np.kron(X,n))
         +w*(np.kron(bdag,b)+np.kron(b,bdag))
         -(mu/2)*(np.kron(n,I)+np.kron(I,n)))
    psi = np.array([0.,-np.sqrt(-w),np.sqrt(-w),lam/np.sqrt(-w)])
    return h,psi,w


def main():
    start = time.perf_counter()
    rng = np.random.default_rng(20261116)
    counts = {1:{},2:{}}
    worst_root = 0.
    for R in (1,2):
        for _ in range(10000):
            ts = rng.normal(size=R)+1j*rng.normal(size=R)
            m,err = count_minima(ts)
            counts[R][m] = counts[R].get(m,0)+1
            worst_root = max(worst_root,err)
            assert m <= R
    three, _ = count_minima([0.,0.,-1.])
    assert three == 3  # Detects a deliberately out-of-scope third-neighbor control.
    cusp_errors = []
    for L in (8,9,10):
        for _ in range(20):
            V = rng.uniform(0.5,5)
            p,g = rng.normal(size=2)+1j*rng.normal(size=2)
            cusp_errors.append(cusp_coefficients(L,V,p,g))
    assert max(cusp_errors) < 1e-12
    u1 = [uniform_u1_check(L) for L in (4,5)]
    worst_projector = 0.
    for _ in range(200):
        V = 10**rng.uniform(-1,1)
        lam = 10**rng.uniform(-2,0)
        h, psi, w = projector_bond(V,lam)
        error = np.max(abs(h-np.outer(psi,psi)))/np.max(abs(h))
        worst_projector = max(worst_projector,float(error))
        assert error < 1e-12
        # Transfer matrix in Bravyi-Gosset convention has a nontrivial Jordan block.
        T = np.array([[psi[1],psi[3]],[-psi[0],-psi[2]]])
        assert abs(np.linalg.det(T)) > 0
        assert T[0,0] == T[1,1] and T[0,1] != 0 and T[1,0] == 0
    rings = []
    V,lam = 1.,0.3
    _,_,w = projector_bond(V,lam)
    for L in range(4,11):
        E,D2,F,W = operators(L)
        # W=-T1+2N, hence H=V*nn+lam*F-w*W when mu=2w.
        H = np.diag(V*E)+lam*F-w*W
        ev = np.linalg.eigvalsh(H)
        vac = np.zeros(1 << L); vac[0]=1
        one = np.zeros(1 << L)
        for j in range(L): one[1 << j]=1/np.sqrt(L)
        res = max(float(np.max(abs(H@vac))),float(np.max(abs(H@one))))
        ground_dim = int(np.count_nonzero(abs(ev)<1e-9))
        assert ground_dim == 2 and ev[0] > -1e-10 and res < 1e-12
        # Every linear combination has Schmidt rank at most two, sampled here.
        psi = (vac+1j*one)/np.sqrt(2)
        singular = np.linalg.svd(psi.reshape(1 << (L//2),-1),compute_uv=False)
        rank = int(np.count_nonzero(singular>1e-10))
        assert rank <= 2
        rings.append({'L':L,'ground_dimension':ground_dim,'gap':float(ev[2]),
                      'ground_state_residual':res,'sample_Schmidt_rank':rank})
    result = {'seed':20261116,'numpy_version':np.__version__,
              'random_bands_per_range':10000,'minima_histograms':counts,
              'range_three_control_minima':three,'largest_stationary_residual':worst_root,
              'cusp_coefficient_cases':len(cusp_errors),'max_cusp_coefficient_error':max(cusp_errors),
              'uniform_U1_checks':u1,
              'projector_cases':200,'max_relative_projector_error':worst_projector,
              'projector_ring_checks':rings,'runtime_seconds':time.perf_counter()-start}
    Path(__file__).with_name('structure_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    main()
