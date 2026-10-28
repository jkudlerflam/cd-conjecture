"""C17 small consistency checks. Python + numpy/scipy, no input or network.

Run OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 checks.py > checks_output.txt
The analytic statements and assumptions are in RESULT.md.
"""
import itertools
import numpy as np
from scipy.linalg import null_space, eigh_tridiagonal
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh


I = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.diag([1, -1]).astype(complex)
paulis = [X, Y, Z]


def tensor(items):
    ans = np.ones((1, 1), dtype=complex)
    for item in items:
        ans = np.kron(ans, item)
    return ans


def local(L, ops):
    return tensor([ops.get(j, I) for j in range(L)])


def coefficient_matrix(mats):
    a = np.array([m.reshape(-1) for m in mats]).T
    return np.concatenate([a.real, a.imag], axis=0)


def algebra_checks():
    L = 4
    bonds = [sum(local(L, {j: a, (j+1) % L: b}) for j in range(L))
             for a in paulis for b in paulis]
    fields = [sum(local(L, {j: a}) for j in range(L)) for a in paulis]
    basis = bonds + fields
    marker = local(L, {0: Z})
    constraints = coefficient_matrix([h @ marker - marker @ h for h in basis])
    ker = null_space(constraints)
    allowed = np.zeros((12, 2))
    allowed[8, 0] = allowed[11, 1] = 1
    residual = np.linalg.norm(ker @ ker.T - allowed @ allowed.T)
    assert ker.shape[1] == 2 and residual < 1e-12
    print('MARKER L=4 nullity=2 (Jzz,Bz), projector error', residual)

    def cross_component(h):
        # Change tensor order from sites (0,1,2,3) to (0,2 ; 1,3).
        a = h.reshape([2]*8).transpose(0, 2, 1, 3, 4, 6, 5, 7).reshape(4, 4, 4, 4)
        even = np.einsum('abcb->ac', a)
        odd = np.einsum('abad->bd', a)
        mat = a.reshape(16, 16)
        return (mat - np.kron(even, np.eye(4))/4
                - np.kron(np.eye(4), odd)/4 + np.trace(mat)*np.eye(16)/16)

    ker = null_space(coefficient_matrix([cross_component(h) for h in basis]))
    allowed = np.zeros((12, 3))
    allowed[9:, :] = np.eye(3)
    residual = np.linalg.norm(ker @ ker.T - allowed @ allowed.T)
    assert ker.shape[1] == 3 and residual < 1e-12
    print('SUBLATTICE L=4 nullity=3 (fields only), projector error', residual)


def classical_checks():
    rng = np.random.default_rng(1701)
    error = 0.0
    for _ in range(40):
        w = rng.normal(size=(2, 2))
        optimum = min(w[0, 0], w[1, 1], (w[0, 1]+w[1, 0])/2)
        for L in range(2, 11):
            words = np.array(list(itertools.product(range(2), repeat=L)))
            energies = np.mean(w[words, np.roll(words, -1, axis=1)], axis=1)
            error = max(error, optimum-float(np.min(energies)))
            assert np.min(energies) >= optimum-1e-12
    print('CLASSICAL 40 random directed binary costs, rings L=2..10, violation', error)


def pair_checks():
    t = 1.0
    print('PAIR relative-coordinate half-line truncated to r=1..200')
    for V, K in [(-3.0, 0.0), (-1.5, np.pi/2), (-1.0, .8*np.pi)]:
        C = 2*t*np.cos(K/2)
        lam = -C/V
        predicted = V+C*C/V
        diagonal = np.zeros(200)
        diagonal[0] = V
        val = eigh_tridiagonal(diagonal, -C*np.ones(199),
                               select='i', select_range=(0, 0), eigvals_only=True)[0]
        err = abs(val-predicted)
        assert abs(lam) < 1 and err < 1e-10
        print(f'  V={V:g} K/pi={K/np.pi:g} lambda={lam:.9f} '
              f'E_bound={val:.12f} error={err:.3g} global_two_atom_min=-4')


def parton_checks():
    charges = np.array([1, 2, 4, 7])
    states = np.array([s for s in itertools.product(range(2), repeat=4)
                       if np.dot(charges, s) == 7])
    rank = np.linalg.matrix_rank(states[1:]-states[0])
    killed = states @ np.array([2, -3, 1, 0])
    assert len(states) == 2 and rank == 1 and np.max(np.abs(killed)) == 0
    print('PARTON states', states.tolist(), 'affine_rank', rank,
          'formal_one_constraint_c=3 completed_torus_c=1')
    print('CONVERSION two bosons K_a=K_p=1: Delta=1.25; '
          'K_a=K_p=0.25: Delta=5.0 (kinematic examples, not realized fluids)')


def spin_hamiltonian(L, model):
    s = np.arange(1 << L, dtype=np.int64)
    rows, cols, vals = [], [], []
    if model == 'Ising':
        diagonal = np.zeros(len(s))
        for j in range(L):
            diagonal -= 1 - 2*((s >> j) & 1)
        rows.append(s)
        cols.append(s)
        vals.append(diagonal)
    for j in range(L):
        k = (j+1) % L
        if model == 'XX':
            keep = ((s >> j) & 1) != ((s >> k) & 1)
            source = s[keep]
        else:
            source = s
        rows.append(source ^ ((1 << j) | (1 << k)))
        cols.append(source)
        vals.append(-np.ones(len(source)))
    return coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                      shape=(1 << L, 1 << L)).tocsr()


def entropy(psi, L, ell):
    # Ordering of the contiguous interval is immaterial to singular values.
    singular = np.linalg.svd(psi.reshape(1 << ell, 1 << (L-ell)), compute_uv=False)
    probabilities = singular**2
    probabilities = probabilities[probabilities > 1e-15]
    return -np.sum(probabilities*np.log(probabilities))


def c_calibrations():
    rng = np.random.default_rng(1702)
    sizes = np.array([6, 8, 10, 12, 14])
    for model, target in [('XX', 1.0), ('Ising', .5)]:
        energies = []
        print('CALIBRATION', model, 'v=2, full spin PBC ED, not a new candidate')
        for L in sizes:
            h = spin_hamiltonian(int(L), model)
            e, vectors = eigsh(h, k=1, which='SA', v0=rng.normal(size=h.shape[0]),
                              tol=1e-12)
            psi = vectors[:, 0]
            expected = -2/np.sin(np.pi/L if model == 'XX' else np.pi/(2*L))
            assert abs(e[0]-expected) < 1e-10
            residual = np.linalg.norm(h @ psi-e[0]*psi)
            ells = np.arange(2, L//2+1)
            entropies = np.array([entropy(psi, int(L), int(ell)) for ell in ells])
            chord = np.log(L/np.pi*np.sin(np.pi*ells/L))
            c_s = 3*np.polyfit(chord, entropies, 1)[0]
            energies.append(e[0])
            print(f'  L={L:2d} E0={e[0]:.12f} c_S={c_s:.8f} residual={residual:.3g}')
        # E/L = e_inf + a/L^2 + b/L^4, a=-pi*c*v/6.
        design = np.array([np.ones(len(sizes)), 1/sizes**2, 1/sizes**4]).T
        fit = np.linalg.lstsq(design, np.array(energies)/sizes, rcond=None)[0]
        c_cas = -6*fit[1]/(np.pi*2)
        print(f'  all-size c_Cas={c_cas:.9f}, expected={target:g}')
        assert abs(c_cas-target) < 0.003
        assert abs(c_s-target) < 0.04


if __name__ == '__main__':
    algebra_checks()
    classical_checks()
    pair_checks()
    parton_checks()
    c_calibrations()
    print('ALL CHECKS PASSED')
