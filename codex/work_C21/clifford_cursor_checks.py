"""Local Clifford cursor: exact Gaussian exchange, but not the desired target.

No history-state or matrix-reflection Hamiltonian is used here. A site
has eight cursor states and one data qubit. Cell cuts are prescribed.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
from pathlib import Path
import json
import numpy as np
import scipy.linalg as la
import scipy.sparse as sp
import scipy.sparse.linalg as sla

DELTA, HOP, GCOUP = 3., 1., .6
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], complex)
Y = np.array([[0, -1j], [1j, 0]], complex)
Z = np.diag([1., -1.]).astype(complex)
A = np.zeros((3, 3), complex)
A[0, 2] = A[2, 0] = 1


def kron(xs):
    out = np.ones((1, 1), complex)
    for x in xs:
        out = np.kron(out, x)
    return out


def source_check(n):
    N, m = 2**n, 2*n
    d = 3*N
    H = DELTA*np.eye(m*d, dtype=complex)
    gamma = []
    for i in range(n):
        for lane, pauli in enumerate((X, Y)):
            pos = 2*i+lane
            op = kron([I2]*i+[pauli]+[I2]*(n-i-1))
            H[pos*d:(pos+1)*d, pos*d:(pos+1)*d] += GCOUP*np.kron(A, op)
            gamma.append(kron([Z]*i+[pauli]+[I2]*(n-i-1)))
    for pos in range(m-1):
        if pos % 2 == 0:
            U = np.eye(d)
        else:
            i = pos//2
            zi = kron([I2]*i+[Z]+[I2]*(n-i-1))
            U = la.block_diag(zi, zi, np.eye(N))
        H[(pos+1)*d:(pos+2)*d, pos*d:(pos+1)*d] = -HOP*U
        H[pos*d:(pos+1)*d, (pos+1)*d:(pos+2)*d] = -HOP*U.conj().T
    P = kron([Z]*n)
    eta = [1j*g@P for g in gamma]
    zn = kron([I2]*(n-1)+[Z])
    K = la.block_diag(1j*zn, 1j*zn, np.eye(N))
    SL = np.zeros((m*d, d), complex); SL[:d] = np.eye(d)
    SR = np.zeros_like(SL); SR[-d:] = K
    ev, U = la.eigh(H)
    distinct = [float(ev[0])]
    for energy in ev[1:]:
        if energy-distinct[-1] > 1e-9:
            distinct.append(float(energy))
    left, right = U.conj().T@SL, U.conj().T@SR
    B = sum((np.kron(right[:, a*N:(a+1)*N], left[:, a*N:(a+1)*N])
             for a in range(3)), np.zeros((len(ev)**2, N*N), complex))
    H2 = -B.conj().T@(B/(ev[:, None]+ev[None, :]).ravel()[:, None])
    scalar = np.trace(H2).real/(N*N)
    pred = scalar*np.eye(N*N, dtype=complex)
    J = np.zeros((m, m))
    for a in range(m):
        for b in range(m):
            op = np.kron(eta[a], gamma[b])
            J[a, b] = np.trace(op@H2).real/(N*N)
            pred += J[a, b]*op
    # A direct resolvent test checks the claimed operator span separately.
    resolvent_errors = []
    for src, majoranas in ((SL, gamma), (SR, eta)):
        R = src.conj().T@la.solve(H+.73*np.eye(len(H)), src, assume_a='her')
        base_xz = np.kron(np.diag([1, 0, 1]), np.eye(N))
        base_y = np.kron(np.diag([0, 1, 0]), np.eye(N))
        fit = np.trace(base_xz@R)/(2*N)*base_xz+np.trace(base_y@R)/N*base_y
        for g in majoranas:
            op = np.kron(A, g)
            fit += np.trace(op@R)/(2*N)*op
        resolvent_errors.append(la.norm(R-fit, 2))
    singular = la.svdvals(J)
    # Reverse of the left filter list appears at the right root.
    gram = -.5*J[::-1, :]
    row = dict(n=n, cursor_sites=n, excitation_vertices=m,
               single_excitation_dimension=len(H), single_excitation_gap=float(ev[0]),
               distinct_accessible_energies=len(distinct),
               scalar=scalar, J=J.tolist(),
               gaussian_residual=la.norm(H2-pred, 2),
               relative_gaussian_residual=la.norm(H2-pred, 2)/singular[0],
               smallest_hopping_singular_value=float(singular[-1]),
               largest_hopping_singular_value=float(singular[0]),
               gram_symmetry_error=la.norm(gram-gram.T, 2),
               gram_least_eigenvalue=float(la.eigvalsh((gram+gram.T)/2)[0]),
               source_resolvent_errors=resolvent_errors)
    assert row['gaussian_residual'] < 1e-12
    assert row['gram_symmetry_error'] < 1e-12
    assert row['gram_least_eigenvalue'] > 0
    assert max(resolvent_errors) < 1e-12
    return row, H2


def local_cell():
    # Carrier states: 0, 0X,0Y,0Z, 1X,1Y,1Z, filler.
    h = np.kron(DELTA*np.diag([0]+[1]*7), I2)
    h[2:8, 2:8] += GCOUP*np.kron(A, X)
    h[8:14, 8:14] += GCOUP*np.kron(A, Y)
    for a in range(3):
        h[2*(1+a):2*(2+a), 2*(4+a):2*(5+a)] -= HOP*I2
        h[2*(4+a):2*(5+a), 2*(1+a):2*(2+a)] -= HOP*I2
    omega = np.array([1, 0, 0, 1], complex)/np.sqrt(2)
    bell = np.column_stack([omega]+[np.kron(p, I2)@omega for p in (X, Y, Z)])
    F4 = [bell.conj().T@np.kron(p, I2)@bell for p in (X, Y, Z)]
    ports = []
    for blocks in (((0, 1, 2, 3), (7, 4, 5, 6)),
                   ((0, 4, 5, 6), (7, 1, 2, 3))):
        port = []
        for f in F4:
            mat = np.zeros((8, 8), complex)
            for block in blocks:
                mat[np.ix_(block, block)] = f
            port.append(np.kron(mat, I2))
        ports.append(port)
    VR = np.eye(16, dtype=complex)
    VR[8:10, 8:10] = 1j*Z
    VR[10:12, 10:12] = 1j*Z
    FL, FR = ports[0], [VR@f@VR.conj().T for f in ports[1]]
    grade = np.kron(np.diag([1, -1, -1, 1, -1, -1, 1, 1]), Z)
    code = np.eye(16, dtype=complex)[:, :2]
    assert la.norm(h@code) == 0
    assert np.count_nonzero(abs(la.eigvalsh(h)) < 1e-12) == 2
    assert la.norm(h@grade-grade@h) == 0
    for fs in (FL, FR):
        for a, f in enumerate(fs):
            assert la.norm(code.conj().T@f@code) < 1e-14
            sign = -1 if a < 2 else 1
            assert la.norm(grade@f-sign*f@grade) < 1e-14
    v = sum((np.kron(r, l) for r, l in zip(FR, FL)), np.zeros((256, 256), complex))
    # Move a lane-1 cursor to the next site's lane 0. Odd syndromes
    # acquire Z on the left data qubit. The vacuum and filler never move.
    transport = np.zeros((256, 256), complex)
    for a in range(3):
        for dl in range(2):
            for dr in range(2):
                source = (2*(4+a)+dl)*16+dr
                target = dl*16+2*(1+a)+dr
                transport[target, source] = (-1)**dl if a < 2 else 1
    hinter = -HOP*(transport+transport.conj().T)
    grade2 = np.kron(grade, grade)
    assert la.norm(hinter@grade2-grade2@hinter) == 0
    assert la.norm(v@grade2-grade2@v) < 1e-12
    return h, v, code, hinter


def chain_ed(n, cells, H2):
    h, v, code, hinter = local_cell()
    sites = n*cells
    dim, N = 16**sites, 2**sites
    H0 = sp.csr_matrix((dim, dim), dtype=complex)
    V = H0.copy()
    prediction = np.zeros((N, N), complex)
    for j in range(sites):
        H0 += sp.kron(sp.eye(16**j), sp.kron(sp.csr_matrix(h), sp.eye(16**(sites-j-1))), format='csr')
    for j in range(sites-1):
        if (j+1) % n:
            H0 += sp.kron(sp.eye(16**j), sp.kron(sp.csr_matrix(hinter), sp.eye(16**(sites-j-2))), format='csr')
        else:
            V += sp.kron(sp.eye(16**j), sp.kron(sp.csr_matrix(v), sp.eye(16**(sites-j-2))), format='csr')
    for j in range(cells-1):
        prediction += np.kron(np.eye(2**(n*j)), np.kron(H2, np.eye(2**(n*(cells-j-2)))))
    G = kron([code]*sites)
    runs = []
    for lam in (.04, .02, .01):
        H = H0+lam*V
        # Dense n=1 two-cell check; sparse n=1 three-cell check.
        if dim < 1000:
            ev, vec = la.eigh(H.toarray(), subset_by_index=(0, N))
        else:
            ev, vec = sla.eigsh(H, k=N+1, which='SA', tol=1e-12,
                               v0=np.random.default_rng(219).normal(size=dim))
            order = np.argsort(ev); ev, vec = ev[order], vec[:, order]
        retained = vec[:, :N]
        retained, _ = la.qr(retained, mode='economic')
        small = retained.conj().T@(H@retained)
        ee, rot = la.eigh(small)
        retained = retained@rot
        polar = la.polar(G.conj().T@retained)[0]
        canonical = (polar*ee)@polar.conj().T
        err = la.norm(canonical-lam*lam*prediction, 2)
        residual = la.norm(H@retained-retained*ee, 2)
        assert residual < 1e-10
        runs.append(dict(lam=lam, canonical_error=err, error_over_lambda3=err/lam**3,
                         eigen_residual=residual, next_gap=float(ev[N]-ee[-1])))
    powers = [np.log(runs[j]['canonical_error']/runs[j+1]['canonical_error'])/np.log(2)
              for j in range(2)]
    assert min(powers) > 2.8 and max(powers) < 3.2
    return dict(n=n, cells=cells, full_dimension=dim, code_dimension=N,
                runs=runs, fitted_pairwise_powers=powers)


def main():
    rows, h2s = [], []
    for n in (1, 2, 3):
        row, h2 = source_check(n)
        rows.append(row); h2s.append(h2)
        print(json.dumps(row), flush=True)
    result = dict(sources=rows, ed=[])
    for n, c in ((1, 2), (1, 3), (2, 2)):
        result['ed'].append(chain_ed(n, c, h2s[n-1]))
        print(json.dumps(result['ed'][-1]), flush=True)
        Path(__file__).with_name('clifford_cursor_results.json').write_text(json.dumps(result, indent=2)+'\n')
    Path(__file__).with_name('clifford_cursor_results.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result['ed'], indent=2))


if __name__ == '__main__':
    main()
