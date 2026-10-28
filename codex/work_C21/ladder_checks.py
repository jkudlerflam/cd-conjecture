"""A constant-alphabet pole ladder and the failure of its source interface.

ED is an exact symmetry-sector calculation for the undecorated control.
Logical data are exact spectators, not truncated virtual degrees of freedom.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
from pathlib import Path
import json
import time
import numpy as np
import scipy.linalg as la
import scipy.sparse as sp
import scipy.sparse.linalg as sla

DELTA, HOP = 3., .4
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], complex)
Y = np.array([[0, -1j], [1j, 0]], complex)
Z = np.diag([1., -1.]).astype(complex)


def kron(xs):
    out = np.ones((1, 1), complex)
    for x in xs:
        out = np.kron(out, x)
    return out


def decorated_source(n):
    m, N = n+1, 2**n
    mix = np.zeros((3, 3), complex)
    mix[0, 2] = mix[2, 0] = 1
    d = 3*N
    H = DELTA*np.eye(m*d, dtype=complex)
    for i in range(m-1):
        xi = kron([I2]*i+[X]+[I2]*(n-i-1))
        U = la.expm(1j*.37*np.kron(mix, xi))
        H[(i+1)*d:(i+2)*d, i*d:(i+1)*d] = -HOP*U
        H[i*d:(i+1)*d, (i+1)*d:(i+2)*d] = -HOP*U.conj().T
    ev, vec = la.eigh(H)
    energies = DELTA-2*HOP*np.cos(np.pi*np.arange(1, m+1)/(m+1))
    weights = 2/(m+1)*np.sin(np.pi*np.arange(1, m+1)/(m+1))**2
    errors = []
    for a in range(m):
        cols = vec[:, a*d:(a+1)*d]
        for root in (0, m-1):
            piece = cols[root*d:(root+1)*d, :]
            errors.append(la.norm(piece@piece.conj().T-weights[a]*np.eye(d), 2))
    grade = np.kron(np.eye(m), np.kron(np.diag([-1, -1, 1]), kron([Z]*n)))
    row = dict(n=n, single_excitation_dimension=m*d, poles=m,
               energies=energies.tolist(), endpoint_weights=weights.tolist(),
               spectral_error=float(np.max(abs(ev-np.repeat(energies, d)))),
               source_weight_error=max(errors), parity_error=la.norm(H@grade-grade@H, 2))
    assert row['source_weight_error'] < 1e-12
    assert row['spectral_error'] < 1e-12
    assert row['parity_error'] < 1e-12
    return row


def local_term_check():
    # Each site is |syndrome=0,X,Y,Z> tensor |data=0,1>.
    mix = np.zeros((3, 3), complex)
    mix[0, 2] = mix[2, 0] = 1
    U = la.expm(1j*.37*np.kron(mix, X))
    T = np.zeros((64, 64), complex)
    for a in range(3):
        for b in range(3):
            for dl in range(2):
                for dr in range(2):
                    for nr in range(2):
                        source = ((a+1)*2+dl)*8+dr
                        target = dl*8+(b+1)*2+nr
                        T[target, source] = U[2*b+nr, 2*a+dr]
    Q = np.kron(np.diag([0, 1, 1, 1]), I2)
    h = DELTA/2*(np.kron(Q, np.eye(8))+np.kron(np.eye(8), Q))-HOP*(T+T.conj().T)
    grade = np.kron(np.diag([1, -1, -1, 1]), Z)
    grade2 = np.kron(grade, grade)
    ev = la.eigvalsh(h)
    assert np.count_nonzero(abs(ev) < 1e-12) == 4
    assert la.norm(h@grade2-grade2@h, 2) < 1e-12
    return dict(d_site=8, bond_dimension=64, local_kernel_dimension=4,
                least_positive_eigenvalue=float(ev[4]),
                parity_error=la.norm(h@grade2-grade2@h, 2))


def matrices(n, cells):
    m, L = n+1, (n+1)*cells
    full = 4**L
    ids = np.arange(full, dtype=np.int64)
    digits = np.array([(ids//(4**i)) % 4 for i in range(L)], dtype=np.int8)
    legal = np.bitwise_xor.reduce(digits, axis=0) == 0
    states = ids[legal]
    digits = digits[:, legal]
    d = len(states)
    lookup = np.full(full, -1, dtype=np.int64)
    lookup[states] = np.arange(d)
    index = np.arange(d)
    hrows, hcols, hvals = [index], [index], [DELTA*np.count_nonzero(digits, axis=0)]
    vrows, vcols, vvals = [], [], []
    bell = np.array([[1, 0, 0, 1], [0, 1, 1, 0],
                     [0, 1, -1, 0], [1, 0, 0, -1]], float)/np.sqrt(2)
    Fs = [bell.T@np.kron(p, I2)@bell for p in (X, Y, Z)]
    v = sum((np.kron(p, p) for p in Fs), np.zeros((16, 16), complex))
    assert np.max(abs(v.imag)) < 1e-14
    for i in range(L-1):
        j = i+1
        a, b = digits[i], digits[j]
        if (i+1) % m:
            mask = (a == 0) ^ (b == 0)
            target = states[mask]+(b[mask].astype(np.int64)-a[mask])*(4**i-4**j)
            dest = lookup[target]
            assert dest.min() >= 0
            hrows.append(dest)
            hcols.append(index[mask])
            hvals.append(np.full(len(dest), -HOP))
        else:
            for channel in (1, 2, 3):
                aa, bb = a ^ channel, b ^ channel
                target = states+(aa.astype(np.int64)-a)*4**i+(bb.astype(np.int64)-b)*4**j
                dest = lookup[target]
                assert dest.min() >= 0
                coef = v[4*aa+bb, 4*a+b].real
                vrows.append(dest)
                vcols.append(index)
                vvals.append(coef)
    H0 = sp.coo_matrix((np.concatenate(hvals), (np.concatenate(hrows), np.concatenate(hcols))),
                      shape=(d, d)).tocsr()
    V = sp.coo_matrix((np.concatenate(vvals), (np.concatenate(vrows), np.concatenate(vcols))),
                     shape=(d, d)).tocsr()
    assert sla.norm(H0-H0.T) == 0 and sla.norm(V-V.T) < 1e-12
    return H0, V


def ed(n, cells):
    started = time.monotonic()
    H0, V = matrices(n, cells)
    m = n+1
    energy = DELTA-2*HOP*np.cos(np.pi*np.arange(1, m+1)/(m+1))
    weight = 2/(m+1)*np.sin(np.pi*np.arange(1, m+1)/(m+1))**2
    a0 = -3*np.sum(weight[:, None]*weight[None, :]/(energy[:, None]+energy[None, :]))
    k2 = (cells-1)*a0
    # Independent resolvent inversion in the complete symmetry sector.
    b = V[:, 0].toarray().ravel()[1:]
    resolvent = sla.spsolve(H0[1:, 1:], b) if H0.shape[0] < 2000 else sla.cg(
        H0[1:, 1:], b, rtol=1e-13, atol=1e-14)[0]
    direct = -b@resolvent
    assert abs(direct-k2) < 1e-11
    rows = []
    seed = np.zeros(H0.shape[0]); seed[0] = 1.
    for lam in (.08, .04, .02):
        H = H0+lam*V
        ev, vec = sla.eigsh(H, k=1, which='SA', tol=1e-12, v0=seed)
        residual = la.norm(H@vec[:, 0]-ev[0]*vec[:, 0])
        err = abs(ev[0]-lam*lam*k2)
        assert residual < 1e-10
        rows.append(dict(lam=lam, energy=float(ev[0]), error=err,
                         error_over_lambda3=err/lam**3, eigen_residual=residual))
    powers = [np.log(rows[i]['error']/rows[i+1]['error'])/np.log(2) for i in range(2)]
    assert min(powers) > 2.9 and max(powers) < 3.1
    return dict(n=n, cells=cells, microscopic_sites=m*cells,
                full_syndrome_dimension=4**(m*cells), sector_dimension=H0.shape[0],
                exact_logical_multiplicity=2**(n*cells),
                k2_scalar=k2, direct_resolvent_error=abs(direct-k2),
                fitted_pairwise_powers=powers, runs=rows,
                seconds=time.monotonic()-started)


def main():
    result = dict(local_term=local_term_check(),
                  decorated_sources=[decorated_source(n) for n in (1, 2, 3)])
    result['ed'] = []
    for n, cells in ((1, 2), (1, 3), (2, 2), (2, 3), (3, 2)):
        row = ed(n, cells)
        result['ed'].append(row)
        print(json.dumps(row), flush=True)
        Path(__file__).with_name('ladder_results.json').write_text(json.dumps(result, indent=2)+'\n')


if __name__ == '__main__':
    main()
