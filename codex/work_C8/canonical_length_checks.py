"""Canonical energy and tail controls, including the physical ring quotient.

These finite tests do not measure a central charge and do not simulate
the large-alphabet C8 bond. Run with single-thread BLAS.
"""
import itertools
import json
from math import factorial
from pathlib import Path

import numpy as np
from scipy.linalg import eigh
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import eigsh


def compositions(total, count, minimum=2):
    if count == 1:
        if total >= minimum:
            yield (total,)
        return
    for first in range(minimum, total - minimum * (count - 1) + 1):
        for tail in compositions(total - first, count - 1, minimum):
            yield (first,) + tail


def atom(ell, b=3, gamma=1.):
    return 0. if ell in (b, b+1) else gamma * (1+abs(ell-b))


def open_matrix(M, L, t, b=3):
    basis = list(compositions(L, M))
    index = {state: i for i, state in enumerate(basis)}
    rows, cols, values = [], [], []
    for col, state in enumerate(basis):
        rows.append(col); cols.append(col)
        values.append(sum(atom(ell, b) for ell in state))
        for j in range(M-1):
            for direction in (-1, 1):
                nxt = list(state)
                nxt[j] += direction; nxt[j+1] -= direction
                if min(nxt[j:j+2]) >= 2:
                    rows.append(index[tuple(nxt)]); cols.append(col)
                    values.append(-t)
    return basis, coo_matrix((values, (rows, cols)),
                             shape=(len(basis), len(basis))).tocsr()


def physical_matrix(M, L, t, b=3):
    # An independent physical basis: subsets of marker positions on a ring.
    basis = []
    for marks in itertools.combinations(range(L), M):
        gaps = tuple((marks[(j+1) % M]-marks[j]) % L for j in range(M))
        if min(gaps) >= 2:
            basis.append(marks)
    index = {state: i for i, state in enumerate(basis)}
    rows, cols, values = [], [], []
    gaps_by_state = []
    for col, marks in enumerate(basis):
        gaps = tuple((marks[(j+1) % M]-marks[j]) % L for j in range(M))
        gaps_by_state.append(gaps)
        rows.append(col); cols.append(col)
        values.append(sum(atom(ell, b) for ell in gaps))
        for j in range(M):
            for direction in (-1, 1):
                nxt = list(marks); nxt[j] = (nxt[j]+direction) % L
                nxt = tuple(sorted(nxt))
                if nxt in index:
                    rows.append(index[nxt]); cols.append(col); values.append(-t)
    H = coo_matrix((values, (rows, cols)),
                   shape=(len(basis), len(basis))).tocsr()
    return basis, np.array(gaps_by_state), H


def low(H, count=1):
    if H.shape[0] <= 500:
        e, U = eigh(H.toarray(), subset_by_index=(0, min(count, H.shape[0])-1))
    else:
        rng = np.random.default_rng(817)
        e, U = eigsh(H, k=count, which='SA', tol=2e-13,
                     v0=rng.normal(size=H.shape[0]))
        order = np.argsort(e); e, U = e[order], U[:, order]
    residual = np.linalg.norm(H @ U - U * e, axis=0).max()
    assert residual < 3e-10, residual
    return e, U, float(residual)


def random_charge_control():
    # A non-stoquastic, charge-conserving four-state atomic control.
    rng = np.random.default_rng(614)
    charges = np.array([0, 1, -1, 2])
    energies = np.array([0., 0., 2., 3.])
    qpair = (charges[:, None]+charges[None, :]).ravel()
    raw = rng.normal(size=(16, 16))+1j*rng.normal(size=(16, 16))
    bond = (raw+raw.conj().T)/2
    bond[qpair[:, None] != qpair[None, :]] = 0
    bond *= .004 / np.linalg.norm(bond, 2)
    v = 2*np.linalg.norm(bond, 2)
    eta = .2; C = 1.
    assert v <= min(eta/(4*C), eta/16)
    rows = []
    for M in (4, 5, 6):
        ground = {}
        sectors = {}
        for N in range(M+1):
            basis = [s for s in itertools.product(range(4), repeat=M)
                     if sum(charges[list(s)]) == N]
            index = {s: i for i, s in enumerate(basis)}
            rr, cc, vv = [], [], []
            for col, state in enumerate(basis):
                rr.append(col); cc.append(col); vv.append(sum(energies[list(state)]))
                for j in range(M-1):
                    ab = 4*state[j]+state[j+1]
                    for cd in np.flatnonzero(bond[:, ab]):
                        nxt = list(state); nxt[j], nxt[j+1] = divmod(int(cd), 4)
                        rr.append(index[tuple(nxt)]); cc.append(col); vv.append(bond[cd, ab])
            H = coo_matrix((vv, (rr, cc)), shape=(len(basis), len(basis))).tocsr()
            e, U, residual = low(H)
            psi = U[:, 0]
            occupation = np.abs(psi)**2
            mean_A = float(occupation @ [sum(energies[list(s)]) for s in basis])
            n0 = float(occupation @ [s.count(0) for s in basis])
            n1 = float(occupation @ [s.count(1) for s in basis])
            assert mean_A <= 2*M*v+1e-11
            assert min(n0-M*(1-N/M-2*C*v), n1-M*(N/M-2*C*v)) >= -1e-11
            ground[N] = float(e[0]); sectors[N] = (n0, n1)
            rows.append(dict(cells=M, excess=N, dimension=len(basis),
                             ground=float(e[0]), atomic_energy=mean_A,
                             code0=n0, code1=n1, residual=residual))
        for N in range(M):
            n0, _ = sectors[N]; _, n1 = sectors[N+1]
            assert ground[N+1]-ground[N] <= 2*M*v/n0+1e-10
            assert ground[N]-ground[N+1] <= 2*M*v/n1+1e-10
            if eta <= N/M and (N+1)/M <= 1-eta:
                assert abs(ground[N+1]-ground[N]) <= 4*v/eta+1e-10
    return dict(local_norm=float(v), rows=rows)


def quotient_check(M, L, t, b=3):
    physical, gaps, Hphysical = physical_matrix(M, L, t, b)
    pindex = {state: i for i, state in enumerate(physical)}
    cells = list(compositions(L, M)); cindex = {s: i for i, s in enumerate(cells)}
    dim = L*len(cells)
    rr, cc, vv = [], [], []
    orbit_r, orbit_c, orbit_v = [], [], []
    for x in range(L):
        for ci, state in enumerate(cells):
            col = x*len(cells)+ci
            marks = [x]
            for gap in state[:-1]: marks.append((marks[-1]+gap) % L)
            marks = tuple(sorted(marks))
            orbit_r.append(col); orbit_c.append(pindex[marks]); orbit_v.append(1/np.sqrt(M))
            rr.append(col); cc.append(col); vv.append(sum(atom(ell, b) for ell in state))
            for j in range(M):
                k = (j+1) % M
                for direction in (-1, 1):
                    nxt = list(state); nxt[j] += direction; nxt[k] -= direction
                    if min(nxt[j], nxt[k]) < 2: continue
                    # Moving the reference marker right lengthens the last cell.
                    xx = (x+direction) % L if j == M-1 else x
                    row = xx*len(cells)+cindex[tuple(nxt)]
                    rr.append(row); cc.append(col); vv.append(-t)
    lift = coo_matrix((vv, (rr, cc)), shape=(dim, dim)).tocsr()
    U = coo_matrix((orbit_v, (orbit_r, orbit_c)),
                   shape=(dim, len(physical))).tocsr()
    diff = lift @ U-U @ Hphysical
    error = float(np.max(np.abs(diff.data))) if diff.nnz else 0.
    gram = U.T @ U-diags(np.ones(len(physical)))
    gram_error = float(np.max(np.abs(gram.data))) if gram.nnz else 0.
    assert max(error, gram_error) < 2e-13
    # A genuinely superposed pinned-marker trial, with arbitrary cell lengths.
    rng = np.random.default_rng(319)
    phi = rng.normal(size=len(cells))+1j*rng.normal(size=len(cells))
    phi /= np.linalg.norm(phi)
    raw = np.zeros(dim, dtype=complex); raw[:len(cells)] = phi
    physical_trial = np.sqrt(M) * (U.T @ raw)
    trial_norm = float(np.vdot(physical_trial, physical_trial).real)
    trial_error = float(abs(np.vdot(physical_trial, Hphysical @ physical_trial)
                            -np.vdot(raw, lift @ raw)))
    assert abs(trial_norm-1) < 1e-13 and trial_error < 2e-12
    return dict(cells=M, length=L, rooted_dimension=dim,
                physical_dimension=len(physical), intertwining_error=error,
                isometry_error=gram_error, pinned_trial_error=trial_error)


def physical_tail_control(M, L, t, eta, b=3):
    v = 4*t; gamma = 1.; C = 1.; N = L-b*M
    assert M >= 1/eta and 2*eta <= N/M <= 1-2*eta
    assert v <= min(eta/(4*C), eta*gamma/16)
    basis, gaps, H = physical_matrix(M, L, t, b)
    e, U, residual = low(H, count=8)
    _, Hopen = open_matrix(M, L, t, b)
    eopen, _, ropen = low(Hopen)
    assert abs(e[0]-eopen[0]) <= v+1e-11
    # Spectral-window norm in rooted coordinates equals the largest
    # eigenvalue of Pi times the mean long-cell occupation times Pi.
    Omega = float(e[-1]-e[0]); assert Omega <= gamma/2
    rows = []
    prev = 1.
    for r in range(2, min(7, int(gaps.max())-b+1)):
        tail_weight = np.mean(gaps >= b+r, axis=1)
        amp = float(np.sqrt(np.dot(np.abs(U[:, 0])**2, tail_weight)))
        gram = U.conj().T @ (tail_weight[:, None]*U)
        window_amp = float(np.sqrt(max(0., eigh(gram, eigvals_only=True)[-1])))
        ground_bound = (2*v/gamma)**(r-1)/factorial(r)
        window_bound = (2*v/gamma)**(r-1)/factorial(r-1)
        recurrence = 2*v/(gamma*r)*prev
        assert amp <= min(ground_bound, recurrence)+2e-12
        assert window_amp <= window_bound+2e-12
        rows.append(dict(cutoff=b+r, ground_tail=amp, ground_bound=ground_bound,
                         recursion_bound=recurrence, window_norm=window_amp,
                         window_bound=window_bound))
        prev = amp
    # Test the min-max cutoff estimate with the same eight low levels.
    cutoff = b+3; r = cutoff-b
    keep = np.flatnonzero(gaps.max(axis=1) < cutoff)
    ecut, _, rcut = low(H[keep][:, keep], count=8)
    delta = np.sqrt(M)*(2*v/gamma)**(r-1)/factorial(r-1)
    error_bound = (Omega*delta**2+M*v*delta)/(1-delta**2)
    shifts = ecut-e
    assert shifts.min() > -2e-10 and shifts.max() <= error_bound+2e-10
    return dict(cells=M, length=L, dimension=len(basis), excess=N,
                eta=eta, hopping=t, local_norm=v, ground=float(e[0]),
                open_ground=float(eopen[0]), window_levels=len(e),
                window_width=Omega, max_residual=max(residual, ropen, rcut),
                tails=rows, cutoff=cutoff, cutoff_dimension=len(keep),
                largest_cutoff_shift=float(shifts.max()),
                cutoff_error_bound=float(error_bound))


def main():
    out = dict(status='finite algebra and spectral controls, no CFT measurement')
    out['random_charge_control'] = random_charge_control()
    out['quotient'] = [quotient_check(M, L, .003) for M, L in ((3, 10), (4, 14))]
    out['physical_tails'] = [physical_tail_control(M, L, .003, eta)
                             for M, L, eta in ((4, 14, .25), (5, 17, .2),
                                               (5, 18, .2), (6, 21, .2))]
    Path('canonical_length_checks.json').write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps(dict(status='passed',
                         charge_cases=len(out['random_charge_control']['rows']),
                         quotient=out['quotient'],
                         physical=[{k: row[k] for k in
                                    ('cells', 'length', 'dimension', 'max_residual',
                                     'largest_cutoff_shift', 'cutoff_error_bound')}
                                   for row in out['physical_tails']]), indent=2))


if __name__ == '__main__':
    main()
