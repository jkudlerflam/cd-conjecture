"""Bounded single-worker verification; no DMRG or parameter scan."""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
from pathlib import Path
from itertools import product, combinations_with_replacement
from collections import defaultdict
import importlib.util
import json
import time
import numpy as np
import scipy.linalg as la
import scipy.sparse as sp
from scipy.sparse.linalg import splu
import mpmath as mp
import model
from moments import proof_identities
from response import coefficients, low_limit, roots_and_response, carrier

HERE = Path(__file__).resolve().parent


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def literal_and_grammar():
    old = load('original_bond9', HERE.parent/'markerless_d9_d10'/'bond9.py')
    p = dict(lam=.23, g=.17, Delta=2.8, D=4., nu=.41, f=0.,
             phiA=.27, phiB=.04, td=.09, psiA=.7, psiB=1.2)
    W, o, bad = model.bond(dict(p, q=1., w=1., bias=0., nu_head=p['nu']))
    W0, o0, bad0 = old.bond(p)
    assert np.array_equal(W, W0) and np.array_equal(o, o0) and np.array_equal(bad, bad0)
    W, o, bad = model.bond(dict(p, q=.13, w=.37, nu_head=.41, hL=.03))
    assert np.max(abs(W-W.conj().T)) == 0
    q2 = np.add.outer(model.CHARGE, model.CHARGE).ravel()
    assert np.max(abs(W*(q2[:, None]-q2[None, :]))) == 0
    isL = lambda x: x in (0, 1)
    isR = lambda x: x in (2, 3)
    def invariants(word):
        return (sum(isL(a) and isR(b) for a, b in zip(word, word[1:])),
                tuple(not isL(a) and not isR(b) for a, b in zip(word, word[1:])))
    count = 0
    for row, col in zip(*np.nonzero(W)):
        if row == col:
            continue
        for x, y in product(range(9), repeat=2):
            before = (x, col//9, col % 9, y)
            after = (x, row//9, row % 9, y)
            assert invariants(before) == invariants(after)
            count += 1
    # Clock and readout from literal diagonal entries, not a chosen logical field.
    n, q = 3, .125
    p = dict(q=q, w=.4, g=0., lam=0., D=4., hL=0., nu_head=0.)
    W, o, _ = model.bond(p)
    ev, v, gamma = model.clock(n, q)
    max_clock_error = max_field_error = 0.
    for bits in product((0, 1), repeat=n):
        words = [(model.IDX['R0'],)+tuple(
            [model.IDX[f'L{a}'] for a in bits[:j]]+[model.IDX['W']]+
            [model.IDX[f'R{a}'] for a in bits[j:]])+(model.IDX['L0'],)
                 for j in range(n+1)]
        H = np.zeros((n+1, n+1), complex)
        for j, word in enumerate(words):
            H[j, j] += sum(o[x] for x in word)
            for site in range(len(word)-1):
                col = word[site]*9+word[site+1]
                for i, target in enumerate(words):
                    if word[:site] == target[:site] and word[site+2:] == target[site+2:]:
                        H[i, j] += W[target[site]*9+target[site+1], col]
        expected = -q*(np.eye(n+1, k=1)+np.eye(n+1, k=-1))
        expected[0, 0] = -1.
        max_clock_error = max(max_clock_error, float(np.max(abs(H-expected))))
        projected = sum(v[j]**2*sum(x == model.IDX['L1'] for x in word)
                        for j, word in enumerate(words))
        max_field_error = max(max_field_error, abs(projected-np.dot(bits, gamma)))
    assert max_clock_error == 0 and max_field_error < 1e-15
    return dict(original_bond_exactly_recovered=True, local_contexts_checked=count,
                dimension=9, hermitian=True, charge_conserved=True,
                clock_error=max_clock_error, field_error=max_field_error,
                clock_positions=n+1, field_profile=gamma.tolist())


def gauge_and_counts():
    old = load('original_scalar9', HERE.parent/'markerless_d9_d10'/'scalar9.py')
    reports = []
    for n in (2, 3):
        phi = np.pi/(2*n)
        HA = np.array(old.Hs(n, 4., 3., phi, .05, np.pi).tolist(), complex)
        HB = np.array(old.Hs(n, 4., 3., 0., .05, np.pi/2).tolist(), complex)
        U = np.diag(np.r_[np.exp(-1j*phi*np.arange(n+1)), np.ones(n+1)])
        perm = list(range(n+1, 2*n+2))+list(range(n+1))
        actual, _ = carrier(n, 3., 1., 4., .05)
        error = max(la.norm((U@HA@U.conj().T)[np.ix_(perm, perm)]-actual),
                    la.norm(HB[np.ix_(perm, perm)]-actual))
        assert error < 1e-14 and abs(U[n, n]+1j) < 1e-14
        reports.append(dict(n=n, gauge_error=float(error), source='X_n,D_0',
                            conversions=list(range(n)), inactive_endpoint=n))
    p = dict(q=.125, w=.5, D=10., nu_head=1.05, Delta=3., g=.001,
             lam=.0001, td=.00001, hL=1e-10, phiA=np.pi/6,
             psiA=np.pi, psiB=np.pi/2)
    W, o, _ = model.bond(p)
    off = W-np.diag(np.diag(W))
    directed = int(np.count_nonzero(off))
    real_basis = sum(int(abs(z.real) > 1e-13*abs(z))+int(abs(z.imag) > 1e-13*abs(z))
                     for z in off[np.triu_indices(81, 1)] if abs(z) > 0)
    reduced = off.copy(); reduced_o = o.copy(); reduced_o[model.IDX['W']] -= 1
    for a in (0, 1):
        i = model.IDX[f'L{a}']*9+model.IDX['W']; reduced[i, i] = 1
    for word in product(range(9), repeat=3):
        raw = sum(o[x] for x in word)+sum(W[word[j]*9+word[(j+1) % 3], word[j]*9+word[(j+1) % 3]].real for j in range(3))
        new = sum(reduced_o[x] for x in word)+sum(reduced[word[j]*9+word[(j+1) % 3], word[j]*9+word[(j+1) % 3]].real for j in range(3))
        assert abs(raw-new) < 1e-14
    assert directed == 82 and real_basis == 45 and np.count_nonzero(reduced_o) == 6
    return dict(gauges=reports, diagonal_rewriting_words=9**3,
                directed_offdiagonal_units=directed, real_Hermitian_offdiagonal_units=real_basis,
                reduced_two_site_diagonals=2, onsite_diagonals=6,
                directed_plus_diagonal_total=90, real_Hermitian_total=53)


def exact_shortest_walks():
    """Integer Neumann coefficients, independently of low_limit()."""
    reports = []
    for n in (2, 3, 4):
        p = n+1; size = p*p
        path = np.eye(p, k=1, dtype=np.int64)+np.eye(p, k=-1, dtype=np.int64)
        K = np.kron(path, np.eye(p, dtype=np.int64))+np.kron(np.eye(p, dtype=np.int64), path)
        powers = [np.eye(size, dtype=np.int64)]
        for _ in range(n):
            powers.append(powers[-1]@K)
        zero = [P[:, 0] for P in powers]; end = [P[:, n*p] for P in powers]
        leading = None
        for degree in range(n+1):
            B = np.zeros((n, n), dtype=np.int64)
            for i, j in product(range(n), repeat=2):
                for r in range(degree+1):
                    for s in range(degree-r+1):
                        t = degree-r-s
                        for coord in (0, 1):
                            maski = np.array([(k//p if coord == 0 else k % p) == i for k in range(size)])
                            maskj = np.array([(k//p if coord == 0 else k % p) == j for k in range(size)])
                            B[i, j] += (zero[r]*maski)@powers[s]@(end[t]*maskj)
                            B[i, j] -= (end[r]*maski)@powers[s]@(zero[t]*maskj)
            if degree < n:
                assert np.count_nonzero(B) == 0
            else:
                leading = B
                assert np.array_equal(B, np.triu(np.ones((n, n), dtype=np.int64), 1)-np.tril(np.ones((n, n), dtype=np.int64), -1))
        reports.append(dict(n=n, exact=True, lower_coefficients_zero=True,
                            coefficient_of_w_to_n_at_2a_equal_1=leading.tolist()))
    return reports


def build_two_cells(n, p):
    cells = model.cell_words(n)
    states = [((model.IDX['R0'],)+a[0]+b[0]+(model.IDX['L0'],),
               (a[1], b[1]), (a[2], b[2]), a[3]+b[3]) for a in cells for b in cells]
    index = {s[0]: i for i, s in enumerate(states)}
    W, o, _ = model.bond(p)
    by_in = defaultdict(list)
    for row, col in zip(*np.nonzero(W)):
        by_in[int(col)].append((int(row), W[row, col]))
    rows, cols, vals = [], [], []
    leaks = 0
    for i, (word, _, _, _) in enumerate(states):
        rows.append(i); cols.append(i); vals.append(sum(o[x] for x in word))
        for site in range(len(word)-1):
            for row, value in by_in[word[site]*9+word[site+1]]:
                target = word[:site]+(row//9, row % 9)+word[site+2:]
                if target not in index:
                    leaks += 1
                else:
                    rows.append(index[target]); cols.append(i); vals.append(value)
    return states, sp.csc_matrix((vals, (rows, cols)), shape=(len(states),)*2), leaks


def two_cell_coefficient():
    reports = []
    for n in (1, 2):
        p = dict(q=.21, w=.6, g=1., lam=1., D=2., Delta=3.,
                 phiA=np.pi/(2*n), phiB=0., psiA=np.pi, psiB=np.pi/2, td=.08)
        states, full, leaks = build_two_cells(n, p)
        _, nog, _ = build_two_cells(n, dict(p, g=0.))
        _, nol, _ = build_two_cells(n, dict(p, lam=0.))
        port, conv = full-nol, full-nog
        H0 = nog-port
        vac = np.array([i for i, s in enumerate(states) if s[1] == ('W', 'W')])
        exc = np.array([i for i, s in enumerate(states) if 'W' not in s[1]])
        ev, v, _ = model.clock(n, p['q'])
        bits = list(product((0, 1), repeat=2*n)); bidx = {b: i for i, b in enumerate(bits)}
        P = np.zeros((len(vac), len(bits)))
        for i, state in enumerate(vac):
            s = states[state]
            P[i, bidx[s[3]]] = v[s[2][0]]*v[s[2][1]]
        code_error = float(np.max(abs(H0[vac][:, vac]@P-2*ev*P)))
        source = port[exc][:, vac]@P
        lu = splu(H0[exc][:, exc]-2*ev*sp.eye(len(exc), format='csc'))
        X = lu.solve(source)
        X1 = lu.solve(conv[exc][:, exc]@X)
        X2 = lu.solve(conv[exc][:, exc]@X1)
        h0, h1, h2 = [-source.conj().T@Y for Y in (X, X1, X2)]
        i0 = bidx[(0,)*(2*n)]
        inds = [bidx[tuple(int(i == j) for i in range(2*n))] for j in range(2*n)]
        h = h2[np.ix_(inds, inds)]-h2[i0, i0]*np.eye(2*n)
        rebuilt = h2[i0, i0]*np.eye(len(bits), dtype=complex)
        for col, bit in enumerate(bits):
            for j in range(2*n):
                if not bit[j]:
                    continue
                mid = list(bit); mid[j] = 0
                for i in range(2*n):
                    if mid[i]:
                        continue
                    end = mid.copy(); end[i] = 1
                    rebuilt[bidx[tuple(end)], col] += h[i, j]*(-1)**(sum(bit[:j])+sum(mid[:i]))
        T, M = coefficients(n, p['D'], p['w'], p['Delta']-ev, p['td'])
        visibility = (v[0]*v[-1])**2
        scale = la.norm(T)*visibility
        terr = float(la.norm(h[:n, n:]-visibility*T)/scale)
        merr = float(la.norm(h[:n, :n]+h[n:, n:]-visibility*M)/scale)
        gaussian = float(la.norm(h2-rebuilt)/scale)
        assert leaks == 0 and code_error < 1e-13
        assert terr < 1e-10 and merr < 1e-8 and gaussian < 1e-8
        assert la.norm(h1) < 1e-14
        assert la.norm(h0-h0[0, 0]*np.eye(len(bits))) < 1e-12*max(la.norm(h0), 1e-99)
        reports.append(dict(n=n, states=len(states), leaks=leaks, code_error=code_error,
                            visibility=visibility, hopping_relative_error=terr,
                            onsite_relative_error=merr, Fock_quadratic_relative_error=gaussian))
    return reports


def finite_delay():
    reports = []
    for n in (2, 3):
        w, q = .5, .0625
        ev, _, gamma = model.clock(n, q); a = 3-ev
        G, B, Q = low_limit(n, w, a)
        F0 = np.diag(gamma/q**2)
        epsilon = .1/max(la.eigvalsh(F0, G)); F = epsilon*F0
        e, X = la.eigh(F, G); e = e[::-1]; X = X[:, ::-1]
        QQ, BB = X.T@Q@X, X.T@B@X
        expected = np.array([-16*(n+1)*sum(QQ[j, l]*BB[j, l]/(e[j]-e[l])
                                         for l in range(n) if l != j) for j in range(n)])
        rows = []
        for D in (1e4, 1e5, 1e6):
            T, M1 = coefficients(n, D, w, a)
            scale = D**(2*n+2)
            responses = np.array([x[1]*D for x in roots_and_response(T*scale, M1*scale/D**n, F)])
            rows.append(dict(D=D, D_times_pair_response=responses.tolist(),
                             response_relative_error=float(la.norm(responses-expected)/la.norm(expected)),
                             delay_relative_error=float(la.norm(scale*D*(T-T.T)/(4*(n+1))-Q)/la.norm(Q)),
                             source_relative_error=float(la.norm(M1*scale/D**n-2j*B)/la.norm(B))))
        for key in ('response_relative_error', 'delay_relative_error', 'source_relative_error'):
            assert rows[-1][key] < 1e-3 and rows[-1][key] < rows[0][key]/20
        T, M = coefficients(n, 20., w, a, td=1e-5)
        Tm, Mm = coefficients(n, 20., w, a, td=-1e-5)
        assert np.array_equal(T, Tm) and np.array_equal(M, -Mm)
        reports.append(dict(n=n, vertices=2*n+2, express_vertices=n+1,
                            a=a, q=q, w=w, epsilon=epsilon, expected=expected.tolist(), rows=rows,
                            hopping_even=True, onsite_odd=True))
    return reports


def mp_blocks(n, D, w, a, td):
    """Independent dense high-precision inverse, not the sparse response code."""
    p, s = n+1, 2*n+2
    H = mp.diag([a]*p+[a+D]*p)
    for j in range(n):
        H[j, j+1] = H[j+1, j] = -w
        H[p+j, p+j+1] = H[p+j+1, p+j] = -1
    H[0, p] = H[p, 0] = -1
    H[-1, n], H[n, -1] = -1j*td, 1j*td
    L = mp.matrix(s*s)
    for i, j, k, l in product(range(s), repeat=4):
        L[i*s+j, k*s+l] = H[i, k]*(j == l)+(i == k)*H[j, l]
    R = L**-1; x = R[:, (s-1)*s]
    U, V = mp.matrix(s*s, n), mp.matrix(s*s, n)
    for j in range(n):
        for b in range(s):
            U[j*s+b, j] = x[j*s+b]
            V[b*s+j, j] = x[b*s+j]
    C, A = U.H*R*V, U.H*R*U+V.H*R*V
    return (mp.matrix([[-2*mp.re(C[i, j]) for j in range(n)] for i in range(n)]),
            mp.matrix([[-2j*mp.im(A[i, j]) for j in range(n)] for i in range(n)]))


def high_precision_roots():
    mp.mp.dps = 70
    n = 2; D = mp.mpf(10000); w = mp.mpf(1)/2; q = mp.mpf(1)/16
    eps = mp.mpf('1e-6'); eta = mp.mpf('1e-6')
    clock = mp.matrix(n+1)
    clock[0, 0] = -1
    for j in range(n):
        clock[j, j+1] = clock[j+1, j] = -q
    ev, v = mp.eigsy(clock); a = 3-ev[0]
    F = mp.diag([eps*sum(v[j, 0]**2 for j in range(i+1, n+1))/q**2 for i in range(n)])
    scale = D**(2*n+2)
    T0, _ = mp_blocks(n, D, w, a, 0); T0 *= scale
    tf = np.array(T0.tolist(), float); ff = np.array(F.tolist(), float)
    seeds = [x[0] for x in roots_and_response(tf, np.zeros((n, n)), ff)]
    reports = []
    for eeta in (eta, eta/2):
        T, M = mp_blocks(n, D, w, a, eeta/D**n)
        T *= scale; M *= scale
        determinant = lambda k: mp.re(mp.det(T*mp.exp(1j*k)+T.T*mp.exp(-1j*k)+M+F))
        roots = [mp.findroot(determinant, (mp.mpf(k)-mp.mpf('1e-5'), mp.mpf(k)+mp.mpf('1e-5')))
                 for k in seeds+[-k for k in seeds]]
        wrap = lambda x: abs((x+mp.pi) % (2*mp.pi)-mp.pi)
        pairs = list(combinations_with_replacement(range(2*n), 2))
        d2 = min(wrap(x-y) for i, x in enumerate(roots) for y in roots[i+1:])
        d4 = min(wrap(roots[a]+roots[b]-roots[c]-roots[d])
                 for i, (a, b) in enumerate(pairs) for c, d in pairs[i+1:])
        sums = [roots[j]+roots[n+j] for j in range(n)]
        assert d2 > 0 and d4 > 0 and abs(sum(sums)) < mp.mpf('1e-50')
        reports.append(dict(eta=str(eeta), roots=list(map(str, roots)),
                            pair_sums=list(map(str, sums)), delta2=str(d2), delta4=str(d4)))
    assert abs(mp.mpf(reports[0]['pair_sums'][0])/mp.mpf(reports[1]['pair_sums'][0])-2) < mp.mpf('1e-5')
    return dict(n=n, D=str(D), q=str(q), w=str(w), epsilon=str(eps),
                precision=mp.mp.dps, status='validation, not an interval certificate', checks=reports)


def main():
    start = time.monotonic()
    out = dict(workers=1, status='bounded algebra checks; all-n proof is in MOMENT_SIGN.md and DERIVATION.md')
    for name, run in [('literal_and_grammar', literal_and_grammar),
                      ('gauge_and_term_counts', gauge_and_counts),
                      ('two_cell_coefficient', two_cell_coefficient),
                      ('exact_moment_identities', lambda: [proof_identities(n) for n in range(2, 7)]),
                      ('exact_shortest_walks', exact_shortest_walks),
                      ('finite_delay', finite_delay), ('high_precision_roots', high_precision_roots)]:
        out[name] = run()
        print(name, 'passed', flush=True)
        (HERE/'checks.json').write_text(json.dumps(out, indent=2)+'\n')
    out['runtime_seconds'] = time.monotonic()-start
    (HERE/'checks.json').write_text(json.dumps(out, indent=2)+'\n')
    print('All checks passed in', round(out['runtime_seconds'], 2), 'seconds.')


if __name__ == '__main__':
    main()
