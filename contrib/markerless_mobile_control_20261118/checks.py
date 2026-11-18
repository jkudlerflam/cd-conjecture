"""Bounded one-worker checks of the literal markerless mobile vacuum sector.

The eigenvectors are only certificate witnesses. Integer Collatz bounds, not
floating eigenvalue residuals, certify the finite comparisons. No thermodynamic
gap or central charge is inferred. All output stays beside this script.
"""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'

from itertools import combinations_with_replacement, product
from fractions import Fraction
from pathlib import Path
import argparse
import json
import sys
import time

import mpmath as mp
import numpy as np
import scipy
import scipy.linalg as la
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'markerless_d9_d10'))
import bond9
import bond10


def selector(b):
    a = lambda k: 2*(1-np.cos(np.pi/(k+1)))
    d = lambda k: (k+1)*a(k)-k*a(k+1)
    return 2-(d(b)+d(b-1))/2


def alpha(N, V):
    if V == 0:
        return 1.
    L = N+V
    return N/L-np.sin(np.pi*N/L)**2/(N*L*np.sin(np.pi/L)**2)


def bath(N, V):
    # Multisets of vacancy gap labels: dimension binomial(N+V-1,V).
    gaps = []
    for labels in combinations_with_replacement(range(N), V):
        g = [0]*N
        for j in labels:
            g[j] += 1
        gaps.append(tuple(g))
    index = {g: i for i, g in enumerate(gaps)}
    rows, cols = [], []
    for i, g in enumerate(gaps):
        for j, count in enumerate(g):
            if not count:
                continue
            v = list(g); v[j] -= 1; v[(j-1) % N] += 1
            k = index[tuple(v)]
            rows += [i, k]; cols += [k, i]
    H = sp.coo_matrix((-np.ones(len(rows), dtype=np.int64), (rows, cols)),
                      shape=(len(gaps),)*2).tocsr()
    contacts = np.array(gaps, dtype=np.int64) == 0
    return H, contacts


def cell(b, nu):
    """b head positions, b+1 contact sources with labels -1,...,b-1."""
    v = np.sqrt(2/(b+1))*np.sin(np.pi*np.arange(1, b+1)/(b+1))
    sources = np.zeros((b+1, b, b))
    sources[0, 0, 0] = nu/2
    sources[-1, -1, -1] = nu/2
    for j in range(b-1):
        sources[j+1, j, j] = sources[j+1, j+1, j+1] = nu/2
        sources[j+1, j, j+1] = sources[j+1, j+1, j] = -1
    a = np.einsum('i,kij,j->k', v, sources, v)
    xi = sources@v-a[:, None]*v
    eps = nu-2*np.cos(np.pi/(b+1))
    assert la.norm(sources.sum(axis=0)@v-eps*v) < 2e-14
    assert la.norm(xi.sum(axis=0)) < 2e-14
    return dict(v=v, sources=sources, a=a, xi=xi, eps=eps)


def embed(mat, i, dims):
    dtype = np.asarray(mat).dtype
    out = sp.csr_matrix(np.array([[1]], dtype=dtype))
    for j, d in enumerate(dims):
        out = sp.kron(out, sp.csr_matrix(mat) if j == i else sp.eye(d, dtype=dtype), format='csr')
    return out


def sources(word, nu):
    N, D = sum(word), int(np.prod(word))
    out = [sp.csr_matrix((D, D)) for _ in range(N)]
    start = 0
    for c, b in enumerate(word):
        for r, mat in enumerate(cell(b, nu)['sources']):
            out[(start+r-1) % N] += embed(mat, c, word)
        start += b
    return out


def integer_hamiltonian(word, V):
    """Exactly 5000 H for t=1, J=1/25, nu=177/100, f=lambda=0."""
    N, D = sum(word), int(np.prod(word))
    Hc, contacts = bath(N, V)
    out = 5000*sp.kron(Hc, sp.eye(D, dtype=np.int64), format='csr')
    # Build integer material sources directly; no rounding of an operator.
    hs = [sp.csr_matrix((D, D), dtype=np.int64) for _ in range(N)]
    start = 0
    for c, b in enumerate(word):
        for r in range(b+1):
            mat = np.zeros((b, b), dtype=np.int64)
            if r == 0:
                mat[0, 0] = 177
            elif r == b:
                mat[-1, -1] = 177
            else:
                mat[r-1, r-1] = mat[r, r] = 177
                mat[r-1, r] = mat[r, r-1] = -200
            hs[(start+r-1) % N] += embed(mat, c, word).astype(np.int64)
        start += b
    for j, h in enumerate(hs):
        out += sp.kron(sp.diags(contacts[:, j].astype(np.int64), dtype=np.int64), h, format='csr')
    out = out.astype(np.int64)
    assert (out-out.T).nnz == 0
    return out


def collatz(B, witness):
    """Exact integer floor/ceiling of every local-energy ratio, grid 1e-10."""
    assert witness.dtype == np.int64 and witness.min() > 0
    row_norm = int(np.max(np.asarray(abs(B).sum(axis=1))))
    assert row_norm*int(witness.max()) < np.iinfo(np.int64).max
    z = B@witness  # exact int64 arithmetic, no overflow by preceding check
    grid = 10**10
    lows = [(int(a)*grid)//(5000*int(b)) for a, b in zip(z, witness)]
    highs = [-((-int(a)*grid)//(5000*int(b))) for a, b in zip(z, witness)]
    lo, hi = min(lows), max(highs)
    return dict(lower_numerator=lo, upper_numerator=hi, denominator=grid,
                lower=lo/grid, upper=hi/grid, width=(hi-lo)/grid,
                min_witness=int(witness.min()), max_witness=int(witness.max()))


WORDS = {'target': [10, 10], 'balanced': [9, 11],
         'one_head': [20], 'three_heads': [6, 7, 7]}


def finite_certificates(verify_only=False):
    saved = np.load(HERE/'witnesses.npz') if verify_only else None
    witnesses, records = {}, {}
    for name, word in WORDS.items():
        B = integer_hamiltonian(word, 2)
        if verify_only:
            w = saved[name]
        else:
            H = B.astype(float)/5000
            # Positive product start and deterministic witness rounding.
            Hc, _ = bath(sum(word), 2)
            _, oc = eigsh(Hc.astype(float), k=1, which='SA', tol=1e-14,
                          v0=np.ones(Hc.shape[0]))
            v = np.ones(1)
            for b in word:
                v = np.kron(v, cell(b, 1.77)['v'])
            v0 = np.kron(abs(oc[:, 0]), v)
            e, vec = eigsh(H, k=1, which='SA', tol=3e-14, maxiter=10000, v0=v0)
            psi = vec[:, 0]*np.sign(vec[0, 0])
            assert psi.min() > 0
            w = np.rint(psi/psi.max()*2**40).astype(np.int64)
        witnesses[name] = w
        records[name] = dict(word=word, dimension=B.shape[0], nnz=B.nnz,
                             **collatz(B, w))
        print(name, records[name], flush=True)
    if not verify_only:
        np.savez_compressed(HERE/'witnesses.npz', **witnesses)
    target = records['target']
    differences = {}
    for name, rec in records.items():
        if name == 'target':
            continue
        lo = rec['lower_numerator']-target['upper_numerator']
        hi = rec['upper_numerator']-target['lower_numerator']
        assert lo > 0
        # nu_mid lies in (1.77, 1.7703), proved separately by rational Taylor bounds.
        mid_lo = lo-240000
        window_lo = mid_lo-5*(len(rec['word'])+2)*37703
        full_lo = mid_lo-400000*len(rec['word'])
        assert min(mid_lo, window_lo, full_lo) > 0
        differences[name] = dict(lower_numerator=lo, upper_numerator=hi,
                                 denominator=10**10,
                                 lower=lo/10**10, upper=hi/10**10,
                                 lower_for_nu_interval_numerator=mid_lo,
                                 lower_for_nu_in_1p77_to_1p7703=mid_lo/10**10,
                                 lower_with_J_window_numerator=window_lo,
                                 lower_with_J_window=window_lo/10**10,
                                 lower_with_port_and_field_numerator=full_lo,
                                 lower_with_port_and_field=full_lo/10**10)
    return dict(N=20, V=2, rho=20/22, J_over_t='1/25', nu='177/100',
                alpha=alpha(20, 2), sectors=records, differences=differences)


def selector_rational_enclosure():
    def atan_bounds(x, terms):
        partial = sum(((-1)**j*x**(2*j+1)/ (2*j+1) for j in range(terms)), Fraction(0))
        next_partial = partial+(-1)**terms*x**(2*terms+1)/(2*terms+1)
        return min(partial, next_partial), max(partial, next_partial)
    a, b = atan_bounds(Fraction(1, 5), 24), atan_bounds(Fraction(1, 239), 8)
    pi_lo, pi_hi = 16*a[0]-4*b[1], 16*a[1]-4*b[0]
    def cos_bounds(x):
        term, partial = Fraction(1), Fraction(1)
        for j in range(1, 16):
            term *= -x*x/((2*j-1)*(2*j))
            partial += term
        other = partial+term*(-x*x)/(31*32)
        return min(partial, other), max(partial, other)
    cs = {}
    for b in (9, 10, 11):
        cs[b] = (cos_bounds(pi_hi/(b+1))[0], cos_bounds(pi_lo/(b+1))[1])
    lower = 2*cs[10][0]+10*cs[9][0]-10*cs[11][1]
    upper = 2*cs[10][1]+10*cs[9][1]-10*cs[11][0]
    assert Fraction(177, 100) < lower < upper < Fraction(17703, 10000)
    return dict(nu_lower=float(lower), nu_upper=float(upper),
                rational_177_over_100_lower=True, rational_17703_over_10000_upper=True,
                method='Machin alternating arctangent bounds; alternating cosine bounds; exact fractions')


def literal_audit():
    p = dict(nu=1.77, lam=0., g=.2, Delta=2.1, D=5., f=0., td=.03, K=2.)
    bond = bond9.full_bond(p)
    records = []
    for word in ([2, 3], [3, 2, 2], [10, 10], [9, 11]):
        words = []
        for heads in product(*(range(b) for b in word)):
            letters = []
            for b, j in zip(word, heads):
                letters += ['L0']*j+['W']+['R0']*(b-j-1)
            words.append(tuple(bond9.IDX[s] for s in letters))
        index = {w: i for i, w in enumerate(words)}
        N, dim = sum(word), len(words)
        actual = [np.zeros((dim, dim)) for _ in range(N)]
        leakage = 0
        for col, w in enumerate(words):
            for j in range(N):
                k = (j+1) % N
                ab = w[j]*9+w[k]
                for cd in np.flatnonzero(abs(bond[:, ab]) > 1e-14):
                    new = list(w); new[j], new[k] = divmod(int(cd), 9)
                    row = index.get(tuple(new))
                    if row is None:
                        leakage += 1
                    else:
                        actual[j][row, col] += bond[cd, ab].real
        expected = sources(word, p['nu'])
        error = max(float(np.max(abs(a-b.toarray()))) for a, b in zip(actual, expected))
        assert error < 1e-14 and leakage == 0
        records.append(dict(word=word, dimension=dim, max_source_error=error, leakage=leakage))
    return records


def response(N, V):
    H, contacts = bath(N, V)
    es, us = la.eigh(H.toarray().astype(float))
    omega = us[:, 0]
    X = (contacts.astype(float)-alpha(N, V))*omega[:, None]
    Z = us[:, 1:].T@X
    G = (Z.T/(es[1:]-es[0]))@Z
    L = N+V
    free = -2*np.sin(np.pi*N/L)/np.sin(np.pi/L)
    assert abs(es[0]-free) < 2e-12
    assert np.max(abs(omega@X)) < 2e-13
    assert np.max(abs(G-np.roll(np.roll(G, 1, axis=0), 1, axis=1))) < 2e-11
    return G, dict(N=N, V=V, rho=N/L, alpha=alpha(N, V),
                   bath_dimension=H.shape[0], bath_gap=float(es[1]-es[0]))


def one_vacancy_response(N):
    eig = np.zeros(N)
    eig[1:] = 1/(4*N*np.sin(np.pi*np.arange(1, N)/N)**2)
    col = np.fft.ifft(eig).real
    return la.circulant(col)


def response_word(word, nu, G):
    N = sum(word)
    sigma = np.zeros(N)
    start, Cclock, energy = 0, 0., 0.
    for b in word:
        info = cell(b, nu)
        js = (start+np.arange(b+1)-1) % N
        np.add.at(sigma, js, info['a'])
        Cclock += float(np.sum(G[np.ix_(js, js)]*(info['xi']@info['xi'].T)))
        energy += info['eps']
        start += b
    Cscalar = float(sigma@G@sigma)
    return dict(frozen_energy=energy, scalar_response=Cscalar,
                clock_response=Cclock, total_response=Cscalar+Cclock)


def coefficient_checks():
    nu = selector(10)
    out = []
    G, meta = response(20, 2)
    pairs = [(WORDS['target'], WORDS[k], k) for k in ('balanced', 'one_head', 'three_heads')]
    for target, other, label in pairs:
        a, b = response_word(target, nu, G), response_word(other, nu, G)
        out.append(dict(label=label, target=target, competitor=other, **meta,
                        first=meta['alpha']*(b['frozen_energy']-a['frozen_energy']),
                        second_scalar=a['scalar_response']-b['scalar_response'],
                        second_clock=a['clock_response']-b['clock_response'],
                        second=a['total_response']-b['total_response'], A=a, B=b))
    for length in (9, 11):
        N = 10*length
        G = one_vacancy_response(N)
        target, other = [10]*length, [length]*10
        a, b = response_word(target, nu, G), response_word(other, nu, G)
        out.append(dict(label=f'weak_head_change_{length}', target=target, competitor=other,
                        N=N, V=1, rho=N/(N+1), alpha=1-1/N,
                        first=(1-1/N)*(b['frozen_energy']-a['frozen_energy']),
                        second_scalar=a['scalar_response']-b['scalar_response'],
                        second_clock=a['clock_response']-b['clock_response'],
                        second=a['total_response']-b['total_response'], A=a, B=b))
    test, _ = response(12, 1)
    assert np.max(abs(test-one_vacancy_response(12))) < 2e-12
    return out


def relative_response_constants():
    nu = selector(10)
    def profile(word):
        # Unwrapped patch, including its TWO DISTINCT exterior contacts.
        out = np.zeros(sum(word)+1)
        start = 0
        for b in word:
            out[start:start+b+1] += cell(b, nu)['a']
            start += b
        return out
    p = np.zeros(10)
    np.add.at(p, (np.arange(11)-1) % 10, cell(10, nu)['a'])
    p_fourier_l1 = float(sum(abs(np.fft.fft(p)/10)))
    out = []
    for label, a, b in (
            ('balanced', [10, 10], [9, 11]),
            ('more_heads', [10]*9, [9]*10),
            ('fewer_heads', [10]*11, [11]*10)):
        d = profile(b)-profile(a)
        scalar = 2*sum(abs(d))*p_fourier_l1+float(d@d)
        clock = sum(float(np.sum(cell(k, nu)['xi']**2)) for k in a+b)
        excess = sum(cell(k, nu)['eps'] for k in b)-sum(cell(k, nu)['eps'] for k in a)
        out.append(dict(label=label, patch_material_length=sum(a),
                        difference_l1=float(sum(abs(d))), difference_l2_sq=float(d@d),
                        target_fourier_l1=p_fourier_l1, scalar_constant=scalar,
                        clock_constant=clock, C_patch=scalar+clock,
                        frozen_excess=excess,
                        conditional_J_over_t_times_Gamma=.29*excess/(scalar+clock)))
    return out


def form_factor_checks():
    nu = selector(10)
    qs = np.linspace(0, 2*np.pi, 2001)
    out = []
    for b in (9, 10, 11):
        info = cell(b, nu)
        phases = np.exp(1j*np.outer(qs, np.arange(b+1)-1))
        direct = np.sum(abs(phases@info['xi'])**2, axis=1)
        D = 4*np.sin(np.pi/(b+1))**2*(b-1)/(b+1)
        formula = (info['eps']**2*np.cos(qs/2)**2+D*np.sin(qs/2)**2
                   -abs(phases@info['a'])**2)
        error = float(np.max(abs(direct-formula)))
        assert error < 5e-15
        cumulative = np.cumsum(info['xi'], axis=0)[:-1]
        U, V = float(np.sum(cumulative**2)), float(np.sum(cumulative.sum(axis=0)**2))
        N = 110
        G = one_vacancy_response(N)
        js = (np.arange(b+1)-1) % N
        R = float(np.sum(G[np.ix_(js, js)]*(info['xi']@info['xi'].T)))
        assert abs(R-(U/N-V/N**2)) < 2e-15
        out.append(dict(b=b, epsilon=info['eps'], derivative_norm=D,
                        form_factor_identity_error=error,
                        source_l1=float(sum(abs(info['a']))),
                        centered_source_sq=float(np.sum(info['xi']**2)),
                        single_vacancy_U=U, single_vacancy_V=V))
    return out


def charge_audit():
    p = dict(lam=.3, g=.2, Delta=2.1, D=5., nu=1.77, f=.05,
             phiA=.2, td=.03, J=.7, t=1.3, mu=0., K=2.)
    q = np.array([bond10.CHARGE[s] for s in bond10.LET])
    q2 = np.add.outer(q, q).ravel()
    h = bond10.bond(p)
    p0 = dict(p, f=0.)
    new = bond10.bond(p0)+p['J']*.4*p['f']/4*np.diag(q2)
    assert np.max(abs(h*(q2[:, None]-q2[None, :]))) == 0
    assert np.max(abs(new*(q2[:, None]-q2[None, :]))) == 0
    for labels in ([1]*9+[0], [int(s in bond9.HEADS) for s in bond10.LET]):
        number = np.add.outer(labels, labels).ravel()
        assert np.max(abs(new*(number[:, None]-number[None, :]))) == 0
    # Field-only commutator with vacancy hopping: |L1,V,L0,V> -> |V,L1,L0,V>.
    # Old contact field differs by J f/2; new charge field is identical.
    def old_field(word):
        value = 0.
        for a, b in zip(word, word[1:]+word[:1]):
            if a != 'V' and b != 'V':
                value += p['J']*p['f']/4*(bond10.CHARGE[a]+bond10.CHARGE[b])
        return value
    left, right = ['L1', 'V', 'L0', 'V'], ['V', 'L1', 'L0', 'V']
    difference = old_field(right)-old_field(left)
    assert abs(difference-p['J']*p['f']/2) < 1e-16
    mp.mp.dps = 100
    b, n = 10, 9
    den0, den1 = mp.mpf(1), mp.mpf(1004)
    for _ in range(1, n+1):
        den0, den1 = den1, 1004*den1-den0
    beta = 1/den1
    ac4 = (mp.mpf(2)/(n+2))**2*mp.sin(mp.pi/(n+2))**4
    f9 = mp.power(2, -718)*ac4*beta**2*mp.power(2, -326-83)
    contact_half = mp.mpf(1)/2-2/mp.pi**2
    count = lambda H: dict(diagonal_entries=int(np.count_nonzero(np.diag(H))),
                           hermitian_offdiagonal_pairs=int(np.count_nonzero(H-np.diag(np.diag(H)))//2))
    return dict(charge_commutator_max=0., old_contact_field_counterexample=difference,
                original_matrix_unit_count=count(h), external_matrix_unit_count=count(new),
                n9_log10_beta=float(mp.log10(beta)), n9_log10_f9=float(mp.log10(f9)),
                half_density_contact=float(contact_half),
                log10_D_over_smallest_external_field_entry=float(mp.log10(4000/(contact_half*f9))))


def literal_single_cell_sources(b, heads, bits, p):
    words = []
    for h, payload, r in product(heads, bits, range(b)):
        word = [f'L{x}' for x in payload[:r]]+[h]+[f'R{x}' for x in payload[r:]]
        words.append(tuple(bond9.IDX[x] for x in word))
    index = {w: i for i, w in enumerate(words)}
    bond = bond9.full_bond(p)
    out = [np.zeros((len(words),)*2, complex) for _ in range(b)]
    for col, word in enumerate(words):
        for j in range(b):
            k = (j+1) % b
            ab = word[j]*9+word[k]
            for cd in np.flatnonzero(abs(bond[:, ab]) > 1e-14):
                changed = list(word); changed[j], changed[k] = divmod(int(cd), 9)
                assert tuple(changed) in index
                out[j][index[tuple(changed)], col] += bond[cd, ab]
    return out


def additional_literal_checks():
    p = dict(nu=selector(10), Delta=4-2*np.cos(np.pi/11), D=0.,
             g=.5, td=.4, lam=0., f=0., K=2., phiA=.27)
    charged = literal_single_cell_sources(3, ['A', 'B', 'XA', 'XB'],
                                          list(product((0, 1), repeat=2)), p)
    vacuum = sources([3], p['nu'])
    differences = []
    for mask in product((0, 1), repeat=3):
        hc = sum((a*h for a, h in zip(mask, charged)), np.zeros_like(charged[0]))
        hv = sum((a*h.toarray() for a, h in zip(mask, vacuum)), np.zeros((3, 3)))
        difference = float(la.eigvalsh(hc)[0]-la.eigvalsh(hv)[0])
        assert difference > -2e-14
        differences.append(difference)
    # Reproduce the referee's position-dependent contact field in a small ring.
    p = dict(p, nu=selector(4), g=0., f=0.)
    Hc, contacts = bath(4, 3)
    def old_energy(bits, f):
        hs = literal_single_cell_sources(4, ['W'], [bits], dict(p, f=f))
        h = 25*sp.kron(Hc, sp.eye(4), format='csr')
        for j in range(4):
            h += sp.kron(sp.diags(contacts[:, j].astype(float)), sp.csr_matrix(hs[j]), format='csr')
        return float(la.eigvalsh(h.toarray())[0])
    base = old_energy((0, 0, 0), 0.)
    shifts = []
    for j in range(3):
        bits = tuple(int(k == j) for k in range(3))
        assert abs(old_energy(bits, 0.)-base) < 1e-12
        shifts.append((old_energy(bits, .001)-base)/.001)
    rho = 4/7
    alpha_infinite = rho-np.sin(np.pi*rho)**2/(np.pi**2*rho)
    return dict(charged_contact_masks=8, charged_test_dimension=48,
                charged_minus_vacuum= differences,
                old_field_shifts_divided_by_f=shifts,
                old_free_contact=alpha(4, 3),
                external_field_shifts_divided_by_f=[alpha_infinite]*3,
                external_field_reason='Exact fixed-Q scalar, using the size-independent density normalization')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify-only', action='store_true')
    args = parser.parse_args()
    started = time.perf_counter()
    results = dict(status='finite integer certificates and finite-volume response; no thermodynamic theorem',
                   selector_enclosure=selector_rational_enclosure(),
                   literal_audit=literal_audit(), finite=finite_certificates(args.verify_only),
                   response=coefficient_checks(), relative_response_constants=relative_response_constants(),
                   form_factors=form_factor_checks(), charge=charge_audit(),
                   additional_literal_checks=additional_literal_checks())
    results['runtime_seconds'] = time.perf_counter()-started
    results['versions'] = dict(numpy=np.__version__, scipy=scipy.__version__)
    if not args.verify_only:
        (HERE/'checks.json').write_text(json.dumps(results, indent=2)+'\n')
    print(json.dumps(results, indent=2))


if __name__ == '__main__':
    main()
