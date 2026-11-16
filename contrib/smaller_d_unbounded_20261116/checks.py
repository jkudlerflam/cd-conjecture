"""Reproducible algebra and finite-matrix checks; no thermodynamic c fit."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')

import itertools
import json
from pathlib import Path
import time

import mpmath as mp
import numpy as np
import scipy.linalg as la

import model as m

HERE = Path(__file__).resolve().parent


def literal_checks():
    results = []
    for n in (1, 2, 3):
        p = m.parameters(n, g=.13, field=.001)
        for mobile in (False, True):
            letters = m.LETTERS if mobile else m.MATERIAL
            h = m.mobile_bond(p, lam=.17, J=.8, t=1.1, mu=.3) if mobile else m.material_bond(p, lam=.17)
            d = len(letters)
            q = np.array([m.charge(x) + m.charge(y) for x in letters for y in letters])
            occ = np.array([(x[0] != 'hole') + (y[0] != 'hole') for x in letters for y in letters])
            marker_number = np.array([int(m.marker(x)) + int(m.marker(y)) for x in letters for y in letters])
            rows, cols = np.nonzero(h)
            outside, inner = 0, 0
            for row, col in zip(rows, cols):
                if row == col:
                    continue
                a, b = letters[col // d], letters[col % d]
                c, e = letters[row // d], letters[row % d]
                outside += (m.arrows(a)[0], m.arrows(b)[1]) != (m.arrows(c)[0], m.arrows(e)[1])
                inner += not (m.allowed(a, b) and m.allowed(c, e))
            r = dict(n=n, dimension=d, nonzero_entries=int(len(rows)),
                     hermiticity=float(la.norm(h - h.conj().T)),
                     charge_commutator=float(np.max(np.abs((q[rows] - q[cols]) * h[rows, cols]))),
                     occupancy_commutator=float(np.max(np.abs((occ[rows] - occ[cols]) * h[rows, cols]))),
                     marker_number_commutator=float(np.max(np.abs((marker_number[rows] - marker_number[cols]) * h[rows, cols]))),
                     outer_arrow_failures=int(outside), internal_grammar_failures=int(inner))
            assert max(r[k] for k in ('hermiticity', 'charge_commutator', 'occupancy_commutator',
                                      'marker_number_commutator', 'outer_arrow_failures',
                                      'internal_grammar_failures')) < 1e-12
            results.append(r)
    return results


def grammar_checks():
    # Enumerate the grammar itself, without using the cell-word generator.
    results = []
    backgrounds = [x for x in m.MATERIAL if not m.marker(x)]
    markers = [x for x in m.MATERIAL if m.marker(x)]
    for n in range(1, 5):
        generated, _ = m.cell_words(n)
        direct = set()
        def visit(prefix, outgoing):
            if len(prefix) == n:
                for end in markers:
                    if m.arrows(end)[0] == outgoing:
                        direct.add(tuple(prefix + [end]))
                return
            for letter in backgrounds:
                if m.arrows(letter)[0] == outgoing:
                    visit(prefix + [letter], m.arrows(letter)[1])
        visit([], 0)
        assert direct == set(generated)
        results.append(dict(n=n, legal_cell_dimension=len(direct), expected=(9*n+5)*2**n))
    # Each legal material pair has one possible color for any intervening
    # nonempty vacancy string. A forbidden pair has none.
    counts = []
    for length in (1, 2, 3):
        for left, right in itertools.product(m.MATERIAL, repeat=2):
            colorings = 0
            for colors in itertools.product((0, 1), repeat=length):
                word = [left] + [('hole', c) for c in colors] + [right]
                colorings += all(m.allowed(a, b) for a, b in zip(word, word[1:]))
            assert colorings == int(m.allowed(left, right))
        counts.append(dict(vacancies=length, material_pairs=400, unique_coloring=True))
    return dict(cells=results, vacancy_coloring=counts)


def mp_path(n, q):
    H = mp.zeros(n+1)
    for j in range(n+1):
        H[j,j] = q*q + (j > 0)
        if j:
            H[j,j-1] = H[j-1,j] = -q
    e, u = mp.eigsy(H)
    v = u[:,0]
    if v[0] < 0:
        v = -v
    return e, v


def selector_checks():
    mp.mp.dps = 160
    results = []
    for q in (mp.mpf(1)/4, mp.mpf(1)/8, mp.mpf(1)/32):
        energy = [mp_path(n, q)[0][0] for n in range(25)]
        lower_ratios, upper_ratios, gaps = [], [], []
        for n, e in enumerate(energy):
            lower = (1-q*q)**2*q**(2*n+2)
            upper = q**(2*n+2)
            assert lower <= e * (1+mp.mpf('1e-110')) and e <= upper * (1+mp.mpf('1e-110'))
            lower_ratios.append(e/lower)
            upper_ratios.append(e/upper)
        for n in range(1, 13):
            mu = -(n+2)*energy[n]
            excess = [energy[k]+mu+(k+1)*energy[n] for k in range(25) if k != n]
            assert min(excess) > 0
            # The analytic argument covers lengths beyond this numerical table.
            assert abs(min(excess)/energy[n+1]-1) < mp.mpf('1e-100')
            gaps.append(mp.nstr(energy[n+1], 18))
        results.append(dict(q=mp.nstr(q, 10), lengths_checked=25, targets_checked=12,
                            minimum_ratio_to_lower=mp.nstr(min(lower_ratios), 18),
                            maximum_ratio_to_upper=mp.nstr(max(upper_ratios), 18),
                            wrong_cell_gaps=gaps))
    return results


def cell_checks():
    out = []
    for n in range(1,5):
        p = m.parameters(n,q=.25)
        words, info, H, C, G, _ = m.cell_matrices(p)
        eigenvalues = la.eigvalsh((H+.13*C).toarray())
        expected = p['E']+p['selector']
        multiplicity = int(np.sum(np.abs(eigenvalues-expected)<1e-10))
        gap = eigenvalues[2**n]-expected
        assert multiplicity == 2**n and gap >= .5-1e-12
        readout = np.array([word.count(('L',1)) for word in words])
        projected = G.T@(readout[:,None]*G)
        gamma = [sum(p['vacuum'][j]**2 for j in range(i+1,n+1)) for i in range(n)]
        target = np.diag([sum(x*y for x,y in zip(bits,gamma))
                          for bits in itertools.product((0,1),repeat=n)])
        error = la.norm(projected-target)
        assert error < 1e-12
        out.append(dict(n=n,cell_dimension=len(words),ground_multiplicity=multiplicity,
                        auxiliary_gap=float(gap),readout_error=float(error)))
    return out


def ring_checks():
    """Independent assembly on full legal no-excited-token rings."""
    letters = [x for x in m.MATERIAL if not m.excited(x)]
    out = []
    for n,L in ((1,2),(1,4),(2,3),(3,4)):
        p = m.parameters(n,q=.25)
        h = m.material_bond(p,lam=0)
        words = [w for w in itertools.product(letters,repeat=L)
                 if all(m.allowed(w[j],w[(j+1)%L]) for j in range(L))]
        index = {w:i for i,w in enumerate(words)}
        H = np.zeros((len(words),len(words)),complex)
        for col,word in enumerate(words):
            for j in range(L):
                k = (j+1)%L
                local = m.INDEX[word[j]]*20+m.INDEX[word[k]]
                for row in np.flatnonzero(h[:,local]):
                    a,b = m.MATERIAL[row//20],m.MATERIAL[row%20]
                    assert not m.excited(a) and not m.excited(b)
                    changed = list(word);changed[j]=a;changed[k]=b
                    H[index[tuple(changed)],col] += h[row,local]
        energies = la.eigvalsh(H)
        expected = -L*p['E']
        multiplicity = int(np.sum(np.abs(energies-expected)<1e-10))
        expected_mult = (n+1)*2**(n*(L//(n+1)))
        assert abs(energies[0]-expected)<1e-11 and multiplicity == expected_mult
        out.append(dict(n=n,L=L,legal_no_excited_dimension=len(words),
                        energy_error=float(abs(energies[0]-expected)),
                        ground_multiplicity=multiplicity,expected_multiplicity=expected_mult))
    return out


def pair_checks():
    results = []
    for n, q, w, v, eta, D in [(1,.2,.7,.6,.35,8.), (2,.2,.7,.6,.35,8.),
                               (2,.13,.45,.3,.2,5.)]:
        p = m.parameters(n, q=q, w=w, v=v, eta=eta, D=D)
        H22, H20, diagnostics = m.pair_coefficient(p)
        matrix, residual = m.gaussian_decomposition(H22, 2*n)
        overlap = p['vacuum'][0]**2*p['vacuum'][-1]**2
        target = -2*overlap*(m.scalar_exchange(n,w,D) + 1j*eta*m.scalar_exchange(n,w*v,D))
        hopping_error = la.norm(matrix[:n,n:]-target)/la.norm(target)
        onsite_error = (la.norm(matrix[:n,:n])+la.norm(matrix[n:,n:]))/la.norm(target)
        zeroth = la.norm(H20-H20[0,0]*np.eye(len(H20)))/la.norm(H20)
        assert residual < 1e-11 and hopping_error < 1e-10 and onsite_error < 1e-8 and zeroth < 1e-12
        results.append(dict(n=n, q=q, w=w, v=v, eta=eta, D=D, **diagnostics,
                            non_gaussian_relative_residual=residual,
                            scalar_exchange_relative_error=float(hopping_error),
                            onsite_relative_residual=float(onsite_error),
                            zeroth_order_nonscalar_relative_residual=float(zeroth)))
    return results


def gram_limit_checks():
    out = []
    for n in range(1,5):
        gamma = m.gram(n,.7)
        previous = 1.
        for D in (100.,1000.,10000.):
            scaled = 2*D**(2*n+2)*m.scalar_exchange(n,.7,D)
            error = la.norm(scaled-gamma)/la.norm(gamma)
            assert error < previous
            previous = error
            out.append(dict(n=n,D=D,relative_error=float(error)))
    return out


def mp_gram(n, w):
    A = mp.zeros(n)
    for i in range(n):
        A[i,i] = 4
        if i:
            A[i,i-1] = A[i-1,i] = -w
    K = mp.matrix(n*n,n*n)
    for i in range(n*n):
        for j in range(n*n):
            K[i,j] = A[i//n,j//n]*(i%n == j%n) + A[i%n,j%n]*(i//n == j//n)
    R = K**-1
    z = R[:,0]
    return mp.matrix([[2*mp.fsum(z[i*n+y]*R[i*n+y,x*n+j]*z[x*n+j]
                                for x in range(n) for y in range(n))
                       for j in range(n)] for i in range(n)])


def band_check(n):
    """High-precision numerical roots, not interval certificates."""
    mp.mp.dps = 100+45*n
    w, u, v, eta = mp.mpf(1)/4, mp.mpf(1)/64, mp.mpf(1)/1024, mp.mpf(1)/100
    q = (w*u)**2
    nu = mp.mpf(10)**(-24*n-40)
    gram1, gram2 = mp_gram(n,w), mp_gram(n,w*v)
    A = mp.matrix([[gram1[i,j]/w**(2*i+2*j) for j in range(n)] for i in range(n)])
    B = mp.matrix([[gram2[i,j]/w**(2*i+2*j) for j in range(n)] for i in range(n)])
    energies, vacuum = mp_path(n,q)
    E = mp.diag([mp.fsum(vacuum[j]**2 for j in range(i+1,n+1))/(q*q*w**(4*i))
                 for i in range(n)])
    L = mp.cholesky(A)
    Li = L**-1
    Bwhite, Ewhite = Li*B*Li.T, Li*E*Li.T
    rs, vecs = mp.eigsy(Bwhite)
    rs = list(reversed(list(rs)))
    vecs = mp.matrix([[vecs[i,n-1-j] for j in range(n)] for i in range(n)])
    theta = [-mp.atan(eta*r) for r in rs]
    slopes = [(vecs[:,j].T*Ewhite*vecs[:,j])[0]/(2*mp.sqrt(1+(eta*rs[j])**2))
              for j in range(n)]
    if n > 1:
        assert max(abs(theta[j+1]/theta[j]) for j in range(n-1)) < mp.mpf(1)/32
        assert max(slopes[j+1]/slopes[j] for j in range(n-1)) < mp.mpf(1)/32
    identity = mp.eye(n)
    Q2, Q1, Q0 = 2*identity+nu*Ewhite, 4*eta*Bwhite, -2*identity+nu*Ewhite
    def determinant(t):
        return mp.det(Q2*t*t+Q1*t+Q0)
    roots = []
    residuals = []
    for j in range(n):
        for sign in (-1,1):
            guess_k = theta[j]+sign*mp.pi/2-sign*nu*slopes[j]
            guess = mp.tan(guess_k/2)
            step = nu*slopes[j]/10
            root = mp.findroot(determinant,(guess-step,guess+step),solver='secant',
                               tol=mp.mpf(10)**(-mp.mp.dps+35),maxsteps=40,verify=False)
            # Normwise spectral residual on the A-whitened Bloch matrix.
            k = 2*mp.atan(root)
            eig = mp.eigsy(-2*identity*mp.cos(k)+2*eta*Bwhite*mp.sin(k)+nu*Ewhite,eigvals_only=True)
            residuals.append(min(abs(x) for x in eig))
            assert abs(k-guess_k) < nu*slopes[j]/100
            roots.append(k)
    def distance(x):
        return abs(x-2*mp.pi*mp.nint(x/(2*mp.pi)))
    delta2 = min(distance(a-b) for i,a in enumerate(roots) for b in roots[:i])
    sums = [roots[i]+roots[j] for i,j in itertools.combinations_with_replacement(range(2*n),2)]
    delta4 = min(distance(a-b) for i,a in enumerate(sums) for b in sums[:i])
    assert delta2 > 0 and delta4 > 0
    assert max(residuals) < delta4*mp.mpf('1e-30')
    return dict(n=n,precision_digits=mp.mp.dps,root_count=len(roots),pair_sums=len(sums),
                q=mp.nstr(q,20),nu=mp.nstr(nu,12),vacuum_energy=mp.nstr(energies[0],24),
                smallest_center=mp.nstr(abs(theta[-1]),24),smallest_slope=mp.nstr(slopes[-1],24),
                delta2=mp.nstr(delta2,24),delta4=mp.nstr(delta4,24),
                max_spectral_residual=mp.nstr(max(residuals),12),
                roots=[mp.nstr(k,mp.mp.dps-25) for k in roots],
                status='high precision numerical; not an interval certificate')


def main():
    started = time.time()
    results = dict(versions=dict(numpy=np.__version__,scipy=__import__('scipy').__version__,mpmath=mp.__version__))
    for name, function in [('literal',literal_checks),('grammar',grammar_checks),
                           ('selector',selector_checks),('cells',cell_checks),('rings',ring_checks),
                           ('pair',pair_checks),('gram_limit',gram_limit_checks)]:
        results[name] = function()
        print(name, 'passed', flush=True)
        (HERE/'checks.json').write_text(json.dumps(results,indent=2)+'\n')
    results['bands'] = []
    for n in range(1,7):
        result = band_check(n)
        results['bands'].append(result)
        print('bands',n,result['root_count'],'roots; delta4',result['delta4'],flush=True)
        (HERE/'checks.json').write_text(json.dumps(results,indent=2)+'\n')
    results['runtime_seconds'] = time.time()-started
    (HERE/'checks.json').write_text(json.dumps(results,indent=2)+'\n')
    print('all checks passed in',round(results['runtime_seconds'],2),'seconds',flush=True)


if __name__ == '__main__':
    main()
