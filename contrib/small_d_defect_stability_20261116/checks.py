"""Independent literal frozen-bond/sector checks; no effective Hamiltonian.

Run with Python + numpy, scipy, mpmath. Writes small checks.json next to
this file. The parent model is imported only for a final entrywise bond
comparison; all enumeration, assembly, and spectra here are independent.
"""
import os
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")

from itertools import product
import importlib.util
import json
from pathlib import Path
import platform
import time

import mpmath as mp
import numpy as np
import scipy
import scipy.linalg as la
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh

HERE = Path(__file__).resolve().parent
LETTERS = ([("L", a) for a in (0, 1)] + [("R", a) for a in (0, 1)]
           + [("v",)]
           + [(rail, r, s) for rail in ("d", "x") for r in (0, 1)
              for s in (-1, 1)]
           + [("e", a) for a in (0, 1)] + [("m",)]
           + [("m", r, s) for r in (0, 1) for s in (-1, 1)])
INDEX = {x: i for i, x in enumerate(LETTERS)}
MARK = np.array([x[0] in ("e", "m") for x in LETTERS])
EXC = np.array([len(x) == 3 for x in LETTERS])
IN = np.array([int(x[0] in ("R", "e")) for x in LETTERS])
OUT = np.array([int(x[0] in ("R", "v", "d", "x")) for x in LETTERS])
CHARGE = np.array([2*x[1] if x[0] in ("L", "R", "e")
                   else x[2] if len(x) == 3 else 0 for x in LETTERS])


def qmatrix(k, q):
    # Construct from edge projectors, not from the parent's tridiagonal.
    mat = np.zeros((k+1, k+1))
    for j in range(k):
        edge = np.zeros(k+1)
        edge[j], edge[j+1] = q, -1
        mat += np.outer(edge, edge)
    mat[k, k] += q*q
    return mat


def config(n, *, field=None, lam=None):
    q, eta = .25, .25
    energies, vec = la.eigh(qmatrix(n, q))
    en = float(energies[0])
    low = (1-q*q)**2*q**(2*n+4)
    port = low/(10*np.sqrt(2+2*eta)) if lam is None else lam
    p = dict(n=n, q=q, w=.4, v=.35, eta=eta, D=3., g=.13,
             Delta=4+en, mu=-(n+2)*en, K=1., lam=port, E=en,
             gamma=float(la.eigvalsh(qmatrix(n+1, q))[0]), gamma_low=low)
    p["field"] = (port**2*p["g"]**2*abs(vec[0, 0]*vec[-1, 0])**2
                  *p["D"]**(-2*n-2)*.01/q**2 if field is None else field)
    p["B"] = np.sqrt(2+2*eta)*abs(port)
    return p


def literal_bond(p):
    """Rebuild all 400x400 entries directly from the displayed rule table."""
    kinetic = np.zeros((400, 400), complex)
    port = np.zeros_like(kinetic)

    def state(pair):
        return 20*INDEX[pair[0]]+INDEX[pair[1]]

    def hop(source, target, amplitude, into=kinetic):
        i, j = state(source), state(target)
        into[j, i] += amplitude
        into[i, j] += np.conjugate(amplitude)

    for a in (0, 1):
        for source, target in (((('v',), ('R', a)), (('L', a), ('v',))),
                               ((('v',), ('e', a)), (('L', a), ('m',)))):
            edge = np.zeros(400)
            edge[state(source)], edge[state(target)] = p['q'], -1
            kinetic += np.outer(edge, edge)
        for r, s in product((0, 1), (-1, 1)):
            parity = (-1)**a if s == -1 else 1
            phase = np.exp(1j*np.pi/(2*p['n'])) if s == -1 else 1
            hop((('d', r, s), ('R', a)), (('L', a), ('d', r, s)),
                -p['w']*(p['v'] if r else 1)*parity)
            hop((('x', r, s), ('R', a)), (('L', a), ('x', r, s)),
                -phase*parity)
            hop((('x', r, s), ('e', a)), (('L', a), ('m', r, s)),
                -phase*parity)
    for r, s in product((0, 1), (-1, 1)):
        for mark in [x for x in LETTERS if x[0] in ('e', 'm')]:
            hop((mark, ('d', r, s)), (mark, ('x', r, s)), -1)
    for r in (0, 1):
        for kind in ('R', 'e'):
            hop((('d', r, -1), (kind, 1)), (('d', r, 1), (kind, 0)), p['g'])
        weight = np.sqrt(p['eta']) if r else 1
        hop((('m',), ('v',)), (('m', r, -1), ('d', r, 1)),
            weight*(-1 if r else -1j), port)
        hop((('m',), ('v',)), (('m', r, 1), ('d', r, -1)), weight, port)
    onsite = np.zeros(20)
    for i, x in enumerate(LETTERS):
        onsite[i] = (p['q']**2*(x == ('m',)) + p['mu']*MARK[i]
                     + p['field']*(x == ('L', 1)))
        if EXC[i]:
            onsite[i] += p['Delta'] + p['D']*(x[0] in ('x', 'm'))
    forbidden = OUT[:, None] != IN[None, :]
    diagonal = ((onsite[:, None]+onsite[None, :])/2+p['K']*forbidden)
    h = kinetic + p['lam']*port + np.diag(diagonal.ravel())
    return h, port, kinetic, onsite


def sector_basis(length, markers, bad):
    """Enumerate literal site letters with fixed markers and bad-link bits."""
    choices = [np.flatnonzero(MARK == bool(markers & (1 << j))).tolist()
               for j in range(length)]
    words = []

    def visit(word):
        j = len(word)
        if j == length:
            if bool(OUT[word[-1]] != IN[word[0]]) == bool(bad & (1 << (length-1))):
                words.append(tuple(word))
            return
        for a in choices[j]:
            if j and bool(OUT[word[-1]] != IN[a]) != bool(bad & (1 << (j-1))):
                continue
            visit(word+[a])
    visit([])
    return words


def assemble(words, local, *, periodic=True, onsite=None, junction=False):
    """Apply the literal two-site rules to site words, with closure assertions."""
    size = len(words)
    index = {word: j for j, word in enumerate(words)}
    columns = [[(int(row), local[row, col]) for row in np.flatnonzero(local[:, col])]
               for col in range(400)]
    rows, cols, vals = [], [], []
    for col, word in enumerate(words):
        length = len(word)
        if onsite is not None:
            rows.append(col); cols.append(col); vals.append(sum(onsite[list(word)]))
        for j in range(length if periodic else length-1):
            right = (j+1) % length
            for row, value in columns[20*word[j]+word[right]]:
                changed = list(word)
                changed[j], changed[right] = divmod(row, 20)
                target = tuple(changed)
                assert target in index, (word, target, j)
                rows.append(index[target]); cols.append(col); vals.append(value)
        if junction and LETTERS[word[0]][0] in ('d', 'x'):
            x = LETTERS[word[0]]
            swapped = ('x' if x[0] == 'd' else 'd', x[1], x[2])
            changed = (INDEX[swapped],) + word[1:]
            assert changed in index
            rows.append(index[changed]); cols.append(col); vals.append(-1.)
    mat = sp.coo_matrix((vals, (rows, cols)), shape=(size, size)).tocsr()
    assert sp.linalg.norm(mat-mat.getH()) < 1e-11
    return mat


def minimum(mat):
    if mat.shape[0] <= 100:
        vals, vec = la.eigh(mat.toarray(), subset_by_index=[0, 0])
        value, vector = float(vals[0]), vec[:, 0]
    else:
        rng = np.random.default_rng(615)
        start = rng.standard_normal(mat.shape[0])
        vals, vec = eigsh(mat, k=1, which='SA', tol=2e-12, maxiter=30000, v0=start)
        value, vector = float(vals[0]), vec[:, 0]
    residual = float(la.norm(mat@vector-value*vector))
    assert residual < 2e-9
    return value, residual


def pattern(lengths):
    pos, mask = 0, 0
    for k in lengths:
        pos += k+1
        mask |= 1 << (pos-1)
    return pos, mask


def defect_counts(length, markers, bad, n):
    """Count complete cells solely from the cyclic marker/cut pattern."""
    positions = [j for j in range(length) if markers & (1 << j)]
    wrong = 0
    if positions:
        for a, z in zip(positions, positions[1:]+[positions[0]+length]):
            # A cell is complete only if its preceding-marker link and
            # all links through its terminating marker are legal.
            if not any(bad & (1 << (j % length)) for j in range(a, z)):
                wrong += z-a-1 != n
    return int(wrong), bad.bit_count()


def selector_checks():
    mp.mp.dps = 110
    reports = []
    for q in (mp.mpf(1)/4, mp.mpf(1)/8, mp.mpf(1)/32):
        energies = []
        gaps = []
        for k in range(19):
            mat = mp.zeros(k+1)
            for j in range(k):
                mat[j, j] += q*q
                mat[j+1, j+1] += 1
                mat[j, j+1] -= q
                mat[j+1, j] -= q
            mat[k, k] += q*q
            spectrum = mp.eigsy(mat, eigvals_only=True)
            e = spectrum[0]
            energies.append(e)
            assert e >= (1-q*q)**2*q**(2*k+2)*(1-mp.mpf('1e-45'))
            assert e <= q**(2*k+2)*(1+mp.mpf('1e-45'))
            if k:
                assert spectrum[1]-e >= 1-2*q
                gaps.append(spectrum[1]-e)
        for n in range(1, 13):
            excess = [(energies[k]+(k-n-1)*energies[n], k)
                      for k in range(19) if k != n]
            val, arg = min(excess)
            assert arg == n+1
            assert abs(val/energies[n+1]-1) < mp.mpf('1e-45')
        reports.append(dict(q=str(q), lengths=19, targets=12,
                            minimum_auxiliary_gap=mp.nstr(min(gaps), 20)))
    return reports


def local_checks():
    spec = importlib.util.spec_from_file_location(
        'parent_model', HERE.parent/'smaller_d_unbounded_20261116'/'model.py')
    parent = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(parent)
    results = []
    for n in (1, 2, 3):
        p = config(n, field=.003, lam=.17)
        h, port, _, _ = literal_bond(p)
        pp = parent.parameters(n, q=p['q'], w=p['w'], v=p['v'], eta=p['eta'],
                               D=p['D'], g=p['g'], K=p['K'], field=p['field'])
        parent_h = parent.material_bond(pp, lam=p['lam'])
        assert LETTERS == parent.MATERIAL
        difference = float(np.max(np.abs(h-parent_h)))
        assert difference < 1e-14
        count = 0
        for row, col in zip(*np.nonzero(h)):
            if row == col:
                continue
            a, b = divmod(col, 20)
            c, d = divmod(row, 20)
            assert MARK[a] == MARK[c] and MARK[b] == MARK[d]
            assert IN[a] == IN[c] and OUT[b] == OUT[d]
            assert OUT[a] == IN[b] and OUT[c] == IN[d]
            assert CHARGE[a]+CHARGE[b] == CHARGE[c]+CHARGE[d]
            count += 1
        norm = float(la.eigvalsh(port)[-1])
        assert abs(norm-np.sqrt(2+2*p['eta'])) < 1e-14
        results.append(dict(n=n, offdiagonal_entries=count,
                            individual_marker_and_adjacent_bad_projectors=True,
                            bond_difference=difference, port_norm=norm))
    return results


def cell_checks():
    results = []
    for k in range(6):
        p = config(2, field=.003)
        _, _, kinetic, onsite = literal_bond(p)
        # Legal cell = start arrow 0, one final marker, no other markers.
        words = [w for w in sector_basis(k+1, 1 << k, 0) if IN[w[0]] == 0]
        assert len(words) == (9*k+5)*2**k
        mat = assemble(words, kinetic, periodic=False, onsite=onsite, junction=True)
        energy, residual = minimum(mat)
        expected = p['mu'] + float(la.eigvalsh(qmatrix(k, p['q']))[0])
        assert abs(energy-expected) < 2e-11
        results.append(dict(k=k, dimension=len(words), ground_energy=energy,
                            expected=expected, residual=residual))
    return results


def all_sector_checks():
    """Exhaust all d^L letters at L=2,3; compare on compatible rings."""
    reports = []
    for n, length in ((1, 2), (2, 3)):
        p = config(n, field=.003)
        h = literal_bond(p)[0]
        _, selected_mask = pattern([n]*(length//(n+1)))
        selected, _ = minimum(assemble(sector_basis(length, selected_mask, 0), h))
        minimum_slack = float('inf')
        minimum_illegal_slack = float('inf')
        total_dim, sectors, max_residual = 0, 0, 0.
        for markers, bad in product(range(1 << length), repeat=2):
            words = sector_basis(length, markers, bad)
            if not words:
                continue
            sectors += 1
            total_dim += len(words)
            energy, residual = minimum(assemble(words, h))
            max_residual = max(max_residual, residual)
            nw, ni = defect_counts(length, markers, bad, n)
            if nw+ni:
                bound = nw*(p['gamma']-5*p['B'])+ni*(p['K']-2*abs(p['mu'])-7*p['B'])
                slack = energy-selected-bound
                assert slack > -2e-10
                minimum_slack = min(minimum_slack, slack)
                if ni:
                    minimum_illegal_slack = min(minimum_illegal_slack, slack)
            elif markers == 0:
                assert energy >= -2e-10 and energy > selected
            else:
                assert abs(energy-selected) < 2e-10
        assert total_dim == 20**length
        reports.append(dict(n=n, L=length, sectors=sectors,
                            total_dimension=total_dim, selected_energy=selected,
                            min_comparison_slack=minimum_slack,
                            min_illegal_slack=minimum_illegal_slack,
                            max_eigenpair_residual=max_residual))
    return reports


def mixed_ring_checks():
    reports = []
    cases = [(1, (1, 1), [(0, 2), (3,), (0, 0, 0, 0)]),
             (1, (1, 1, 1), [(1, 0, 2), (1, 3), (0, 4), (5,)]),
             (2, (2, 2), [(1, 3), (0, 4), (5,), (1, 1, 1)]),
             (3, (3,), [(0, 2), (1, 1)])]
    for n, selected_pattern, competitors in cases:
        p = config(n)  # Nonzero field with the parent's lambda^2 scale.
        h = literal_bond(p)[0]
        length, selected_mask = pattern(selected_pattern)
        selected_words = sector_basis(length, selected_mask, 0)
        selected, selected_residual = minimum(assemble(selected_words, h))
        for ks in competitors:
            other_length, mask = pattern(ks)
            assert length == other_length
            words = sector_basis(length, mask, 0)
            energy, residual = minimum(assemble(words, h))
            wrong = [k for k in ks if k != n]
            delta = sum(float(la.eigvalsh(qmatrix(k, p['q']))[0])
                        +(k-n-1)*p['E'] for k in wrong)
            bound = delta-5*len(wrong)*p['B']
            assert energy-selected >= bound-2e-10
            reports.append(dict(n=n, L=length, pattern=list(ks),
                                dimension=len(words), selected_dimension=len(selected_words),
                                field=p['field'], lam=p['lam'],
                                energy_difference=energy-selected,
                                stronger_comparison_bound=bound,
                                slack=energy-selected-bound,
                                residual=max(residual, selected_residual)))
    return reports


def boundary_checks():
    """Finite instances of the almost-additive step, with resolved dressing.

    Large diagnostic lambda here is outside the sufficient defect schedule.
    Boundary almost-additivity itself holds at every lambda.
    """
    p = config(1, field=.003, lam=.08)
    _, port, kinetic, onsite = literal_bond(p)
    h = kinetic+p['lam']*port
    energies = {}
    for cells in (1, 2, 3):
        length, mask = pattern([1]*cells)
        words = sector_basis(length, mask, 0)
        mat = assemble(words, h, periodic=False, onsite=onsite, junction=True)
        energies[cells] = minimum(mat)[0]
    violations = [abs(energies[a+b]-energies[a]-energies[b])
                  for a, b in ((1, 1), (1, 2))]
    assert max(violations) <= p['B']+1e-10
    return dict(lam=p['lam'], field=p['field'], B=p['B'],
                open_energies=energies, joining_energy_changes=violations,
                status='Boundary lemma diagnostic; outside defect certificate')


def main():
    start = time.time()
    results = dict(status='numerical corroboration of the analytical theorem; no CFT claim',
                   versions=dict(python=platform.python_version(), numpy=np.__version__,
                                 scipy=scipy.__version__, mpmath=mp.__version__))
    for key, function in [('selector', selector_checks), ('local', local_checks),
                          ('cells', cell_checks), ('all_small_sectors', all_sector_checks),
                          ('mixed_rings', mixed_ring_checks), ('boundary', boundary_checks)]:
        results[key] = function()
        print(key, 'passed', flush=True)
    results['runtime_seconds'] = time.time()-start
    (HERE/'checks.json').write_text(json.dumps(results, indent=2)+'\n')
    print('all checks passed in', round(results['runtime_seconds'], 2), 'seconds')


if __name__ == '__main__':
    main()
