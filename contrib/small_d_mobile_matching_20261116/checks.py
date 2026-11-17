"""Literal mobile/gap checks and charged-rail pin diagnostics.

No finite-volume perturbation denominator is used as a uniform bound.
All saved eigenvalues are numerical corroboration, not certificates.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
from itertools import combinations, product
from pathlib import Path
import importlib.util
import json
import time

import numpy as np
import scipy
import scipy.linalg as la
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh

HERE = Path(__file__).resolve().parent


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


frozen = load('independent_frozen', HERE.parent/'small_d_defect_stability_20261116'/'checks.py')
parent = load('parent_mobile', HERE.parent/'smaller_d_unbounded_20261116'/'model.py')


def compositions(v, n):
    if n == 1:
        return [(v,)]
    return [(j,)+tail for j in range(v+1) for tail in compositions(v-j, n-1)]


def gap_matrices(n, v, phi=0.):
    gaps = compositions(v, n)
    index = {g: i for i, g in enumerate(gaps)}
    rows, cols, values = [], [], []
    for col, g in enumerate(gaps):
        for j in range(n):
            if g[j]:
                moved = list(g)
                moved[j] -= 1
                moved[(j-1) % n] += 1
                row = index[tuple(moved)]
                rows.extend((row, col)); cols.extend((col, row))
                values.extend((-np.exp(-1j*phi), -np.exp(1j*phi)))
    h = sp.coo_matrix((values, (rows, cols)), shape=(len(gaps), len(gaps))).tocsr()
    contacts = [sp.diags([float(g[j] == 0) for g in gaps], format='csr') for j in range(n)]
    return gaps, h, contacts


def free_energy(n, v):
    length = n+v
    return -2*np.sin(np.pi*n/length)/np.sin(np.pi/length)


def alpha(n, v):
    length = n+v
    rho = n/length
    return rho-np.sin(np.pi*rho)**2/(n*length*np.sin(np.pi/length)**2)


def ground(h, vector=False):
    if h.shape[0] < 100:
        energies, states = la.eigh(h.toarray(), subset_by_index=[0, 0])
    else:
        start = np.random.default_rng(771).standard_normal(h.shape[0])
        energies, states = eigsh(h, k=1, which='SA', tol=2e-11, maxiter=30000, v0=start)
    energy, state = float(energies[0]), states[:, 0]
    residual = float(la.norm(h@state-energy*state))
    assert residual < 3e-8
    return (energy, state, residual) if vector else energy


def material_words(ks, cyclic=False, vacuum_zero=False):
    n, mask = frozen.pattern(ks)
    masks = {mask}
    if cyclic:
        masks = {((mask << j) | (mask >> (n-j))) & ((1 << n)-1) for j in range(n)}
    words = sorted(set(w for m in masks for w in frozen.sector_basis(n, m, 0)))
    if vacuum_zero:
        allowed = {frozen.INDEX[x] for x in [('L', 0), ('R', 0), ('v',), ('e', 0), ('m',)]}
        words = [w for w in words if all(x in allowed for x in w)]
    return words


def material_sources(words, h):
    n = len(words[0])
    size = len(words)
    index = {w: i for i, w in enumerate(words)}
    columns = [[(int(row), h[row, col]) for row in np.flatnonzero(h[:, col])]
               for col in range(400)]
    sources = []
    for j in range(n):
        rows, cols, values = [], [], []
        for col, w in enumerate(words):
            right = (j+1) % n
            for row, value in columns[20*w[j]+w[right]]:
                changed = list(w)
                changed[j], changed[right] = divmod(row, 20)
                target = tuple(changed)
                assert target in index, (w, target, j)
                rows.append(index[target]); cols.append(col); values.append(value)
        sources.append(sp.coo_matrix((values, (rows, cols)), shape=(size, size)).tocsr())
    return sources


def mobile_gap(sources, v, J, phi=0.):
    n, d = len(sources), sources[0].shape[0]
    gaps, charge, contacts = gap_matrices(n, v, phi)
    h = sp.kron(charge, sp.eye(d), format='csr')
    for a, source in zip(contacts, sources):
        h += J*sp.kron(a, source, format='csr')
    return gaps, h


def literal_mobile(p, J=.08, t=1., mu=.17):
    material = frozen.literal_bond(p)[0]
    letters = frozen.LETTERS+[('hole', 0), ('hole', 1)]
    h = np.zeros((484, 484), complex)
    ids = [22*i+j for i in range(20) for j in range(20)]
    h[np.ix_(ids, ids)] = J*material
    for i in range(20):
        source = (20+int(frozen.IN[i]))*22+i
        target = i*22+20+int(frozen.OUT[i])
        h[source, target] = h[target, source] = -t
    ins = np.r_[frozen.IN, 0, 1]
    outs = np.r_[frozen.OUT, 0, 1]
    for i, j in product(range(22), repeat=2):
        if i >= 20 or j >= 20:
            h[22*i+j, 22*i+j] += J*p['K']*(outs[i] != ins[j])
        h[22*i+j, 22*i+j] -= mu*(int(i < 20)+int(j < 20))/2
    pp = parent.parameters(p['n'], q=p['q'], w=p['w'], v=p['v'], eta=p['eta'],
                           D=p['D'], g=p['g'], K=p['K'], field=p['field'])
    expected = parent.mobile_bond(pp, lam=p['lam'], J=J, t=t, mu=mu)
    assert np.max(np.abs(h-expected)) < 1e-13
    return letters, h


def physical_basis(material, v):
    n, length = len(material[0]), len(material[0])+v
    words = set()
    for occupied in combinations(range(length), n):
        for mat in material:
            word = [-1]*length
            for pos, x in zip(occupied, mat):
                word[pos] = x
            for j, pos in enumerate(occupied):
                following = occupied[(j+1) % n]
                cursor = (pos+1) % length
                while cursor != following:
                    word[cursor] = 20+int(frozen.OUT[mat[j]])
                    cursor = (cursor+1) % length
            words.add(tuple(word))
    return sorted(words)


def physical_matrix(words, h):
    index = {w: i for i, w in enumerate(words)}
    columns = [[(int(row), h[row, col]) for row in np.flatnonzero(h[:, col])]
               for col in range(484)]
    rows, cols, values = [], [], []
    for col, w in enumerate(words):
        for j in range(len(w)):
            right = (j+1) % len(w)
            for row, value in columns[22*w[j]+w[right]]:
                changed = list(w)
                changed[j], changed[right] = divmod(row, 22)
                assert tuple(changed) in index
                rows.append(index[tuple(changed)]); cols.append(col); values.append(value)
    return sp.coo_matrix((values, (rows, cols)), shape=(len(words), len(words))).tocsr()


def rotation_subspace(gaps, words, theta):
    basis = [(g, w) for g in gaps for w in words]
    index = {x: j for j, x in enumerate(basis)}
    visited = set()
    rows, cols, values = [], [], []
    count = 0
    for first in range(len(basis)):
        if first in visited:
            continue
        orbit, current = [], first
        while current not in orbit:
            orbit.append(current); visited.add(current)
            g, w = basis[current]
            current = index[(g[1:]+g[:1], w[1:]+w[:1])]
        size = len(orbit)
        if abs(np.exp(1j*theta*size)-1) > 1e-10:
            continue
        for j, row in enumerate(orbit):
            rows.append(row); cols.append(count); values.append(np.exp(-1j*theta*j)/np.sqrt(size))
        count += 1
    return sp.coo_matrix((values, (rows, cols)), shape=(len(basis), count)).tocsr()


def representation_checks():
    reports = []
    for n, vacancies in ((1, 1), (2, 1)):
        p = frozen.config(n, field=.003, lam=.031)
        J, mu = .08, .17
        _, local = literal_mobile(p, J=J, mu=mu)
        words = material_words([n], cyclic=True)
        sources = material_sources(words, frozen.literal_bond(p)[0])
        phys = physical_basis(words, vacancies)
        microscopic = physical_matrix(phys, local)
        literal_spectrum = la.eigvalsh(microscopic.toarray())
        occupied, length = n+1, n+1+vacancies
        all_spectrum = []
        for m in range(length):
            momentum = 2*np.pi*m/length
            gaps, rooted = mobile_gap(sources, vacancies, J, phi=momentum/occupied)
            rotation = rotation_subspace(gaps, words, momentum*length/occupied)
            compressed = (rotation.getH()@rooted@rotation).toarray()
            all_spectrum.extend(la.eigvalsh(compressed)-mu*occupied)
        error = float(np.max(np.abs(literal_spectrum-np.sort(all_spectrum))))
        assert error < 2e-11, error
        reports.append(dict(n=n, N=occupied, V=vacancies, physical_dimension=len(phys),
                            full_spectrum_error=error, onsite_contact_halves_retained=True))
    return reports


def scalar_checks():
    reports = []
    for n, v in ((3, 2), (3, 3), (4, 4), (5, 5), (6, 6)):
        gaps, h, contacts = gap_matrices(n, v)
        e0, state, _ = ground(h, vector=True)
        assert abs(e0-free_energy(n, v)) < 1e-11
        actual_a = float(np.real(np.vdot(state, contacts[0]@state)))
        assert abs(actual_a-alpha(n, v)) < 1e-11
        selected = np.array([i for i, g in enumerate(gaps) if g[0] >= 1])
        old = compositions(v-1, n)
        reorder = [gaps.index((g[0]+1,)+g[1:]) for g in old]
        compressed = h[reorder][:, reorder]
        reduced = gap_matrices(n, v-1)[1]
        error = float(sp.linalg.norm(compressed-reduced))
        assert error < 1e-13
        pressure = free_energy(n, v-1)-e0
        u = pressure/3
        energy, psi, residual = ground(h+u*contacts[0], vector=True)
        actual_contact = float(np.real(np.vdot(psi, contacts[0]@psi)))
        lower = u*pressure**2/(pressure**2+64)
        assert energy-e0 >= lower-1e-11
        reports.append(dict(N=n, V=v, dimension=len(gaps), alpha=actual_a,
                            compression_error=error, pressure_cost=pressure, u=u,
                            impurity_cost=energy-e0, lower_bound=lower,
                            contact=actual_contact, residual=residual))
    return reports


def replacement_checks():
    """Check the core pin inequality for every binary local contact pattern."""
    reports = []
    for n in (1, 2, 3):
        p = frozen.config(n, field=.003, lam=0.)
        p['D'], p['g'] = p['q']**2, .5  # Worst permitted extra mass and conversion.
        _, _, kinetic, onsite = frozen.literal_bond(p)
        words = material_words([n])
        index = {w: i for i, w in enumerate(words)}
        b = n+1
        pinned = np.array([i for i, w in enumerate(words) if any(frozen.EXC[x] for x in w)])
        kappa = 1-p['g']-p['q']**2
        minimum = float('inf')
        for contact in product((0, 1), repeat=b+1):
            h = np.zeros((len(words), len(words)), complex)
            for col, w in enumerate(words):
                h[col, col] += sum((contact[j]+contact[j+1])*onsite[w[j]]/2 for j in range(b))
                for j in range(b-1):
                    localcol = 20*w[j]+w[j+1]
                    for row in np.flatnonzero(kinetic[:, localcol]):
                        changed = list(w)
                        changed[j], changed[j+1] = divmod(int(row), 20)
                        h[index[tuple(changed)], col] += contact[j+1]*kinetic[row, localcol]
                x = frozen.LETTERS[w[0]]
                if x[0] in ('d', 'x'):
                    other = ('x' if x[0] == 'd' else 'd', x[1], x[2])
                    changed = (frozen.INDEX[other],)+w[1:]
                    h[index[changed], col] -= contact[0]
            replacement_diagonal, head_contacts = [], []
            for i in pinned:
                w = words[i]
                pos = next(j for j in range(b) if frozen.EXC[w[j]])
                if pos == n:
                    target = [('L', 0)]*n+[('m',)]
                else:
                    target = [('L', 0)]*pos+[('v',)]+[('R', 0)]*(n-pos-1)+[('e', 0)]
                target = tuple(frozen.INDEX[x] for x in target)
                replacement_diagonal.append(h[index[target], index[target]].real)
                head_contacts.append(contact[pos+1])
            difference = h[np.ix_(pinned, pinned)]-np.diag(replacement_diagonal)-kappa*np.diag(head_contacts)
            smallest = float(la.eigvalsh(difference, subset_by_index=[0, 0])[0])
            minimum = min(minimum, smallest)
            assert smallest >= -2e-12
        reports.append(dict(n=n, charged_dimension=len(pinned), contact_patterns=2**(b+1),
                            D=p['D'], g=p['g'], kappa=kappa, minimum_operator_slack=minimum))
    return reports


def channel_checks():
    """Exact Kraus identities with a quantum head and a vacancy-number cutoff."""
    n, vmax, phi = 3, 3, .37
    gaps = [g for v in range(vmax+1) for g in compositions(v, n)]
    index = {g: i for i, g in enumerate(gaps)}
    # A three-position head is sufficient: all identities tensor with its rails/bits.
    d = n
    hgap = sp.block_diag([gap_matrices(n, v, phi)[1] for v in range(vmax+1)], format='csr')
    h = sp.kron(hgap, sp.eye(d), format='csr')
    diag = np.array([float(g[j] == 0) for g in gaps for j in range(d)])
    P = sp.diags(diag, format='csr')
    Q = sp.eye(len(diag), format='csr')-P
    rows, cols = [], []
    for i, g in enumerate(gaps):
        for j in range(d):
            if g[j]:
                lower = list(g); lower[j] -= 1
                rows.append(index[tuple(lower)]*d+j); cols.append(i*d+j)
    lower = sp.coo_matrix((np.ones(len(rows)), (rows, cols)), shape=h.shape).tocsr()
    number = sp.diags([sum(g) for g in gaps for _ in range(d)], dtype=float, format='csr')
    errors = dict(completeness=float(sp.linalg.norm(P@P+lower.getH()@lower-sp.eye(h.shape[0]))),
                  kinetic_identity=float(sp.linalg.norm(P@h@P+lower.getH()@h@lower-h+P@h@Q+Q@h@P)),
                  tilt_identity=float(sp.linalg.norm(P@number@P+lower.getH()@number@lower-number+Q)))
    assert max(errors.values()) < 1e-12
    return dict(N=n, total_vacancy_cutoff=vmax, quantum_head_positions=d,
                dimension=h.shape[0], phase=phi, **errors)


def pin_checks():
    reports = []
    # Full payload and both charged buses are retained in these small checks.
    for n, ks, lam in ((1, [1, 1], 0.), (1, [1, 1], 1e-5), (2, [2], 0.)):
        p = frozen.config(n, field=.003, lam=lam)
        J, s = .0005, .6
        words = material_words(ks)
        sources = material_sources(words, frozen.literal_bond(p)[0])
        N, d = len(words[0]), len(words)
        b, R = n+1, n+2
        M = p['Delta']+p['D']+abs(p['mu'])+p['field']+4+p['B']
        legal_pairs = np.array([20*i+j for i, j in product(range(20), repeat=2)
                                if frozen.OUT[i] == frozen.IN[j]])
        local = frozen.literal_bond(p)[0][np.ix_(legal_pairs, legal_pairs)]
        actual_norm = float(np.max(np.abs(la.eigvalsh(local))))
        assert actual_norm <= M
        assert J*R*M <= 1/32
        pinned = np.array([i for i, w in enumerate(words) if any(frozen.EXC[x] for x in w[:b])])
        head = {i: next(j for j in range(b) if frozen.EXC[words[i][j]]) for i in pinned}
        pin_sources = [h[pinned][:, pinned] for h in sources]
        # A free-bath norm enclosure excludes all but a few vacancy numbers.
        trial_v = N
        trial_upper = free_energy(N, trial_v)+s*trial_v+J*M*N
        vmax = int(np.ceil((trial_upper+2*N+J*M*N)/s))+1
        candidates = [v for v in range(vmax+1)
                      if free_energy(N, v)+s*v-J*M*N <= trial_upper+1e-12]
        full_data, pin_data = {}, {}
        for v in candidates:
            gaps, hfull = mobile_gap(sources, v, J)
            e, psi, residual = ground(hfull, vector=True)
            full_data[v] = (e, psi, residual)
            _, hpin = mobile_gap(pin_sources, v, J)
            ep, psip, rp = ground(hpin, vector=True)
            head_contact = np.array([float(g[head[i]] == 0) for g in gaps for i in pinned])
            a = float(np.dot(abs(psip)**2, head_contact))
            pin_data[v] = (ep, a, rp)
        vf = min(candidates, key=lambda v: full_data[v][0]+s*v)
        vp = min(candidates, key=lambda v: pin_data[v][0]+s*v)
        full = full_data[vf][0]+s*vf
        pin = pin_data[vp][0]+s*vp
        dex = (1-abs(p['g'])-p['q']**2)/1024
        assert pin-full >= J*dex-2e-10
        assert pin_data[vp][1] >= 1/1024
        canonical_cost = pin_data[vf][0]-full_data[vf][0]
        assert canonical_cost >= J*dex-2e-10
        full_psi = full_data[vf][1].reshape(-1, d)
        occupation = float(np.sum(abs(full_psi[:, pinned])**2))
        if lam:
            assert occupation <= (2*p['B']/dex)**2+1e-10
        reports.append(dict(n=n, N=N, material_dimension=d, pinned_dimension=len(pinned),
                            lam=lam, J=J, s=s, M_bound=M, actual_legal_bond_norm=actual_norm, JRM=J*R*M,
                            vacancy_candidates=candidates, excluded_v_above=vmax,
                            full_minimizer_V=vf, pin_minimizer_V=vp,
                            tilted_pin_cost=pin-full, canonical_pin_cost=canonical_cost,
                            theorem_cost=J*dex, conditioned_head_contact=pin_data[vp][1],
                            charged_occupation=occupation,
                            occupation_bound=(2*p['B']/dex)**2,
                            max_residual=max(x[2] for x in list(full_data.values())+list(pin_data.values()))))
    return reports


def balanced_checks():
    reports = []
    for n in (1, 2):
        p = frozen.config(n, field=.003, lam=0.)
        good = material_words([n, n], vacuum_zero=True)
        bad = material_words([n-1, n+1], vacuum_zero=True)
        h = frozen.literal_bond(p)[0]
        good_sources, bad_sources = material_sources(good, h), material_sources(bad, h)
        N = 2*(n+1)
        V = N
        a = alpha(N, V)
        delta = (la.eigvalsh(frozen.qmatrix(n-1, p['q']))[0]
                 +la.eigvalsh(frozen.qmatrix(n+1, p['q']))[0]-2*p['E'])
        for J in (.04, .02, .01, .005):
            eg = ground(mobile_gap(good_sources, V, J)[1])
            eb = ground(mobile_gap(bad_sources, V, J)[1])
            leading = J*a*delta
            reports.append(dict(n=n, N=N, V=V, J=J,
                                selected_clock_dimension=len(good), wrong_clock_dimension=len(bad),
                                alpha=a, exact_cost=eb-eg, first_order_cost=leading,
                                relative_remainder=eb-eg-leading,
                                remainder_over_J_squared=(eb-eg-leading)/J**2,
                                status='finite canonical zero-bit vacuum block; no uniform remainder theorem'))
    return reports


def main():
    start = time.time()
    out = dict(status='finite numerical corroboration; legal sectors only; no CFT claim',
               versions=dict(numpy=np.__version__, scipy=scipy.__version__))
    for name, function in [('representation', representation_checks), ('scalar', scalar_checks),
                           ('replacement', replacement_checks), ('controlled_channel', channel_checks),
                           ('charged_pin', pin_checks), ('balanced_pair', balanced_checks)]:
        out[name] = function()
        print(name, 'passed', flush=True)
        (HERE/'checks.json').write_text(json.dumps(out, indent=2)+'\n')
    out['runtime_seconds'] = time.time()-start
    (HERE/'checks.json').write_text(json.dumps(out, indent=2)+'\n')
    print('all checks passed in', round(out['runtime_seconds'], 2), 'seconds')


if __name__ == '__main__':
    main()
