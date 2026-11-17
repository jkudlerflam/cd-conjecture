"""Independent vacuum-clock matrices, coherent-channel identities, and audit.

Numerics corroborate the proof; no thermodynamic threshold or c is fitted.
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


mobile = load('previous_mobile', HERE.parent/'small_d_mobile_matching_20261116'/'checks.py')
parent = mobile.parent


def clock(n, q=.25):
    edges = []
    for j in range(n):
        vector = np.zeros(n+1)
        vector[j], vector[j+1] = q, -1
        edges.append(np.outer(vector, vector))
    Q = sum(edges)
    u = q**np.arange(n+1, dtype=float)
    u /= la.norm(u)
    terminal = np.zeros_like(Q); terminal[-1, -1] = 1
    d = 1+q*q-2*q*np.cos(np.pi/(n+1))
    c = 1+q*q+2*q*np.cos(np.pi/(n+1))
    energy = float(la.eigvalsh(Q+q*q*terminal)[0])
    return dict(edges=edges, Q=Q, u=u, P=np.outer(u, u), T=terminal,
                d=d, c=c, E=energy, mu=-(n+2)*energy,
                tau=q*q*u[-1]**2, epsilon=q*q*abs(u[-1])*np.sqrt(1-u[-1]**2))


def gap_basis(n, vacancies):
    out = []
    for bars in combinations(range(vacancies+n-1), n-1):
        endpoints = (-1,)+bars+(vacancies+n-1,)
        out.append(tuple(endpoints[j+1]-endpoints[j]-1 for j in range(n)))
    return out


def bath(n, vacancies, phase=0.):
    gaps = gap_basis(n, vacancies)
    index = {g: j for j, g in enumerate(gaps)}
    rows, cols, values = [], [], []
    for col, g in enumerate(gaps):
        for i in range(n):
            if g[i] == 0:
                continue
            changed = list(g); changed[i] -= 1; changed[(i-1) % n] += 1
            row = index[tuple(changed)]
            rows += [row, col]; cols += [col, row]
            values += [-np.exp(-1j*phase), -np.exp(1j*phase)]
    h = sp.coo_matrix((values, (rows, cols)), shape=(len(gaps), len(gaps))).tocsr()
    contacts = [sp.diags([float(g[i] == 0) for g in gaps], format='csr') for i in range(n)]
    return gaps, h, contacts


def embed(matrix, position, dims):
    result = sp.csr_matrix([[1.]])
    for j, d in enumerate(dims):
        result = sp.kron(result, sp.csr_matrix(matrix) if j == position else sp.eye(d), format='csr')
    return result


def sources(ks, q, mu, terminal=True):
    dims = [k+1 for k in ks]
    size = int(np.prod(dims))
    N = sum(dims)
    h = [sp.csr_matrix((size, size)) for _ in range(N)]
    offset = 0
    for cell, k in enumerate(ks):
        info = clock(k, q) if k else None
        edges = info['edges'] if k else []
        for j, edge in enumerate(edges):
            h[offset+j] += embed(edge, cell, dims)
        marker = offset+k
        T = np.zeros((k+1, k+1)); T[-1, -1] = 1
        for bond in ((marker-1) % N, marker):
            h[bond] += mu/2*sp.eye(size, format='csr')
            if terminal:
                h[bond] += q*q/2*embed(T, cell, dims)
        offset += k+1
    return h


def assemble(hsources, vacancies, J, phase=0.):
    gaps, hcharge, contacts = bath(len(hsources), vacancies, phase)
    d = hsources[0].shape[0]
    base = sp.kron(hcharge, sp.eye(d), format='csr')
    material = sum((sp.kron(a, h, format='csr') for a, h in zip(contacts, hsources)),
                   sp.csr_matrix(base.shape))
    return gaps, base+J*material


def minimum(h):
    if h.shape[0] <= 100:
        energies, vectors = la.eigh(h.toarray(), subset_by_index=[0, 0])
    else:
        start = np.random.default_rng(1747).normal(size=h.shape[0])
        energies, vectors = eigsh(h, k=1, which='SA', tol=1e-12, maxiter=30000, v0=start)
    e, psi = float(energies[0]), vectors[:, 0]
    residual = float(la.norm(h@psi-e*psi))
    assert residual < 2e-8
    return e, psi, residual


def source_checks():
    reports = []
    for n in (1, 2):
        info = clock(n)
        ks, dims = [n, n], [n+1, n+1]
        analytic = sources(ks, .25, info['mu'])
        p = mobile.frozen.config(n, field=0., lam=0.)
        literal_words = mobile.material_words(ks, vacuum_zero=True)
        literal = mobile.material_sources(literal_words, mobile.frozen.literal_bond(p)[0])
        index = {w: i for i, w in enumerate(literal_words)}
        mapping = []
        for heads in product(range(n+1), repeat=2):
            word = []
            for head in heads:
                letters = ([('L', 0)]*n+[('m',)] if head == n else
                           [('L', 0)]*head+[('v',)]+[('R', 0)]*(n-head-1)+[('e', 0)])
                word += [mobile.frozen.INDEX[x] for x in letters]
            mapping.append(index[tuple(word)])
        errors = [float(sp.linalg.norm(a-b[mapping][:, mapping])) for a, b in zip(analytic, literal)]
        assert max(errors) < 1e-13
        gaps, base = assemble(sources(ks, .25, info['mu'], terminal=False), 2*(n+1), .0125)
        P = sp.kron(sp.eye(len(gaps)), embed(info['P'], 0, dims), format='csr')
        ff_commutator = float(sp.linalg.norm(base@P-P@base))
        full = assemble(analytic, 2*(n+1), .0125)[1]
        contacts = bath(2*(n+1), 2*(n+1))[2]
        C = (contacts[n-1]+contacts[n])/2
        local_comm = info['T']@info['P']-info['P']@info['T']
        predicted = .0125*.25**2*sp.kron(C, embed(local_comm, 0, dims), format='csr')
        leakage_error = float(sp.linalg.norm(full@P-P@full-predicted))
        assert ff_commutator < 1e-12 and leakage_error < 1e-12
        reports.append(dict(n=n, source_difference=max(errors),
                            ff_commutator=ff_commutator, terminal_commutator_identity_error=leakage_error,
                            terminal_leakage_per_J=info['epsilon'], d=info['d'], c=info['c']))
    return reports


def channel_checks():
    # Two nonorthogonal gradients: this is not a classical head measurement.
    n, N, vmax, phase = 2, 3, 3, .31
    info = clock(n)
    gaps = [g for v in range(vmax+1) for g in gap_basis(N, v)]
    index = {g: i for i, g in enumerate(gaps)}
    hg = sp.block_diag([bath(N, v, phase)[1] for v in range(vmax+1)], format='csr')
    H = sp.kron(hg, sp.eye(n+1), format='csr')
    G = sp.diags([sum(g) for g in gaps for _ in range(n+1)], dtype=float, format='csr')
    kraus, kinetic_removed = [], sp.csr_matrix(H.shape)
    contact_energy = sp.csr_matrix(H.shape)
    for j, edge in enumerate(info['edges']):
        A = sp.diags([float(g[j] == 0) for g in gaps], format='csr')
        notA = sp.eye(len(gaps), format='csr')-A
        rows, cols = [], []
        for col, g in enumerate(gaps):
            if g[j]:
                lower = list(g); lower[j] -= 1
                rows.append(index[tuple(lower)]); cols.append(col)
        lowering = sp.coo_matrix((np.ones(len(rows)), (rows, cols)), shape=hg.shape).tocsr()
        sqrt_edge = edge/np.sqrt((1+.25**2)*info['c'])
        kraus += [sp.kron(A, sqrt_edge, format='csr'), sp.kron(lowering, sqrt_edge, format='csr')]
        kinetic_removed += sp.kron(A@hg@notA+notA@hg@A, edge/info['c'], format='csr')
        contact_energy += sp.kron(A, edge, format='csr')
    vals, vec = la.eigh(np.eye(n+1)-info['Q']/info['c'])
    remainder = (vec*np.sqrt(np.maximum(vals, 0)))@vec.T
    kraus.append(sp.kron(sp.eye(len(gaps)), remainder, format='csr'))
    completeness = sum((k.getH()@k for k in kraus), sp.csr_matrix(H.shape))
    transformed_H = sum((k.getH()@H@k for k in kraus), sp.csr_matrix(H.shape))
    transformed_G = sum((k.getH()@G@k for k in kraus), sp.csr_matrix(H.shape))
    totalQ = sp.kron(sp.eye(len(gaps)), info['Q'], format='csr')
    errors = dict(completeness=float(sp.linalg.norm(completeness-sp.eye(H.shape[0]))),
                  kinetic=float(sp.linalg.norm(transformed_H-H+kinetic_removed)),
                  vacancy=float(sp.linalg.norm(transformed_G-G+(totalQ-contact_energy)/info['c'])))
    assert max(errors.values()) < 2e-13
    return dict(n=n, N=N, V_cutoff=vmax, phase=phase,
                edge_noncommutativity=float(la.norm(info['edges'][0]@info['edges'][1]-info['edges'][1]@info['edges'][0])),
                **errors)


def finite_clocks():
    reports = []
    for n in (1, 2):
        q, s = .25, .6
        info = clock(n, q)
        ks, dims = [n, n], [n+1, n+1]
        N, V = 2*(n+1), 2*(n+1)
        hs = sources(ks, q, info['mu'])
        Qbasis = la.null_space(info['u'][None, :])
        material_embedding = sp.kron(Qbasis, sp.eye(n+1), format='csr')
        gaps = gap_basis(N, V)
        pin_embedding = sp.kron(sp.eye(len(gaps)), material_embedding, format='csr')
        Pperp = sp.kron(sp.eye(len(gaps)), embed(np.eye(n+1)-info['P'], 0, dims), format='csr')
        for J in (1/160, 1/80, 1/32, .1, .3, 1.):
            h = assemble(hs, V, J)[1]
            e, psi, residual = minimum(h)
            ep, _, rp = minimum(pin_embedding.getH()@h@pin_embedding)
            occupation = float(np.real(np.vdot(psi, Pperp@psi)))
            r = info['c']+q*q+abs(info['mu'])/2
            allowance = .5*info['d']/(8*info['c']*r)
            lower = J*(info['d']/257-info['tau'])
            occ_bound = (65*info['tau']+8*info['c']*J*r)/info['d']
            # Explicit free-sea lower bounds exclude all unseen V when possible.
            upper = ep+s*V
            vmax = int(np.ceil((upper+2*N+J*abs(info['mu'])*N)/s))+1
            candidates = [v for v in range(vmax+1)
                          if mobile.free_energy(N, v)+s*v-J*abs(info['mu'])*N <= upper+1e-11]
            certified_number = candidates == [V]
            if J <= allowance and certified_number:
                assert ep-e >= lower-2e-10
                assert occupation <= occ_bound+2e-10
            assert ep-e > 0
            reports.append(dict(n=n, N=N, V=V, J_over_t=J, q=q,
                                dimension=h.shape[0], pin_cost=ep-e, pin_cost_over_J=(ep-e)/J,
                                ff_excitation_weight=occupation, theorem_pin_lower=lower,
                                theorem_weight_upper=occ_bound,
                                sufficient_J_over_t=allowance,
                                tilted_minimum_certified=certified_number,
                                possible_tilted_V=candidates,
                                residual=max(residual, rp),
                                scope='exact zero-bit vacuum sector; fixed finite N,V'))
    return reports


def simplified_variants():
    reports = []
    for n, D in product((1, 2), (0., 3.)):
        p = parent.parameters(n, q=.25, w=.4, v=.35, eta=.25, D=D, g=0., field=0.)
        p['Delta'] = 4.
        _, _, H, C, G, _ = parent.cell_matrices(p)
        spectrum = la.eigvalsh((H+.13*C).toarray())
        expected = p['E']+p['selector']
        assert np.max(np.abs(spectrum[:2**n]-expected)) < 2e-12
        gap = float(spectrum[2**n]-expected)
        assert gap >= .5
        local = parent.material_bond(p, lam=.017)
        failures = 0
        for row, col in zip(*np.nonzero(local)):
            if row == col:
                continue
            a, b = parent.MATERIAL[col//20], parent.MATERIAL[col % 20]
            c, d = parent.MATERIAL[row//20], parent.MATERIAL[row % 20]
            failures += ((parent.marker(a), parent.marker(b)) != (parent.marker(c), parent.marker(d)))
            failures += ((parent.arrows(a)[0], parent.arrows(b)[1]) != (parent.arrows(c)[0], parent.arrows(d)[1]))
            failures += not (parent.allowed(a, b) and parent.allowed(c, d))
        assert failures == 0
        h22, _, dimensions = parent.pair_coefficient(p)
        _, gaussian_error = parent.gaussian_decomposition(h22, 2*n)
        assert gaussian_error < 1e-10
        reports.append(dict(n=n, Delta=4., D=D, ground_energy=expected,
                            auxiliary_gap=gap, grammar_failures=failures,
                            gaussian_pair_relative_error=gaussian_error, **dimensions))
    return reports


def full_material_checks():
    """Both buses, every payload, readout, and nonzero ports are retained."""
    reports = []
    n, N, V, J, f, lam = 1, 4, 4, 1/80, .001, 1e-5
    info = clock(n)
    cell_words, descriptions = parent.cell_words(n)
    dc = len(cell_words)
    P0, Pvac = np.zeros((dc, dc)), np.zeros((dc, dc))
    Qcolumns = []
    for bits in product((0, 1), repeat=n):
        indices = [next(i for i, (data, token, pos) in enumerate(descriptions)
                        if data == bits and pos == j and token in (('v',), ('m',))) for j in range(n+1)]
        P0[np.ix_(indices, indices)] = info['P']
        Pvac[indices, indices] = 1
        column = np.zeros(dc)
        column[indices] = la.null_space(info['u'][None, :])[:, 0]
        Qcolumns.append(column)
    ordered = [tuple(parent.INDEX[z] for z in x+y) for x in cell_words for y in cell_words]
    words = sorted(ordered)
    order_index = {w: i for i, w in enumerate(ordered)}
    permutation = [order_index[w] for w in words]
    pin_isometry = sp.kron(np.array(Qcolumns).T, sp.eye(dc), format='csr')[permutation]
    projector = sp.kron(Pvac-P0, sp.eye(dc), format='csr')[permutation][:, permutation]
    for D in (0., 25.):
        p = parent.parameters(n, q=.25, w=.4, v=.35, eta=.25, D=D, g=.13, K=1., field=f)
        p['Delta'] = 4.
        hs = mobile.material_sources(words, parent.material_bond(p, lam=lam))
        gaps, h = mobile.mobile_gap(hs, V, J)
        U = sp.kron(sp.eye(len(gaps)), pin_isometry, format='csr')
        e, psi, residual = minimum(h)
        ep, _, rp = minimum(U.getH()@h@U)
        P = sp.kron(sp.eye(len(gaps)), projector, format='csr')
        weight = float(np.real(np.vdot(psi, P@psi)))
        B = np.sqrt(2+2*p['eta'])*lam
        r = info['c']+.25**2+abs(info['mu'])/2+n*f
        bound = (65*(info['tau']+2*B)+8*info['c']*J*(r+2*B))/info['d']
        pin_lower = J*(info['d']/257-info['tau'])
        assert weight <= bound and ep-e >= pin_lower
        s = .6
        upper = ep+s*V
        negative = abs(info['mu'])+B
        vmax = int(np.ceil((upper+2*N+J*negative*N)/s))+1
        candidates = [v for v in range(vmax+1)
                      if mobile.free_energy(N, v)+s*v-J*negative*N <= upper+1e-11]
        assert candidates == [V]
        reports.append(dict(n=n, N=N, V=V, J_over_t=J, Delta=4., D=D,
                            f=f, lam=lam, dimension=h.shape[0],
                            clock_pin_cost=ep-e, theorem_pin_lower=pin_lower,
                            clock_excitation_weight=weight, theorem_weight_upper=bound,
                            tilted_minimum_certified=True, possible_tilted_V=candidates,
                            residual=max(residual, rp)))
    return reports


def constants():
    rows = []
    for n in range(1, 7):
        info = clock(n)
        r = info['c']+.25**2+abs(info['mu'])/2
        rows.append(dict(n=n, q=.25, tau=info['tau'], terminal_leakage=info['epsilon'],
                         pin_per_J=info['d']/257-info['tau'],
                         J_over_t_allowance=.5*info['d']/(8*info['c']*r),
                         frozen_B_allowance=(1-.25**2)**2*.25**(2*n+4)/10,
                         selector_mu=info['mu']))
    assert min(row['pin_per_J'] for row in rows) >= 2/4369-1e-15
    return rows


def main():
    start = time.time()
    out = dict(status='numerical corroboration; no uniform threshold or CFT fit',
               versions=dict(numpy=np.__version__, scipy=scipy.__version__))
    for name, function in [('literal_sources', source_checks), ('coherent_channel', channel_checks),
                           ('finite_clocks', finite_clocks), ('full_material', full_material_checks),
                           ('simplified_variants', simplified_variants),
                           ('parameter_constants', constants)]:
        out[name] = function()
        print(name, 'passed', flush=True)
        (HERE/'checks.json').write_text(json.dumps(out, indent=2)+'\n')
    out['runtime_seconds'] = time.time()-start
    (HERE/'checks.json').write_text(json.dumps(out, indent=2)+'\n')
    print('all checks passed in', round(out['runtime_seconds'], 2), 'seconds')


if __name__ == '__main__':
    main()
