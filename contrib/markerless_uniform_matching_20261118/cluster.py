"""Bounded literal-cell matching for the modified physical d=9 bond.

All matrices retain the complete allowed clock/head spaces. The only sector
restriction used for port matching is even charged-head number, which is exact.
"""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
from pathlib import Path
from itertools import product, combinations
from collections import defaultdict
import importlib.util
import numpy as np
import scipy.linalg as la
import scipy.sparse as sp
from scipy.sparse.linalg import splu

BASE = Path(__file__).resolve().parents[1]/'markerless_all_n_20261118'
spec = importlib.util.spec_from_file_location('all_n_literal_model', BASE/'model.py')
model = importlib.util.module_from_spec(spec)
spec.loader.exec_module(model)


def words(n, cells):
    local = model.cell_words(n)
    out = []
    for parts in product(local, repeat=cells):
        heads = tuple(a[1] for a in parts)
        if sum(h != 'W' for h in heads) % 2:
            continue
        word = (model.IDX['R0'],)+sum((a[0] for a in parts), ())+ (model.IDX['L0'],)
        out.append((word, heads, tuple(a[2] for a in parts), sum((a[3] for a in parts), ())))
    return out


def literal(n, cells, params, states=None, active_ports=None):
    states = words(n, cells) if states is None else states
    index = {s[0]: i for i, s in enumerate(states)}
    W, onsite, _ = model.bond(params)
    by_in = defaultdict(list)
    for row, col in zip(*np.nonzero(W)):
        by_in[int(col)].append((int(row), W[row, col]))
    rows, cols, vals = [], [], []
    leaks = 0
    source = (model.IDX['W'], model.IDX['W'])
    targets = {(model.IDX['XA'], model.IDX['B']), (model.IDX['XB'], model.IDX['A'])}
    for col, state in enumerate(states):
        word = state[0]
        rows.append(col); cols.append(col); vals.append(sum(onsite[a] for a in word))
        for site in range(len(word)-1):
            before = word[site:site+2]
            for rowpair, value in by_in[9*before[0]+before[1]]:
                after = (rowpair//9, rowpair % 9)
                isport = (before == source and after in targets) or (after == source and before in targets)
                if isport and active_ports is not None:
                    edge = site//(n+1)-1
                    if edge not in active_ports:
                        continue
                target = word[:site]+after+word[site+2:]
                if target not in index:
                    leaks += 1
                else:
                    rows.append(index[target]); cols.append(col); vals.append(value)
    H = sp.csc_matrix((vals, (rows, cols)), shape=(len(states),)*2)
    assert leaks == 0
    return states, H


class Cluster:
    def __init__(self, n, cells, params):
        self.n, self.cells, self.params = n, cells, dict(params)
        self.states, H = literal(n, cells, dict(params, lam=0., hL=0.))
        ev, v, gamma = model.clock(n, params['q'])
        self.clock_energy, self.v, self.gamma = ev, v, gamma
        self.bits = list(product((0, 1), repeat=n*cells))
        self.bidx = {b: i for i, b in enumerate(self.bits)}
        P = sp.lil_matrix((len(self.states), len(self.bits)), dtype=complex)
        for row, state in enumerate(self.states):
            if all(h == 'W' for h in state[1]):
                P[row, self.bidx[state[3]]] = np.prod([v[j] for j in state[2]])
        self.P = P.tocsc()
        self.H0 = H-cells*ev*sp.eye(len(self.states), format='csc')
        assert la.norm((self.H0@self.P).toarray()) < 1e-12
        self.PP = (self.P@self.P.conj().T).tocsc()
        self.lu = splu(self.H0+self.PP)
        self.ports = []
        for e in range(cells-1):
            _, HV = literal(n, cells, dict(params, lam=1., hL=0.), self.states, {e})
            self.ports.append(HV-H)
        self.V = sum(self.ports, sp.csc_matrix(H.shape))
        field = np.array([sum(a == model.IDX['L1'] for a in s[0]) for s in self.states])
        self.F = sp.diags(field, format='csc', dtype=float)

    def R(self, A):
        if sp.issparse(A):
            A = A.toarray()
        return self.lu.solve(A)-self.P@(self.P.conj().T@A)

    def coefficients(self, active_ports=None, field_scale=0.):
        V = self.V if active_ports is None else sum((self.ports[e] for e in active_ports), sp.csc_matrix(self.H0.shape))
        P = self.P
        VP = (V@P).toarray()
        X = self.R(VP)
        B, A = VP.conj().T@X, X.conj().T@X
        Y = V@X
        K4 = -Y.conj().T@self.R(Y)+(A@B+B@A)/2
        FP = (self.F@P).toarray()
        Z = self.R(FP)
        FPP = (P.conj().T@FP)
        mixed = X.conj().T@(self.F@X)+Z.conj().T@(V@X)+X.conj().T@(V@Z)-(A@FPP+FPP@A)/2
        field2 = -FP.conj().T@Z
        return dict(K2=-B+field_scale*FPP, port4=K4, mixed=mixed,
                    field2=field2, K4=K4+field_scale*mixed+field_scale**2*field2)

    def conversion_series(self, degree=4):
        assert self.params.get('g', 0.) == 0.
        _, HG = literal(self.n, self.cells, dict(self.params, g=1., lam=0., hL=0.), self.states)
        C = HG-(self.H0+self.cells*self.clock_energy*sp.eye(len(self.states), format='csc'))
        VP = (self.V@self.P).toarray()
        X = self.R(VP)
        out = [-VP.conj().T@X]
        for order in range(1, degree+1):
            X = -self.R(C@X)
            out.append(-VP.conj().T@X)
        return out

    def bloch_fourth(self, field_scale=0.):
        """Independent wave-operator recursion and metric correction."""
        P = self.P.toarray()
        F = field_scale*self.F
        omega = [P]
        energy = [np.zeros((P.shape[1], P.shape[1]), complex)]
        for r in range(1, 5):
            source = self.V@omega[r-1]
            if r >= 2:
                source += F@omega[r-2]
            energy.append(P.conj().T@source)
            for j in range(1, r):
                source -= omega[j]@energy[r-j]
            omega.append(-self.R(source))
        metric2 = omega[1].conj().T@omega[1]
        canonical4 = energy[4]+(metric2@energy[2]-energy[2]@metric2)/2
        return energy[2], canonical4


def annihilate(bit, j):
    if not bit[j]:
        return None
    sign = (-1)**sum(bit[:j]); result = list(bit); result[j] = 0
    return tuple(result), sign


def create(bit, i):
    if bit[i]:
        return None
    sign = (-1)**sum(bit[:i]); result = list(bit); result[i] = 1
    return tuple(result), sign


def fock_extract(H, modes):
    """Normal order relative to the empty logical vacuum; no fitted constants."""
    bits = list(product((0, 1), repeat=modes)); idx = {b: i for i, b in enumerate(bits)}
    zero = idx[(0,)*modes]
    one = [idx[tuple(int(j == i) for j in range(modes))] for i in range(modes)]
    scalar = H[zero, zero]
    h = H[np.ix_(one, one)]-scalar*np.eye(modes)
    quadratic = scalar*np.eye(len(bits), dtype=complex)
    for col, bit in enumerate(bits):
        for j in range(modes):
            mid = annihilate(bit, j)
            if mid is None:
                continue
            for i in range(modes):
                end = create(mid[0], i)
                if end is not None:
                    quadratic[idx[end[0]], col] += h[i, j]*mid[1]*end[1]
    pairs = list(combinations(range(modes), 2))
    two = [idx[tuple(int(r in pair) for r in range(modes))] for pair in pairs]
    U = (H-quadratic)[np.ix_(two, two)]
    quartic = np.zeros_like(H)
    for col, bit in enumerate(bits):
        for b, (k, l) in enumerate(pairs):
            first = annihilate(bit, k)
            second = None if first is None else annihilate(first[0], l)
            if second is None:
                continue
            for a, (i, j) in enumerate(pairs):
                third = create(second[0], j)
                fourth = None if third is None else create(third[0], i)
                if fourth is not None:
                    quartic[idx[fourth[0]], col] += U[a, b]*first[1]*second[1]*third[1]*fourth[1]
    return dict(scalar=scalar, quadratic=h, quartic=U, quartic_operator=quartic,
                pair_indices=pairs, higher=H-quadratic-quartic)


def clock_field(n, q):
    ev, v, gamma = model.clock(n, q)
    H = -q*(np.eye(n+1, k=1)+np.eye(n+1, k=-1)); H[0, 0] = -1
    R = la.inv(H-ev*np.eye(n+1)+np.outer(v, v))-np.outer(v, v)
    sources = np.column_stack([v*(np.arange(n+1) > i) for i in range(n)])
    C = sources.T@R@sources
    cumulative = np.cumsum(v*v)
    conductance = q*v[:-1]*v[1:]
    currents = np.array([[-cumulative[min(i, j)]*(1-cumulative[max(i, j)])
                          for j in range(n)] for i in range(n)])
    current_C = (currents/conductance)@currents.T
    visibility_derivative = -2*(v[0]*v[-1])**2*((R@sources)[0, :]/v[0]+(R@sources)[-1, :]/v[-1])
    return dict(ev=ev, v=v, gamma=gamma, C=C, current_C=current_C,
                visibility_derivative=visibility_derivative)
