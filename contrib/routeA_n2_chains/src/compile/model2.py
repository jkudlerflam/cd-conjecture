"""Lane version of the dimer-port model: cells with an express lane E and
a data lane U attached only at a central switch site M.

layout 'uniform':  P D^n P            (lane E only; data on E)  == model.py
layout 'switch' :  P L^m M R^m P      (E on every site; U + data on R sites;
                                       E(M) -- U(R_1) is the only E/U link)
Every lane slot is a dimer (beta, alpha) with states none/T/H/D as in
model.py; conversions act between a U-dimer (or the E-dimer in 'uniform')
and the data orbital on the same site.
"""
import itertools
import numpy as np
import scipy.sparse as sp
from model import Chain, NONE, CARRIERS, DEFAULT


class LaneChain(Chain):
    def __init__(self, layout, m, M, K=2, params=None, max_data=None, one_per_cell=False):
        self.layout, self.m, self.M, self.K = layout, m, M, K
        self.max_data, self.one_per_cell = max_data, one_per_cell
        self.p = dict(DEFAULT)
        if params:
            self.p.update(params)
        if layout == 'uniform':
            kinds = ['P'] + ['D'] * m + ['P']
        elif layout == 'switch':
            kinds = ['P'] + ['L'] * m + ['M'] + ['R'] * m + ['P']
        else:
            raise ValueError(layout)
        self.kinds = kinds
        self.n = m
        self.dimers = []     # (cell, site, lane, beta, alpha)
        self.datas = []      # (cell, site, lane_of_conversion, mode)
        self.site_dimer = {}
        mo = 0
        for j in range(M):
            for s, k in enumerate(kinds):
                lanes = ['E'] + (['U'] if k == 'R' else [])
                for ln in lanes:
                    self.site_dimer[(j, s, ln)] = len(self.dimers)
                    self.dimers.append((j, s, ln, mo, mo + 1))
                    mo += 2
                if k == 'D':
                    self.datas.append((j, s, 'E', mo)); mo += 1
                if k == 'R':
                    self.datas.append((j, s, 'U', mo)); mo += 1
        self.nmodes = mo
        self.data_modes = [d[3] for d in self.datas]
        # intra-cell hop edges between dimers on adjacent sites
        self.edges = []
        for j in range(M):
            for s in range(len(kinds) - 1):
                self.edges.append((self.site_dimer[(j, s, 'E')], self.site_dimer[(j, s + 1, 'E')]))
                if kinds[s] == 'R' and kinds[s + 1] == 'R':
                    self.edges.append((self.site_dimer[(j, s, 'U')], self.site_dimer[(j, s + 1, 'U')]))
                if kinds[s] == 'M' and kinds[s + 1] == 'R':
                    self.edges.append((self.site_dimer[(j, s, 'E')], self.site_dimer[(j, s + 1, 'U')]))
        last = len(kinds) - 1
        self.junctions = [(self.site_dimer[(j, last, 'E')], self.site_dimer[(j + 1, 0, 'E')])
                          for j in range(M - 1)]
        self._basis2()

    def _basis2(self):
        nd = len(self.data_modes)
        states = []
        Dn = len(self.dimers)
        for k in range(self.K + 1):
            for subset in itertools.combinations(range(Dn), k):
                if self.one_per_cell and k > 0:
                    cells = [self.dimers[i][0] for i in subset]
                    if k != self.M or len(set(cells)) != self.M:
                        continue
                for types in itertools.product(CARRIERS, repeat=k):
                    carr = dict(zip(subset, types))
                    base = 0
                    for i, d in enumerate(self.dimers):
                        b, a = carr.get(i, NONE)
                        base |= (b << d[3]) | (a << d[4])
                    for bits in range(1 << nd):
                        if self.max_data is not None and bin(bits).count('1') > self.max_data:
                            continue
                        st = base
                        for q, dm in enumerate(self.data_modes):
                            if (bits >> q) & 1:
                                st |= 1 << dm
                        states.append(st)
        self.states = np.array(sorted(states), dtype=object)
        self.index = {st: i for i, st in enumerate(self.states)}
        self.dim = len(self.states)
        self.number = np.array([bin(st).count('1') for st in self.states])
        code = []
        for i, st in enumerate(self.states):
            if all(((st >> d[3]) & 1, (st >> d[4]) & 1) == NONE for d in self.dimers):
                code.append(i)
        self.code = np.array(code)

    def lane_parts(self):
        P = {}
        D = self.dimers
        P['ET'] = self.build([(1.0, [], [(d[3], 0), (d[4], 1)]) for d in D])
        P['EH'] = self.build([(1.0, [], [(d[3], 0), (d[4], 0)]) for d in D])
        P['ED'] = self.build([(1.0, [], [(d[3], 1), (d[4], 1)]) for d in D])
        P['MU'] = self.build([(1.0, [], [(d[3], 1)]) for d in self.datas])
        AT, AH, AD, CA, CB = [], [], [], [], []
        for (i1, i2) in self.edges:
            b, a = D[i1][3], D[i1][4]
            b2, a2 = D[i2][3], D[i2][4]
            AT.append((1.0, [(b, 1), (a, 0), (a2, 1), (b2, 0)], []))
            AH.append((1.0, [(b, 1), (b2, 0)], [(a, 0), (a2, 0)]))
            AD.append((1.0, [(a2, 1), (a, 0)], [(b, 1), (b2, 1)]))
        for (j, s, ln, c) in self.datas:
            d = D[self.site_dimer[(j, s, ln)]]
            b, a = d[3], d[4]
            CA.append((1.0, [(c, 1), (a, 0)], [(b, 0)]))
            CB.append((1.0, [(b, 1), (c, 0)], [(a, 1)]))
        for k, L in (('AT', AT), ('AH', AH), ('AD', AD), ('CA', CA), ('CB', CB)):
            P[k] = self.build(L)
        hop, dens = [], []
        for (iR, iL) in self.junctions:
            R, L = D[iR], D[iL]
            for mr in (R[3], R[4]):
                for ml in (L[3], L[4]):
                    hop.append((0.5, [(mr, 1), (ml, 0)], []))
            NR = self._Nlist((None, None, R[3], R[4]))
            NL = self._Nlist((None, None, L[3], L[4]))
            for cr, opr in NR:
                for cl, opl in NL:
                    dens.append((cr * cl, opr + opl, []))
        P['BH'] = self.build(hop)
        P['BN'] = self.build(dens)
        return P


def segment_parts(ch):
    """Like lane_parts, but the carrier hops are split by segment so that the two
    express arms can be given different hopping and phases:
    'L' = E-lane edges from P_L to M, 'R' = E-lane edges from M to P_R,
    'U' = data-lane edges, 'S' = the switch edge E(M)-U(R_1)."""
    P = ch.lane_parts()
    D = ch.dimers
    seg = {k: {'AT': [], 'AH': [], 'AD': []} for k in 'LRUS'}
    mid = ch.kinds.index('M')
    for (i1, i2) in ch.edges:
        s1, l1 = D[i1][1], D[i1][2]
        s2, l2 = D[i2][1], D[i2][2]
        if l1 == 'E' and l2 == 'E':
            k = 'L' if s2 <= mid else 'R'
        elif l1 == 'U' and l2 == 'U':
            k = 'U'
        else:
            k = 'S'
        b, a = D[i1][3], D[i1][4]
        b2, a2 = D[i2][3], D[i2][4]
        seg[k]['AT'].append((1.0, [(b, 1), (a, 0), (a2, 1), (b2, 0)], []))
        seg[k]['AH'].append((1.0, [(b, 1), (b2, 0)], [(a, 0), (a2, 0)]))
        seg[k]['AD'].append((1.0, [(a2, 1), (a, 0)], [(b, 1), (b2, 1)]))
    for k in 'LRUS':
        for t in ('AT', 'AH', 'AD'):
            P[t + '_' + k] = ch.build(seg[k][t])
    return P


def assemble_seg(P, p):
    """p may contain per-segment multipliers xL, xR (E arms) and a phase
    phiT_L on carrier-T hops in the left arm only."""
    off = 0
    for k in 'LRUS':
        x = p.get('x' + k, 1.0)
        phT = p['phiT'] + (p.get('phiT_L', 0.0) if k == 'L' else 0.0)
        off = off + x * (-p['wT'] * np.exp(1j * phT) * P['AT_' + k]
                         - p['wH'] * P['AH_' + k] - p['wD'] * P['AD_' + k])
    off = off - p['tauA'] * P['CA'] - p['tauB'] * P['CB']
    H0 = (off + off.getH() + p['eT'] * P['ET'] + p['eH'] * P['EH']
          + p['eD'] * P['ED'] + p['muc'] * P['MU'])
    V = P['BH'] + P['BH'].getH() + p['kappa'] * P['BN']
    return H0.tocsr(), V.tocsr()
