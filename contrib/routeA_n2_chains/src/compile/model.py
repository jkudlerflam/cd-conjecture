"""Fermionic dimer-port cell with mobile carriers (route A compile attempt).

Every site carries a "dimer" of two fermion orbitals (beta, alpha).
Interior sites also carry one data orbital c (the logical mode).
Dimer states: none=(1,0) [code], T=(0,1), H=(0,0), D=(1,1).
Ports are dimer sites without data.  Cells: P D^n P, cells concatenated.

All operators are built from physical fermion operators with a
cell-major, site-major Jordan-Wigner order (beta, alpha, c per site),
so total fermion number and parity are exact by construction.

The Hilbert space is truncated to configurations with at most K
non-'none' dimers (carriers).  H0 conserves the carrier number and the
bond V changes it by 0 or +-2, so the truncation is exact for every
perturbative coefficient through order lambda^(K+1) (K even).
"""
import itertools
import numpy as np
import scipy.sparse as sp

NONE, TT, HH, DD = (1, 0), (0, 1), (0, 0), (1, 1)
CARRIERS = (TT, HH, DD)

DEFAULT = dict(eT=1.0, eH=1.0, eD=1.0, wT=0.3, wH=0.3, wD=0.3,
               phiT=0.0, tauA=0.3, tauB=0.3, kappa=1.0, muc=0.0)


class Chain:
    def __init__(self, n, M, K=2, params=None, max_data=None, one_per_cell=False):
        self.n, self.M, self.K = n, M, K
        self.max_data = max_data
        self.one_per_cell = one_per_cell
        self.p = dict(DEFAULT)
        if params:
            self.p.update(params)
        # sites
        self.sites = []          # (cell, x, beta_mode, alpha_mode, data_mode or None)
        m = 0
        for j in range(M):
            for x in range(n + 2):
                interior = 1 <= x <= n
                s = (j, x, m, m + 1, m + 2 if interior else None)
                m += 3 if interior else 2
                self.sites.append(s)
        self.nmodes = m
        self.data_modes = [s[4] for s in self.sites if s[4] is not None]
        self._basis()

    # ---------------------------------------------------------------- basis
    def _basis(self):
        S = len(self.sites)
        states = []
        nd = len(self.data_modes)
        for k in range(self.K + 1):
            for subset in itertools.combinations(range(S), k):
                if self.one_per_cell and k > 0:
                    cells = [self.sites[i][0] for i in subset]
                    if k != self.M or len(set(cells)) != self.M:
                        continue
                for types in itertools.product(CARRIERS, repeat=k):
                    base = 0
                    carr = dict(zip(subset, types))
                    for i, s in enumerate(self.sites):
                        b, a = carr.get(i, NONE)
                        base |= (b << s[2]) | (a << s[3])
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
            ok = True
            for s in self.sites:
                if ((st >> s[2]) & 1, (st >> s[3]) & 1) != NONE:
                    ok = False
                    break
            if ok:
                code.append(i)
        self.code = np.array(code)

    # ------------------------------------------------------------ operators
    @staticmethod
    def apply(ops, conds, st):
        for (mo, val) in conds:
            if ((st >> mo) & 1) != val:
                return None
        sign = 1
        for (mo, dag) in reversed(ops):
            occ = (st >> mo) & 1
            if occ == dag:
                return None
            if bin(st & ((1 << mo) - 1)).count('1') & 1:
                sign = -sign
            st ^= (1 << mo)
        return sign, st

    def build(self, terms):
        rows, cols, vals = [], [], []
        for i, st in enumerate(self.states):
            for coef, ops, conds in terms:
                r = self.apply(ops, conds, st)
                if r is None:
                    continue
                sgn, new = r
                jdx = self.index.get(new)
                if jdx is None:
                    continue
                rows.append(jdx)
                cols.append(i)
                vals.append(coef * sgn)
        return sp.csr_matrix((np.array(vals, dtype=complex), (rows, cols)),
                             shape=(self.dim, self.dim))

    # ----------------------------------------------------------- the terms
    def cell_terms(self):
        """One-directional off-diagonal terms (h.c. added later) and diagonal."""
        p = self.p
        off, diag = [], []
        for i, s in enumerate(self.sites):
            j, x, b, a, c = s
            # on-site carrier energies
            diag += [(p['eT'], [], [(b, 0), (a, 1)]),
                     (p['eH'], [], [(b, 0), (a, 0)]),
                     (p['eD'], [], [(b, 1), (a, 1)])]
            if c is not None:
                if p['muc']:
                    diag.append((p['muc'], [], [(c, 1)]))
                # conversion A: T + c0 -> H + c1   (deposit)
                off.append((-p['tauA'], [(c, 1), (a, 0)], [(b, 0)]))
                # conversion B: T + c1 -> D + c0   (pickup)
                off.append((-p['tauB'], [(b, 1), (c, 0)], [(a, 1)]))
            # intra-cell hops to the right neighbour (not across junctions)
            if x < self.n + 1:
                s2 = self.sites[i + 1]
                b2, a2 = s2[2], s2[3]
                off.append((-p['wT'] * np.exp(1j * p['phiT']),
                            [(b, 1), (a, 0), (a2, 1), (b2, 0)], []))
                off.append((-p['wH'], [(b, 1), (b2, 0)], [(a, 0), (a2, 0)]))
                off.append((-p['wD'], [(a2, 1), (a, 0)], [(b, 1), (b2, 1)]))
        return off, diag

    def junction_terms(self):
        """lambda-independent bond V, split as (hop one-direction, NN hermitian)."""
        hop, dens = [], []
        r2 = 1 / np.sqrt(2)
        for j in range(self.M - 1):
            R = self.sites[j * (self.n + 2) + self.n + 1]
            L = self.sites[(j + 1) * (self.n + 2)]
            for mr in (R[2], R[3]):
                for ml in (L[2], L[3]):
                    hop.append((0.5, [(mr, 1), (ml, 0)], []))
            NR = self._Nlist(R)
            NL = self._Nlist(L)
            for cr, opr in NR:
                for cl, opl in NL:
                    dens.append((self.p['kappa'] * cr * cl, opr + opl, []))
        return hop, dens

    @staticmethod
    def _Nlist(s):
        b, a = s[2], s[3]
        return [(0.5, [(b, 1), (b, 0)]), (0.5, [(b, 1), (a, 0)]),
                (0.5, [(a, 1), (b, 0)]), (0.5, [(a, 1), (a, 0)]), (-0.5, [])]

    def hamiltonians(self):
        off, diag = self.cell_terms()
        A = self.build(off)
        H0 = A + A.getH() + self.build(diag)
        hop, dens = self.junction_terms()
        B = self.build(hop)
        V = B + B.getH() + self.build(dens)
        return H0.tocsr(), V.tocsr()

    # ------------------------------------------------------- logical helpers
    def code_op(self, ops):
        """Matrix of a physical fermion operator restricted to the code basis."""
        nc = len(self.code)
        out = np.zeros((nc, nc), complex)
        pos = {int(c): q for q, c in enumerate(self.code)}
        for q, ci in enumerate(self.code):
            r = self.apply(ops, [], self.states[ci])
            if r is None:
                continue
            sgn, new = r
            k = self.index.get(new)
            if k is not None and k in pos:
                out[pos[k], q] += sgn
        return out


def parts(ch):
    """Unit-coefficient pieces so that H0 and V can be re-assembled quickly.

    H0 = eT*ET + eH*EH + eD*ED + muc*MU
         + [-wT e^{i phiT} AT - wH AH - wD AD - tauA CA - tauB CB + h.c.]
    V  = (BH + BH^dag) + kappa*BN
    """
    P = {}
    sites = ch.sites
    P['ET'] = ch.build([(1.0, [], [(s[2], 0), (s[3], 1)]) for s in sites])
    P['EH'] = ch.build([(1.0, [], [(s[2], 0), (s[3], 0)]) for s in sites])
    P['ED'] = ch.build([(1.0, [], [(s[2], 1), (s[3], 1)]) for s in sites])
    P['MU'] = ch.build([(1.0, [], [(s[4], 1)]) for s in sites if s[4] is not None])
    AT, AH, AD, CA, CB = [], [], [], [], []
    for i, s in enumerate(sites):
        j, x, b, a, c = s
        if c is not None:
            CA.append((1.0, [(c, 1), (a, 0)], [(b, 0)]))
            CB.append((1.0, [(b, 1), (c, 0)], [(a, 1)]))
        if x < ch.n + 1:
            b2, a2 = sites[i + 1][2], sites[i + 1][3]
            AT.append((1.0, [(b, 1), (a, 0), (a2, 1), (b2, 0)], []))
            AH.append((1.0, [(b, 1), (b2, 0)], [(a, 0), (a2, 0)]))
            AD.append((1.0, [(a2, 1), (a, 0)], [(b, 1), (b2, 1)]))
    for k, L in (('AT', AT), ('AH', AH), ('AD', AD), ('CA', CA), ('CB', CB)):
        P[k] = ch.build(L)
    hop, dens = [], []
    r = ch.junction_terms.__func__  # reuse geometry with kappa=1
    saved = ch.p['kappa']
    ch.p['kappa'] = 1.0
    hop, dens = ch.junction_terms()
    ch.p['kappa'] = saved
    P['BH'] = ch.build(hop)
    P['BN'] = ch.build(dens)
    return P


def assemble(P, p):
    off = (-p['wT'] * np.exp(1j * p['phiT']) * P['AT']
           - p['wH'] * np.exp(1j * p.get('phiH', 0.0)) * P['AH']
           - p['wD'] * np.exp(1j * p.get('phiD', 0.0)) * P['AD']
           - p['tauA'] * P['CA'] - p['tauB'] * P['CB'])
    H0 = (off + off.getH() + p['eT'] * P['ET'] + p['eH'] * P['EH']
          + p['eD'] * P['ED'] + p['muc'] * P['MU'])
    V = P['BH'] + P['BH'].getH() + p['kappa'] * P['BN']
    return H0.tocsr(), V.tocsr()
