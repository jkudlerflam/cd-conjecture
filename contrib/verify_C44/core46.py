"""Independent literal build of C39's 46-letter initialized core.

Written from the prose of codex/work_C39/PROOFS.md Sections 1-2 and
codex/work_C32/PROOFS.md Section 1 only. Nothing is imported from any
C32, C39, C44 or C45 script.

Letters are tuples (kind, bus, species, bit):
  ('L',None,None,a), ('R',None,None,a)            data letters, a = 0,1
  ('d',0,'v',a)                                  vacuum token on the data rail
  ('A',None,'e',None), ('B',None,'e',None)       empty markers
  ('A',0,'v',None),    ('B',0,'v',None)          vacuum token on a marker
  ('d',r,s,a), ('x',r,s,a), ('A',r,s,None), ('B',r,s,None)
                                                 bus r = 1,2, species s = '0','+','-'
"""
import itertools
import numpy as np
import scipy.sparse as sp

SPECIES = ('0', '+', '-')
QCH = {'v': 0, '0': 0, '+': 1, '-': -1, 'e': 0, None: 0}


def alphabet():
    lets = []
    for a in (0, 1):
        lets += [('L', None, None, a), ('R', None, None, a), ('d', 0, 'v', a)]
    lets += [('A', None, 'e', None), ('B', None, 'e', None),
             ('A', 0, 'v', None), ('B', 0, 'v', None)]
    for r in (1, 2):
        for s in SPECIES:
            for a in (0, 1):
                lets += [('d', r, s, a), ('x', r, s, a)]
            lets += [('A', r, s, None), ('B', r, s, None)]
    assert len(lets) == 46 and len(set(lets)) == 46
    return lets


def arrows(l):
    """(incoming, outgoing) arrow of a letter; None marks a boundary side."""
    k, bus, s, a = l
    if k == 'L':
        return ('L', 'L')
    if k == 'R':
        return ('R', 'R')
    if k in ('d', 'x'):
        return ('L', 'R')
    if k == 'A':
        return (None, 'L' if s == 'e' else 'R')
    if k == 'B':
        return ('R' if s == 'e' else 'L', None)
    raise ValueError(l)


def charge(l):
    k, bus, s, a = l
    return (a or 0) + QCH[s]


def excited(l):
    return l[1] in (1, 2)


def high(l):
    return excited(l) and l[0] in ('A', 'B', 'x')


class Params:
    def __init__(self, n, a=4.0, D=0.2, g=0.4, q=0.8, w=0.8, v=0.65, eta=0.6, K=10.0):
        self.n, self.a, self.D, self.g, self.q = n, a, D, g, q
        self.w = {1: w, 2: w * v}
        self.eta = eta
        self.K = K
        self.phx = np.exp(1j * np.pi / (2 * n))   # neutral express phase


def onsite(l, P):
    """Onsite block O=(a-3)Q+D Q_high+g C acting on one letter: list of (l', amp)."""
    out = []
    if excited(l):
        out.append((l, (P.a - 3.0) + (P.D if high(l) else 0.0)))
    k, r, s, a = l
    if k == 'd' and r in (1, 2):
        # g|d(+,0)><d(0,1)| + g|d(-,1)><d(0,0)| + h.c.
        conv = {(  '0', 1): ('+', 0), ('+', 0): ('0', 1),
                ('0', 0): ('-', 1), ('-', 1): ('0', 0)}
        if (s, a) in conv:
            s2, a2 = conv[(s, a)]
            out.append((('d', r, s2, a2), P.g))
    return out


def kinetic(x, y, P):
    """Excited off-diagonal moves on an ordered pair, list of ((x',y'), amp).

    Rightward amplitudes: left link -1, switch -1, rail -w_r*(-1)^a[s=0],
    express and express-to-B -1*(e^{i pi/2n}(-1)^a)[s=0]. Reverse moves carry
    the complex conjugate amplitude.
    """
    out = []
    for r in (1, 2):
        for s in SPECIES:
            for a in (0, 1):
                z = (-1) ** a if s == '0' else 1
                ph = P.phx * z if s == '0' else 1.0
                pairs = [
                    ((('A', r, s, None), ('R', None, None, a)), (('A', None, 'e', None), ('d', r, s, a)), -1.0),
                    ((('A', None, 'e', None), ('d', r, s, a)), (('A', None, 'e', None), ('x', r, s, a)), -1.0),
                    ((('x', r, s, a), ('B', None, 'e', None)), (('L', None, None, a), ('B', r, s, None)), -1.0 * ph),
                ]
                for b in (0, 1):
                    pairs.append(((('d', r, s, a), ('R', None, None, b)), (('L', None, None, a), ('d', r, s, b)), -P.w[r] * z))
                    pairs.append(((('x', r, s, a), ('R', None, None, b)), (('L', None, None, a), ('x', r, s, b)), -1.0 * ph))
                for ini, fin, amp in pairs:
                    if (x, y) == ini:
                        out.append((fin, amp))
                    if (x, y) == fin:
                        out.append((ini, np.conj(amp)))
    return out


def vacuum_edges():
    """Forward vacuum edges (initial pair, final pair)."""
    E = []
    for a in (0, 1):
        E.append(((('A', 0, 'v', None), ('R', None, None, a)), (('A', None, 'e', None), ('d', 0, 'v', a))))
        E.append(((('d', 0, 'v', a), ('B', None, 'e', None)), (('L', None, None, a), ('B', 0, 'v', None))))
        for b in (0, 1):
            E.append(((('d', 0, 'v', a), ('R', None, None, b)), (('L', None, None, a), ('d', 0, 'v', b))))
    return E


VAC = vacuum_edges()


def bond_action(x, y, P):
    """Repeated positive bond h_c on |x y>: dict {(x',y'): amp}."""
    out = {}

    def add(key, amp):
        out[key] = out.get(key, 0.0) + amp
    # Q allocation: Q on the left letter, 2Q on the right letter
    d = (1.0 if excited(x) else 0.0) + (2.0 if excited(y) else 0.0)
    # grammar penalties
    ax, ay = arrows(x), arrows(y)
    if ax[1] is None or ay[0] is None or ax[1] != ay[0]:
        d += P.K
    if d:
        add((x, y), d)
    # onsite halves
    for l2, amp in onsite(x, P):
        add((l2, y), 0.5 * amp)
    for l2, amp in onsite(y, P):
        add((x, l2), 0.5 * amp)
    for fin, amp in kinetic(x, y, P):
        add(fin, amp)
    # vacuum rank-one edges (q|i>-|f>)(q<i|-<f|)
    q = P.q
    for ini, fin in VAC:
        if (x, y) == ini:
            add(ini, q * q)
            add(fin, -q)
        if (x, y) == fin:
            add(fin, 1.0)
            add(ini, -q)
    return out


def left_boundary(l, P):
    out = {}
    if excited(l):
        out[l] = 2.0
    for l2, amp in onsite(l, P):
        out[l2] = out.get(l2, 0.0) + 0.5 * amp
    if l[0] != 'A':
        out[l] = out.get(l, 0.0) + P.K
    return out


def right_boundary(l, P):
    out = {}
    if excited(l):
        out[l] = 1.0
    for l2, amp in onsite(l, P):
        out[l2] = out.get(l2, 0.0) + 0.5 * amp
    if l[0] != 'B':
        out[l] = out.get(l, 0.0) + P.K
    return out


def full_bond_matrix(P):
    """The 2116 x 2116 repeated bond, for positivity and Hermiticity checks."""
    lets = alphabet()
    idx = {l: i for i, l in enumerate(lets)}
    M = np.zeros((46 * 46, 46 * 46), complex)
    for x in lets:
        for y in lets:
            col = idx[x] * 46 + idx[y]
            for (x2, y2), amp in bond_action(x, y, P).items():
                M[idx[x2] * 46 + idx[y2], col] += amp
    return M


def legal_words(n):
    lets = alphabet()
    k = n + 2
    As = [l for l in lets if l[0] == 'A']
    Bs = [l for l in lets if l[0] == 'B']
    Ds = [l for l in lets if l[0] in ('L', 'R', 'd', 'x')]
    words = []
    for w0 in As:
        for mid in itertools.product(Ds, repeat=n):
            for wl in Bs:
                w = (w0,) + mid + (wl,)
                if all(arrows(w[i])[1] == arrows(w[i + 1])[0] for i in range(k - 1)):
                    words.append(w)
    return words


class Core:
    """Core cell operators in the legal sector of one initialized word."""

    def __init__(self, P):
        self.P = P
        n = P.n
        self.k = n + 2
        self.words = legal_words(n)
        self.index = {w: i for i, w in enumerate(self.words)}
        self.dim = len(self.words)
        assert self.dim == (13 * n + 14) * 2 ** n
        self.bonds = [self._site_op(j) for j in range(self.k - 1)]
        self.aL = self._boundary(0, left_boundary)
        self.bR = self._boundary(self.k - 1, right_boundary)
        self.Qch = np.array([sum(charge(l) for l in w) for w in self.words])

    def _site_op(self, j):
        rows, cols, vals = [], [], []
        for c, w in enumerate(self.words):
            for (x2, y2), amp in bond_action(w[j], w[j + 1], self.P).items():
                w2 = w[:j] + (x2, y2) + w[j + 2:]
                if abs(amp) == 0:
                    continue
                if w2 not in self.index:
                    raise RuntimeError('legal sector not invariant')
                rows.append(self.index[w2]); cols.append(c); vals.append(amp)
        return sp.csr_matrix((vals, (rows, cols)), shape=(self.dim, self.dim), dtype=complex)

    def _boundary(self, j, fn):
        rows, cols, vals = [], [], []
        for c, w in enumerate(self.words):
            for l2, amp in fn(w[j], self.P).items():
                w2 = w[:j] + (l2,) + w[j + 1:]
                rows.append(self.index[w2]); cols.append(c); vals.append(amp)
        return sp.csr_matrix((vals, (rows, cols)), shape=(self.dim, self.dim), dtype=complex)

    def H(self, Dmask=None):
        """Cell core Hamiltonian; Dmask[j]=0 removes occupied bond j entirely."""
        if Dmask is None:
            Dmask = [1] * (self.k - 1)
        H = self.aL + self.bR
        for j, dj in enumerate(Dmask):
            if dj:
                H = H + self.bonds[j]
        return H

    def token(self, w):
        for pos, l in enumerate(w):
            if l[1] is not None:
                return pos, l
        raise RuntimeError

    def sector(self, bus):
        """Indices of words whose token belongs to `bus` (0 = vacuum)."""
        return np.array([i for i, w in enumerate(self.words) if self.token(w)[1][1] == bus])

    def data_bits(self, w):
        bits = []
        for l in w[1:-1]:
            bits.append(l[3])
        return tuple(bits)

    def code(self):
        """Kernel basis G|b>, b in {0,1}^n, from the rules (checked against H)."""
        n, q = self.P.n, self.P.q
        S = sum(q ** (2 * l) for l in range(n + 2))
        G = np.zeros((self.dim, 2 ** n), complex)
        for i, w in enumerate(self.words):
            pos, l = self.token(w)
            if l[1] != 0:
                continue
            b = self.data_bits(w)
            col = int(''.join(map(str, b)), 2) if n else 0
            G[i, col] = q ** pos / np.sqrt(S)
        return G

    def source(self, end, bus, sigma):
        """Endpoint source F_(bus,sigma) at end 'A' (site 0) or 'B' (site k-1)."""
        j = 0 if end == 'A' else self.k - 1
        vac = (end, 0, 'v', None)
        exc = {s: (end, bus, s, None) for s in SPECIES}
        if sigma == '0':
            pairs = [(exc['0'], vac), (vac, exc['0'])]           # |0><v| + |v><0|
        elif sigma == '+':
            pairs = [(exc['+'], vac), (vac, exc['-'])]           # |+><v| + |v><-|
        else:
            pairs = [(vac, exc['+']), (exc['-'], vac)]           # F_- = F_+^dagger
        rows, cols, vals = [], [], []
        for c, w in enumerate(self.words):
            for out_l, in_l in pairs:
                if w[j] == in_l:
                    w2 = w[:j] + (out_l,) + w[j + 1:]
                    rows.append(self.index[w2]); cols.append(c); vals.append(1.0)
        return sp.csr_matrix((vals, (rows, cols)), shape=(self.dim, self.dim), dtype=complex)


PORT = {1: (1.0, np.pi / 2), 2: (None, np.pi)}   # b_r, phi_r ; b_2 = sqrt(eta)


def port_terms(P, r):
    b = 1.0 if r == 1 else np.sqrt(P.eta)
    phi = PORT[r][1]
    return b, [('0', '0', 1.0), ('+', '-', np.exp(1j * phi)), ('-', '+', np.exp(-1j * phi))]
