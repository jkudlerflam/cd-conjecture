"""Referee rebuild of the C84 11-letter bond from codex/work_C84/PROOF.md Section 1.

Written from the text of PROOF Section 1 (allowed-pair table, move table,
diagonal). It does not import bond11.py, parent_shift.py or bond18.py.
Letters are strings. Also provides the referee's legal-word enumeration,
decoder, and a generic invariant-closure builder (open chain with a frozen
marker on the left, or a ring).
"""
import itertools
import numpy as np
import scipy.sparse as sp

LETTERS = ['L0', 'L1', 'R0', 'R1', 'DA', 'DB', 'XA', 'XB', 'E', 'MA', 'MB']
IX = {l: i for i, l in enumerate(LETTERS)}
d = len(LETTERS)
MARKERS = ['E', 'MA', 'MB']
DATA = [l for l in LETTERS if l not in MARKERS]


def charge(l):
    if l in ('L0', 'L1', 'R0', 'R1'):
        return 2 * int(l[1])
    if l == 'E':
        return 0
    return -1 if l[1] == 'A' else 1


def excited(l):
    return l[0] in 'DX' or l in ('MA', 'MB')


GROUP1 = ('L0', 'L1', 'R0', 'DA', 'DB', 'XA', 'XB', 'E', 'MA', 'MB')
GROUP2 = ('L0', 'L1', 'R0', 'DA', 'DB', 'XA', 'XB')
GROUP3 = ('R0', 'R1', 'E')


def allowed(l, r):
    """PROOF Section 1 complete allowed ordered-pair table."""
    if l == 'L0' or l in MARKERS:
        return r in GROUP1
    if l == 'L1':
        return r in GROUP2
    return r in GROUP3


def default_params(p, **kw):
    """p = number of data sites (n = p - 1 active bits)."""
    q = dict(p=p, D=4., Delta=4 - 2 * np.cos(np.pi / (p + 2)), g=0.04, td=1e-3,
             phiA=np.pi / (2 * p), phiB=0., psiA=np.pi, psiB=np.pi / 2,
             K=np.inf, musel=0., field=0.)
    q.update(kw)
    return q


def z(s, a):
    return (-1) ** a if s == 'A' else 1


def move_rules(p):
    """(in_left, in_right, out_left, out_right, amplitude) for |out><in|; h.c. added later."""
    R = []
    ph = {'A': p['phiA'], 'B': p['phiB']}
    ps = {'A': p['psiA'], 'B': p['psiB']}
    for a in (0, 1):
        R.append(('R0', f'R{a}', f'L{a}', 'R0', -1.))
        for s in 'AB':
            R.append((f'D{s}', f'R{a}', f'L{a}', f'D{s}', -z(s, a)))
            R.append((f'X{s}', f'R{a}', f'L{a}', f'X{s}', -np.exp(1j * ph[s]) * z(s, a)))
    R.append(('R0', 'E', 'L0', 'E', -1.))
    for s in 'AB':
        R.append((f'X{s}', 'E', 'L0', f'M{s}', -np.exp(1j * ph[s])))
        R.append((f'D{s}', 'E', 'L0', f'M{s}', -p['td'] * np.exp(1j * ps[s])))
    for m in MARKERS:
        for s in 'AB':
            R.append((m, f'D{s}', m, f'X{s}', -1.))
    R.append(('DA', 'R1', 'DB', 'R0', p['g']))
    return R


def port_rules(p):
    return [('E', 'R0', 'MA', 'DB', np.exp(-1j * np.pi / 4)),
            ('E', 'R0', 'MB', 'DA', np.exp(1j * np.pi / 4))]


def onsite(p):
    o = np.zeros(d)
    for l in LETTERS:
        i = IX[l]
        if excited(l):
            o[i] += p['Delta']
        if l[0] == 'X' or l in ('MA', 'MB'):
            o[i] += p['D']
        if l in MARKERS:
            o[i] += p['musel']
        if l in ('L0', 'L1', 'R0', 'R1'):
            o[i] += p['field'] * int(l[1])
    return o


def two_site(rules, letters=LETTERS):
    ix = {l: i for i, l in enumerate(letters)}
    dd = len(letters)
    h = np.zeros((dd * dd, dd * dd), complex)
    for il, ir, ol, orr, amp in rules:
        i = ix[il] * dd + ix[ir]
        j = ix[ol] * dd + ix[orr]
        h[j, i] += amp
        h[i, j] += np.conj(amp)
    return h


def bond(p, lam=1.):
    pen = np.diag([float(p['K']) * (not allowed(l, r)) for l in LETTERS for r in LETTERS])
    o = np.diag(onsite(p))
    return (two_site(move_rules(p)) + pen + lam * two_site(port_rules(p))
            + 0.5 * (np.kron(o, np.eye(d)) + np.kron(np.eye(d), o)))


# ==== legal words
def legal_words(p):
    """All zero-penalty words of p data letters plus terminal marker, after any marker."""
    out = []
    for w in itertools.product(DATA, repeat=p):
        for m in MARKERS:
            word = ('E',) + w + (m,)
            if all(allowed(a, b) for a, b in zip(word, word[1:])):
                out.append(w + (m,))
    return out


def decode(word):
    """Referee's reading: (j, species, bits) with n = p-1 bits."""
    p = len(word) - 1
    j = next((i for i, l in enumerate(word[:p]) if l[0] != 'L'), p)
    if j < p:
        bits = [int(l[1]) for l in word[:j]] + [int(l[1]) for l in word[j + 1:p]]
        tok = word[j]
        sp_ = 'v' if tok == 'R0' else tok
    else:
        bits = [int(l[1]) for l in word[:p - 1]]
        sp_ = 'v' if word[p] == 'E' else word[p]
    return j, sp_, tuple(bits)


def encode_word(j, sp_, bits):
    """PROOF eqs. (3)-(4)."""
    n = len(bits)
    p = n + 1
    if j == p:
        return tuple(f'L{a}' for a in bits) + ('L0', 'E' if sp_ == 'v' else sp_)
    tok = 'R0' if sp_ == 'v' else sp_
    return (tuple(f'L{a}' for a in bits[:j]) + (tok,) + tuple(f'R{a}' for a in bits[j:]) + ('E',))


# ==== closures
def rule_table(rules):
    T = {}
    for il, ir, ol, orr, amp in rules:
        T.setdefault((il, ir), {}).setdefault((ol, orr), 0)
        T[(il, ir)][(ol, orr)] += amp
        T.setdefault((ol, orr), {}).setdefault((il, ir), 0)
        T[(ol, orr)][(il, ir)] += np.conj(amp)
    return {k: [(a, b, v) for (a, b), v in w.items() if v != 0] for k, w in T.items()}


def closure(seeds, p, bonds, left=None, port_bonds=(), project=False, cap=600000,
            rules=None, prules=None, allowed_fn=None, onsite_fn=None):
    """Invariant closure of seed words under the literal rules.
    bonds: list of (i,j) position pairs; (-1,0) couples the frozen letter `left`
    to site 0. port_bonds: bonds where the port acts. If project, transitions into
    words with a bad pair are dropped (K = infinity). Returns states, H0, V."""
    allowed_fn = allowed_fn or allowed
    mv = rule_table(rules if rules is not None else move_rules(p))
    pt = rule_table(prules if prules is not None else port_rules(p))
    o = onsite_fn(p) if onsite_fn else dict(zip(LETTERS, onsite(p)))
    states = [tuple(s) for s in seeds]
    ix = {s: i for i, s in enumerate(states)}
    hr, hc, hv, vr, vc, vv = [], [], [], [], [], []
    cur = 0

    def bad_count(w):
        return sum(not allowed_fn(left if i == -1 else w[i], w[j]) for (i, j) in bonds)

    while cur < len(states):
        w = states[cur]

        def emit(z_, val, isport=False):
            z_ = tuple(z_)
            if project and bad_count(z_) > 0:
                return
            if z_ not in ix:
                ix[z_] = len(states)
                states.append(z_)
                if len(states) > cap:
                    raise RuntimeError('cap')
            if isport:
                vr.append(ix[z_]); vc.append(cur); vv.append(val)
            else:
                hr.append(ix[z_]); hc.append(cur); hv.append(val)
        diag = sum(o.get(l, 0.) for l in w)
        for (i, j) in bonds:
            a = left if i == -1 else w[i]
            b = w[j]
            for a2, b2, val in mv.get((a, b), []):
                if i == -1 and a2 != a:
                    continue
                z_ = list(w)
                if i != -1:
                    z_[i] = a2
                z_[j] = b2
                emit(z_, val)
            if (i, j) in port_bonds:
                for a2, b2, val in pt.get((a, b), []):
                    z_ = list(w); z_[i] = a2; z_[j] = b2; emit(z_, val, True)
        if not project and np.isfinite(p['K']):
            diag += p['K'] * bad_count(w)
        if diag != 0:
            hr.append(cur); hc.append(cur); hv.append(diag)
        cur += 1
    N = len(states)
    H = sp.coo_matrix((hv, (hr, hc)), shape=(N, N)).tocsc()
    V = sp.coo_matrix((vv, (vr, vc)), shape=(N, N)).tocsc()
    return states, H, V


def cell_bonds(p, cells):
    """Open chain: frozen marker 'E' at -1, then cells of p data + marker."""
    L = cells * (p + 1)
    bonds = [(-1, 0)] + [(i, i + 1) for i in range(L - 1)]
    ports = tuple((c * (p + 1) + p, c * (p + 1) + p + 1) for c in range(cells - 1))
    return bonds, ports


def nbad(word, left='E'):
    w = (left,) + tuple(word)
    return sum(not allowed(a, b) for a, b in zip(w, w[1:]))
