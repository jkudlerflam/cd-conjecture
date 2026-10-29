"""Referee rebuild of the C80 12-letter bond from PROOF.md Section 1 only.

Written from the text of codex/work_C80/PROOF.md Section 1. It does not
import bond_shift.py or bond18.py. Letters are strings.
Also contains a generic invariant-closure builder for open chains with
an optional frozen letter on the left, and the referee's own map from
legal words to (token position, species, bits).
"""
import itertools
import numpy as np
import scipy.sparse as sp

LETTERS = ['L0', 'L1', 'R0', 'R1', 'DA', 'DB', 'XA', 'XB', 'E0', 'E1', 'MA', 'MB']
IX = {l: i for i, l in enumerate(LETTERS)}
d = len(LETTERS)
MARKERS = ['E0', 'E1', 'MA', 'MB']
DATA = [l for l in LETTERS if l not in MARKERS]


def charge(l):
    if l[0] in 'LRE':
        return 2 * int(l[1])
    return -1 if l[1] == 'A' else 1


def allowed(l, r):
    """PROOF Section 1 allowed-pair table."""
    if l[0] == 'L' or l in MARKERS:
        return r in ('L0', 'L1', 'R0', 'DA', 'DB', 'XA', 'XB', 'E0', 'MA', 'MB')
    return r in ('R0', 'R1', 'E0', 'E1')


def default_params(n, **kw):
    p = dict(n=n, D=4., Delta=4 - 2 * np.cos(np.pi / (n + 2)), g=0.04, td=1e-3,
             phiA=np.pi / (2 * n), phiB=0., psiA=np.pi, psiB=np.pi / 2,
             K=np.inf, musel=0., field=0.)
    p.update(kw)
    return p


def z(s, a):
    return (-1) ** a if s == 'A' else 1


def move_rules(p):
    """(in_left, in_right, out_left, out_right, amplitude): amplitude |out><in|, h.c. added later."""
    R = []
    ph = {'A': p['phiA'], 'B': p['phiB']}
    ps = {'A': p['psiA'], 'B': p['psiB']}
    for a in (0, 1):
        R.append(('R0', f'R{a}', f'L{a}', 'R0', -1.))
        R.append(('R0', f'E{a}', f'L{a}', 'E0', -1.))
        for s in 'AB':
            R.append((f'D{s}', f'R{a}', f'L{a}', f'D{s}', -z(s, a)))
            R.append((f'X{s}', f'R{a}', f'L{a}', f'X{s}', -np.exp(1j * ph[s]) * z(s, a)))
            R.append((f'X{s}', f'E{a}', f'L{a}', f'M{s}', -np.exp(1j * ph[s]) * z(s, a)))
            R.append((f'D{s}', f'E{a}', f'L{a}', f'M{s}', -p['td'] * np.exp(1j * ps[s]) * z(s, a)))
    for m in MARKERS:
        for s in 'AB':
            R.append((m, f'D{s}', m, f'X{s}', -1.))
    R.append(('DA', 'R1', 'DB', 'R0', p['g']))
    R.append(('DA', 'E1', 'DB', 'E0', p['g']))
    return R


def port_rules(p):
    return [('E0', 'R0', 'MA', 'DB', np.exp(-1j * np.pi / 4)),
            ('E0', 'R0', 'MB', 'DA', np.exp(1j * np.pi / 4))]


def onsite(p):
    o = np.zeros(d)
    for l in LETTERS:
        i = IX[l]
        if l[0] in 'DX' or l in ('MA', 'MB'):
            o[i] += p['Delta']
        if l[0] == 'X' or l in ('MA', 'MB'):
            o[i] += p['D']
        if l in MARKERS:
            o[i] += p['musel']
        if l[0] in 'LRE':
            o[i] += p['field'] * int(l[1])
    return o


def two_site(rules):
    h = np.zeros((d * d, d * d), complex)
    for il, ir, ol, orr, amp in rules:
        i = IX[il] * d + IX[ir]
        j = IX[ol] * d + IX[orr]
        h[j, i] += amp
        h[i, j] += np.conj(amp)
    return h


def bond(p, lam=1.):
    pen = np.diag([float(p['K']) * (not allowed(l, r)) for l in LETTERS for r in LETTERS])
    o = np.diag(onsite(p))
    return (two_site(move_rules(p)) + pen + lam * two_site(port_rules(p))
            + 0.5 * (np.kron(o, np.eye(d)) + np.kron(np.eye(d), o)))


# ==== words
def legal_words(n):
    """All words of n data letters plus terminal marker, after any marker,
    with no bad pair. Brute force over the 8 data letters and 4 markers."""
    out = []
    for w in itertools.product(DATA, repeat=n):
        for m in MARKERS:
            word = ('E0',) + w + (m,)
            if all(allowed(a, b) for a, b in zip(word, word[1:])):
                out.append(w + (m,))
    return out


def decode(word):
    """Referee's reading of a legal word: (j, species, bits).
    species 'v' for the vacuum, 'DA','DB','XA','XB' or 'MA','MB'."""
    n = len(word) - 1
    j = next((i for i, l in enumerate(word[:n]) if l[0] != 'L'), n)
    bits = [int(l[1]) for l in word[:j]]
    if j < n:
        tok = word[j]
        sp_ = 'v' if tok == 'R0' else tok
        bits += [int(l[1]) for l in word[j + 1:n]] + [int(word[n][1])]
    else:
        sp_ = 'v' if word[n] == 'E0' else word[n]
    return j, sp_, tuple(bits)


def encode_word(j, sp_, bits):
    """PROOF eq. (2)-(3): token at j displaces bit j to the right."""
    n = len(bits)
    if j == n:
        return tuple(f'L{a}' for a in bits) + ('E0' if sp_ == 'v' else sp_,)
    tok = 'R0' if sp_ == 'v' else sp_
    tail = [f'R{a}' for a in bits[j:n - 1]]
    return tuple(f'L{a}' for a in bits[:j]) + (tok,) + tuple(tail) + (f'E{bits[n - 1]}',)


# ==== closures
def rule_table(rules):
    T = {}
    for il, ir, ol, orr, amp in rules:
        T.setdefault((il, ir), {}).setdefault((ol, orr), 0)
        T[(il, ir)][(ol, orr)] += amp
        T.setdefault((ol, orr), {}).setdefault((il, ir), 0)
        T[(ol, orr)][(il, ir)] += np.conj(amp)
    return {k: [(a, b, v) for (a, b), v in w.items() if v != 0] for k, w in T.items()}


def nbad(word, left=None):
    w = ((left,) if left is not None else ()) + tuple(word)
    return sum(not allowed(a, b) for a, b in zip(w, w[1:]))


def closure(seeds, p, bonds, left=None, port_bonds=(), project=False, cap=400000,
            rules=None, prules=None, allowed_fn=None):
    """Invariant closure of seed words. bonds: list of (i,j) position pairs,
    (-1,0) means the frozen letter `left` with site 0. port_bonds: bonds on
    which the port acts. If project, transitions into words with a bad pair
    are dropped (K = infinity). Returns states, H0, V (scipy csc)."""
    allowed_fn = allowed_fn or allowed
    mv = rule_table(rules if rules is not None else move_rules(p))
    pt = rule_table(prules if prules is not None else port_rules(p))
    o = onsite(p)
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
        diag = sum(o[IX[l]] for l in w)
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


def cell_bonds(n, cells):
    """Open chain: frozen marker 'E0' at -1, then cells of n data + marker."""
    L = cells * (n + 1)
    bonds = [(-1, 0)] + [(i, i + 1) for i in range(L - 1)]
    ports = tuple((c * (n + 1) + n, c * (n + 1) + n + 1) for c in range(cells - 1))
    return bonds, ports


def excited(l):
    return l[0] in 'DX' or l in ('MA', 'MB')
