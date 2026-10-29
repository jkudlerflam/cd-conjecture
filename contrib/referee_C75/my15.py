"""Independent rebuild of the C75 15-letter bond from PROOF.md Section 1.

Written by the referee from the text only. It does not import bond15.py,
bond17.py or bond18.py. Letters are named by strings.
"""
import itertools
import numpy as np
import scipy.sparse as sp

LETTERS = (['L0', 'L1', 'R0', 'R1']
           + [f'D{s}{a}' for s in 'AB' for a in (0, 1)]
           + [f'X{s}{a}' for s in 'AB' for a in (0, 1)]
           + ['M0', 'MA', 'MB'])
IX = {l: i for i, l in enumerate(LETTERS)}
d = len(LETTERS)
assert d == 15
MARKERS = ['M0', 'MA', 'MB']


def kind(l):
    return l[0] if l[0] != 'M' else l


def bit(l):
    return int(l[-1]) if l[0] in 'LRDX' else None


def charge(l):
    if l[0] in 'LR':
        return 2 * bit(l)
    if l[0] in 'DX':
        return 2 * bit(l) + (-1 if l[1] == 'A' else 1)
    return {'M0': 0, 'MA': -1, 'MB': 1}[l]


def allowed(l, r):
    """PROOF Section 1 table: markers and L may precede anything;
    R, D, X may only precede R_0, R_1, M_0."""
    if l[0] == 'M' or l[0] == 'L':
        return True
    return r in ('R0', 'R1', 'M0')


def default_params(n, **kw):
    p = dict(n=n, D=10., Delta=4 - 2 * np.cos(np.pi / (n + 2)), g=0.03, td=1e-3,
             phiA=np.pi / (2 * n), phiB=0., psiA=np.pi, psiB=np.pi / 2,
             portA=np.exp(-1j * np.pi / 4), portB=np.exp(1j * np.pi / 4),
             K=np.inf, musel=0., field=0.)
    p.update(kw)
    return p


def z(s, a):
    return (-1) ** a if s == 'A' else 1


def move_rules(p):
    """List of (in_left, in_right, out_left, out_right, amplitude), without h.c."""
    R = []
    ph = {'A': p['phiA'], 'B': p['phiB']}
    ps = {'A': p['psiA'], 'B': p['psiB']}
    for a, b in itertools.product((0, 1), repeat=2):
        R.append((f'R{a}', f'R{b}', f'L{a}', f'R{b}', -1.))
        for s in 'AB':
            R.append((f'D{s}{a}', f'R{b}', f'L{a}', f'D{s}{b}', -z(s, a)))
            R.append((f'X{s}{a}', f'R{b}', f'L{a}', f'X{s}{b}', -np.exp(1j * ph[s]) * z(s, a)))
    for a in (0, 1):
        R.append((f'R{a}', 'M0', f'L{a}', 'M0', -1.))
        for s in 'AB':
            R.append((f'X{s}{a}', 'M0', f'L{a}', f'M{s}', -np.exp(1j * ph[s]) * z(s, a)))
            R.append((f'D{s}{a}', 'M0', f'L{a}', f'M{s}', -p['td'] * np.exp(1j * ps[s]) * z(s, a)))
            for m in MARKERS:
                R.append((m, f'D{s}{a}', m, f'X{s}{a}', -1.))
    return R


def port_rules(p):
    R = []
    for a in (0, 1):
        R.append(('M0', f'R{a}', 'MA', f'DB{a}', p['portA']))
        R.append(('M0', f'R{a}', 'MB', f'DA{a}', p['portB']))
    return R


def onsite(p):
    o = np.zeros((d, d), complex)
    for l in LETTERS:
        i = IX[l]
        if l[0] in 'DX' or l in ('MA', 'MB'):
            o[i, i] += p['Delta']
        if l[0] == 'X' or l in ('MA', 'MB'):
            o[i, i] += p['D']
        if l[0] == 'M':
            o[i, i] += p['musel']
        if bit(l) is not None:
            o[i, i] += p['field'] * bit(l)
    i, j = IX['DB0'], IX['DA1']
    o[i, j] += p['g']
    o[j, i] += p['g']
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
    K = p['K']
    pen = np.diag([float(K) * (not allowed(l, r)) for l in LETTERS for r in LETTERS])
    o = onsite(p)
    return (two_site(move_rules(p)) + pen + lam * two_site(port_rules(p))
            + 0.5 * (np.kron(o, np.eye(d)) + np.kron(np.eye(d), o)))


# Many-site builders
def rule_table(rules):
    """dict (l,r) -> list of (l',r',amp), h.c. included, duplicates summed."""
    T = {}
    for il, ir, ol, orr, amp in rules:
        T.setdefault((il, ir), {}).setdefault((ol, orr), 0)
        T[(il, ir)][(ol, orr)] += amp
        T.setdefault((ol, orr), {}).setdefault((il, ir), 0)
        T[(ol, orr)][(il, ir)] += np.conj(amp)
    return {k: [(a, b, v) for (a, b), v in w.items()] for k, w in T.items()}


def nbad(word, left=None):
    w = ((left,) if left is not None else ()) + tuple(word)
    return sum(not allowed(a, b) for a, b in zip(w, w[1:]))


def closure(seeds, p, bonds, project=False, left=None, with_port_at=None, cap=400000):
    """Invariant closure of seed words under onsite terms, the move rules on
    the listed bonds (pairs of positions), and the port on bond with_port_at.
    left is an optional frozen letter to the left of position 0 (bond (-1,0)
    is then included in bonds as (-1,0)). If project, transitions into words
    with a bad pair are dropped and the penalty is not added (K = infinity).
    Returns states, H0 (sparse), V (sparse)."""
    mv = rule_table(move_rules(p))
    pt = rule_table(port_rules(p))
    o = onsite(p)
    oc = {l: [(LETTERS[j], o[j, IX[l]]) for j in np.flatnonzero(o[:, IX[l]])] for l in LETTERS}
    states = [tuple(s) for s in seeds]
    ix = {s: i for i, s in enumerate(states)}
    hr, hc, hv, vr, vc, vv = [], [], [], [], [], []
    cur = 0

    def get(w):
        if w[0] == -1:
            return left
        return w

    while cur < len(states):
        w = states[cur]

        def emit(z_, val, isport=False):
            z_ = tuple(z_)
            if project and nbad(z_, left) > 0:
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
        for pos, l in enumerate(w):
            for l2, val in oc[l]:
                z_ = list(w); z_[pos] = l2; emit(z_, val)
        for (i, j) in bonds:
            a = left if i == -1 else w[i]
            b = w[j]
            for a2, b2, val in mv.get((a, b), []):
                if i == -1:
                    assert a2 == a
                z_ = list(w)
                if i != -1:
                    z_[i] = a2
                z_[j] = b2
                emit(z_, val)
            if (i, j) == with_port_at:
                for a2, b2, val in pt.get((a, b), []):
                    z_ = list(w); z_[i] = a2; z_[j] = b2; emit(z_, val, True)
        if not project and np.isfinite(p['K']):
            nb = sum(not allowed(left if i == -1 else w[i], w[j]) for (i, j) in bonds)
            if nb:
                hr.append(cur); hc.append(cur); hv.append(p['K'] * nb)
        cur += 1
    N = len(states)
    H = sp.coo_matrix((hv, (hr, hc)), shape=(N, N)).tocsc()
    V = sp.coo_matrix((vv, (vr, vc)), shape=(N, N)).tocsc() if vv else sp.csc_matrix((N, N))
    return states, H, V
