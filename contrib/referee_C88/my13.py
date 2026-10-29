"""Referee rebuild of the C88 thirteen-letter mobile bond from codex/work_C88/PROOF.md
Section 1 (o, I table, move table, diagonal u and eq. (1)) alone.
It does not import bond13.py, bond11.py or any Astra file. Letters are strings.
Also provides the pin counter F (PROOF eq. (2)), ring helpers, the grammar
counters C, A, Z (PROOF Section 4) and a generic invariant-closure builder."""
import numpy as np
import scipy.sparse as sp

MAT = ['L0', 'L1', 'R0', 'R1', 'DA', 'DB', 'XA', 'XB', 'E', 'MA', 'MB']
LET = MAT + ['V0', 'V1']
MARKERS = ('E', 'MA', 'MB')
IXM = {l: i for i, l in enumerate(MAT)}
IX = {l: i for i, l in enumerate(LET)}

# PROOF Section 1 table
O = {'L0': 0, 'L1': 0, 'R0': 1, 'R1': 1, 'DA': 1, 'DB': 1, 'XA': 1, 'XB': 1,
     'E': 0, 'MA': 0, 'MB': 0, 'V0': 0, 'V1': 1}
I = {'L0': (0,), 'L1': (0,), 'R0': (0, 1), 'R1': (1,), 'DA': (0,), 'DB': (0,),
     'XA': (0,), 'XB': (0,), 'E': (0, 1), 'MA': (0,), 'MB': (0,), 'V0': (0,), 'V1': (1,)}


def is_vac(l):
    return l[0] == 'V'


def allowed(a, b):
    return O[a] in I[b]


def charge(l):
    if l[0] in 'LR':
        return 2 * int(l[1])
    if l in ('E', 'V0', 'V1', 'V'):
        return 0
    return -1 if l[1] == 'A' else 1


def params(p, **kw):
    q = dict(p=p, D=4., Delta=4 - 2 * np.cos(np.pi / (p + 2)), g=.07, td=.003,
             phiA=np.pi / (2 * p), phiB=0., psiA=np.pi, psiB=np.pi / 2, theta=-np.pi / 4,
             K=100., musel=0., field=0.)
    q.update(kw)
    return q


def z(s, a):
    return (-1) ** a if s == 'A' else 1


def material_rules(p, lam):
    """(in_l, in_r, out_l, out_r, amp) meaning amp |out><in|; h.c. added."""
    ph = {'A': p['phiA'], 'B': p['phiB']}
    ps = {'A': p['psiA'], 'B': p['psiB']}
    R = []
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
    R.append(('E', 'R0', 'MA', 'DB', lam * np.exp(1j * p['theta'])))
    R.append(('E', 'R0', 'MB', 'DA', lam * np.exp(-1j * p['theta'])))
    return R


def u_onsite(p):
    o = {}
    for l in MAT:
        v = 0.
        if l[0] in 'DX' or l in ('MA', 'MB'):
            v += p['Delta']
        if l[0] == 'X' or l in ('MA', 'MB'):
            v += p['D']
        if l in MARKERS:
            v += p['musel']
        if l[0] in 'LR':
            v += p['field'] * int(l[1])
        o[l] = v
    return o


def b11(p, lam):
    d = 11
    h = np.zeros((d * d, d * d), complex)
    for il, ir, ol, orr, amp in material_rules(p, lam):
        i = IXM[il] * d + IXM[ir]; j = IXM[ol] * d + IXM[orr]
        h[j, i] += amp; h[i, j] += np.conj(amp)
    o = u_onsite(p)
    for a in MAT:
        for b in MAT:
            k = IXM[a] * d + IXM[b]
            h[k, k] += (o[a] + o[b]) / 2
    return h


def h11(p, lam):
    pen = np.diag([p['K'] * (not allowed(a, b)) for a in MAT for b in MAT])
    return b11(p, lam) + pen


def bond(p, lam, J=1., t=1., mu=0., coloured=True):
    """PROOF eq. (1). Returns (h, letters). Uncoloured control uses one vacancy 'V'
    with no penalty on vacancy pairs (material pairs still penalised)."""
    ls = LET if coloured else MAT + ['V']
    dd = len(ls); ix = {l: i for i, l in enumerate(ls)}
    h = np.zeros((dd * dd, dd * dd), complex)
    occ = [ix[a] * dd + ix[b] for a in MAT for b in MAT]
    h[np.ix_(occ, occ)] = J * h11(p, lam)
    for a in MAT:
        if coloured:
            for q in I[a]:
                src = ix[f'V{q}'] * dd + ix[a]; dst = ix[a] * dd + ix[f'V{O[a]}']
                h[dst, src] += -t; h[src, dst] += -t
        else:
            src = ix['V'] * dd + ix[a]; dst = ix[a] * dd + ix['V']
            h[dst, src] += -t; h[src, dst] += -t
    for a in ls:
        for b in ls:
            k = ix[a] * dd + ix[b]
            if coloured and (is_vac(a) or is_vac(b)) and not allowed(a, b):
                h[k, k] += J * p['K']
            h[k, k] -= mu * ((not is_vac(a)) + (not is_vac(b))) / 2
    return h, ls


# ---- ring helpers
def squeeze(w):
    return [l for l in w if not is_vac(l)]


def F(w):
    m = squeeze(w)
    n = len(m)
    if n == 0:
        return 0
    return sum(m[i] == 'L1' and m[(i + 1) % n] in MARKERS for i in range(n))


def nbad(w, coloured=True):
    L = len(w)
    if coloured:
        return sum(not allowed(w[i], w[(i + 1) % L]) for i in range(L))
    return sum((not is_vac(w[i])) and (not is_vac(w[(i + 1) % L])) and not allowed(w[i], w[(i + 1) % L])
               for i in range(L))


def nmatbad(w):
    m = squeeze(w); n = len(m)
    return sum(not allowed(m[i], m[(i + 1) % n]) for i in range(n)) if n else 0


cnum = {a: len(I[a]) - 1 for a in MAT}


def counters(w):
    """PROOF Section 4 C, A, Z on a cyclic physical word."""
    L = len(w)
    C = sum(cnum[w[i]] for i in range(L) if not is_vac(w[i]) and is_vac(w[(i + 1) % L]))
    A = sum(1 for i in range(L) if O[w[i - 1]] == 1 and w[i] == 'R0' and w[(i + 1) % L] in ('R0', 'R1', 'E'))
    Zc = sum(1 for i in range(L) if O[w[i - 1]] == 1 and w[i] == 'E' and w[(i + 1) % L] == 'R0')
    return C, A, Zc


def transitions(h):
    n = h.shape[0]
    T = {}
    for col in range(n):
        rows = np.flatnonzero(abs(h[:, col]) > 0)
        T[col] = [(int(r), h[r, col]) for r in rows if r != col]
    return T


def ring_closure(seed, h, ls, coloured=True, project=False, cap=600000):
    """Invariant closure of a ring seed under the bond h with the penalty diagonal
    stripped out (h should be built at K=0). Returns states, H (no penalty), pen
    (penalty count per state, to be multiplied by J K), matbad."""
    dd = len(ls); ix = {l: i for i, l in enumerate(ls)}
    T = transitions(h); diag = np.real(np.diag(h))
    L = len(seed)
    states = [tuple(seed)]; idx = {states[0]: 0}
    r, c, v = [], [], []
    cur = 0
    while cur < len(states):
        w = states[cur]
        dg = 0.
        for i in range(L):
            a, b = w[i], w[(i + 1) % L]
            col = ix[a] * dd + ix[b]
            dg += diag[col]
            for row, val in T[col]:
                zz = list(w); zz[i] = ls[row // dd]; zz[(i + 1) % L] = ls[row % dd]; zz = tuple(zz)
                if project and nbad(zz, coloured):
                    continue
                if zz not in idx:
                    idx[zz] = len(states); states.append(zz)
                    if len(states) > cap:
                        raise RuntimeError('cap')
                r.append(idx[zz]); c.append(cur); v.append(val)
        if dg != 0:
            r.append(cur); c.append(cur); v.append(dg)
        cur += 1
    N = len(states)
    H = sp.coo_matrix((v, (r, c)), shape=(N, N)).tocsc()
    pen = np.array([nbad(w, coloured) for w in states], float)
    mb = np.array([nmatbad(w) for w in states], float)
    return states, H, pen, mb


def colour(seed):
    """Give each vacancy the colour of the preceding (cyclic) letter's outgoing state."""
    s = list(seed)
    L = len(s)
    for _ in range(L):
        for i in range(L):
            if s[i] in ('V', 'V0', 'V1'):
                prev = s[i - 1]
                if prev in ('V', 'V0', 'V1') and prev == 'V':
                    continue
                s[i] = f'V{O[prev]}'
    s = tuple(s)
    assert nbad(s) == 0, s
    return s


def encode(j, species, bits):
    """C84 cell word: p = len(bits)+1 data sites then a marker. Head j in 0..p-1,
    or j = p for the terminal all-L word. species 'v' is the vacuum."""
    p = len(bits) + 1
    if j == p:
        return tuple(f'L{a}' for a in bits) + ('L0', 'E' if species == 'v' else species)
    tok = 'R0' if species == 'v' else species
    return tuple(f'L{a}' for a in bits[:j]) + (tok,) + tuple(f'R{a}' for a in bits[j:]) + ('E',)
