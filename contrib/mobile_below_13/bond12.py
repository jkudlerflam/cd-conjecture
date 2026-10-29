"""d = 12 mobile bond: the eleven C84/C88 material letters plus ONE uncoloured
vacancy V.  V carries automaton state 0 only: o(V)=0, I(V)={0}.

Differences from C88's bond13 (V_0, V_1):
  * V_1 is deleted.  A vacancy may sit only on a bond whose left letter has
    outgoing state 0 (L_a, E, M_s, V) and whose right letter accepts state 0
    (everything except R_1).
  * Vacancy hopping is the swap  a V <-> V a  for a in S = {L_0, L_1, E, M_A, M_B},
    amplitude -t, bit independent.
  * Vacancy-context rail switch  -J(|V X_s><V D_s| + h.c.), the same amplitude as
    the marker-context switch, so a vacancy in an inter-cell gap does not cut the
    neighbouring cell's folded carrier.
  * Contact energy w on the ordered pair (V, R_0) (a vacancy immediately left
    of the walker).  w<0 makes the walker-mediated gap-gap interaction
    repulsive; see vacancy_scales.py and README.
  * Material on-site energies are genuine on-site terms (o(V)=0), not contact
    terms, so a vacancy next to a marker or a token changes no diagonal energy.
The material block is literally C88's material_bond (pinless C84 grammar).
All float parameters are diagnostic.  Python numpy only.
"""
import numpy as np

MAT = [('L', 0), ('L', 1), ('R', 0), ('R', 1), ('d', 'A'), ('d', 'B'),
       ('x', 'A'), ('x', 'B'), ('e', 0), ('m', 'A'), ('m', 'B')]
V = ('V', 0)
LET = MAT + [V]
IX = {a: i for i, a in enumerate(LET)}
NM = 11
D12 = 12
S_HOP = [('L', 0), ('L', 1), ('e', 0), ('m', 'A'), ('m', 'B')]


def outgoing(a):
    if a[0] == 'V':
        return 0
    return int(a[0] in ('R', 'd', 'x'))


def incoming(a):
    if a == ('R', 0) or a == ('e', 0):
        return (0, 1)
    if a == ('R', 1):
        return (1,)
    return (0,)


def allowed(a, b):
    return outgoing(a) in incoming(b)


def charge(a):
    if a[0] == 'V':
        return 0
    if a[0] in ('L', 'R', 'e'):
        return 2 * a[1]
    return -1 if a[1] == 'A' else 1


def marker(a):
    return a[0] in ('e', 'm')


def excited(a):
    return a[0] in ('d', 'x', 'm')


def params(n=3, **kw):
    """n counts physical data sites p (n-1 active bits), as in C88's bond13."""
    p = dict(n=n, D=4., Delta=4 - 2 * np.cos(np.pi / (n + 2)), g=.07,
             td=.003, K=100., musel=0., field=0., flux=np.pi / (2 * n), theta=-np.pi / 4, w=0.)
    p.update(kw)
    return p


def onsite(p):
    o = []
    for a in MAT:
        v = p['musel'] * marker(a)
        if a[0] in ('d', 'x', 'm'):
            v += p['Delta'] + p['D'] * (a[0] in ('x', 'm'))
        if a[0] in ('L', 'R', 'e'):
            v += p['field'] * a[1]
        o.append(v)
    return o + [0.]


def material_hop(p, lam):
    """Off-diagonal material terms on the 11x11 block (C88 material_bond without
    its diagonal).  Returns a 121x121 matrix."""
    d = 11
    h = np.zeros((d * d, d * d), complex)

    def add(out, inn, z):
        i = IX[out[0]] * d + IX[out[1]]
        j = IX[inn[0]] * d + IX[inn[1]]
        h[i, j] += z
        h[j, i] += np.conj(z)
    for a in (0, 1):
        la, ra, ea = ('L', a), ('R', a), ('e', a)
        add((la, ('R', 0)), (('R', 0), ra), -1.)
        if a == 0:
            add((la, ('e', 0)), (('R', 0), ea), -1.)
        for s in ('A', 'B'):
            z = (-1) ** a if s == 'A' else 1
            phi = p['flux'] if s == 'A' else 0.
            psi = p['n'] * p['flux'] + np.pi / 2 if s == 'A' else np.pi / 2
            add((la, ('d', s)), (('d', s), ra), -z)
            add((la, ('x', s)), (('x', s), ra), -z * np.exp(1j * phi))
            if a == 0:
                add((la, ('m', s)), (('x', s), ea), -z * np.exp(1j * phi))
                add((la, ('m', s)), (('d', s), ea), -p['td'] * z * np.exp(1j * psi))
    for m in [a for a in MAT if marker(a)]:
        for s in ('A', 'B'):
            add((m, ('x', s)), (m, ('d', s)), -1.)
    add((('d', 'B'), ('R', 0)), (('d', 'A'), ('R', 1)), p['g'])
    add((('m', 'A'), ('d', 'B')), (('e', 0), ('R', 0)), lam * np.exp(1j * p['theta']))
    add((('m', 'B'), ('d', 'A')), (('e', 0), ('R', 0)), lam * np.exp(-1j * p['theta']))
    return h


def bond(p, lam=.1, J=1., t=1., mu=0., vswitch=True, transparent=False):
    """144x144 bond.  transparent=True gives the C80/C75-F5 control: an
    uncoloured vacancy allowed on every bond and hopping across every letter."""
    d = D12
    h = np.zeros((d * d, d * d), complex)
    occ = [d * i + j for i in range(NM) for j in range(NM)]
    h[np.ix_(occ, occ)] = J * material_hop(p, lam)
    vi = IX[V]
    if vswitch and not transparent:
        for s in ('A', 'B'):
            i = vi * d + IX[('x', s)]
            j = vi * d + IX[('d', s)]
            h[i, j] += -J
            h[j, i] += -J
    hop = MAT if transparent else S_HOP
    for a in hop:
        i = IX[a]
        h[i * d + vi, vi * d + i] += -t
        h[vi * d + i, i * d + vi] += -t
    o = onsite(p)
    if not transparent:  # vacancy-walker contact energy on the ordered pair (V, R_0)
        k = vi * d + IX[('R', 0)]
        h[k, k] += J * p.get('w', 0.)
    for i, a in enumerate(LET):
        for j, b in enumerate(LET):
            k = d * i + j
            h[k, k] += J * (o[i] + o[j]) / 2
            if transparent:
                bad = (i < NM and j < NM) and not allowed(a, b)
            else:
                bad = not allowed(a, b)
            h[k, k] += J * p['K'] * bad
            h[k, k] -= mu * ((i < NM) + (j < NM)) / 2
    return h
