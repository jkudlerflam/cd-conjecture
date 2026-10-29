"""Independent literal C75 bond and its one-state vacancy lift.

Only resolved diagnostic parameters should be passed to this float builder.
The certified hierarchy is specified by integers in PROOF.md.
"""
import numpy as np

LET = [(kind, s, a) for kind, species in [('L', ['']), ('R', ['']),
       ('D', ['A', 'B']), ('X', ['A', 'B'])] for s in species for a in (0, 1)]
LET += [('M', s, None) for s in ('0', 'A', 'B')]
IX = {x: i for i, x in enumerate(LET)}


def charge(x):
    kind, s, a = x
    return 2*(a or 0) + {'': 0, '0': 0, 'A': -1, 'B': 1}[s]


def allowed(x, y):
    return x[0] in ('M', 'L') or y[0] == 'R' or y == ('M', '0', None)


def params(n=1, **kw):
    p = dict(n=n, D=4., Delta=4-2*np.cos(np.pi/(n+2)), g=.07,
             td=.003, K=30., musel=1., field=.017,
             phiA=np.pi/(2*n), phiB=0., psiA=np.pi, psiB=np.pi/2)
    p.update(kw)
    return p


def parts(p):
    d = 15
    moves = np.zeros((d*d, d*d), complex)
    port = np.zeros_like(moves)
    onsite = np.zeros((d, d), complex)

    def add(mat, out, inn, amp):
        i = IX[out[0]]*d+IX[out[1]]
        j = IX[inn[0]]*d+IX[inn[1]]
        mat[i, j] += amp
        mat[j, i] += np.conj(amp)

    m0 = ('M', '0', None)
    for a in (0, 1):
        la, ra = ('L', '', a), ('R', '', a)
        for b in (0, 1):
            rb = ('R', '', b)
            add(moves, (la, rb), (ra, rb), -1.)
            for s in ('A', 'B'):
                z = (-1)**a if s == 'A' else 1
                for k, amp in [('D', -z), ('X', -z*np.exp(1j*p['phi'+s]))]:
                    add(moves, (la, (k, s, b)), ((k, s, a), rb), amp)
        add(moves, (la, m0), (ra, m0), -1.)
        for s in ('A', 'B'):
            z = (-1)**a if s == 'A' else 1
            ms = ('M', s, None)
            add(moves, (la, ms), (('X', s, a), m0), -z*np.exp(1j*p['phi'+s]))
            add(moves, (la, ms), (('D', s, a), m0), -p['td']*z*np.exp(1j*p['psi'+s]))
            for m in [x for x in LET if x[0] == 'M']:
                add(moves, (m, ('X', s, a)), (m, ('D', s, a)), -1.)
        add(port, (('M', 'A', None), ('D', 'B', a)), (m0, ra), np.exp(-1j*np.pi/4))
        add(port, (('M', 'B', None), ('D', 'A', a)), (m0, ra), np.exp(1j*np.pi/4))
    for i, (k, s, a) in enumerate(LET):
        onsite[i, i] = p['field']*(a or 0)+p['musel']*(k == 'M')
        if s in ('A', 'B'):
            onsite[i, i] += p['Delta']+p['D']*(k in ('M', 'X'))
    i, j = IX['D', 'A', 1], IX['D', 'B', 0]
    onsite[i, j] = onsite[j, i] = p['g']
    penalty = np.diag([p['K']*(not allowed(x, y)) for x in LET for y in LET])
    return moves, penalty, onsite, port


def occupied_bond(p, lam):
    a, b, o, v = parts(p)
    return a+b+lam*v+(np.kron(o, np.eye(15))+np.kron(np.eye(15), o))/2


def lift(h, t=1., J=.01, mu=0.):
    out = np.zeros((256, 256), complex)
    ix = np.array([16*a+b for a in range(1, 16) for b in range(1, 16)])
    out[np.ix_(ix, ix)] = J*h
    for a in range(1, 16):
        out[a, 16*a] = out[16*a, a] = -t
    number = np.r_[0, np.ones(15)]
    out -= mu/2*np.diag(np.add.outer(number, number).ravel())
    return out
