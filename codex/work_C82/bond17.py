"""C82 two-color vacancy repair of the C75 bond. Resolved floats only."""
import numpy as np

MAT = [(k, s, a) for k, ss in [('L', ['']), ('R', ['']),
       ('D', ['A', 'B']), ('X', ['A', 'B'])] for s in ss for a in (0, 1)]
MAT += [('M', s, None) for s in ('0', 'A', 'B')]
LET = MAT + [('V', '0', None), ('V', '1', None)]
IX = {a: i for i, a in enumerate(LET)}


def outgoing(a):
    if a[0] == 'V':
        return int(a[1])
    return int(a[0] in ('R', 'D', 'X'))


def incoming(a):
    if a[0] == 'V':
        return (int(a[1]),)
    return (0, 1) if a[0] == 'R' or a == ('M', '0', None) else (0,)


def allowed(a, b):
    return outgoing(a) in incoming(b)


def params(n=3, **kw):
    p = dict(n=n, D=4., Delta=4-2*np.cos(np.pi/(n+2)), g=.07,
             td=.003, K=100., musel=1., field=.017,
             phiA=np.pi/(2*n), phiB=0., psiA=np.pi, psiB=np.pi/2)
    p.update(kw)
    return p


def material_bond(p, lam):
    d = 15
    h = np.zeros((d*d, d*d), complex)
    o = np.zeros((d, d), complex)

    def add(out, inn, amp):
        i = IX[out[0]]*d+IX[out[1]]
        j = IX[inn[0]]*d+IX[inn[1]]
        h[i, j] += amp
        h[j, i] += np.conj(amp)

    m0 = ('M', '0', None)
    for a in (0, 1):
        la, ra = ('L', '', a), ('R', '', a)
        for b in (0, 1):
            rb = ('R', '', b)
            add((la, rb), (ra, rb), -1.)
            for s in ('A', 'B'):
                z = (-1)**a if s == 'A' else 1
                for k, amp in [('D', -z), ('X', -z*np.exp(1j*p['phi'+s]))]:
                    add((la, (k, s, b)), ((k, s, a), rb), amp)
        add((la, m0), (ra, m0), -1.)
        for s in ('A', 'B'):
            z = (-1)**a if s == 'A' else 1
            ms = ('M', s, None)
            add((la, ms), (('X', s, a), m0), -z*np.exp(1j*p['phi'+s]))
            add((la, ms), (('D', s, a), m0), -p['td']*z*np.exp(1j*p['psi'+s]))
            for m in MAT[-3:]:
                add((m, ('X', s, a)), (m, ('D', s, a)), -1.)
        add((('M', 'A', None), ('D', 'B', a)), (m0, ra), lam*np.exp(-1j*np.pi/4))
        add((('M', 'B', None), ('D', 'A', a)), (m0, ra), lam*np.exp(1j*np.pi/4))
    for i, (k, s, a) in enumerate(MAT):
        o[i, i] = p['field']*(a or 0)+p['musel']*(k == 'M')
        if s in ('A', 'B'):
            o[i, i] += p['Delta']+p['D']*(k in ('M', 'X'))
    i, j = IX['D', 'A', 1], IX['D', 'B', 0]
    o[i, j] = o[j, i] = p['g']
    h += (np.kron(o, np.eye(d))+np.kron(np.eye(d), o))/2
    h += np.diag([p['K']*(not allowed(a, b)) for a in MAT for b in MAT])
    return h


def bond(p, lam=.03, J=1., t=1., mu=0.):
    d = 17
    h = np.zeros((d*d, d*d), complex)
    occ = np.array([d*a+b for a in range(15) for b in range(15)])
    h[np.ix_(occ, occ)] = J*material_bond(p, lam)
    for ai, a in enumerate(MAT):
        for q in incoming(a):
            i, j = (15+q)*d+ai, ai*d+15+outgoing(a)
            h[i, j] = h[j, i] = -t
    for i, a in enumerate(LET):
        for j, b in enumerate(LET):
            if i >= 15 or j >= 15:
                h[d*i+j, d*i+j] += J*p['K']*(not allowed(a, b))
            h[d*i+j, d*i+j] -= mu*((i < 15)+(j < 15))/2
    return h
