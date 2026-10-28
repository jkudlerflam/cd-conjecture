"""Referee rebuild of the 576 x 576 vacancy-lifted bond, from the prose only.

The 23-letter bond is written from the term list of
contrib/explicit_candidate/README.md section 1.2 and the lift from
contrib/counterexample_campaign/COUNTEREXAMPLE_D24.md section 1.  Nothing is
imported from bond23.py or bond24.py.  The result is compared entry by entry
with contrib/counterexample_campaign/doping/bond24_fixture.npz, and the
claimed structural properties are checked.

Run:  python3 -B rebuild_bond.py      (numpy only, about one second)
Output: rebuild_bond.json beside this script.
"""
import json
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
FIXTURE = HERE.parent / 'counterexample_campaign' / 'doping' / 'bond24_fixture.npz'

# Alphabet in the order of the README table (section 1.1).
LET = ([('L', None, a) for a in (0, 1)] + [('R', None, a) for a in (0, 1)]
       + [('Xd', s, a) for s in 'v0+-' for a in (0, 1)]
       + [('Xx', s, a) for s in '0+-' for a in (0, 1)]
       + [('M0', None, None)] + [('M', s, None) for s in 'v0+-'])
assert len(LET) == 23
IX = {l: i for i, l in enumerate(LET)}
Q_S = {'v': 0, '0': 0, '+': 1, '-': -1}
EXC = '0+-'
ARROWS = {'L': ('L', 'L'), 'R': ('R', 'R'), 'Xd': ('L', 'R'), 'Xx': ('L', 'R'),
          'M0': ('R', 'L'), 'M': ('L', 'L')}


def zfac(s, a):
    return (-1) ** a if s == '0' else 1


def h23(p):
    D = 23
    H = np.zeros((D * D, D * D), complex)

    def pair(l, r):
        return IX[l] * D + IX[r]

    def hop(out, inn, c):          # c |out><inn| + h.c.
        H[pair(*out), pair(*inn)] += c
        H[pair(*inn), pair(*out)] += np.conj(c)

    for a in (0, 1):
        for b in (0, 1):
            hop((('L', None, a), ('Xd', 'v', b)), (('Xd', 'v', a), ('R', None, b)), -p['t'])
            for s in EXC:
                hop((('L', None, a), ('Xd', s, b)), (('Xd', s, a), ('R', None, b)),
                    -p['t'] * np.exp(1j * p['ph'][s]) * zfac(s, a))
                hop((('L', None, a), ('Xx', s, b)), (('Xx', s, a), ('R', None, b)),
                    -p['tx'] * np.exp(1j * p['phx'][s]) * zfac(s, a))
        hop((('L', None, a), ('M', 'v', None)), (('Xd', 'v', a), ('M0', None, None)), -p['t'])
        for s in EXC:
            hop((('L', None, a), ('M', s, None)), (('Xx', s, a), ('M0', None, None)),
                -p['tx'] * np.exp(1j * p['phx'][s]) * zfac(s, a))
            hop((('L', None, a), ('M', s, None)), (('Xd', s, a), ('M0', None, None)),
                -p['td'] * np.exp(1j * p['phd'][s]) * zfac(s, a))
            for m in [('M0', None, None)] + [('M', q, None) for q in 'v0+-']:
                hop((m, ('Xx', s, a)), (m, ('Xd', s, a)),
                    -p['tsw'] * np.exp(1j * p['phsw'][s]))
    # penalty on arrow mismatches
    for l in LET:
        for r in LET:
            if ARROWS[l[0]][1] != ARROWS[r[0]][0]:
                H[pair(l, r), pair(l, r)] += p['K']
    # port lambda V (already Hermitian as written)
    FM = {f: np.zeros((D, D)) for f in '0+-'}
    FX = {f: np.zeros((D, D)) for f in '0+-'}
    FM['0'][IX[('M', '0', None)], IX[('M', 'v', None)]] = 1
    FM['0'][IX[('M', 'v', None)], IX[('M', '0', None)]] = 1
    FM['+'][IX[('M', '+', None)], IX[('M', 'v', None)]] = 1
    FM['+'][IX[('M', 'v', None)], IX[('M', '-', None)]] = 1
    for a in (0, 1):
        FX['0'][IX[('Xd', '0', a)], IX[('Xd', 'v', a)]] = 1
        FX['0'][IX[('Xd', 'v', a)], IX[('Xd', '0', a)]] = 1
        FX['+'][IX[('Xd', '+', a)], IX[('Xd', 'v', a)]] = 1
        FX['+'][IX[('Xd', 'v', a)], IX[('Xd', '-', a)]] = 1
    FM['-'] = FM['+'].T.copy()
    FX['-'] = FX['+'].T.copy()
    H += p['lam'] * (np.kron(FM['0'], FX['0'])
                     + np.exp(1j * p['phi']) * np.kron(FM['+'], FX['-'])
                     + np.exp(-1j * p['phi']) * np.kron(FM['-'], FX['+']))
    # on-site operator o, half on each factor
    o = np.zeros((D, D), complex)
    dl = {'0': 0.0, '+': p['delta'], '-': -p['delta']}
    for i, (k, s, a) in enumerate(LET):
        if a == 1:
            o[i, i] += p['eps']
        if s is not None and s in EXC:
            o[i, i] += p['Delta'] + dl[s]
            if k == 'Xx':
                o[i, i] += p['Dx']
            if k == 'M':
                o[i, i] += p['DM']
        if k in ('M0', 'M'):
            o[i, i] += p['musel']
    for (lo, li, g) in ((('Xd', '+', 0), ('Xd', '0', 1), p['gp']),
                        (('Xd', '-', 1), ('Xd', '0', 0), p['gm'])):
        o[IX[lo], IX[li]] += g
        o[IX[li], IX[lo]] += np.conj(g)
    I = np.eye(D)
    H += 0.5 * (np.kron(o, I) + np.kron(I, o))
    return H


def h24(h, t, J, mu):
    d = 24
    out = np.zeros((d * d, d * d), complex)
    occ = np.array([a * d + b for a in range(1, d) for b in range(1, d)])
    out[np.ix_(occ, occ)] = J * h
    for a in range(1, d):
        out[0 * d + a, a * d + 0] += -t       # |vac,a><a,vac|
        out[a * d + 0, 0 * d + a] += -t       # |a,vac><vac,a|
    nvec = np.r_[0.0, np.ones(23)]
    out -= 0.5 * mu * np.diag(np.add.outer(nvec, nvec).reshape(-1))
    return out


def comm_err(h, v):
    tot = np.add.outer(v, v).reshape(-1)
    return float(np.abs(h * (tot[:, None] - tot[None, :])).max())


def main():
    # fixture parameters as recorded in bond24_checks.json
    p = dict(t=1.0, tx=1.0, tsw=1.0, td=0.013, Delta=2.5, delta=0.0, Dx=7.0, DM=8.0,
             gp=0.031, gm=0.027, eps=0.017, musel=1.1, K=10.0, lam=0.023, phi=np.pi / 2,
             ph={'0': .21, '+': -.18, '-': .32}, phx={'0': .44, '+': .07, '-': -.11},
             phsw={'0': .12, '+': -.13, '-': .19}, phd={'0': .24, '+': .39, '-': -.26})
    H23 = h23(p)
    H24 = h24(H23, t=1.3, J=0.0027, mu=0.0)
    fix = np.load(FIXTURE)
    diff = np.abs(H24 - fix['h'])
    nvec = np.r_[0, np.ones(23)]
    q = np.r_[0, [(a or 0) + (Q_S[s] if s else 0) for (_, s, a) in LET]]
    mk = np.r_[0, [1 if k in ('M0', 'M') else 0 for (k, _, _) in LET]]
    res = dict(
        fixture_max_abs_difference=float(diff.max()),
        fixture_entries_differing_above_1em15=int((diff > 1e-15).sum()),
        nonzero_rebuilt=int(np.count_nonzero(np.abs(H24) > 0)),
        nonzero_fixture=int(np.count_nonzero(fix['h'])),
        hermiticity=float(np.abs(H24 - H24.conj().T).max()),
        occupied_number_commutator=comm_err(H24, nvec),
        internal_charge_commutator=comm_err(H24, q),
        marker_number_commutator=comm_err(H24, mk),
        fixture_charge_vector_matches=bool(np.array_equal(q, fix['charge'])),
        fixture_marker_vector_matches=bool(np.array_equal(mk, fix['marker'])),
    )
    # with mu != 0 the lifted bond is still Hermitian and conserving
    H24mu = h24(H23, t=1.3, J=0.0027, mu=0.37)
    res['with_mu_hermiticity'] = float(np.abs(H24mu - H24mu.conj().T).max())
    res['with_mu_number_commutator'] = comm_err(H24mu, nvec)
    # squeezed-word check: vacancy swap preserves label order, h23 preserves
    # marker positions (every off-diagonal entry keeps each site's marker bit)
    off = np.argwhere(np.abs(H23 - np.diag(np.diag(H23))) > 0)
    keep = all(mk[1 + i // 23] == mk[1 + j // 23] and mk[1 + i % 23] == mk[1 + j % 23]
               for i, j in off)
    res['h23_preserves_marker_positions_sitewise'] = bool(keep)
    (HERE / 'rebuild_bond.json').write_text(json.dumps(res, indent=2) + '\n')
    print(json.dumps(res, indent=2))


if __name__ == '__main__':
    main()
