"""(6) Referee rebuild of the C83 14-letter coloured-vacancy bond from C83
PROOF.md eq. (3), on top of my12.py (C80 PROOF Section 1), compared entry by
entry with C83's bond14.py. Also: Hermiticity, the three charges (internal
U(1), occupied number, marker number), and the uncoloured 13-letter variant.
Argument: path to a scratch copy of codex/work_C83."""
import sys, json, os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import numpy as np
import my12

W = sys.argv[1]
sys.path.insert(0, W)
import bond14  # object under test

LET = my12.LETTERS + ['V0', 'V1']
# C83 PROOF Section 1 table
O = {**{l: 0 for l in ['L0', 'L1', 'E0', 'E1', 'MA', 'MB']}, **{l: 1 for l in ['R0', 'R1', 'DA', 'DB', 'XA', 'XB']},
     'V0': 0, 'V1': 1}
I = {'L0': (0,), 'L1': (0,), 'R0': (0, 1), 'R1': (1,), 'DA': (0,), 'DB': (0,), 'XA': (0,), 'XB': (0,),
     'E0': (0, 1), 'E1': (1,), 'MA': (0,), 'MB': (0,), 'V0': (0,), 'V1': (1,)}
ext_allowed = lambda a, b: O[a] in I[b]
grammar_ok = all(ext_allowed(a, b) == my12.allowed(a, b) for a in my12.LETTERS for b in my12.LETTERS)


def my14(pm, lam, J, t, mu, coloured=True):
    ls = LET if coloured else my12.LETTERS + ['V']
    d = len(ls); ix = {l: i for i, l in enumerate(ls)}
    h = np.zeros((d * d, d * d), complex)
    occ = [d * i + j for i in range(12) for j in range(12)]
    h[np.ix_(occ, occ)] = J * my12.bond(pm, lam)
    for a in my12.LETTERS:
        if coloured:
            for q in I[a]:
                src = ix[f'V{q}'] * d + ix[a]; dst = ix[a] * d + ix[f'V{O[a]}']
                h[dst, src] += -t; h[src, dst] += -t
        else:
            src = ix['V'] * d + ix[a]; dst = ix[a] * d + ix['V']
            h[dst, src] += -t; h[src, dst] += -t
    for a in ls:
        for b in ls:
            k = ix[a] * d + ix[b]
            if coloured and (a[0] == 'V' or b[0] == 'V') and not ext_allowed(a, b):
                h[k, k] += J * pm['K']
            h[k, k] -= mu * ((a[0] != 'V') + (b[0] != 'V')) / 2
    return h, ls


out = dict(grammar_matches_C80_table=grammar_ok, tests=[])
name = {('L', 0): 'L0', ('L', 1): 'L1', ('R', 0): 'R0', ('R', 1): 'R1', ('d', 'A'): 'DA', ('d', 'B'): 'DB',
        ('x', 'A'): 'XA', ('x', 'B'): 'XB', ('e', 0): 'E0', ('e', 1): 'E1', ('m', 'A'): 'MA', ('m', 'B'): 'MB',
        ('V', 0): 'V0', ('V', 1): 'V1'}
for n, lam, J, t, mu, K in [(2, .1, 1., 1., 0., 100.), (13, .37, .3, 1.7, .23, 1e4), (3, .05, 2., .4, -.5, 37.)]:
    flux = np.pi / (2 * n) + .011
    pc = bond14.params(n, D=7., g=.071, td=.003, K=K, musel=.4, field=.17, flux=flux)
    pm = my12.default_params(n, D=7., Delta=pc['Delta'], g=.071, td=.003, phiA=flux, phiB=0.,
                             psiA=n * flux + np.pi / 2, psiB=np.pi / 2, K=K, musel=.4, field=.17)
    for col in (True, False):
        theirs = bond14.bond(pc, lam=lam, J=J, t=t, mu=mu, colored=col)
        mine, ls = my14(pm, lam, J, t, mu, col)
        d = len(ls)
        tl = bond14.LET if col else bond14.MAT + [('V', 0)]
        perm = [ls.index(name[l]) if (col or l[0] != 'V') else ls.index('V') for l in tl]
        P2 = np.array([perm[i] * d + perm[j] for i in range(d) for j in range(d)])
        mp = mine[np.ix_(P2, P2)]
        q = np.array([my12.charge(l) if l[0] != 'V' else 0 for l in ls]); qq = (q[:, None] + q[None, :]).ravel()
        no = np.array([l[0] != 'V' for l in ls]).astype(int); nn = (no[:, None] + no[None, :]).ravel()
        mk = np.array([l in my12.MARKERS for l in ls]).astype(int); mm = (mk[:, None] + mk[None, :]).ravel()
        rec = dict(n=n, coloured=col, dim=d * d, max_entry_difference=float(np.max(abs(mp - theirs))),
                   hermiticity=float(np.max(abs(mine - mine.conj().T))),
                   charge_comm=float(np.max(abs(mine * (qq[:, None] - qq[None, :])))),
                   occupied_comm=float(np.max(abs(mine * (nn[:, None] - nn[None, :])))),
                   marker_number_comm=float(np.max(abs(mine * (mm[:, None] - mm[None, :])))))
        print(rec); out['tests'].append(rec)
json.dump(out, open('check_bond14.json', 'w'), indent=1)
