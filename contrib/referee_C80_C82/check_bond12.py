"""(1) Entry-by-entry comparison of the referee rebuild (my12.py) with C80's
bond_shift.bond(p, size=12), plus Hermiticity, charge and marker checks.
Argument: path to a scratch copy of codex/work_C80."""
import sys, json, os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import numpy as np
import my12

W = sys.argv[1]
sys.path.insert(0, W)
import bond_shift, bond18  # objects under test

out = []
lab = bond_shift.alphabet(12)
name = {('L', 0): 'L0', ('L', 1): 'L1', ('R', 0): 'R0', ('R', 1): 'R1', ('d', 'A'): 'DA', ('d', 'B'): 'DB',
        ('x', 'A'): 'XA', ('x', 'B'): 'XB', ('e', 0): 'E0', ('e', 1): 'E1', ('m', 'A'): 'MA', ('m', 'B'): 'MB'}
perm = [my12.IX[name[l]] for l in lab]
P2 = np.array([perm[i] * 12 + perm[j] for i in range(12) for j in range(12)])
for n, lam, K in [(2, 1., 100.), (2, 0.37, 1e4), (12, 1., 2.0 ** 20), (12, 0.21, 37.)]:
    flux = np.pi / (2 * n) + 0.013  # generic flux, psi_A tied as in bond_shift
    pc = bond18.params(n, g=.071, td=.003, flux=flux, theta=-np.pi / 4, field=.17, musel=.4, K=K, D=7.)
    theirs = bond_shift.bond(pc, 12, lam=lam)
    pm = my12.default_params(n, D=7., Delta=pc['Delta'], g=.071, td=.003, phiA=flux, phiB=0.,
                             psiA=n * flux + np.pi / 2, psiB=np.pi / 2, K=K, musel=.4, field=.17)
    mine = my12.bond(pm, lam=lam)
    mine_p = mine[np.ix_(P2, P2)]
    q = np.array([my12.charge(l) for l in my12.LETTERS]); qq = (q[:, None] + q[None, :]).ravel()
    mk = np.array([l in my12.MARKERS for l in my12.LETTERS]).astype(int)
    m1 = np.kron(np.diag(mk), np.eye(12)); m2 = np.kron(np.eye(12), np.diag(mk))
    rec = dict(n=n, lam=lam, K=K, max_entry_difference=float(np.max(abs(mine_p - theirs))),
               nonzeros_mine=int(np.count_nonzero(mine)), nonzeros_theirs=int(np.count_nonzero(theirs)),
               hermiticity=float(np.max(abs(mine - mine.conj().T))),
               charge_commutator=float(np.max(abs(mine * (qq[:, None] - qq[None, :])))),
               marker_commutator_site1=float(np.max(abs(mine @ m1 - m1 @ mine))),
               marker_commutator_site2=float(np.max(abs(mine @ m2 - m2 @ mine))))
    print(rec); out.append(rec)
# The cited nonzero count (205) at a generic point with small K excludes nothing; report ours.
json.dump(out, open('check_bond12.json', 'w'), indent=1)
