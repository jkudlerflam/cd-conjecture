"""(1) C84 bond: referee rebuild (my11.py, from PROOF Section 1) against bond11.py,
entry by entry, at generic points. Hermiticity, U(1) charge, marker projectors.
Compression identity H_11 = J^dag H_12 J + K Pi_(L1,marker) against the previous
referee's independent C80 rebuild (contrib/referee_C80_C82/my12.py).
Inter-cell factorisation: on a (marker, data) bond the zero-port operator is
identity on the marker (F3 bookkeeping), and the pin sits on the intra-cell bond.
Usage: python3 check_bond11.py <scratch copy of codex/work_C84>"""
import sys, os, json
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(1, os.path.join(HERE, '..', 'referee_C80_C82'))
import my11, my12
W = sys.argv[1]
sys.path.insert(0, W)
import bond11  # object under test

name = {('L', 0): 'L0', ('L', 1): 'L1', ('R', 0): 'R0', ('R', 1): 'R1', ('d', 'A'): 'DA', ('d', 'B'): 'DB',
        ('x', 'A'): 'XA', ('x', 'B'): 'XB', ('e', 0): 'E', ('m', 'A'): 'MA', ('m', 'B'): 'MB'}
d = my11.d
out = dict(tests=[])
theirs_letters = [name[l] for l in bond11.alphabet()]
assert theirs_letters == my11.LETTERS, theirs_letters
out['allowed_table_matches_bond11'] = all(
    my11.allowed(name[a], name[b]) == bond11.allowed(a, b) for a in bond11.alphabet() for b in bond11.alphabet())
for p_, lam, K, flux, mus, fld in [(2, .1, 100., np.pi / 4 + .013, .4, .17), (12, .37, 2 ** 20, np.pi / 24, .9, .021),
                                   (3, .21, 37., .37, -.3, .5), (5, 1., 1e3, .11, 0., 0.)]:
    pt = dict(n=p_, D=7., Delta=4 - 2 * np.cos(np.pi / (p_ + 2)), g=.071, td=.0031, flux=flux, theta=-np.pi / 4,
              K=K, musel=mus, field=fld)
    pm = my11.default_params(p_, D=7., g=.071, td=.0031, phiA=flux, phiB=0., psiA=p_ * flux + np.pi / 2,
                             psiB=np.pi / 2, K=K, musel=mus, field=fld)
    theirs = bond11.bond(pt, lam=lam)
    mine = my11.bond(pm, lam)
    q = np.array([my11.charge(l) for l in my11.LETTERS]); qq = np.add.outer(q, q).ravel()
    mk = np.array([l in my11.MARKERS for l in my11.LETTERS]).astype(int)
    mL = np.kron(mk, np.ones(d)); mR = np.kron(np.ones(d), mk)
    # compression of the previous referee's C80 rebuild
    pm12 = my12.default_params(p_, D=7., Delta=pm['Delta'], g=.071, td=.0031, phiA=flux, phiB=0.,
                               psiA=p_ * flux + np.pi / 2, psiB=np.pi / 2, K=K, musel=mus, field=fld)
    h12 = my12.bond(pm12, lam)
    keep = [my12.IX['E0' if l == 'E' else l] for l in my11.LETTERS]
    sel = np.array([12 * i + j for i in keep for j in keep])
    comp = h12[np.ix_(sel, sel)]
    pin = np.diag([K * (l == 'L1' and r in my11.MARKERS) for l in my11.LETTERS for r in my11.LETTERS])
    rec = dict(p=p_, lam=lam, K=K,
               max_entry_difference_vs_bond11=float(np.max(abs(mine - theirs))),
               nonzero_entries_offdiag_and_diag=int(np.sum(abs(mine) > 1e-15)),
               hermiticity=float(np.max(abs(mine - mine.conj().T))),
               charge_commutator=float(np.max(abs(mine * (qq[:, None] - qq[None, :])))),
               marker_left_commutator=float(np.max(abs(mine * (mL[:, None] - mL[None, :])))),
               marker_right_commutator=float(np.max(abs(mine * (mR[:, None] - mR[None, :])))),
               compression_identity_error=float(np.max(abs(mine - (comp + pin)))))
    # F3 bookkeeping: zero-port bond with a marker on the left is identity on the marker
    h0 = my11.bond(pm, 0.).reshape(d, d, d, d)  # [l_out, r_out, l_in, r_in]
    o = my11.onsite(pm)
    dev = 0.
    for m in my11.MARKERS:
        for m2 in my11.MARKERS:
            blk = h0[my11.IX[m], :, my11.IX[m2], :].copy()
            if m == m2:
                blk -= np.eye(d) * o[my11.IX[m]] / 2
                ref = h0[my11.IX['E'], :, my11.IX['E'], :] - np.eye(d) * o[my11.IX['E']] / 2
                dev = max(dev, float(np.max(abs(blk - ref))))
            else:
                dev = max(dev, float(np.max(abs(blk))))
    rec['marker_left_bond_identity_on_marker_dev'] = dev
    print(rec, flush=True)
    out['tests'].append(rec)
# generic-point nonzero count as C84 states (166): count at a generic point with field, musel nonzero
out['note'] = ('C84 states 166 nonzero entries at a generic diagnostic point; our count includes the diagonal '
               'penalties and onsite terms at the same kind of point.')
json.dump(out, open(os.path.join(HERE, 'check_bond11.json'), 'w'), indent=1)
