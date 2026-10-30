"""Check of the charged MPO and block-sparse DMRG against exact diagonalization
(L = 4 cells, fixed payload and length numbers). python3 -B edcheck.py -> edcheck.json"""
import json
import numpy as np
from model2 import build
from bsdmrg import local_charges, mpo_charged, finite_dmrg, QSHIFT

qs = local_charges()
out = []
for name in ('M0', 'M1', 'M2', 'C1'):
    h2, h1, info, ops, T = build(name)
    W, delta, err, viol = mpo_charged(h2, h1, 8, qs)
    L, d = 4, 8
    H = np.zeros((d ** L, d ** L), complex)
    for i in range(L):
        H += np.kron(np.kron(np.eye(d ** i), h1), np.eye(d ** (L - i - 1)))
    for i in range(L - 1):
        H += np.kron(np.kron(np.eye(d ** i), h2), np.eye(d ** (L - i - 2)))
    Q = np.array([sum(qs[x] for x in np.unravel_index(k, [d] * L)) for k in range(d ** L)])
    for Nl in (0, 1, 2):
        Qt = 4 * QSHIFT + Nl
        sel = np.where(Q == Qt)[0]
        e = np.linalg.eigvalsh(H[np.ix_(sel, sel)])[0]
        o = finite_dmrg(W, delta, L, Qt, [256] * 4, qs)
        rec = dict(model=name, Nl=Nl, E_ED=float(e), E_DMRG=o['E'], diff=o['E'] - float(e), charge_violation=viol, mpo_err=err)
        print(rec, flush=True)
        out.append(rec)
json.dump(out, open('edcheck.json', 'w'), indent=1)
