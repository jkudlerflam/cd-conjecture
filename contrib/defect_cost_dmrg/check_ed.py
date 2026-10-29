"""DMRG against independent sparse ED (C89 construction, open chain) and the exact free energy."""
import json
from model import selector, ed_open, free_energy_exact
from dmrg_run import ground

rows = []
for b, cfgs, V, J in [(3, [([3, 3], [0, 0]), ([2, 4], [0, 0]), ([3, 3], [1, 0]), ([6], [0])], 5, 0.04),
                      (3, [([3, 3, 3], [0, 0, 0]), ([2, 4, 3], [0, 0, 0]), ([3, 3, 3], [0, 1, 0])], 4, 0.3),
                      (4, [([4, 4], [0, 0]), ([3, 5], [0, 0])], 6, 0.1)]:
    mu = selector(b)
    for lengths, pins in cfgs:
        e, e0, contact, dim = ed_open(lengths, pins, V, J, mu)
        r = ground(lengths, pins, V, J, mu, G=V, chis=(16, 64))
        r0 = ground(lengths, pins, V, 0.0, mu, G=V, chis=(16, 64))
        ex0 = free_energy_exact(sum(lengths), V)
        row = dict(lengths=lengths, pins=pins, V=V, J=J, dim=dim, E_ed=e, E_dmrg=r['E'], dE=r['E'] - e,
                   E0_ed=e0, E0_dmrg=r0['E'], E0_exact=ex0,
                   contact_err=max(abs(a - c) for a, c in zip(r0['contact'], contact)))
        print(row, flush=True)
        rows.append(row)
        assert abs(row['dE']) < 1e-9 and abs(e0 - ex0) < 1e-9 and abs(r0['E'] - ex0) < 1e-9 and row['contact_err'] < 1e-8
json.dump(rows, open('check_ed.json', 'w'), indent=1)
print('ALL CHECKS PASS')
