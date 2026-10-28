"""Exhaustive arrow repair, quantum-sector energies, and C8 penalty checks."""
from dataclasses import replace
import itertools
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh

from certified_periodic_bond import CertifiedPeriodicBond
from periodic_bond_checks import actual_columns
from periodic_clean_bond import PeriodicBond


def arrow_checks():
    cases = 0
    for L in range(2, 8):
        for state in itertools.product(range(4), repeat=L):
            # Local selector here is the integer 2*incoming+outgoing.
            defects = [(state[j] % 2) != (state[(j+1) % L] // 2) for j in range(L)]
            fixed = list(state)
            for j, defect in enumerate(defects):
                if defect: fixed[(j+1) % L] ^= 2
            assert all(fixed[j] % 2 == fixed[(j+1) % L] // 2 for j in range(L))
            assert sum(x != y for x, y in zip(state, fixed)) == sum(defects)
            cases += 1
    return cases


def quantum_checks():
    # Site labels contain incoming, outgoing, and an independent spectator bit.
    labels = list(itertools.product(range(2), repeat=3))
    index = {s: i for i, s in enumerate(labels)}
    rng = np.random.default_rng(9421)
    h = np.zeros((64, 64), complex)
    groups = {}
    for a, left in enumerate(labels):
        for b, right in enumerate(labels):
            key = left[0], right[1], left[1] ^ right[0]
            groups.setdefault(key, []).append(8*a+b)
    for members in groups.values():
        z = rng.normal(size=(len(members), len(members)))+1j*rng.normal(size=(len(members), len(members)))
        h[np.ix_(members, members)] = (z+z.conj().T)/2
    h /= np.linalg.norm(h, 2)
    B = 1.; K = 4.5
    rows = []
    for L in (3, 4):
        bitwords = list(itertools.product(range(2), repeat=L))
        codewords = list(itertools.product(bitwords, bitwords))
        sector_matrices = []
        for f in bitwords:
            basis = [tuple(index[(out[(j-1) % L] ^ f[(j-1) % L], out[j], spec[j])]
                           for j in range(L)) for out, spec in codewords]
            lookup = {s: i for i, s in enumerate(basis)}
            H = np.zeros((len(basis), len(basis)), complex)
            for col, state in enumerate(basis):
                for j in range(L):
                    nxtj = (j+1) % L
                    source = 8*state[j]+state[nxtj]
                    for target in np.flatnonzero(h[:, source]):
                        state2 = list(state); state2[j], state2[nxtj] = divmod(int(target), 8)
                        H[lookup[tuple(state2)], col] += h[target, source]
            assert np.max(np.abs(H-H.conj().T)) < 1e-13
            sector_matrices.append(H)
        elegal = float(eigh(sector_matrices[0], eigvals_only=True, subset_by_index=(0, 0))[0])
        for f, H in zip(bitwords, sector_matrices):
            k = sum(f)
            e = float(eigh(H, eigvals_only=True, subset_by_index=(0, 0))[0])+K*k
            change = eigh(H-sector_matrices[0], eigvals_only=True)
            change_norm = float(max(abs(change[0]), abs(change[-1])))
            assert change_norm <= 4*B*k+2e-12
            assert e >= elegal+(K-4*B)*k-2e-12
            rows.append(dict(length=L, defects=k, pattern=f, dimension=len(H),
                             energy=e, legal_energy=elegal, repair_norm=change_norm,
                             repair_bound=4*B*k, sector_lower=elegal+(K-4*B)*k))
    return dict(bond_norm=B, penalty=K, rows=rows)


def oracle_checks():
    H = CertifiedPeriodicBond(q=4, epsilon_max=1.)
    original = PeriodicBond(q=4, flat_vacancies=True)
    rng = np.random.default_rng(7103)
    legal_error = 0.; illegal_relative_error = 0.; columns = 0
    for _ in range(512):
        a = H.decode(int(rng.integers(H.d))); b = H.decode(int(rng.integers(H.d)))
        for sa, sb in itertools.product(range(4), repeat=2):
            aa, bb = replace(a, selector=sa), replace(b, selector=sb)
            key = H.encode(aa), H.encode(bb)
            new = H.column(*key, perturbation=.0005)
            old = original.column(*key, perturbation=.0005)
            for target in new.keys() | old.keys():
                diff = new.get(target, 0)-old.get(target, 0)
                expected = H.Ksel-original.Ksel if target == key and not H.legal(sa, sb) else 0.
                if H.legal(sa, sb): legal_error = max(legal_error, abs(diff))
                else: illegal_relative_error = max(illegal_relative_error, abs(diff-expected)/(1+abs(expected)))
            columns += 1
    assert legal_error == 0 and illegal_relative_error < 2e-15
    full = actual_columns(H)
    return dict(q=H.q, dimension=H.d, columns_compared=columns,
                legal_error=float(legal_error), illegal_relative_error=float(illegal_relative_error),
                old_penalty=H.original_selector_penalty, new_penalty=H.Ksel,
                norm_without_penalty_bound=H.repair_norm_bound,
                defect_gap_bound=H.selector_defect_gap,
                full_space_checks=full)


def main():
    out = dict(status='passed, selector-sector controls only',
               arrow_words=arrow_checks(), quantum=quantum_checks(), oracle=oracle_checks())
    Path('selector_repair_checks.json').write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps(dict(status='passed', arrow_words=out['arrow_words'],
                         quantum_sectors=len(out['quantum']['rows']), oracle=out['oracle']), indent=2))


if __name__ == '__main__':
    main()
