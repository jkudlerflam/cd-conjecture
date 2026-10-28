"""Independent finite-matrix checks of the pinned-marker pressure comparison.

Requires numpy. Run from any directory; prints JSON and writes
checks.json beside this script. No data from C44/C45 are imported.
"""

from itertools import product
from math import comb
from pathlib import Path
import json

import numpy as np


def logsumexp(values):
    largest = np.max(values)
    return largest + np.log(np.sum(np.exp(values-largest)))


def random_bond(d, q, rng):
    """Generic complex Hermitian bond conserving total marker count."""
    markers = np.array([int(a < q) + int(b < q)
                        for a, b in product(range(d), repeat=2)])
    z = rng.normal(size=(d*d, d*d)) + 1j*rng.normal(size=(d*d, d*d))
    z[markers[:, None] != markers[None, :]] = 0
    h = (z + z.conj().T)/2
    return h / np.max(np.abs(np.linalg.eigvalsh(h)))


def matrix_elements(h, d):
    return {
        (a, b): [(c, e, h[c*d+e, a*d+b])
                 for c, e in product(range(d), repeat=2)
                 if abs(h[c*d+e, a*d+b]) > 0]
        for a, b in product(range(d), repeat=2)
    }


def physical_matrix(d, length, entries):
    words = list(product(range(d), repeat=length))
    index = {w: i for i, w in enumerate(words)}
    h = np.zeros((len(words), len(words)), complex)
    for col, w in enumerate(words):
        for j in range(length):
            k = (j+1) % length
            for c, e, value in entries[w[j], w[k]]:
                out = list(w)
                out[j], out[k] = c, e
                h[index[tuple(out)], col] += value
    return words, h


def cells(word, q):
    starts = [j for j, a in enumerate(word) if a < q]
    assert starts[0] == 0
    ends = starts[1:] + [len(word)]
    return tuple(word[j:k] for j, k in zip(starts, ends))


def rooted_matrix(d, q, length, entries):
    """All cell lengths, periodic boundary, root seam moves retained.

    Unlike compression, a move of the root marker reanchors the word at
    that same marker. Its cell lists change only in the first/last cells.
    """
    words = [w for w in product(range(d), repeat=length) if w[0] < q]
    index = {w: i for i, w in enumerate(words)}
    h = np.zeros((len(words), len(words)), complex)
    seam_entries = 0
    for col, w in enumerate(words):
        before = cells(w, q)
        for j in range(length):
            k = (j+1) % length
            for c, e, value in entries[w[j], w[k]]:
                out = list(w)
                out[j], out[k] = c, e
                if out[0] >= q:
                    # Marker-number conservation identifies the moved root.
                    root = 1 if j == 0 else length-1
                    assert out[root] < q and (j == 0 or k == 0)
                    out = out[root:] + out[:root]
                    after = cells(tuple(out), q)
                    assert len(after) == len(before)
                    assert after[1:-1] == before[1:-1]
                    seam_entries += 1
                h[index[tuple(out)], col] += value
    return words, h, seam_entries


def main():
    rng = np.random.default_rng(20261028)
    cases = []
    for d, q, lengths in [(2, 1, range(3, 8)), (3, 1, range(3, 6)),
                           (3, 2, range(3, 6)), (4, 2, range(3, 5))]:
        bond = random_bond(d, q, rng)
        entries = matrix_elements(bond, d)
        for length in lengths:
            words, physical = physical_matrix(d, length, entries)
            pinned_index = [j for j, w in enumerate(words) if w[0] < q]
            pin = physical[np.ix_(pinned_index, pinned_index)]
            root_words, root, seam_entries = rooted_matrix(d, q, length, entries)
            assert root_words == [words[j] for j in pinned_index]
            herm_error = float(np.linalg.norm(root-root.conj().T))
            assert herm_error < 1e-12

            counts = {}
            for w in root_words:
                cell_list = cells(w, q)
                assert tuple(a for c in cell_list for a in c) == w
                count = len(cell_list)
                counts[count] = counts.get(count, 0)+1
            for count in range(1, length+1):
                expected = comb(length-1, count-1)*q**count*(d-q)**(length-count)
                assert counts[count] == expected
            assert sum(counts.values()) == q*d**(length-1)

            labels = np.array([len(cells(w, q)) for w in root_words])
            assert np.max(np.abs(root[labels[:, None] != labels[None, :]])) < 1e-14
            seam_norm = float(np.max(np.abs(np.linalg.eigvalsh(root-pin))))
            seam_bound = 2*d**4
            assert seam_norm <= seam_bound+1e-12
            spectra = [np.linalg.eigvalsh(a) for a in (physical, pin, root)]
            rows = []
            for beta in (0.2, 1.7, 5.0):
                log_z, log_pin, log_root = [float(logsumexp(-beta*ev))
                                          for ev in spectra]
                pin_residual = log_z-log_pin-np.log(d/q)
                full_residual = log_z-log_root-np.log(d/q)
                assert abs(pin_residual) <= 4*beta+1e-11
                # Use the measured seam norm as an independent sharper check.
                assert abs(full_residual) <= beta*(4+seam_norm)+1e-11
                rows.append(dict(beta=beta, pin_log_residual=float(pin_residual),
                                 rooted_log_residual=float(full_residual)))
            cases.append(dict(d=d, marker_rank=q, length=length,
                              physical_dimension=len(words),
                              pinned_dimension=len(root_words),
                              root_hermiticity_error=herm_error,
                              root_seam_entries=seam_entries,
                              root_seam_norm=seam_norm, trace_checks=rows))
    result = dict(status="all checks passed", cases=len(cases),
                  thermal_checks=sum(len(c["trace_checks"]) for c in cases),
                  scope="Generic complex marker-conserving bonds; no CFT inference",
                  results=cases)
    output = Path(__file__).with_name("checks.json")
    output.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({k: v for k, v in result.items() if k != "results"}, indent=2))


if __name__ == "__main__":
    main()
