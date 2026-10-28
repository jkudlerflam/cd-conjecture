#!/usr/bin/env python3
"""Full local permutations and adversarial helper extraction checks."""

from itertools import product
import json
from pathlib import Path
import random

from parameter_oracle import (Site, alphabet, extract_until_end, extraction_word,
                              inverse_op, local, run_word, species)


def dimensions(p):
    u = p * p + p + 1
    q, logical = 2 * u, (p + 1) * u
    m = max(1, (q - 1).bit_length(), 2 * (logical - 1).bit_length())
    return q, m


def local_permutations(p):
    q, m = dimensions(p)
    symbols = alphabet(q, m)
    sites = [Site(a, b, c) for a in symbols for b in (0, 1) for c in (0, 1)]
    ops = ('PARAM', 'END', 'COPY_DIGIT', 'COPY_END', 'SHIFT', 'SHIFT_DAG')
    for op in ops:
        images = set()
        for site in sites:
            target = local(op, site, p, m)
            assert local(inverse_op(op), target, p, m) == site
            assert species(target.sigma) == species(site.sigma)
            images.add(target)
        assert len(images) == len(sites)
    assert p < 2 ** m
    return {'p': p, 'q': q, 'index_count_m': m, 'Sigma_dimension': len(symbols),
            'with_two_neutral_tracks': len(sites),
            'all_local_gate_columns_checked': len(ops) * len(sites)}


def expected_output(initial, p, m):
    out = list(initial)
    n = len(out)
    for j in range(m):
        at = n - 1 - j
        x = out[at]
        out[at] = Site(x.sigma, x.digit ^ ((p >> j) & 1),
                       x.end ^ int(j == m - 1))
    return tuple(out)


def exhaustive_small():
    # Exhaust every untrusted Sigma word and every value of BOTH target tracks.
    # The right-end H helper is pinned; other H/T symbols can be malformed.
    p, q, m, n = 2, 2, 2, 3
    symbols = alphabet(q, m)
    word = extraction_word(n, m)
    cases = 0
    for unknown in product(symbols, repeat=n - 1):
        sigma = unknown + (('H', 0, 0),)
        for tracks in product(range(4), repeat=n):
            initial = tuple(Site(x, t & 1, t >> 1) for x, t in zip(sigma, tracks))
            final = run_word(initial, word, p, m)
            assert final == expected_output(initial, p, m)
            assert run_word(final, word, p, m, reverse=True) == initial
            cases += 1
    return {'p': p, 'q': q, 'm': m, 'sites': n,
            'complete_promised_helper_basis_cases': cases,
            'arbitrary_target_tracks_included': True,
            'all_untrusted_payloads_restored': True}


def arbitrary_payload_walks():
    rng = random.Random(240317242)
    cases, vertices = 0, 0
    largest = None
    for p in (2, 3, 5, 7, 11):
        q, m = dimensions(p)
        symbols = alphabet(q, m)
        for extra in (0, 1, 5):
            n = m + extra
            for trial in range(20):
                initial = tuple([Site(rng.choice(symbols)) for _ in range(n - 1)]
                                + [Site(('H', 0, 0))])
                final, trajectory = extract_until_end(initial, p, m)
                assert trajectory == extraction_word(n, m)
                assert final == expected_output(initial, p, m)
                assert sum(a.end for a in final) == 1
                decoded = sum(final[n - 1 - j].digit << j for j in range(m))
                assert decoded == p
                assert run_word(final, trajectory, p, m, reverse=True) == initial
                cases += 1
                vertices += len(trajectory) + 1
                largest = {'p': p, 'sites': n, 'digits': m,
                           'elementary_gates': len(trajectory)}
    return {'malformed_payload_walks': cases, 'gate_vertices_checked': vertices,
            'last_example': largest, 'one_stopping_oracle_for_all_sizes': True}


def main():
    report = {
        'status': 'passed',
        'scope': 'Exact parameter oracle and extraction word; no compiled control table',
        'local_permutations': [local_permutations(p) for p in (2, 3, 5, 7, 11)],
        'complete_small_basis': exhaustive_small(),
        'uniform_algorithm': arbitrary_payload_walks(),
    }
    Path(__file__).with_name('parameter_results.json').write_text(
        json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
