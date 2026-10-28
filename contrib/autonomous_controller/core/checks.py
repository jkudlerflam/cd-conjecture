#!/usr/bin/env python3
"""Full local boundary checks, bulk samples, and actual initialized histories."""
from collections import defaultdict
from itertools import product
import json
from pathlib import Path
import random

from core import *


def local_columns():
    for work in range(WORK_DIM):
        assert WORK.f(WORK.f(work)) == work
        assert WORK.seed(WORK.seed(work)) == work
        a = work_decode(work)[0].static
        assert work_decode(WORK.f(work))[0].static == a
        assert work_decode(WORK.seed(work))[0].static == a
    accepted = [work for work in range(WORK_DIM) if initial_allowed(work)]
    assert len(accepted) == 10  # Five static marks, one free qubit each.
    checked = 0
    for column in (left_column, right_column):
        for source in range(LOCAL_DIM):
            for target, coefficient in column(source).items():
                assert column(target).get(source, 0) == coefficient.conjugate()
            checked += 1
    rng = random.Random(240317242)
    sources = {(rng.randrange(LOCAL_DIM), rng.randrange(LOCAL_DIM)) for _ in range(8000)}
    # Cover every clock-pair label with one-head states, including each of
    # the controls and both fresh/updated phases. Quantum bits are varied.
    for ca, cb in product(ALPHABET, repeat=2):
        for control in WORK.machine.controls:
            symbol_a = WORK.machine.symbols[rng.randrange(WORK.machine.B)]
            symbol_b = WORK.machine.symbols[rng.randrange(WORK.machine.B)]
            phase = rng.randrange(2)
            a = WORK.encode(symbol_a, rng.randrange(2), control.name, phase)
            b = WORK.encode(symbol_b, rng.randrange(2))
            if rng.randrange(2):
                a, b = b, a
            sources.add((encode(ca, a), encode(cb, b)))
    entries = 0
    for source in sources:
        for target, coefficient in bond_column(*source).items():
            assert bond_column(*target).get(source, 0) == coefficient.conjugate()
            assert (4 * coefficient.real).is_integer() and (4 * coefficient.imag).is_integer()
            entries += 1
    return dict(work_dimension=WORK_DIM, core_local_dimension=LOCAL_DIM,
                both_boundary_hermiticity_columns_exhausted=checked,
                full_work_F_and_seed_involution_columns=WORK_DIM,
                initialized_local_states=10,
                bulk_hermiticity_source_columns_sampled=len(sources),
                bulk_sample_nonzero_entries=entries,
                all_checked_entries_exact_complex_quarter_integers=True)


def clock_checks():
    results = []
    for n in range(1, 5):
        legal = set()
        transitions = 0
        for word in product(ALPHABET, repeat=n):
            defects = clock_defects(word)
            if not any(defects):
                legal.add(word)
            for target, label, site in clock_forward(word):
                assert clock_defects(target) == defects
                transitions += 1
        words, steps = clock_path(n)
        assert set(words) == legal and len(words) == len(legal)
        assert sum(label == "finish" for label, site in steps) == 2**n
        results.append(dict(n=n, full_clock_words=15**n, legal_vertices=len(legal),
                            all_directed_clock_transitions_checked=transitions))
    # The grammar is locally unique at every tested admissible length;
    # this is graph traversal through the five allowed static marks.
    for n in range(4, 30):
        paths = [("LEFT",)]
        for _ in range(n - 1):
            paths = [word + (right,) for word in paths
                     for left, right in STATIC_PAIRS if left == word[-1]]
        accepted = [word for word in paths if word[-1] == "RIGHT"]
        assert accepted == [("LEFT", "OUT_L") + ("BULK",) * (n-4) + ("OUT_R", "RIGHT")]
    return dict(all_clock_sectors=results, static_grammar_lengths=[4, 29])


def history(n, inputs):
    """Return one exact unnormalized history from a sparse quantum input."""
    words, steps = clock_path(n)
    current = {initial_work_word(n, bits): amplitude for bits, amplitude in inputs.items()}
    result = {}
    for t, clock_word in enumerate(words):
        for work_word, amplitude in current.items():
            result[tuple(encode(c, w) for c, w in zip(clock_word, work_word))] = amplitude
        if t == len(steps):
            break
        label, site = steps[t]
        updated = defaultdict(complex)
        for work_word, amplitude in current.items():
            for target, coefficient in clock_work_column(label, site, work_word).items():
                updated[target] += coefficient * amplitude
        current = cleaned(updated)
    return result


def apply_hamiltonian(n, vector):
    result = defaultdict(complex)
    for source, amplitude in vector.items():
        for target, coefficient in left_column(source[0]).items():
            result[(target,) + source[1:]] += coefficient * amplitude
        for target, coefficient in right_column(source[-1]).items():
            result[source[:-1] + (target,)] += coefficient * amplitude
        for site in range(n - 1):
            for pair, coefficient in bond_column(source[site], source[site+1]).items():
                result[source[:site] + pair + source[site+2:]] += coefficient * amplitude
    return cleaned(result)


def histories():
    # At n=4 the clock is too short to reach the demonstration's readout.
    # This case verifies every allowed input and the initialization theorem,
    # without attributing a readout weight to it.
    n = 4
    columns = []
    for bits in product((0, 1), repeat=n):
        column = history(n, {bits: 1+0j})
        assert not apply_hamiltonian(n, column)
        columns.append(column)
    C = len(clock_path(n)[0])
    for j, left in enumerate(columns):
        for k, right in enumerate(columns):
            overlap = sum(coefficient.conjugate() * right.get(word, 0)
                          for word, coefficient in left.items())
            assert overlap == (C if j == k else 0)
    # n=6 reaches genuine coherent gate execution. Use a superposed input
    # with exact dyadic complex coefficients, then test the complete global
    # Hamiltonian residual, not only the propagation recurrence.
    n = 6
    inputs = {tuple((x >> (n-1-j)) & 1 for j in range(n)): a
              for x, a in ((0, .5), (63, .5j), (43, -.5), (21, -.5j))}
    coherent = history(n, inputs)
    residual = apply_hamiltonian(n, coherent)
    assert not residual, max(map(abs, residual.values()), default=0)
    quantum_branched = len(coherent) > len(inputs) * len(clock_path(n)[0])
    assert quantum_branched
    n = 9
    input_bits = (1, 1, 1, 1, 0, 1, 1, 1, 1)
    complete = history(n, {input_bits: 1+0j})
    assert not apply_hamiltonian(n, complete)
    clock_words, clock_steps = clock_path(n)
    assert len(complete) > len(clock_words)
    completed_sweeps = [t + 1 for t, (label, site) in enumerate(clock_steps) if label == "finish"]
    cleanup_clock = clock_words[completed_sweeps[WORK.machine.cycle_steps(n) - 1]]
    cleanup_work = initial_work_word(n, input_bits)
    cleanup_work = (WORK.seed(cleanup_work[0]),) + cleanup_work[1:]
    expected_cleanup = tuple(encode(c, w) for c, w in zip(cleanup_clock, cleanup_work))
    cleanup_slice = {word: amplitude for word, amplitude in complete.items()
                     if tuple(decode(index)[0] for index in word) == cleanup_clock}
    assert cleanup_slice == {expected_cleanup: 1+0j}
    active_clocks = [set(), set()]
    for word in complete:
        for which, endpoint in enumerate((0, n - 1)):
            clock, work = decode(word[endpoint])
            symbol, record, control, phase = work_decode(work)
            if symbol.cursor and control is None:
                active_clocks[which].add(tuple(decode(index)[0] for index in word))
    assert list(map(len, active_clocks)) == [658, 656]
    return dict(n4_all_input_basis_histories=16, n4_clock_vertices=C,
                n4_complete_global_H_residual=0,
                n4_history_Gram_equals_C_times_identity=True,
                n6_coherent_clock_vertices=len(clock_path(6)[0]),
                n6_input_basis_support=len(inputs), n6_history_support=len(coherent),
                n6_contains_genuine_coherent_gate_branching=quantum_branched,
                n6_complete_global_H_residual=0,
                n9_input="|111101111>; zero-based site4=0 activates CV; classical input blank",
                n9_clock_vertices=len(clock_path(9)[0]),
                n9_history_support=len(complete),
                n9_genuine_coherent_gate_branching=True,
                n9_exact_quantum_and_classical_cleanup_at_first_complete_cycle=True,
                n9_complete_controller_cycles_within_clock=2**9 // WORK.machine.cycle_steps(9),
                n9_active_endpoint_clock_vertices=list(map(len, active_clocks)),
                n9_complete_global_H_residual=0,
                scope="Explicit zero-mode residuals; completeness and positive gap follow from the stated path-and-pins proof")


def main():
    result = dict(scope="Explicit initialized open-chain demonstration core; no high-channel or c-d violation claim",
                  local=local_columns(), clock=clock_checks(), initialized_histories=histories(),
                  all_passed=True)
    Path(__file__).with_name("results.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
