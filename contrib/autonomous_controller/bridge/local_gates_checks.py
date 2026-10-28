#!/usr/bin/env python3
"""Full local F/seed basis and representative complete G block checks."""

from itertools import product
import json
from pathlib import Path

from local_gates import WORK, WORK_DIM, decode, encode, f_column, g_column, seed_column


def check_g_pair(left, right):
    source = (left, right)
    column = g_column(*source)
    assert abs(sum(abs(x) ** 2 for x in column.values()) - 1) < 1e-15
    twice = {}
    for target, coefficient in column.items():
        a, ra, qa, pa = decode(source[0])
        b, rb, qb, pb = decode(source[1])
        ta, _, _, _ = decode(target[0])
        tb, _, _, _ = decode(target[1])
        assert (a, b) == (ta, tb)
        reverse = g_column(*target)
        assert reverse.get(source, 0) == coefficient.conjugate()
        for final, second in reverse.items():
            twice[final] = twice.get(final, 0) + coefficient * second
    assert {out: coefficient for out, coefficient in twice.items() if coefficient} == {source: 1+0j}


def main():
    machine = WORK.machine
    initial_allowed = 0
    seed_changed = 0
    for index in range(WORK_DIM):
        assert encode(*decode(index)) == index
        assert WORK.f(WORK.f(index)) == index
        assert WORK.seed(WORK.seed(index)) == index
        assert f_column(index) == {WORK.f(index): 1+0j}
        assert seed_column(index) == {WORK.seed(index): 1+0j}
        a, record, control, phase = decode(index)
        b, target_record, _, _ = decode(WORK.f(index))
        assert a.static == b.static and record == target_record
        seeded_symbol, seeded_record, _, _ = decode(WORK.seed(index))
        assert a == seeded_symbol and record == seeded_record
        if a.static != "LEFT":
            assert WORK.seed(index) == index
        initial_allowed += int(WORK.initial_allowed(index))
        seed_changed += int(WORK.seed(index) != index)
    assert initial_allowed == 10 and seed_changed == 64

    # Every control, phase, orientation, quantum basis input, and every
    # left neutral symbol. Three right neutral records include nonblank
    # boundary/frontier/program/cursor combinations. Neutral records are
    # spectators in G; its unitary head/oracle block is independent of them.
    anchors = (machine.symbols[0], machine.symbols[31], machine.symbols[-1])
    one_head_cases = 0
    for q, phase, on_left, a, b, ra, rb in product(
        machine.controls, range(2), (False, True), machine.symbols, anchors, range(2), range(2)
    ):
        left = encode(a, ra, q.name if on_left else None, phase if on_left else 0)
        right = encode(b, rb, None if on_left else q.name, 0 if on_left else phase)
        check_g_pair(left, right)
        one_head_cases += 1

    zero_head_cases = 0
    for a, b, ra, rb in product(machine.symbols, machine.symbols, range(2), range(2)):
        pair = encode(a, ra), encode(b, rb)
        assert g_column(*pair) == {pair: 1+0j}
        zero_head_cases += 1
    two_head_cases = 0
    for qa, qb, pa, pb, ra, rb in product(machine.controls, machine.controls,
                                        range(2), range(2), range(2), range(2)):
        pair = encode(anchors[1], ra, qa.name, pa), encode(anchors[-1], rb, qb.name, pb)
        assert g_column(*pair) == {pair: 1+0j}
        two_head_cases += 1

    # Execute one complete six-site controller cycle by literal sparse
    # local basis columns. Start with all quantum records one, which gives
    # a nontrivial coherent superposition under the first CV gate. This is
    # separate from the vector-only and collision-only simulators.
    n = 6
    configuration = machine.initialize(n)
    initial = tuple(encode(a, 1, configuration.control if j == 0 else None)
                    for j, a in enumerate(configuration.tape))
    state = {initial: 1+0j}
    maximum_support = 1
    local_gates = 0
    for _ in range(machine.cycle_steps(n)):
        machine.step(configuration)
        for left in range(n - 1):
            out = {}
            for basis, amplitude in state.items():
                for (a, b), coefficient in g_column(basis[left], basis[left+1]).items():
                    target = basis[:left] + (a, b) + basis[left+2:]
                    out[target] = out.get(target, 0) + amplitude * coefficient
            state = {basis: amplitude for basis, amplitude in out.items() if amplitude}
            maximum_support = max(maximum_support, len(state))
            local_gates += 1
        for site in reversed(range(n)):
            out = {}
            for basis, amplitude in state.items():
                target = basis[:site] + (WORK.f(basis[site]),) + basis[site+1:]
                out[target] = out.get(target, 0) + amplitude
            state = {basis: amplitude for basis, amplitude in out.items() if amplitude}
            local_gates += 1
        for basis in state:
            heads = []
            for j, index in enumerate(basis):
                symbol, record, control, phase = decode(index)
                assert symbol == configuration.tape[j]
                if control is not None:
                    heads.append((j, control, phase))
            assert heads == [(configuration.position, configuration.control, 0)]
    assert state == {initial: 1+0j}
    assert maximum_support > 1

    result = {
        "all_passed": True,
        "work_dimension": WORK_DIM,
        "two_site_matrix_dimension": WORK_DIM ** 2,
        "all_F_seed_encode_decode_basis_states_checked": WORK_DIM,
        "initially_allowed_local_states": initial_allowed,
        "seed_changed_local_states": seed_changed,
        "one_head_G_columns_and_involution_checks": one_head_cases,
        "all_zero_head_G_columns_checked": zero_head_cases,
        "two_head_control_phase_quantum_cases_checked": two_head_cases,
        "neutral_symbol_note": "All left symbols and three right anchors for one-head blocks; neutral factors are spectators by construction",
        "literal_sparse_local_gate_cycle": {
            "sites": n,
            "local_gate_applications": local_gates,
            "maximum_quantum_basis_support": maximum_support,
            "exact_final_return": True,
        },
        "scope": "Explicit full-space F/G/LEFT-seed gates; no global clock-sector theorem supplied here",
    }
    Path(__file__).with_name("local_gates_results.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
