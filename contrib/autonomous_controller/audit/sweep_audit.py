#!/usr/bin/env python3
"""Independent C8-clock refinement and journal-machine adversarial audit.

The readout input is supplied, so this is not an initialization proof.  Every
actual C8 path vertex, including counter/initialization idles, is inspected.
"""

from __future__ import annotations

import importlib.util
from itertools import product
import json
from pathlib import Path
import random
import sys
import types

import numpy as np


ROOT = Path(__file__).resolve().parents[3]


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    sys.modules[name] = result
    spec.loader.exec_module(result)
    return result


# Use the existing literal clock rules/path without importing its unrelated
# spectral-check SciPy dependency. No definitions below full_check are needed.
clock = types.ModuleType('audit_outer_clock')
clock_source = (ROOT / 'codex/work_C8/sweeping_clock_checks.py').read_text()
clock_source = clock_source[clock_source.index('ALPHABET='):clock_source.index('\ndef full_check(')]
exec(compile(clock_source, 'existing_C8_clock_rules', 'exec'), clock.__dict__)
readout = module('audit_readout', ROOT / 'contrib/autonomous_controller/readout/transport.py')
ram = module('audit_generator_machine', ROOT / 'contrib/autonomous_controller/generator/machine.py')
sys.path.insert(0, str(ROOT / 'contrib/autonomous_controller/bridge'))
integrated = module('audit_integrated_machine', ROOT / 'contrib/autonomous_controller/bridge/integrated_readout.py')
sys.path.pop(0)


def active(tape, endpoint):
    x = tape[endpoint]
    return bool(x.ready and x.control is None)


def refined_history(n, left, right):
    machine = readout.Readout()
    tape = machine.initialize(n, left, right)
    words, steps = clock.path(n)
    result = [tape]
    for label, position in steps:
        work = list(tape)
        if label == 'out':
            work[position], work[position + 1] = machine.g(
                work[position], work[position + 1])
        elif label == 'back':
            work[position + 1] = machine.f(work[position + 1])
        elif label == 'finish':
            work[0] = machine.f(work[0])
        # The initial clock sweep is identity on the supplied work input.
        # All carry/turn/return/begin/bounce/activation vertices are idle here.
        tape = tuple(work)
        result.append(tape)
    assert len(result) == len(words)
    return result


def full_clock_checks():
    reports = []
    for n in range(4, 11):
        geometries = [(1, n - 2)]
        if n >= 6:
            geometries.append((2, n - 3))
        for left, right in geometries:
            history = refined_history(n, left, right)
            counts = [0, 0]
            for tape in history:
                assert sum(x.control is not None for x in tape) == 1
                for which, endpoint, source in ((0, 0, left), (1, n - 1, right)):
                    if active(tape, endpoint):
                        assert tape[endpoint].record == source
                        counts[which] += 1
            assert min(counts) > 0
            reports.append({'sites': n, 'output_positions': [left, right],
                            'actual_clock_vertices': len(history),
                            'full_sweeps': 2 ** n, 'active_vertices': counts,
                            'exact_weights': [f'{c}/{len(history)}' for c in counts]})
    return reports


def move_code(code, records):
    """Independent integer-permutation implementation; bit0 is physical site0."""
    n = len(records)
    out = np.zeros_like(code)
    for old in range(len(code)):
        new = sum(((old >> records[j]) & 1) << j for j in range(n))
        out[new] = code[old]
    return out


def apply_creation(code, site):
    out = np.zeros_like(code)
    for index in range(len(code)):
        if not (index >> site) & 1:
            out[index | (1 << site)] = code[index]
    return out


def coherent_clock_compression():
    n, left, right = 6, 1, 4
    rng = np.random.default_rng(17242240)
    raw = rng.normal(size=(2 ** n, 5)) + 1j * rng.normal(size=(2 ** n, 5))
    code, _ = np.linalg.qr(raw)
    history = refined_history(n, left, right)
    reports = []
    for endpoint, source in ((0, left), (n - 1, right)):
        target = code.conj().T @ apply_creation(code, source)
        total = np.zeros_like(target)
        largest, active_count = 0.0, 0
        for tape in history:
            if not active(tape, endpoint):
                continue
            evolved = move_code(code, [a.record for a in tape])
            value = evolved.conj().T @ apply_creation(evolved, endpoint)
            largest = max(largest, float(np.max(np.abs(value - target))))
            total += value
            active_count += 1
        weight = active_count / len(history)
        error = float(np.max(np.abs(total / len(history) - weight * target)))
        assert error < 1e-13 and largest < 1e-13
        reports.append({'endpoint': endpoint, 'arbitrary_entangled_code_dimension': 5,
                        'active_vertices': active_count, 'clock_vertices': len(history),
                        'compression_error': error, 'worst_active_vertex_error': largest})
    return reports


def word_machine_rollback():
    # Deliberately combine nested calls, returned arrays, aliasing, overwritten
    # cells/registers, branches, negative values, and output before HALT.
    source = '''
def helper(a, x):
    old = a[0]
    a[0] = x
    a.append(old)
    return a

def main(p):
    a = [p, -3]
    b = a
    for j in range(3):
        a = helper(a, j - p)
        emit(a[0], a[len(a) - 1])
    if p % 2 == 0:
        b[1] = b[0] * p
    else:
        b[1] = b[0] + p
    for j in range(2, -1, -1):
        emit(j, b[j])
    return b
'''
    program = ram.Compiler(source).payload
    cases, steps = 0, 0
    for p in range(8):
        vm = ram.Machine(program, 'main', [p], journal=True)
        snapshots = [vm.state()]
        while not vm.halted:
            vm.step()
            snapshots.append(vm.state())
        steps += vm.steps
        while vm.journal:
            vm.undo()
            snapshots.pop()
            assert vm.state() == snapshots[-1]
        assert len(snapshots) == 1 and vm.steps == 0
        cases += 1
    return {'adversarial_program_runs': cases, 'every_instruction_inverse_checks': steps,
            'covers': ['nested calls/RET', 'array aliasing', 'NEW/APPEND/SET',
                       'register overwrite', 'both branch outcomes',
                       'negative integers', 'EMIT', 'halt reversal'],
            'scope': 'Journalled word-machine semantic states, not local machine transitions'}


def initialization_prefix():
    # Four head-overlay types suffice to exercise the exact full-space rule:
    # ordinary <-> chosen fresh start, while every other head type is fixed.
    seed = (1, 0, 2, 3)
    assert all(seed[seed[x]] == x for x in range(4))
    reports = []
    for n in range(2, 7):
        accepted = []
        for initial in product(range(4), repeat=n):
            work = list(initial)
            passed = True
            for i in range(n):
                # The ordinary-head pin is checked while I is at site i.
                passed &= work[i] == 0
                # Decorate IZ -> LI at its departed site, controlled LEFT.
                if i == 0:
                    work[0] = seed[work[0]]
            if passed:
                accepted.append(initial)
                assert work == [1] + [0] * (n - 1)
        assert accepted == [(0,) * n]
        reports.append({'sites': n, 'complete_initial_head_words': 4 ** n,
                        'accepted_words': 1, 'final_head_position': 0})
    return reports


def initialized_integrated_clock():
    """Start no-head/blank; I sweep creates LEFT head; actual C8 then runs it."""
    machine = integrated.IntegratedReadout()
    reports = []
    for n in range(4, 11):
        reference = machine.initialize(n)
        tape = list(reference.tape)
        position, control, phase = None, None, 0
        words, steps = clock.path(n)
        sweeps = 0
        active_counts = [0, 0]
        oracle_events, reference_events = [], []
        records = list(range(n))
        for label, at in steps:
            if label == 'initialize' and at == 0:
                # This is the LEFT-conditioned B on the departing I site.
                assert control is None
                position, control, phase = 0, 'F:ADVANCE', 0
            elif label == 'out':
                q = machine.control_by_name[control]
                if q.direction == 1:
                    forward = phase == 0 and position == at
                    reverse = phase == 1 and position == at + 1
                    target = at + 1 if forward else at
                else:
                    forward = phase == 0 and position == at + 1
                    reverse = phase == 1 and position == at
                    target = at if forward else at + 1
                if forward or reverse:
                    if q.oracle is not None:
                        oracle_events.append((q.oracle, at, q.adjoint ^ reverse))
                        if q.oracle == 'SWAP':
                            records[at], records[at + 1] = records[at + 1], records[at]
                    position, phase = target, 1 - phase
            elif label in ('back', 'finish'):
                site = at + 1 if label == 'back' else 0
                if position == site:
                    assert phase == 1
                    source = machine.encode(control, tape[site])
                    assert source in machine.partial
                    control, tape[site] = machine.decode(machine.collision[source])
                    phase = 0
                if label == 'finish':
                    event = machine.step(reference)
                    if event is not None:
                        reference_events.append(event)
                    assert (position, control, tuple(tape)) == reference.snapshot()
                    assert phase == 0
                    sweeps += 1
            # Every actual path vertex, including all idle clock operations.
            for which, endpoint, source in ((0, 0, 1), (1, n - 1, n - 2)):
                if tape[endpoint].cursor and position != endpoint:
                    assert records[endpoint] == source
                    assert control.startswith('RO:')
                    active_counts[which] += 1
        assert sweeps == 2 ** n
        assert oracle_events == reference_events
        if n >= 9:
            assert 2 ** n >= machine.cycle_steps(n) and min(active_counts) > 0
        reports.append({'sites': n, 'clock_vertices': len(words),
                        'sweeps_checked_against_macro_machine': sweeps,
                        'actual_circuit_oracle_events': len(oracle_events),
                        'complete_work_cycles': sweeps // machine.cycle_steps(n),
                        'active_endpoint_vertices': active_counts,
                        'active_weights': [f'{x}/{len(words)}' for x in active_counts]})
    return reports


def hamiltonian_columns():
    core = module('audit_initialized_core', ROOT / 'contrib/autonomous_controller/core/core.py')
    rng = random.Random(194400)
    assert core.LOCAL_DIM == 194400 and core.WORK_DIM == 12960
    counts = {mark: 0 for mark in ('LEFT', 'BULK', 'RIGHT', 'OUT_L', 'OUT_R')}
    for work in range(core.WORK_DIM):
        if core.initial_allowed(work):
            counts[core.work_decode(work)[0].static] += 1
    assert set(counts.values()) == {2}
    # Local checks intentionally sample, rather than claiming exhaustive
    # enumeration of the 37,791,360,000 full two-site input columns.
    columns, entries, complex_entries = 0, 0, 0
    one_indices = [core.encode(c, rng.randrange(core.WORK_DIM))
                   for c in core.ALPHABET for _ in range(40)]
    for function in (core.onsite_column, core.left_column, core.right_column):
        for source in one_indices:
            for target, value in function(source).items():
                reverse = function(target).get(source, 0)
                assert abs(value - reverse.conjugate()) < 1e-14
                entries += 1
            columns += 1
    for before, after, label in core.RULES:
        for clocks in (before, after):
            for _ in range(40):
                source = tuple(core.encode(c, rng.randrange(core.WORK_DIM)) for c in clocks)
                column = core.bond_column(*source)
                for target, value in column.items():
                    reverse = core.bond_column(*target).get(source, 0)
                    assert abs(value - reverse.conjugate()) < 1e-14
                    entries += 1
                    complex_entries += int(value.imag != 0)
                columns += 1
    # Random selection overwhelmingly finds two-head identity blocks, so add
    # every nontrivial one-head oracle control/direction/phase explicitly.
    for q in core.WORK.machine.controls:
        for phase in (0, 1):
            for head_left in (False, True):
                neutral = rng.choice(core.WORK.machine.symbols)
                head = core.WORK.encode(neutral, 1, q.name, phase)
                ordinary = core.WORK.encode(rng.choice(core.WORK.machine.symbols), 0)
                works = (head, ordinary) if head_left else (ordinary, head)
                for clocks in (('D0', '0'), ('A0', 'D0')):
                    source = tuple(core.encode(c, w) for c, w in zip(clocks, works))
                    for target, value in core.bond_column(*source).items():
                        reverse = core.bond_column(*target).get(source, 0)
                        assert abs(value - reverse.conjugate()) < 1e-14
                        entries += 1
                        complex_entries += int(value.imag != 0)
                    columns += 1
    assert complex_entries > 0
    return {'local_dimension': core.LOCAL_DIM, 'work_dimension': core.WORK_DIM,
            'complete_work_initial_pin_columns': core.WORK_DIM,
            'allowed_qubit_states_per_static_mark': counts,
            'sampled_Hamiltonian_columns': columns,
            'Hermitian_reverse_entries_checked': entries,
            'nonreal_entries_checked': complex_entries,
            'all_two_site_columns_enumerated': False}


def main():
    report = {'status': 'passed',
              'scope': 'Independent supplied-readout clock refinement and word-machine audit',
              'full_outer_clock': full_clock_checks(),
              'coherent_compression': coherent_clock_compression(),
              'journal_machine': word_machine_rollback(),
              'left_seed_initialization_prefix': initialization_prefix(),
              'initialized_integrated_C8_clock': initialized_integrated_clock(),
              'full_Hamiltonian_interface': hamiltonian_columns()}
    Path(__file__).with_name('sweep_audit_results.json').write_text(
        json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
