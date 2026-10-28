#!/usr/bin/env python3
"""Whole-space gates, promised-helper basis, schedules and cleanup checks."""

from collections import Counter
from dataclasses import replace
from itertools import product
import json
from pathlib import Path
import random

from walk import WRITE, ParameterWalk, Site, species
from parameter_oracle import extraction_word, inverse_op


def dimensions(p):
    u = p*p+p+1
    q, nlogical = 2*u, (p+1)*u
    return q, max((q-1).bit_length(), 2*(nlogical-1).bit_length())


def all_sites(machine):
    for mark, digit, end in machine.symbols:
        for sigma in machine.alphabet:
            yield Site(mark, sigma, digit, end)
            for control in machine.controls:
                for phase in range(2):
                    yield Site(mark, sigma, digit, end, control, phase)


def charge(sigma_pair):
    return Counter(s for sigma in sigma_pair if (s := species(sigma)) is not None)


def local_checks(p, q, m):
    machine = ParameterWalk(p, q, m)
    assert sorted(machine.collision) == list(range(96))
    assert all(machine.inverse[b] == a for a, b in enumerate(machine.collision))
    for a, b in enumerate(machine.collision):
        assert machine.decode(a)[1] == machine.decode(b)[1]
    count = 0
    for a in all_sites(machine):
        b = machine.f(a)
        assert machine.f(b) == a
        assert a.mark == b.mark
        assert species(a.sigma) == species(b.sigma)
        count += 1
    assert count == machine.dimension
    return {"p": p, "q": q, "m": m, "Sigma_dimension": len(machine.alphabet),
            "complete_F_columns_checked": count,
            "full_collision_inverse_and_all_species_charges": True}


def movement_checks():
    # G ignores neutral marks/digits/end; exhaust every Sigma pair and
    # every one-head selector, including its inactive orientations/phases.
    machine = ParameterWalk(2, 2, 2)
    count = 0
    for sigma_a, sigma_b in product(machine.alphabet, repeat=2):
        for q, phase, on_left in product(machine.controls, range(2), (False, True)):
            a = Site("LEFT", sigma_a, 1, 0, q if on_left else None, phase if on_left else 0)
            b = Site("RIGHT", sigma_b, 0, 1, q if not on_left else None, phase if not on_left else 0)
            c, d = machine.g(a, b)
            assert machine.g(c, d) == (a, b)
            assert (c.mark, c.digit, c.end) == (a.mark, a.digit, a.end)
            assert (d.mark, d.digit, d.end) == (b.mark, b.digit, b.end)
            assert charge((c.sigma, d.sigma)) == charge((a.sigma, b.sigma))
            count += 1
    # All head selectors with zero or two heads are identity, independently
    # of payload by the explicit top-level branch in g.
    tags = [(None, 0)] + [(q, ph) for q in machine.controls for ph in range(2)]
    inactive = 0
    for (qa, pa), (qb, pb) in product(tags, repeat=2):
        if (qa is None) != (qb is None):
            continue
        a = Site("BULK", ("T", 1, 1, 1), 1, 1, qa, pa)
        b = Site("RIGHT", ("H", 0, 1), 0, 1, qb, pb)
        assert machine.g(a, b) == (a, b)
        inactive += 1
    return {"all_Sigma_pairs_times_one_head_selectors": count,
            "zero_or_two_head_selectors": inactive,
            "neutral_factors_unchanged_by_definition_and_checked": True}


def schedule(tape):
    return tuple((j, a.control, a.phase) for j, a in enumerate(tape) if a.control is not None)


def run_case(machine, payload, micros=False):
    initial = machine.initialize(payload)
    tape, trace = initial, [schedule(initial)]
    midpoint = None
    events = []
    for tick in range(1, 4*machine.m+1):
        at, control, _ = schedule(tape)[0]
        destination = at + machine.direction[control]
        movement = machine.oracle[control]
        if movement == "PUSH":
            events.extend([("SHIFT", at), ("SWAP", min(at, destination))])
        elif movement == "POP":
            events.extend([("SWAP", min(at, destination)), ("SHIFT_DAG", destination)])
        if machine.write_forward(replace(tape[destination], control=control)):
            events.extend((op, destination) for op in WRITE)
        if micros:
            actual = tape
            for actual, _ in machine.micro_sweep(tape):
                assert sum(a.control is not None for a in actual) == 1
            next_tape = machine.step(tape)
            assert actual == next_tape
        else:
            next_tape = machine.step(tape)
        tape = next_tape
        if machine.write_reverse(tape[destination]):
            events.extend((op, destination) for op in reversed(WRITE))
        trace.append(schedule(tape))
        if tick == 2*machine.m:
            midpoint = tape
            assert events == extraction_word(len(tape), machine.m)
            assert tuple(a.sigma for a in tape) == tuple(payload)
            assert tape[-1].control == "R:BACK"
            for j, a in enumerate(tape):
                digit_index = len(tape)-1-j
                assert a.digit == ((machine.p >> digit_index) & 1 if digit_index < machine.m else 0)
                assert a.end == int(digit_index == machine.m-1)
            assert sum(a.digit << j for j, a in enumerate(reversed(tape))) == machine.p
    assert midpoint is not None
    assert tape == initial
    forward_word = extraction_word(len(tape), machine.m)
    inverse_word = [(inverse_op(op), at) for op, at in reversed(forward_word)]
    assert events == forward_word + inverse_word
    return trace


def exhaustive_payloads():
    machine = ParameterWalk(2, 2, 2)
    reference = None
    count = 0
    for sigma in product(machine.alphabet, repeat=2):
        trace = run_case(machine, sigma + (("H", 0, 0),), micros=True)
        if reference is None:
            reference = trace
        assert trace == reference
        count += 1
    return {"q": 2, "m": 2, "n": 3,
            "all_promised_helper_Sigma_basis_cases": count,
            "all_micro_sweeps_match_direct_step": True,
            "primitive_event_word_matches_original_and_exact_inverse": True,
            "same_schedule_for_all_untrusted_payloads": True}


def larger_cases():
    rng = random.Random(28102026)
    cases, sweeps = 0, 0
    records = []
    for p in (2, 3, 5, 7, 11, 17, 31, 127):
        q, m = dimensions(p)
        machine = ParameterWalk(p, q, m)
        for extra in (0, 1, 7):
            n = m + extra
            reference = None
            for trial in range(8):
                payload = tuple(rng.choice(machine.alphabet) for _ in range(n-1)) + (("H", 0, 0),)
                trace = run_case(machine, payload, micros=trial == 0)
                if reference is None:
                    reference = trace
                assert trace == reference
                cases += 1
                sweeps += 4*m
        records.append({"p": p, "q": q, "m": m, "midpoint_sweeps": 2*m,
                        "full_cycle_sweeps": 4*m})
    # Separately cover the small m=2 boundary and non-prime bit strings.
    for m in range(2, 12):
        for p in (0, 1, 2**m-1):
            machine = ParameterWalk(p, 2, m)
            run_case(machine, [("T", m-1, 1, 1)]*(m-1)+[("H", 0, 0)], micros=True)
            cases += 1
            sweeps += 4*m
    return {"payload_cycles": cases, "sweeps": sweeps, "family_examples": records,
            "loaded_T_and_malformed_H_included": True,
            "n_equals_m_and_m_equals_2_boundaries_checked": True}


def main():
    machine = ParameterWalk(2, 2, 2)
    manifest = machine.manifest()
    path = Path(__file__).with_name("transition_table.json")
    assert json.loads(path.read_text()) == json.loads(json.dumps(manifest))
    report = {
        "status": "passed",
        "scope": "Finite parameter reader with initialized helper and exact reverse cleanup; no composed generator or counterexample",
        "neutral_collision_entries": len(machine.collision),
        "prescribed_neutral_entries": len(machine.partial),
        "neutral_dimension_B": machine.B,
        "controls_K": machine.K,
        "local_factor_times_Sigma": machine.local_factor,
        "on_site": [local_checks(2, 2, 2)] + [local_checks(p, *dimensions(p)) for p in (2, 3)],
        "two_site_movement": movement_checks(),
        "complete_payload_basis": exhaustive_payloads(),
        "larger_adversarial_cases": larger_cases(),
    }
    Path(__file__).with_name("results.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
