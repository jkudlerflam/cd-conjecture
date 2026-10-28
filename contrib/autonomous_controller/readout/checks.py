"""Whole-space involutions and every-vertex endpoint-operator checks."""
from itertools import product
import json
from pathlib import Path
import numpy as np
from transport import Readout, Site


def active(tape, endpoint):
    return bool(tape[endpoint].ready and tape[endpoint].control is None)


def check_vertex(tape, left, right):
    assert sum(a.control is not None for a in tape) == 1
    for endpoint, output in ((0, left), (len(tape)-1, right)):
        if active(tape, endpoint):
            assert tape[endpoint].record == output


def orbit(machine, n, left, right):
    initial = tape = machine.initialize(n, left, right)
    states = []
    for _ in range(machine.cycle_steps(n, left, right)):
        for updated, operation in machine.micro_sweep(tape):
            states.append(tape)
            check_vertex(tape, left, right)
            tape = updated
    assert tape == initial
    assert all(any(active(t, e) for t in states) for e in (0, n-1))
    return states


def raised(v, site, n):
    result = np.zeros_like(v)
    for index in range(2**n):
        if not (index >> site) & 1:
            result[index | (1 << site)] = v[index]
    return result


def permuted(v, records):
    n = len(records)
    result = np.zeros_like(v)
    for old in range(2**n):
        new = sum(((old >> source) & 1) << j for j, source in enumerate(records))
        result[new] = v[old]
    return result


def main():
    m = Readout()
    assert len(m.partial) == 20
    table_file = Path(__file__).with_name("transition_table.json")
    assert json.loads(table_file.read_text()) == json.loads(json.dumps(m.manifest()))
    local = []
    for mark, ready in m.symbols:
        for record in range(2):
            local.append(Site(mark, ready, record))
            local.extend(Site(mark, ready, record, q, phase)
                         for q in m.controls for phase in range(2))
    for a in local:
        assert m.f(m.f(a)) == a
        assert (m.f(a).mark, m.f(a).record) == (a.mark, a.record)
    for a, b in product(local, repeat=2):
        x, y = m.g(a, b)
        assert m.g(x, y) == (a, b)
        assert (x.mark, y.mark, x.ready, y.ready) == (a.mark, b.mark, a.ready, b.ready)
        assert sorted((x.record, y.record)) == sorted((a.record, b.record))
    cases = vertices = 0
    for n in range(4, 16):
        for left in range(1, n-2):
            for right in range(left+1, n-1):
                states = orbit(m, n, left, right)
                cases += 1
                vertices += len(states)

    # Arbitrary coherent input code: not restricted to computational-basis
    # payloads or product auxiliaries. Testing the pullbacks checks the
    # operator identity and the convention used by the record permutation.
    n, left, right = 6, 1, 4
    rng = np.random.default_rng(20261028)
    v, _ = np.linalg.qr(rng.normal(size=(2**n, 4)) + 1j*rng.normal(size=(2**n, 4)))
    states = orbit(m, n, left, right)
    projections = []
    for endpoint, output in ((0, left), (n-1, right)):
        actual = np.zeros((4, 4), complex)
        unflagged = np.zeros_like(actual)
        count = 0
        target = v.conj().T @ raised(v, output, n)
        worst = 0.0
        for tape in states:
            vt = permuted(v, [a.record for a in tape])
            matrix = vt.conj().T @ raised(vt, endpoint, n)
            if active(tape, endpoint):
                actual += matrix
                count += 1
                worst = max(worst, float(np.max(np.abs(matrix-target))))
            if tape[endpoint].control is None:
                unflagged += matrix
        weight = count / len(states)
        error = float(np.max(np.abs(actual/len(states)-weight*target)))
        negative_error = float(np.max(np.abs(unflagged/len(states)-weight*target)))
        assert error < 1e-13 and worst < 1e-13 and negative_error > 1e-3
        projections.append(dict(endpoint=endpoint, input_output_site=output,
                                active_vertices=count, sampled_vertices=len(states),
                                uniform_sample_weight=weight, operator_error=error,
                                worst_active_vertex_error=worst,
                                missing_ready_negative_control_error=negative_error))
    report = dict(status="passed", scope="Separate readout primitive; no high-channel preparation or global core-sector proof",
                  neutral_alphabet_B=m.B, controls_K=m.K,
                  head_overlay_factor=(1+2*m.K)*m.B,
                  collision_entries=len(m.table), prescribed_entries=len(m.partial),
                  toy_quantum_dimension=2, one_site_states_checked=len(local),
                  two_site_states_checked=len(local)**2,
                  valid_geometries=cases, elementary_vertices_checked=vertices,
                  cycle_head_moves="4(n-1)", coherent_operator_checks=projections,
                  visibility_note="Weights are for the listed uniform G/F samples, not the full C8 history; safety holds at every vertex and under any additional clock-idle weighting")
    Path(__file__).with_name("results.json").write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
