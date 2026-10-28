"""End-to-end checks of the explicit microscopic bond assembly fixture."""
from itertools import product
from pathlib import Path
import argparse
import json
import math
import time
import numpy as np
from bond import MicroscopicBond, one_fermion_core, plus


def dense_from_columns(dimension, column):
    out = np.zeros((dimension, dimension), complex)
    for j in range(dimension):
        for i, value in column(j).items():
            out[i, j] = value
    return out


def core_checks(bond):
    core = bond.core
    k = core.length
    local = dense_from_columns(4, core.bulk)
    h = np.zeros((2**k, 2**k), complex)
    for j in range(k-1):
        h += np.kron(np.kron(np.eye(2**j), local), np.eye(2**(k-j-2)))
    frame = np.zeros((2**k, 2))
    frame[0, 0] = 1
    for j in range(k):
        frame[2**j, 1] = 1 / math.sqrt(k)
    e = np.linalg.eigvalsh(h)
    assert np.linalg.norm(h @ frame) < 1e-12
    assert sum(abs(e) < 1e-12) == 2 and e[2] > 0
    raising = np.array([[0, 0], [1, 0]], complex)
    left = np.kron(raising, np.eye(2**(k-1)))
    right = np.kron(np.eye(2**(k-1)), raising)
    port_error = max(np.linalg.norm(frame.T @ a @ frame - raising/math.sqrt(k))
                     for a in (left, right))
    assert port_error < 1e-12
    return frame, dict(dimension=2**k, code_dimension=2,
                       gap=float(e[2]), port_error=float(port_error))


def wrapper_states_frame(bond, core_frame, holes_count):
    k, length = bond.k, bond.k + holes_count
    states, rows = [], []
    charge_states = bond.charge_states[holes_count-1]
    chis = bond.chis[holes_count-1]
    tags = [0] + [1]*(k-2) + [2]
    for holes, amplitude in zip(charge_states, chis):
        occupied = [j for j in range(length) if j not in holes]
        for cbits in product(range(2), repeat=k):
            cindex = sum(bit*2**(k-j-1) for j, bit in enumerate(cbits))
            for fbits in product(range(2), repeat=k):
                filt = np.prod([plus(j*bond.theta)[bit] for j, bit in enumerate(fbits)])
                letters = iter(bond.wrap_encode(t, c, f) for t, c, f in zip(tags, cbits, fbits))
                state = tuple((0 if j < occupied[0] else 2 if j > occupied[-1] else 1)
                              if j in holes else next(letters) for j in range(length))
                states.append(state)
                rows.append(amplitude*filt*core_frame[cindex])
    return states, np.asarray(rows)


def wrapper_checks(bond, core_frame):
    output = []
    for holes in (1, 2):
        states, frame = wrapper_states_frame(bond, core_frame, holes)
        index = {s: i for i, s in enumerate(states)}
        h = np.zeros((len(states), len(states)), complex)
        for col, word in enumerate(states):
            for j in range(len(word)-1):
                for pair, value in bond.wrapper_column(*word[j:j+2]).items():
                    target = word[:j] + pair + word[j+2:]
                    h[index[target], col] += value
        energy = bond.charge_energies[holes-1][0] - 4*bond.U
        residual = float(np.linalg.norm(h @ frame - energy*frame))
        orthogonality = float(np.linalg.norm(frame.conj().T @ frame - np.eye(2)))
        eigenvalues = np.linalg.eigvalsh(h)
        # The filter is stiff and the physical gap is tiny; record absolute
        # residuals and finite-block spectra rather than claiming exact arithmetic.
        assert residual < 1e-9 and orthogonality < 1e-12
        assert max(abs(eigenvalues[:2] - energy)) < 1e-9
        assert eigenvalues[2] - energy > 1e-9
        port_errors = []
        for end, tag in ((0, 0), (-1, 2)):
            op_frame = np.zeros_like(frame)
            for col, state in enumerate(states):
                t, c, f = bond.occupied(state[end]) if state[end] >= 3 else (-1, -1, -1)
                if t == tag and c == 0:
                    target = list(state)
                    target[end] = bond.wrap_encode(t, 1, f)
                    op_frame[index[tuple(target)]] += frame[col]
            projected = frame.conj().T @ op_frame
            want = np.array([[0, 0], [1, 0]]) * bond.r[holes-1] / math.sqrt(bond.k)
            port_errors.append(float(np.linalg.norm(projected-want)))
        assert max(port_errors) < 1e-12
        output.append(dict(length=bond.k+holes, sector_dimension=len(states),
                           ground_energy=float(energy), code_residual=residual,
                           gap_in_selected_occupation_sector=float(eigenvalues[2]-energy),
                           port_error=max(port_errors)))
    return output


def equalizer_checks(bond):
    errors = []
    for r, length in zip(bond.r, (bond.b, bond.b+1)):
        u = plus(-(length-1)*bond.phi)
        errors.extend([abs(r*np.vdot(u, bond.EL @ u)-1),
                       abs(r*np.vdot(u.conj(), bond.ER @ u.conj())-1)])
    assert max(errors) < 1e-12
    # Product-state overlap of the local helix conversion, directly over sites.
    before_lengths, after_lengths = (bond.b, bond.b+1), (bond.b+1, bond.b)
    def phases(lengths):
        return [(j*bond.phi, (j-ell+1)*bond.phi) for ell in lengths for j in range(ell)]
    overlap = np.exp(-1j*bond.b*bond.phi)
    for j, (old, new) in enumerate(zip(phases(before_lengths), phases(after_lengths))):
        for theta_old, theta_new in zip(old, new):
            if j == bond.b:
                theta_old += bond.b*bond.phi
            overlap *= np.vdot(plus(theta_new), plus(theta_old))
    error = abs(overlap - bond.beta_helix)
    assert error < 1e-12
    return dict(compensated_port_error=float(max(errors)), helix_transport_error=float(error))


def full_bond_checks(bond, export):
    d, total = bond.dimension, bond.dimension**2
    charges = np.array([bond.site_charge(i) for i in range(d)])
    markers = np.array([bond.decode(i)[0] == 0 for i in range(d)], dtype=int)
    rows, columns, values = [], [], []
    max_decomposition_error = 0.0
    for col in range(total):
        out = bond.column(col)
        for row, val in out.items():
            rows.append(row)
            columns.append(col)
            values.append(val)
        if col % 127 == 0:
            split = {}
            for component in ("base", "hopping", "move", "field"):
                for row, val in bond.column(col, component).items():
                    split[row] = split.get(row, 0) + val
            max_decomposition_error = max(max_decomposition_error,
                max((abs(out.get(row, 0)-split.get(row, 0)) for row in out.keys() | split.keys()), default=0))
    rows, columns, values = np.array(rows, dtype=np.int64), np.array(columns, dtype=np.int64), np.array(values)
    keys = rows*total + columns
    adjoint_keys = columns*total + rows
    order, adjoint_order = np.argsort(keys), np.argsort(adjoint_keys)
    assert np.array_equal(keys[order], adjoint_keys[adjoint_order])
    hermiticity_error = float(np.max(np.abs(values[order] - values[adjoint_order].conj())))
    qout = charges[rows//d] + charges[rows % d]
    qin = charges[columns//d] + charges[columns % d]
    species_error = int(np.max(np.abs(qout-qin)))
    marker_error = int(np.max(np.abs(markers[rows//d] + markers[rows % d]
                                     - markers[columns//d] - markers[columns % d])))
    assert hermiticity_error < 1e-10 and species_error == marker_error == 0
    assert max_decomposition_error < 1e-10
    one_norm = float(np.max(np.bincount(columns, weights=np.abs(values), minlength=total)))
    if export:
        np.savez_compressed(export, row=rows, column=columns, value=values,
                            shape=np.array([total, total]), local_dimension=d)
    return dict(local_dimension=d, matrix_shape=[total, total],
                every_column_evaluated=True, stored_entries=len(values),
                hermiticity_error=hermiticity_error, species_charge_error=species_error,
                marker_number_error=marker_error, decomposition_error=max_decomposition_error,
                operator_norm_upper_bound=one_norm,
                data_file=str(export.name) if export else None)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--export", action="store_true", help="Save the entire sparse two-site fixture matrix")
    args = parser.parse_args()
    started = time.perf_counter()
    directory = Path(__file__).resolve().parent
    bond = MicroscopicBond(one_fermion_core(), s=1e-7)
    frame, core = core_checks(bond)
    result = dict(manifest=bond.manifest(), core=core,
                  selected_wrapper_sectors=wrapper_checks(bond, frame),
                  equalizer=equalizer_checks(bond),
                  full_bond=full_bond_checks(bond, directory / "fixture_bond.npz" if args.export else None))
    result["runtime_seconds"] = time.perf_counter()-started
    (directory / "results.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
