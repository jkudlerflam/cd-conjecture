"""Independent literal tensor-product audit, L <= 8 and 30 s wall maximum."""

import json
import os
from pathlib import Path
import signal

for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
             "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS"):
    os.environ[name] = "1"

import numpy as np


HERE = Path(__file__).resolve().parent
X = np.array([[0, 1], [1, 0]])
Y = np.array([[0, -1j], [1j, 0]])
Z = np.diag([1, -1])
I = np.eye(2)
b = (X + 1j * Y) / 2  # annihilates a down spin
n = (I - Z) / 2
A, B, HZ = 19 / 32, 13 / 32, 1783 / 1000
TOL = 1.e-9


def tensor(factors):
    product = 1.
    for factor in factors:
        product = np.kron(product, factor)
    return product


def embed(length, site_ops):
    return tensor([site_ops.get(j, I) for j in range(length)])


def bonds(length, periodic):
    return [(j, (j + 1) % length) for j in range(length if periodic else length - 1)]


def literal(length, periodic):
    h = -HZ * sum(embed(length, {j: Z}) for j in range(length))
    for j, k in bonds(length, periodic):
        h = h + A * embed(length, {j: X, k: Y}) + B * embed(length, {j: Y, k: X})
        h = h + embed(length, {j: Z, k: Z})
    return h


def transformed(length, periodic, frame):
    links = bonds(length, periodic)
    # Derive the diagonal Pauli terms independently as occupation projectors.
    h = (len(links) - HZ * length) * np.eye(2 ** length, dtype=complex)
    for j in range(length):
        degree = sum(j in pair for pair in links)
        h += (2 * HZ - 2 * degree) * embed(length, {j: n})
    for j, k in links:
        h += 4 * embed(length, {j: n, k: n})
        pair = embed(length, {j: b, k: b})
        exchange = embed(length, {j: b, k: b.conj().T})
        if frame == "uniform":
            h += pair + pair.conj().T + 1j * (A - B) * (exchange - exchange.conj().T)
        else:
            seam = -1 if k == 0 and length % 4 == 2 else 1
            h += seam * ((-1) ** j * (pair + pair.conj().T)
                         + (A - B) * (exchange + exchange.conj().T))
    return h


def rotation(length, angles):
    return tensor([np.exp(-1j * angle * np.array([1., -1.]) / 2) for angle in angles])


def in_frame(operator, phases):
    return phases.conj()[:, None] * operator * phases[None, :]


def permutation(length, destinations):
    basis = np.arange(2 ** length).reshape((2,) * length)
    # Rows of this sparse permutation would be new states and columns old ones.
    return basis.transpose(np.argsort(destinations)).ravel()


def staggered(length, terms, distance):
    out = np.zeros((2 ** length, 2 ** length), dtype=complex)
    for j in range(length):
        for coefficient, left, right in terms:
            out += (-1) ** j * coefficient * embed(length, {j: left, (j + distance) % length: right})
    return out / np.sqrt(length)


def maxabs(x):
    return float(np.max(np.abs(x)))


def literal_audit():
    out = {}
    for size, periodic in [(6, False), (6, True), (8, True)]:
        h = literal(size, periodic)
        for frame, angles in [("uniform", [np.pi / 4] * size),
                              ("real", np.arange(size) * np.pi / 2)]:
            h_rot = in_frame(h, rotation(size, angles))
            err = maxabs(h_rot - transformed(size, periodic, frame))
            assert err < 1.e-12
            out[f"L{size}_{'PBC' if periodic else 'OBC'}_{frame}_max_error"] = err
    return out


def ring():
    size = 8
    h = literal(size, True)
    # Direct tensor-order indices: old spin j moves to j+1 or -j.
    translate = permutation(size, [(j + 1) % size for j in range(size)])
    invert = permutation(size, [(-j) % size for j in range(size)])
    sign = tensor([np.array([1, (-1) ** j]) for j in range(size)])
    parity = tensor([np.array([1, -1])] * size)
    # A=S K in the literal frame. In the uniform frame Theta=I K.
    assert maxabs(sign[:, None] * h.conj() * sign[None, :] - h) < 1.e-12
    assert maxabs(sign[translate] - parity * sign) == 0
    uniform = rotation(size, [np.pi / 4] * size)
    h_uni = in_frame(h, uniform)
    assert maxabs(h_uni.conj()[np.ix_(invert, invert)] - h_uni) < 1.e-12

    energies, vecs = np.linalg.eigh(h)
    e0 = energies[0]
    assert energies[1] - e0 > TOL
    ground = vecs[:, 0]
    momentum_phase = np.vdot(ground, ground[translate])
    k0 = int(round(np.angle(momentum_phase) * size / (2 * np.pi))) % size
    assert abs(abs(momentum_phase) - 1) < TOL
    odd = np.abs(ground[parity == -1]) ** 2
    p0 = int(round(float(odd.sum())))

    # Fourier-project columns of the raw eigensystem. This works for arbitrary
    # rotations among degenerate eigenstates, including A-linked partners.
    minima = []
    for p in (0, 1):
        sector = parity == (-1) ** p
        shifted = vecs.copy()
        traces = []
        for step in range(size):
            traces.append((vecs[sector].conj() * shifted[sector]).sum(axis=0))
            shifted = shifted[translate]
        weights = np.fft.fft(traces, axis=0).real / size
        for k, row in enumerate(weights):
            minima.append(dict(p=p, k=k, E=float(energies[row > TOL].min())))

    probes = {"kappa2": staggered(size, [(1, X, Y), (-1, Y, X)], 2),
              "B2": staggered(size, [(1, X, X), (-1, Y, Y)], 2),
              "P2": staggered(size, [(1, b, b)], 2),
              "P3": staggered(size, [(1, b, b)], 3)}
    rows = {}
    amps = {}
    for label, op in probes.items():
        u = op @ ground
        amps[label] = vecs.conj().T @ u
        weights = np.abs(amps[label]) ** 2
        active = (energies - e0 > TOL) & (weights > 1.e-11)
        pole = energies[active][0]
        rows[label] = dict(first_gap=float(pole - e0),
                           first_weight=float(weights[np.abs(energies - pole) < TOL].sum()),
                           S_connected=float(np.vdot(u, u).real - abs(np.vdot(ground, u)) ** 2))
        assert maxabs(op[np.ix_(translate, translate)] + op) < 1.e-12

    # Relative phases of the two Hermitian probes are observable in the residue matrix.
    first = np.abs(energies - e0 - rows["kappa2"]["first_gap"]) < TOL
    cross = np.vdot(amps["B2"][first], amps["kappa2"][first])
    ratio = abs(cross) ** 2 / (rows["kappa2"]["first_weight"] * rows["B2"]["first_weight"])
    assert abs(ratio - 1) < 1.e-11 and abs(cross.real) < 1.e-11

    result = dict(L=size, E0=float(e0), p0=p0, k0=k0, sector_minima=minima, q_pi=rows,
                  cross_B2_kappa2=[float(cross.real), float(cross.imag)],
                  normalized_squared_B2_kappa2_coherence=float(ratio))
    reference = json.loads((HERE.parent / "evidence/ring_spectra.json").read_text())["rings"][0]
    assert reference["L"] == size
    assert reference["ground_multiplet"] == [[p0, k0, 0]]
    errors = [abs(e0 - reference["E0"])]
    for row, other in zip(minima, reference["lowest_energy_by_sector"]):
        assert (row["p"], row["k"]) == (other["p"], other["k"])
        errors.append(abs(row["E"] - other["E"]))
    other_pi = reference["momenta"][size // 2]["probes"]
    for name, row in rows.items():
        errors.extend(abs(value - other_pi[name][key]) for key, value in row.items())
    other_cross = min((s for s in reference["six_lowest_neutral_poles"] if s["m"] == size // 2),
                      key=lambda s: abs(s["gap"] - rows["kappa2"]["first_gap"]))
    errors.append(abs(cross - complex(*other_cross["cross_B2_kappa2"])))
    assert max(errors) < 1.e-10
    result["max_discrepancy_against_evidence"] = float(max(errors))
    return result


def main():
    signal.alarm(30)
    out = dict(method="Dense unblocked Pauli tensor products, NumPy, L<=8, single thread",
               frames=literal_audit(), ring=ring())
    print(json.dumps(out, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
