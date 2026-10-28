"""Finite checks of the qubit gap reduction and additive-charge classification.

Python 3 + NumPy only. Run from any directory; writes checks.json beside this file.
These checks validate conventions and examples, not thermodynamic proofs.
"""
import os
for variable in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[variable] = "1"

import json
from pathlib import Path
import time
import numpy as np

PAULI = np.array([[[0, 1], [1, 0]], [[0, -1j], [1j, 0]],
                  [[1, 0], [0, -1]]], dtype=complex)


def cross(n):
    x, y, z = n
    return np.array([[0, -z, y], [z, 0, -x], [-y, x, 0.]])


def rotation(axis, angle):
    axis = np.asarray(axis, float)
    axis /= np.linalg.norm(axis)
    a = cross(axis)
    return np.eye(3) + np.sin(angle)*a + (1-np.cos(angle))*(a@a)


def axis_certificate(J, B, n, tol=1e-10):
    """Floating-point diagnostic of a supplied axis, not a rigorous oracle."""
    n = np.asarray(n, float) / np.linalg.norm(n)
    lam = float(n @ J @ n)
    p = np.eye(3)-np.outer(n, n)
    transverse_norm = float(np.linalg.norm(p @ J @ p, 2))
    residual = float(max(np.linalg.norm(J@n-lam*n),
                         np.linalg.norm(J.T@n-lam*n)))
    passed = residual <= tol and lam <= -transverse_norm+tol
    return dict(axis=n.tolist(), longitudinal_exchange=lam,
                transverse_norm=transverse_norm, axis_residual=residual,
                passes=bool(passed),
                gap_lower_bound=float(2*abs(B@n)) if passed else None)


def ring_hamiltonian(J, B, length):
    """Build H=sum sigma_i^T J sigma_(i+1)+B.sigma_i, bit i = site i."""
    dim = 1 << length
    h = np.zeros((dim, dim), complex)
    states = np.arange(dim)

    def add_word(coefficient, word):
        if coefficient == 0:
            return
        target = states.copy()
        phase = np.ones(dim, complex)
        for site, axis in word:
            bit = (states >> site) & 1
            if axis == 0:
                target ^= 1 << site
            elif axis == 1:
                target ^= 1 << site
                phase *= 1j*(1-2*bit)
            else:
                phase *= 1-2*bit
        h[target, states] += coefficient*phase

    for i in range(length):
        for a in range(3):
            add_word(B[a], [(i, a)])
            for b in range(3):
                add_word(J[a, b], [(i, a), ((i+1) % length, b)])
    assert np.max(abs(h-h.conj().T)) < 1e-12
    return h


def charge_matrix(J, B, length):
    """Local Pauli coefficient equations for [sum n_i.sigma_i,H]=0."""
    m = np.zeros((12*length, 3*length))
    for site in range(length):
        for axis in range(3):
            omega = cross(np.eye(3)[axis])
            m[9*site:9*site+9, 3*site+axis] += (omega@J).ravel()
            following = (site+1) % length
            m[9*site:9*site+9, 3*following+axis] -= (J@omega).ravel()
            m[9*length+3*site:9*length+3*site+3,
              3*site+axis] = omega@B
    return m


def onsite_charge(n):
    length = len(n)
    # A separate Kronecker construction to check the local linear equations.
    q = np.zeros((1 << length, 1 << length), complex)
    for site in range(length):
        term = np.ones((1, 1), complex)
        for j in reversed(range(length)):
            term = np.kron(term, np.einsum("a,aij->ij", n[j], PAULI)
                           if j == site else np.eye(2))
        q += term
    return q


def main():
    start = time.monotonic()
    rng = np.random.default_rng(20261028)
    result = {"seed": 20261028, "gap_examples": [], "charge_examples": []}
    ez = np.array([0., 0., 1.])
    r = rotation([1., 2., -1.], .73)
    examples = [
        ("ferromagnetic_Ising_mixed", np.diag([0., 0., -1.]),
         np.array([.91, .17, .13]), ez),
        ("interacting_chiral_mixed", np.array([[-.32, .21, 0.],
          [-.07, .18, 0.], [0., 0., -1.]]), np.array([.41, -.23, -.19]), ez),
        ("pure_field_saturates", np.zeros((3, 3)), .27*ez, ez),
        ("dominance_equality", np.diag([1., -1., -1.]),
         np.array([.24, 0., .11]), ez),
    ]
    base = examples[1]
    examples.append(("rotated_interacting_chiral", r@base[1]@r.T,
                     r@base[2], r@ez))
    for label, J, B, n in examples:
        cert = axis_certificate(J, B, n)
        assert cert["passes"] and cert["gap_lower_bound"] > 0
        gaps = []
        for length in (4, 6, 8):
            eig = np.linalg.eigvalsh(ring_hamiltonian(J, B, length))
            gap = float(eig[1]-eig[0])
            assert gap >= cert["gap_lower_bound"]-1e-10
            gaps.append({"L": length, "gap": gap})
        result["gap_examples"].append(dict(name=label, J=J.tolist(), B=B.tolist(),
                                            certificate=cert, spectra=gaps))

    # Reject two tempting but invalid extensions of the theorem.
    af = axis_certificate(np.diag([0., 0., 1.]), 2*ez, ez)
    tilted_J = np.array([[.1, 0., .08], [0., .2, 0.], [0., 0., -1.]])
    off_axis = axis_certificate(tilted_J, .2*ez, ez)
    assert not af["passes"] and not off_axis["passes"]
    result["rejected_certificates"] = {"antiferromagnet": af,
                                        "mixed_axis_exchange": off_axis}

    # Exact finite-chain TFIM energy compared to an independent spin matrix.
    critical = []
    for length in (4, 6, 8):
        eig = np.linalg.eigvalsh(ring_hamiltonian(np.diag([0., 0., -1.]),
                                                np.array([-1., 0., 0.]), length))
        expected = -2/np.sin(np.pi/(2*length))
        assert abs(eig[0]-expected) < 1e-10
        critical.append({"L": length, "E0": float(eig[0]),
                         "exact_E0": float(expected),
                         "gap": float(eig[1]-eig[0])})
    result["critical_Ising_ED"] = critical
    result["critical_Ising_exact_Casimir"] = [
        {"L": length, "c_estimator": float(-3*length/np.pi *
          (-2/np.sin(np.pi/(2*length))+4*length/np.pi))}
        for length in (16, 32, 64, 128, 256)]

    # The other gap is ABOVE two exact Neel vacua, not E1-E0.
    staggered = []
    for a, b, offset in ((.8, .6, .25), (.3, 1.4, .1), (-.7, -.4, .3)):
        boundary = np.sqrt(a*a+b*b/4)
        for length in (4, 6, 8):
            c = boundary+offset
            eig = np.linalg.eigvalsh(ring_hamiltonian(np.diag([a, -a, c]),
                                                    b*ez, length))
            degeneracy = int(np.count_nonzero(abs(eig+c*length) < 1e-9))
            assert degeneracy == 2
            gap = float(eig[2]-eig[0])
            assert abs(gap-4*offset) < 1e-9
            edge_eig = np.linalg.eigvalsh(ring_hamiltonian(
                np.diag([a, -a, boundary]), b*ez, length))
            edge_degeneracy = int(np.count_nonzero(abs(edge_eig+boundary*length) < 1e-9))
            assert edge_degeneracy == length+1
            staggered.append(dict(a=a, b=b, c=float(c), L=length,
                                  ground_degeneracy=degeneracy,
                                  gap_above_vacua=gap, exact_gap=4*offset,
                                  boundary_degeneracy=edge_degeneracy))
    result["staggered_exact_wing"] = staggered

    rz = rotation([0., 0., 1.], 2*np.pi/5)
    charge_cases = [
        ("generic_XYZ", np.diag([1., .7, .4]), np.zeros(3), 6, 0),
        ("XXZ_uniform", np.diag([1., 1., .4]), .3*ez, 6, 1),
        ("staggered_even", np.diag([1., -1., .4]), .3*ez, 6, 1),
        ("staggered_odd", np.diag([1., -1., .4]), .3*ez, 5, 0),
        ("Heisenberg", np.eye(3), np.zeros(3), 6, 3),
        ("proper_twist_compatible", rz, np.zeros(3), 5, 3),
        ("proper_twist_incompatible", rz, np.zeros(3), 6, 1),
        ("improper_twist_compatible", -rz, np.zeros(3), 5, 3),
        ("classical_Ising", np.diag([0., 0., -1.]), .3*ez, 6, 6),
        ("transverse_Ising", np.diag([0., 0., -1.]), np.array([.3, 0., 0.]), 6, 0),
        ("oriented_rank_one", np.outer([1., 0., 0.], [0., 1., 0.]), np.zeros(3), 6, 0),
        ("onsite_only", np.zeros((3, 3)), .3*ez, 6, 6),
        ("zero_H", np.zeros((3, 3)), np.zeros(3), 6, 18),
        ("rank_two_staggered", np.diag([1., -1., 0.]), .3*ez, 6, 1),
    ]
    for label, J, B, length, expected in charge_cases:
        m = charge_matrix(J, B, length)
        _, singular, vh = np.linalg.svd(m, full_matrices=False)
        rank = int(np.count_nonzero(singular > 1e-9))
        nullity = 3*length-rank
        assert nullity == expected, (label, nullity, expected)
        mr = charge_matrix(r@J@r.T, r@B, length)
        assert 3*length-np.linalg.matrix_rank(mr, tol=1e-9) == expected
        max_commutator = 0.
        if nullity:
            h = ring_hamiltonian(J, B, length)
            for vec in vh[rank:]:
                q = onsite_charge(vec.reshape(length, 3))
                max_commutator = max(max_commutator,
                                    float(np.linalg.norm(q@h-h@q)))
            assert max_commutator < 1e-9
        result["charge_examples"].append(dict(name=label, L=length,
                                              nullity=nullity,
                                              max_commutator=max_commutator))
    # Independent coefficient-vs-full-matrix commutator norm identity, also
    # for nonconserved vectors. Pauli orthogonality gives the normalization.
    max_norm_error = 0.
    for _ in range(8):
        J, B, n = rng.normal(size=(3, 3)), rng.normal(size=3), rng.normal(size=(4, 3))
        h, q = ring_hamiltonian(J, B, 4), onsite_charge(n)
        local = np.linalg.norm(charge_matrix(J, B, 4)@n.ravel())
        full = np.linalg.norm(q@h-h@q)/(2*np.sqrt(2**4))
        max_norm_error = max(max_norm_error, float(abs(full-local)))
    assert max_norm_error < 1e-10
    result["commutator_norm_identity_max_error"] = max_norm_error
    result["elapsed_seconds"] = time.monotonic()-start
    out = Path(__file__).with_name("checks.json")
    out.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({"gap_examples": len(examples), "spectra": 3*len(examples),
                      "charge_examples": len(charge_cases),
                      "commutator_norm_identity_max_error": max_norm_error,
                      "elapsed_seconds": result["elapsed_seconds"],
                      "output": str(out)}, indent=2))


if __name__ == "__main__":
    main()
