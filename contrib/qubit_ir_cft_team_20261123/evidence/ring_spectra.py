"""Exact target-spin spectral weights on periodic L=8,10,12 rings (NumPy)."""
import json
import os
from pathlib import Path
import signal
import time

for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
             "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[name] = "1"

import numpy as np


HERE = Path(__file__).resolve().parent
SIZES = (8, 10, 12)
A, B, FIELD = 19 / 32, 13 / 32, 1783 / 1000
DEGENERACY = 1.e-9
NAMES = ("n", "kappa2", "B2", "P2", "P3", "b")


def translation_orbits(length):
    size = 1 << length
    shift = np.full(size, -1)
    parent = np.full(size, -1)
    orbits = []
    for s in range(size):
        if parent[s] != -1:
            continue
        t, sites = s, []
        while parent[t] == -1:
            parent[t] = len(orbits)
            shift[t] = len(sites)
            sites.append(t)
            t = ((t << 1) % size) | (t >> (length - 1))
        orbits.append(np.array(sites))
    return orbits, parent, shift


def solve(length):
    """Basis vector for orbit s is sum_r exp(-ikr)|T^r s>/sqrt(period)."""
    orbits, parent, shift = translation_orbits(length)
    sizes = np.array([len(o) for o in orbits])
    parities = np.array([int(o[0]).bit_count() % 2 for o in orbits])
    sectors = {}
    max_residual = 0.
    for p in (0, 1):
        for k in range(length):
            ids = np.flatnonzero((parities == p) & ((sizes * k) % length == 0))
            starts = np.array([orbits[i][0] for i in ids])
            columns = np.arange(len(ids))
            local_id = np.full(len(orbits), -1)
            local_id[ids] = columns
            z = np.array([1 - 2 * ((starts >> j) & 1) for j in range(length)])
            diag = (z * np.roll(z, -1, axis=0)).sum(axis=0) - FIELD * z.sum(axis=0)
            h = np.diag(diag.astype(complex))
            phase = 2 * np.pi * k / length
            for j in range(length):
                b = (j + 1) % length
                dest = starts ^ (1 << j) ^ (1 << b)
                rows = local_id[parent[dest]]
                ok = rows >= 0
                amp = 1j * (A * z[b] + B * z[j])
                factor = np.sqrt(sizes[ids] / sizes[parent[dest]])
                h[rows[ok], columns[ok]] += (amp * factor * np.exp(1j * phase * shift[dest]))[ok]
            assert np.max(np.abs(h - h.conj().T)) < 1.e-12
            e, v = np.linalg.eigh(h)
            residual = np.linalg.norm(h @ v - v * e, axis=0).max()
            max_residual = max(max_residual, float(residual))
            assert residual < 1.e-10
            dest = np.concatenate([orbits[i] for i in ids])
            transform = np.exp(-1j * phase * shift[dest]) / np.sqrt(sizes[parent[dest]])
            sectors[p, k] = (e, v, dest, transform, sizes[ids])
    assert sum(len(s[0]) for s in sectors.values()) == 1 << length
    return sectors, max_residual


def local_probes(ground, length):
    """Original-frame local annihilators, down-spin density, chirality and pair quadrature."""
    sites = np.arange(1 << length)
    n = (sites & 1)
    n2, n3 = (sites >> 2) & 1, (sites >> 3) & 1
    z, z2 = 1 - 2 * n, 1 - 2 * n2
    probes = {}
    for name, mask, amp in (("n", 0, n),
                            ("kappa2", 5, 1j * (z2 - z)),
                            ("B2", 5, 1 + z * z2),
                            ("P2", 5, n * n2),
                            ("P3", 9, n * n3),
                            ("b", 1, n)):
        out = np.empty_like(ground)
        out[sites ^ mask] = amp * ground
        probes[name] = out
    return probes


def project(sector, vector):
    _, v, dest, phase, sizes = sector
    starts = np.cumsum(sizes) - sizes
    basis = np.add.reduceat(phase.conj() * vector[dest], starts)
    return v.conj().T @ basis


def calc(length, sectors, max_residual):
    ground_energy = min(s[0][0] for s in sectors.values())
    ground_indices = [(p, k, i) for (p, k), s in sectors.items()
                      for i in np.flatnonzero(s[0] < ground_energy + DEGENERACY)]
    # List the entire GS multiplet and use its normalized projector in all observables.
    expectation = dict.fromkeys(NAMES, 0j)
    sources = []
    for p, k, i in ground_indices:
        _, eigvec, dest, phase, sizes = sectors[p, k]
        g = np.zeros(1 << length, dtype=complex)
        g[dest] = phase * np.repeat(eigvec[:, i], sizes)
        probes = local_probes(g, length)
        for name in NAMES:
            expectation[name] += np.vdot(g, probes[name]) / len(ground_indices)
        sources.append((p, k, probes))
    momenta, poles = [], []
    # At fixed (q,E), sum M_a^* M_b over all degenerate targets and average all GSs.
    for m in range(length):
        entries = []
        for p, k, probes in sources:
            for target_p in (0, 1):
                sector = sectors[target_p, (k + m) % length]
                gaps = sector[0] - ground_energy
                coeff = np.zeros((len(gaps), len(NAMES)), dtype=complex)
                for a, name in enumerate(NAMES):
                    if target_p == (p ^ int(name == "b")):
                        coeff[:, a] = np.sqrt(length / len(sources)) * project(sector, probes[name])
                entries.extend(zip(gaps, coeff))
        entries.sort(key=lambda entry: entry[0])
        groups = []
        for energy, vec in entries:
            if not groups or energy - groups[-1][0] > DEGENERACY:
                groups.append((max(float(energy), 0.), np.zeros((len(NAMES), len(NAMES)), complex)))
            groups[-1][1][:] += np.outer(vec.conj(), vec)
        elastic = sum((w.diagonal().real for e, w in groups if e < DEGENERACY), np.zeros(len(NAMES)))
        total = sum(w.diagonal().real for e, w in groups)
        inelastic = total - elastic
        connected = total - (length * np.abs(list(expectation.values())) ** 2 if m == 0 else 0)
        mode = {}
        for a, name in enumerate(NAMES):
            active = [(e, w[a, a].real) for e, w in groups if e > DEGENERACY and w[a, a].real > 1.e-11]
            mode[name] = dict(S_connected=float(connected[a]),
                              elastic_connected=float(connected[a] - inelastic[a]),
                              S_inelastic=float(inelastic[a]),
                              first_gap=active[0][0] if active else None,
                              first_weight=float(active[0][1]) if active else 0.)
        momenta.append(dict(m=m, probes=mode))
        if m in (1, length - 1, length // 2):
            selected = [(e, w) for e, w in groups
                        if e > DEGENERACY and np.max(w.diagonal().real[:-1]) > 1.e-11][:6]
            for e, w in selected:
                poles.append(dict(m=m, gap=e, weight={n: float(w[a, a].real) for a, n in enumerate(NAMES[:-1])},
                                  cross_n_kappa2=[float(x) for x in (w[0, 1].real, w[0, 1].imag)],
                                  cross_B2_kappa2=[float(x) for x in (w[2, 1].real, w[2, 1].imag)]))
    # Full momentum sum, including both the elastic and inelastic connected parts.
    for name in NAMES:
        local_variance = sum(np.linalg.norm(probes[name]) ** 2 for p, k, probes in sources) / len(sources)
        local_variance -= abs(expectation[name]) ** 2
        assert abs(sum(r["probes"][name]["S_connected"] for r in momenta) / length - local_variance) < 1.e-10
    return dict(L=length, E0=float(ground_energy), ground_multiplet=[[p, k, int(i)] for p, k, i in ground_indices],
                max_eigenpair_residual=max_residual,
                lowest_energy_by_sector=[dict(p=p, k=k, E=float(s[0][0])) for (p, k), s in sectors.items()],
                momenta=momenta, six_lowest_neutral_poles=poles)


def main():
    def timeout(_signal, _frame):
        raise TimeoutError("90-second hard cap")
    signal.signal(signal.SIGALRM, timeout)
    signal.alarm(90)
    t0 = time.monotonic()
    rings = []
    for length in SIZES:
        rings.append(calc(length, *solve(length)))
    result = dict(parameters=dict(sizes=SIZES, a=A, b=B, field=FIELD, boundary="periodic spin",
                                  momentum="T|g>=exp(i kg)|g>, q=2 pi m/L, O(q)=sum_j exp(-i q j) O_j/sqrt(L)",
                                  energy_group_tolerance=DEGENERACY, weight_threshold=1.e-11,
                                  time_limit_seconds=90, threads=1),
                  runtime_seconds=time.monotonic() - t0, rings=rings)
    (HERE / "ring_spectra.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(f"Checked {SIZES} in {result['runtime_seconds']:.3f} s")
    signal.alarm(0)


if __name__ == "__main__":
    main()
