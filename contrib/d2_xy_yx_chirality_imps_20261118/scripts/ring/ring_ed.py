#!/usr/bin/env python3
"""ring_ed.py
Written by OpenAI Codex (gpt-6.2-astra) to a spec by Claude Opus 5.5, 2026-11-24, self-checked against dense ED.

Usage: python ring_ed.py L=20 delta=0.1875 h=1.783 nev=6 ms=all check=1 vecs=0

Site j is bit j, with 1 denoting down. ms selects the spectra to print.
Unlisted momenta still contribute their lowest eigenvalue to the global minimum.
"""
import os
import sys
from time import perf_counter

# ARPACK's vector operations can slow down with many BLAS threads.
for key in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
            "BLIS_NUM_THREADS", "OMP_NUM_THREADS"):
    os.environ.setdefault(key, "1")

import numpy as np
from scipy.linalg import eigh
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigsh


def translate(x, L):
    return ((x << 1) & ((1 << L) - 1)) | (x >> (L - 1))


def pop(x):
    table = np.array([bin(i).count("1") for i in range(256)], dtype=np.uint8)
    octets = np.ascontiguousarray(x, dtype=np.uint32).view(np.uint8)
    return table[octets.reshape(-1, 4)].sum(axis=1, dtype=np.int32)


class Orbits:
    def __init__(self, L):
        self.L = L
        chunks = []
        for start in range(0, 1 << L, 1 << 18):
            x = np.arange(start, min(start + (1 << 18), 1 << L), dtype=np.int32)
            smallest = x.copy()
            y = x
            for _ in range(1, L):
                y = translate(y, L)
                np.minimum(smallest, y, out=smallest)
            chunks.append(x[x == smallest])

        self.rep = np.concatenate(chunks)
        self.down = pop(self.rep)
        self.period = np.full(self.rep.size, L, dtype=np.int32)

        # If x = T^s rep[i], address[x] = L*i + s, with 0 <= s < period[i].
        self.address = np.empty(1 << L, dtype=np.int32)
        base = L * np.arange(self.rep.size, dtype=np.int32)
        y = self.rep
        for s in range(L):
            live = s < self.period
            self.address[y[live]] = base[live] + s
            y = translate(y, L)
            self.period[(y == self.rep) & live] = s + 1


class Bonds:
    """Cache flipped representatives and their shifts across all momenta."""

    def __init__(self, orb, p, delta, h):
        self.orb = orb
        L = orb.L
        self.ids = np.flatnonzero((orb.down & 1) == p)
        x = orb.rep[self.ids]
        self.diag = (L - 2 * pop(x ^ translate(x, L))
                     - h * (L - 2 * orb.down[self.ids]))

        shape = (L, x.size)
        self.dst = np.empty(shape, dtype=np.int32)
        self.shift = np.empty(shape, dtype=np.uint8)
        self.amp = np.empty(shape, dtype=np.float64)

        # H[flipped,x] / i, indexed by 2*bit_j + bit_{j+1}.
        weight = np.array([1., -delta, delta, -1.])
        norm = np.sqrt(orb.period)
        for j in range(L):
            n = (j + 1) % L
            flip = x ^ ((1 << j) | (1 << n))
            dst, s = np.divmod(orb.address[flip], L)
            bits = 2 * ((x >> j) & 1) + ((x >> n) & 1)
            self.dst[j] = dst
            self.shift[j] = s
            self.amp[j] = weight[bits] * norm[self.ids] / norm[dst]

    def block(self, m):
        orb = self.orb
        L = orb.L
        take = np.flatnonzero((m * orb.period[self.ids]) % L == 0)
        ids = self.ids[take]
        size = ids.size
        index = np.full(orb.rep.size, -1, dtype=np.int32)
        index[ids] = np.arange(size, dtype=np.int32)

        cols = np.empty((size, L + 1), dtype=np.int32)
        vals = np.empty((size, L + 1), dtype=np.complex128)
        cols[:, 0] = index[ids]
        vals[:, 0] = self.diag[take]

        # Row = source, so store the conjugate of H[destination, source].
        phase = -1j * np.exp((-2j * np.pi / L) * ((m * np.arange(L)) % L))
        for j in range(L):
            col = index[self.dst[j, take]]
            cols[:, j + 1] = np.maximum(col, 0)
            v = phase[self.shift[j, take]] * self.amp[j, take]
            v[col < 0] = 0
            vals[:, j + 1] = v

        ptr = np.arange(size + 1, dtype=np.int32) * (L + 1)
        H = csr_matrix((vals.ravel(), cols.ravel(), ptr), shape=(size, size))
        H.sum_duplicates()
        H.eliminate_zeros()
        return H, ids


def lowest(H, count, vector=False):
    n = H.shape[0]
    count = min(count, n)
    if n == 0:
        return np.empty(0)
    if n < 2000 or count >= n - 1:
        subset = None if count == n else (0, count - 1)
        return eigh(H.toarray(), eigvals_only=not vector, subset_by_index=subset,
                    check_finite=False)

    result = eigsh(H, k=count, which="SA", tol=1e-13,
                   return_eigenvectors=vector)
    if vector:
        e, v = result
        order = np.argsort(e)
        return e[order], v[:, order]
    return np.sort(result)


def direct_hamiltonian(L, delta, h):
    X = np.array([[0, 1], [1, 0]])
    Y = np.array([[0, -1j], [1j, 0]])
    Z = np.diag([1, -1])
    I = np.eye(2)

    def term(operators):
        A = np.ones((1, 1))
        for j in reversed(range(L)):
            A = np.kron(A, operators.get(j, I))
        return A

    H = np.zeros((1 << L, 1 << L), dtype=complex)
    for j in range(L):
        n = (j + 1) % L
        H += ((1 + delta) / 2) * term({j: X, n: Y})
        H += ((1 - delta) / 2) * term({j: Y, n: X})
        H += term({j: Z, n: Z})
        H -= h * term({j: Z})
    return H


def self_check(L, delta, h):
    H = direct_hamiltonian(L, delta, h)
    E = eigh(H, eigvals_only=True)
    x = np.arange(E.size)

    # Apply T directly by permuting tensor axes in order (L-1,...,0).
    # For a vector v, (T v)[x] = v[invT[x]].
    invT = x.reshape((2,) * L).transpose((*range(1, L), 0)).ravel()
    T = np.argsort(invT)
    P = 1 - 2 * (pop(x) & 1)

    # Obtain simultaneous full-space eigenvectors even at energy degeneracies.
    # Q = 2P + (z T + conjugate(z) T†)/2 has a distinct real eigenvalue
    # for each (m,p). In H + scale*Q the sector energy intervals are disjoint.
    z = np.exp(-0.5j * np.pi / L)
    tags = np.cos((2 * np.pi * np.arange(L) - np.pi / 2) / L)
    gap = np.min(np.diff(np.sort(np.r_[2 + tags, -2 + tags])))
    scale = (np.ptp(E) + 1) / gap
    joint = H.copy()
    joint[x, x] += 2 * scale * P
    joint[x, invT] += (scale / 2) * z
    joint[x, T] += (scale / 2) * z.conjugate()
    _, V = eigh(joint, driver="evd", overwrite_a=True)

    tv = V[invT]
    phase = np.sum(V.conj() * tv, axis=0)
    mvec = np.rint(np.angle(phase) * L / (2 * np.pi)).astype(int) % L
    pv = P[:, None] * V
    pvec = (np.sum(V.conj() * pv, axis=0).real < 0).astype(int)

    def maxnorm(A):
        return np.linalg.norm(A, axis=0).max()

    terr = maxnorm(tv - V * np.exp(2j * np.pi * mvec / L))
    perr = maxnorm(pv - V * (1 - 2 * pvec))
    HV = H @ V
    energy = np.sum(V.conj() * HV, axis=0).real
    herr = maxnorm(HV - V * energy)

    orb = Orbits(L)
    spectra = []
    sector_err = 0.
    for p in range(2):
        bonds = Bonds(orb, p, delta, h)
        for m in range(L):
            B, _ = bonds.block(m)
            e = lowest(B, B.shape[0])
            ref = np.sort(energy[(mvec == m) & (pvec == p)])
            if e.size != ref.size:
                raise AssertionError(f"L={L} m={m} p={p}: {e.size} != {ref.size}")
            if e.size:
                sector_err = max(sector_err, np.max(np.abs(e - ref)))
            spectra.append(e)

    union_err = np.max(np.abs(np.sort(np.concatenate(spectra)) - E))
    print(f"check L={L} union={union_err:.3e} sector={sector_err:.3e} "
          f"T_res={terr:.3e} P_res={perr:.3e} H_res={herr:.3e}", flush=True)
    if not max(union_err, sector_err, terr, perr, herr) <= 1e-10:
        raise AssertionError("dense self-check failed")


def observables(orb, ids, m, v):
    L = orb.L
    ids_x, s = np.divmod(orb.address, L)
    coeff = np.zeros(orb.rep.size, dtype=complex)
    coeff[ids] = v / np.sqrt(orb.period[ids])
    psi = coeff[ids_x] * np.exp((-2j * np.pi / L) * ((m * s) % L))
    prob = np.abs(psi) ** 2
    if abs(prob.sum() - 1) > 1e-10:
        raise AssertionError("ground vector normalization failed")

    x = np.arange(prob.size, dtype=np.int32)
    z0 = 1 - 2 * (x & 1)
    mean = prob @ z0
    print(f"<Z_0>={mean:.12f}")
    for r in range(1, L // 2 + 1):
        zr = 1 - 2 * ((x >> r) & 1)
        connected = prob @ (z0 * zr) - mean * (prob @ zr)
        print(f"r={r} <Z_0 Z_r>_connected={connected:.12f}")


def run(L, delta, h, nev, ms, vecs):
    """With vecs=True, return ((m, parity), representatives, ground vector)."""
    orb = Orbits(L)
    best = np.inf
    for p in range(2):
        bonds = Bonds(orb, p, delta, h)
        for m in range(L):
            start = perf_counter()
            H, ids = bonds.block(m)
            e = lowest(H, nev if m in ms else 1)
            if m in ms:
                values = " ".join(f"{v:.12f}" for v in e)
                print(f"L={L} m={m} parity={p} dim={ids.size} "
                      f"time={perf_counter()-start:.3f}s E=[{values}]", flush=True)
            if e.size and e[0] < best:
                best = e[0]
                winner = (m, p)
                if vecs:
                    best_H, best_ids = H, ids

    m, p = winner
    print(f"global_ground={best:.12f} L={L} m={m} parity={p}", flush=True)
    if vecs:
        _, V = lowest(best_H, 1, vector=True)
        v = V[:, 0]
        if L <= 20:
            observables(orb, best_ids, m, v)
        else:
            print("Skipped full-basis observables for L > 20.")
        return winner, orb.rep[best_ids], v


def main(args):
    opts = dict(L="20", delta="0.1875", h="1.783", nev="6",
                ms="all", check="0", vecs="0")
    for arg in args:
        key, sep, val = arg.partition("=")
        if not sep or key not in opts:
            raise SystemExit(f"Unknown argument {arg!r}. Use key=value.")
        opts[key] = val

    L, nev = int(opts["L"]), int(opts["nev"])
    delta, h = float(opts["delta"]), float(opts["h"])
    if not 2 <= L <= 30 or nev < 1 or not np.isfinite([delta, h]).all():
        raise SystemExit("Require 2 <= L <= 30, nev >= 1 and finite delta, h.")
    if any(opts[key] not in ("0", "1") for key in ("vecs", "check")):
        raise SystemExit("check and vecs take 0 or 1.")
    ms = set(range(L)) if opts["ms"] == "all" else set(map(int, opts["ms"].split(",")))
    if not ms or not ms <= set(range(L)):
        raise SystemExit("ms must be all or momenta in 0..L-1.")

    print("# T shifts site j to j+1. T|m> = exp(+2*pi*i*m/L)|m>. "
          "parity = N_down mod 2. Times cover block assembly and solve.", flush=True)
    if opts["check"] == "1":
        for size in (8, 10):
            self_check(size, delta, h)
    return run(L, delta, h, nev, ms, opts["vecs"] == "1")


if __name__ == "__main__":
    main(sys.argv[1:])
