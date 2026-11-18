"""Bounded, one-worker quartet matching in the complete leading blockade model.

No MPS or production calculations. All sectors include the blockade and both
generated attractions. The reported coefficients are finite-volume matrix
elements in a declared spectral convention, not continuum cosine coefficients.
"""
from __future__ import annotations

import os
for _key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
             "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_key] = "1"

import argparse
import itertools
import json
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import scipy.linalg as la
import scipy.sparse as sp
from threadpoolctl import threadpool_limits


@dataclass
class Sector:
    L: int
    N: int
    momentum: int = 0
    blockade: bool = True

    def __post_init__(self):
        self.mask = (1 << self.L) - 1
        self.reps, self.periods, self.lookup = [], [], {}
        for xs in itertools.combinations(range(self.L), self.N):
            state = sum(1 << j for j in xs)
            if not self.allowed(state) or state in self.lookup:
                continue
            orbit, s = [], state
            while s not in orbit:
                orbit.append(s)
                s = self.translate(s)
            # Input combinations are not sorted as integers. Pick a reproducible
            # canonical representative, independently of enumeration order.
            rep = min(orbit)
            orbit, s = [], rep
            while s not in orbit:
                orbit.append(s)
                s = self.translate(s)
            idx = len(self.reps)
            valid = (self.momentum * len(orbit)) % self.L == 0
            for shift, s in enumerate(orbit):
                self.lookup[s] = (idx if valid else -1, shift)
            if valid:
                self.reps.append(rep)
                self.periods.append(len(orbit))
        self.periods = np.array(self.periods)

    def translate(self, s):
        return ((s << 1) & self.mask) | (s >> (self.L - 1))

    def allowed(self, s):
        return not self.blockade or not (s & self.translate(s))

    def element(self, target, source_period):
        row, shift = self.lookup[target]
        if row < 0:
            return None, 0.0
        phase = np.exp(2j * np.pi * self.momentum * shift / self.L)
        return row, np.sqrt(source_period / self.periods[row]) * phase

    def matrix(self, transitions, target=None):
        target = self if target is None else target
        rows, cols, data = [], [], []
        for col, s in enumerate(self.reps):
            for out, amp in transitions(s):
                row, weight = target.element(out, self.periods[col])
                if row is not None:
                    rows.append(row)
                    cols.append(col)
                    data.append(amp * weight)
        return sp.coo_matrix((data, (rows, cols)),
                             shape=(len(target.reps), len(self.reps))).tocsr()

    def hamiltonian(self, kappa=1.0, real_hop=0.0, V2=-2.0, V3=-2/3,
                    t2=-1.0, V1=0.0):
        w = real_hop + 1j * kappa
        def transitions(s):
            diag = sum(V1 * ((s >> j) & 1) * ((s >> ((j+1) % self.L)) & 1)
                       + V2 * ((s >> j) & 1) * ((s >> ((j+2) % self.L)) & 1)
                       + V3 * ((s >> j) & 1) * ((s >> ((j+3) % self.L)) & 1)
                       for j in range(self.L))
            yield s, diag
            for j in range(self.L):
                if not (s >> j) & 1:
                    continue
                for d, amp in ((1, w), (-1, w.conjugate()), (2, t2), (-2, t2)):
                    if amp == 0:
                        continue
                    dest = (j + d) % self.L
                    if (s >> dest) & 1:
                        continue
                    out = s ^ (1 << j) ^ (1 << dest)
                    if self.allowed(out):
                        yield out, amp
        return self.matrix(transitions)

    def pair(self, target, separation=3):
        def transitions(s):
            for j in range(self.L):
                bits = (1 << j) | (1 << ((j + separation) % self.L))
                if s & bits:
                    continue
                out = s | bits
                if target.allowed(out):
                    yield out, 1.0
        return self.matrix(transitions, target)

    def sigma(self):
        assert self.N % 2 == 0
        def transitions(s):
            xs = [j for j in range(self.L) if (s >> j) & 1]
            out = sum(1 << ((-j) % self.L) for j in xs)
            yield out, (-1) ** (self.N // 2 + sum(xs))
        return self.matrix(transitions)

    def sigma_basis(self, sign=1):
        S = self.sigma().tocsc()
        assert np.max(np.abs((S @ S - sp.eye(S.shape[0])).toarray()), initial=0) < 1e-12
        rows, cols, data, seen = [], [], [], set()
        col = 0
        for j in range(S.shape[0]):
            if j in seen:
                continue
            start, stop = S.indptr[j:j+2]
            assert stop - start == 1
            i, coeff = int(S.indices[start]), S.data[start]
            seen.update((i, j))
            if i == j:
                if abs(coeff-sign) < 1e-10:
                    rows.append(j); cols.append(col); data.append(1.0)
                    col += 1
            else:
                rows.extend((j, i)); cols.extend((col, col))
                data.extend((1/np.sqrt(2), sign*coeff/np.sqrt(2)))
                col += 1
        return sp.coo_matrix((data, (rows, cols)), shape=(S.shape[0], col)).tocsr()


def ground(H):
    vals, vecs = la.eigh(H, subset_by_index=(0, min(1, len(H)-1)))
    v = vecs[:, 0]
    # A reproducible phase is useful, although only |A| is invariant under
    # independent phases of the external dressed eigenstates.
    pivot = int(np.argmax(np.abs(v)))
    v = v * np.exp(-1j * np.angle(v[pivot]))
    return float(vals[0]), v, float(vals[1]-vals[0]) if len(vals)>1 else None


def match(L, N=0, kappa=1.0, sigma_sign=None):
    """Exact second-order coupling of two chosen degenerate dressed states.

    H0 contains no chemical potential. Tune mu=(E_(N+4)-E_N)/4 only to
    define an on-shell finite-volume matching experiment. This is not a
    prescription for an L-dependent microscopic Hamiltonian in a CFT claim.
    """
    if sigma_sign is None:
        sigma_sign = (-1) ** (N//2)
    momentum = ((N // 2) * (L // 2)) % L
    sectors = [Sector(L, n, momentum) for n in (N, N+2, N+4)]
    bases = [s.sigma_basis(sigma_sign) for s in sectors]
    hs = [(U.conj().T @ s.hamiltonian(kappa) @ U).toarray()
          for s, U in zip(sectors, bases)]
    err_herm = max(float(np.max(np.abs(H-H.conj().T))) for H in hs)
    assert err_herm < 1e-11
    ei, vi, gap_i = ground(hs[0])
    ef, vf, gap_f = ground(hs[2])
    mu = (ef - ei)/4
    bar = (ei+ef)/2
    Bs = [(bases[i+1].conj().T @ sectors[i].pair(sectors[i+1]) @ bases[i]).toarray()
          for i in range(2)]
    rhs = Bs[0] @ vi
    den = hs[1] - bar*np.eye(len(hs[1]))
    sol = la.solve(den, rhs, assume_a="her")
    A = -np.vdot(vf, Bs[1] @ sol)
    mid_evals, mid_vecs = la.eigh(hs[1])
    terms = -(vf.conj() @ Bs[1] @ mid_vecs) * (mid_vecs.conj().T @ rhs)/(mid_evals-bar)
    connected = None
    if N == 0:
        # D^dagger|vac> = (H2-2mu)^-1 B3^dagger|vac>. Extend the exact
        # two-particle solution as a projected pair creator on all sectors.
        D2 = np.zeros_like(Bs[1])
        reconstructed = np.zeros_like(sol)
        pair_coefficients = {}
        for r in range(2, L//2+1):
            p0 = (bases[1].conj().T @ sectors[0].pair(sectors[1], r) @ bases[0]).toarray()[:, 0]
            norm = np.vdot(p0, p0).real
            if norm < 1e-14:
                continue
            coeff = np.vdot(p0, sol)/norm
            pair_coefficients[str(r)] = [float(coeff.real), float(coeff.imag)]
            reconstructed += coeff*p0
            D2 += coeff*(bases[2].conj().T @ sectors[1].pair(sectors[2], r) @ bases[1]).toarray()
        pair_product = D2 @ sol
        R = (hs[2]-ef*np.eye(len(hs[2]))) @ pair_product - 2*Bs[1]@sol
        Ac = np.vdot(vf, R)/2
        connected = dict(identity_error=float(abs(A-Ac)),
                         pair_reconstruction_error=float(la.norm(reconstructed-sol)),
                         connected_source_norm=float(la.norm(R)),
                         pair_coefficients=pair_coefficients)
    return dict(L=L, N=N, kappa=kappa, momentum_pi=2*momentum/L,
                sigma=sigma_sign, dims=[len(H) for H in hs],
                initial_energy=ei, final_energy=ef, matching_mu=mu,
                initial_gap_in_sector=gap_i, final_gap_in_sector=gap_f,
                intermediate_min_gap=float(np.min(mid_evals-bar)),
                intermediate_pole_distance=float(np.min(np.abs(mid_evals-bar))),
                amplitude_real=float(A.real), amplitude_imag=float(A.imag),
                amplitude_abs=float(abs(A)),
                spectral_terms_abs_sum=float(np.sum(np.abs(terms))),
                spectral_sum_error=float(abs(A-np.sum(terms))),
                connected=connected,
                resolvent_residual=float(la.norm(den@sol-rhs)),
                hermiticity_error=err_herm)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    start = time.perf_counter()
    with threadpool_limits(limits=1):
        rows = [match(14, 0, 2)]
        rows += [match(L) for L in (12, 16, 20, 24)]
        rows += [match(L, 2) for L in (16, 20, 24)]
        rows += [match(20, 0, k) for k in (0.5, 2.0)]
    out = dict(description=__doc__, seconds=time.perf_counter()-start, rows=rows)
    print(json.dumps(out, indent=2))
    if args.output:
        args.output.write_text(json.dumps(out, indent=2)+"\n")


if __name__ == "__main__":
    main()
