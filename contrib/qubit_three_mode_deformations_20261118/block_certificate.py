"""Small open-block positivity tests for the complete m=3 leading model.

The block sum is an identity on every sufficiently long constrained ring;
these are local-operator certificates, not finite-size phase estimates.
One BLAS worker. No background calculation or DMRG.
"""
from __future__ import annotations

import os
for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[key] = "1"

from fractions import Fraction
from itertools import combinations
from pathlib import Path
import json
import time

import numpy as np
from scipy import linalg
from threadpoolctl import threadpool_limits


def words(length: int, number: int):
    return [sum(1 << j for j in occ) for occ in combinations(range(length), number)
            if all(b-a not in (1, 2, 4) for a, b in combinations(occ, 2))]


def integer_operators(length: int, number: int):
    """Open sum h3, T1 and the conditioned distance-eight potential.

    h3 = |01-10><01-10| at sites j,j+3, then compressed to the
    complete open no-1,2,4 space. The compression keeps blocked-hop
    diagonal entries. T1 moves one physical site right.
    """
    states = words(length, number)
    lookup = {s:i for i, s in enumerate(states)}
    size = len(states)
    base = np.zeros((size, size), dtype=np.int64)
    hop = np.zeros_like(base)
    potential = np.zeros(size, dtype=np.int64)
    for col, state in enumerate(states):
        for j in range(length-3):
            k = j+3
            nj, nk = (state >> j) & 1, (state >> k) & 1
            if nj != nk:
                base[col, col] += 1
                dest = state ^ (1 << j) ^ (1 << k)
                if dest in lookup:
                    base[lookup[dest], col] -= 1
        for j in range(length-1):
            if (state >> j) & 1 and not (state >> (j+1)) & 1:
                dest = state ^ (1 << j) ^ (1 << (j+1))
                if dest in lookup:
                    hop[lookup[dest], col] += 1
        for j in range(length-8):
            potential[col] += (((state >> j) & 1)*((state >> (j+8)) & 1)
                               *(1-((state >> (j+3)) & 1))
                               *(1-((state >> (j+5)) & 1)))
    return states, base, hop, potential


def block(length, number, tau, mu):
    states, base, hop, potential = integer_operators(length, number)
    matrix = (base.astype(complex)/(length-3)
              +(tau*hop+np.conj(tau)*hop.T)/(length-1)
              +np.diag((3/16)*potential/(length-8)
                       +float(mu)*(1/3-number/length)))
    return states, matrix


def threshold(length, tau):
    report = []
    for number in range(length//3+1):
        states, matrix = block(length, number, tau, 0)
        minimum = float(linalg.eigvalsh(matrix, subset_by_index=[0, 0])[0])
        holes = 1/3-number/length
        needed = max(0., -minimum/holes) if holes else None
        report.append(dict(N=number, dimension=len(states), lowest_at_mu_zero=minimum,
                           sufficient_mu=needed))
    return report


def exact_ldl(matrix):
    """Fraction LDL, including zero-pivot checks; returns an exact PSD proof."""
    a = [[Fraction(x) for x in row] for row in matrix]
    pivots = []
    for k in range(len(a)):
        pivot = a[k][k]
        if pivot < 0:
            raise AssertionError(("negative pivot", k, str(pivot)))
        column = [a[i][k] for i in range(k+1, len(a))]
        if pivot == 0:
            if any(column):
                raise AssertionError(("zero pivot with nonzero column", k))
        else:
            for ii, x in enumerate(column, k+1):
                if x:
                    for jj in range(ii, len(a)):
                        if a[jj][k]:
                            a[jj][ii] -= x*a[jj][k]/pivot
                            a[ii][jj] = a[jj][ii]
        pivots.append(pivot)
    positive = [p for p in pivots if p]
    return dict(dimension=len(a), rank=len(positive), zero_pivots=len(a)-len(positive),
                minimum_positive_pivot=str(min(positive)) if positive else None)


def exact_block(length, number, tau_real, tau_imag, mu):
    """Integer scaling of the realification, suitable for exact LDL."""
    _, base, hop, potential = integer_operators(length, number)
    rr, ii = Fraction(tau_real), Fraction(tau_imag)
    mu = Fraction(mu)
    real = [[Fraction(int(base[i, j]), length-3)
             +rr*Fraction(int(hop[i, j]+hop[j, i]), length-1)
             for j in range(len(base))] for i in range(len(base))]
    imag = [[ii*Fraction(int(hop[i, j]-hop[j, i]), length-1)
             for j in range(len(base))] for i in range(len(base))]
    for i in range(len(base)):
        real[i][i] += Fraction(3*int(potential[i]), 16*(length-8))+mu*(Fraction(1,3)-Fraction(number,length))
    if not ii:
        return real
    return [real[i]+[-x for x in imag[i]] for i in range(len(base))]+[
        imag[i]+real[i] for i in range(len(base))]


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--length", type=int, default=12)
    parser.add_argument("--gamma", type=float, default=.1)
    parser.add_argument("--exact", action="store_true")
    args = parser.parse_args()
    start = time.monotonic()
    with threadpool_limits(1):
        reports = {}
        for name, phase in (("real_positive", 1), ("real_negative", -1),
                            ("imag_positive", 1j)):
            reports[name] = threshold(args.length, args.gamma*phase)
        if args.exact:
            certificates = {}
            for rr, ii in (("1/10", "0"), ("-1/10", "0"), ("0", "1/10")):
                key = f"{rr}+i({ii})"
                certificates[key] = [exact_ldl(exact_block(args.length, number, rr, ii, "7/100"))
                                     for number in range(args.length//3+1)]
            reports["exact_certificates"] = certificates
    reports["runtime_seconds"] = time.monotonic()-start
    print(json.dumps(reports, indent=2))
