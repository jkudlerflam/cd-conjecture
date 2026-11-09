"""Blockaded (Fibonacci-neighbourhood) effective model of the parity-class qubit chain.

Full model (C152 eq. 1.2, open chain, Jordan-Wigner fermions):
    H = sum_j [ w c+_{j+1} c_j + h.c. - Dl (c_j c_{j+1} + h.c.) + V n_j n_{j+1} - mu n_j ]
with w = Jx + Jy + 2iD, Dl = Jx - Jy, V = 4 Jz, mu = 4 Jz + 2 hz.

At large V the low-energy space P has no two adjacent particles. Pairing only acts virtually.
eff_rules(): exact second-order Schrieffer-Wolff of the pairing inside P, computed by enumeration
(fermion signs included), plus the first-order hop w and -mu n inside P.
eff_analytic(): the closed form derived in PROOF.md section 3:
    t2 = -Dl^2/(2V) NNN hop, V2 = 2 t2 at distance 2, V3 = (2/3) t2 at distance 3,
    on-site shift +3 Dl^2/V relative to the vacuum (bulk), vacuum energy -(L-1) Dl^2/V.
Sites are bits: site j <-> bit j. Fermion order: c+_{j1} c+_{j2} ... |0>, j1 < j2 < ...
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from itertools import combinations


def popc(x):
    return bin(x).count("1")


def configs(L, N=None):
    out = []
    for s in range(1 << L):
        if s & (s >> 1):
            continue
        if N is not None and popc(s) != N:
            continue
        out.append(s)
    return out


def c_dag(s, j):
    """Apply c+_j to basis state s. Returns (sign, s') or None."""
    if (s >> j) & 1:
        return None
    sign = -1 if popc(s & ((1 << j) - 1)) % 2 else 1
    return sign, s | (1 << j)


def c_ann(s, j):
    if not (s >> j) & 1:
        return None
    sign = -1 if popc(s & ((1 << j) - 1)) % 2 else 1
    return sign, s & ~(1 << j)


def apply_seq(s, ops):
    """ops applied right to left: list of ('d'|'a', site), rightmost first in the list order given reversed."""
    sign = 1
    for kind, j in reversed(ops):
        r = c_dag(s, j) if kind == 'd' else c_ann(s, j)
        if r is None:
            return None
        sg, s = r
        sign *= sg
    return sign, s


def nadj(s):
    return popc(s & (s >> 1))


def eff_rules(L, Dl, V, w, mu, N=None):
    """Second-order SW for pairing inside P (open chain), exact by enumeration. Returns (basis, H dense or sparse)."""
    basis = configs(L, N) if N is not None else configs(L)
    idx = {s: i for i, s in enumerate(basis)}
    H = {}

    def add(i, k, v):
        H[(i, k)] = H.get((i, k), 0) + v

    for i, s in enumerate(basis):
        # -mu n
        add(i, i, -mu * popc(s))
        # first-order hop w c+_{j+1} c_j + h.c. within P
        for j in range(L - 1):
            r = apply_seq(s, [('d', j + 1), ('a', j)])
            if r is not None and r[1] in idx and nadj(r[1]) == 0:
                add(idx[r[1]], i, w * r[0])
            r = apply_seq(s, [('d', j), ('a', j + 1)])
            if r is not None and r[1] in idx and nadj(r[1]) == 0:
                add(idx[r[1]], i, np.conj(w) * r[0])
        # second order pairing: H_pair = -Dl sum_b (c_b c_{b+1} + c+_{b+1} c+_b)
        # creation step: -Dl c+_{b+1} c+_b
        for b in range(L - 1):
            r = apply_seq(s, [('d', b + 1), ('d', b)])
            if r is None:
                continue
            sg1, m = r
            Em = V * nadj(m)  # energy above P (P states have nadj = 0)
            amp1 = -Dl * sg1
            for b2 in range(L - 1):
                r2 = apply_seq(m, [('a', b2), ('a', b2 + 1)])
                if r2 is None:
                    continue
                sg2, f = r2
                if nadj(f) != 0 or f not in idx:
                    continue
                add(idx[f], i, -(amp1 * (-Dl * sg2)) / Em)
        # annihilation first is impossible inside P (needs an adjacent pair)
    n = len(basis)
    rows = [k[0] for k in H]; cols = [k[1] for k in H]; vals = [H[k] for k in H]
    M = sp.csr_matrix((vals, (rows, cols)), shape=(n, n), dtype=complex)
    return basis, M


def eff_analytic(L, Dl, V, w, mu, N=None, vac=True):
    """Closed-form second-order model (bulk rules; open-chain edges treated exactly by the same bond counting)."""
    basis = configs(L, N) if N is not None else configs(L)
    idx = {s: i for i, s in enumerate(basis)}
    t2 = -Dl ** 2 / (2 * V)
    rows, cols, vals = [], [], []
    for i, s in enumerate(basis):
        # diagonal: sum over bonds of -Dl^2/E(bond), E = V * (number of adjacent bonds after creating the pair)
        d = -mu * popc(s)
        for b in range(L - 1):
            if (s >> b) & 1 or (s >> (b + 1)) & 1:
                continue
            m = s | (1 << b) | (1 << (b + 1))
            d += -Dl ** 2 / (V * nadj(m))
        rows.append(i); cols.append(i); vals.append(d)
        for j in range(L - 1):
            for (a, bsite, amp) in ((j + 1, j, w), (j, j + 1, np.conj(w))):
                if (s >> bsite) & 1 and not (s >> a) & 1:
                    f = (s & ~(1 << bsite)) | (1 << a)
                    if nadj(f) == 0:
                        rows.append(idx[f]); cols.append(i); vals.append(amp)
        for j in range(L - 2):
            for (a, bsite) in ((j + 2, j), (j, j + 2)):
                if (s >> bsite) & 1 and not (s >> a) & 1:
                    f = (s & ~(1 << bsite)) | (1 << a)
                    if nadj(f) == 0:
                        rows.append(idx[f]); cols.append(i); vals.append(t2)
    n = len(basis)
    return basis, sp.csr_matrix((vals, (rows, cols)), shape=(n, n), dtype=complex)


def full_open(L, Dl, V, w, mu):
    """Full 2^L fermion Hamiltonian, open chain."""
    n = 1 << L
    rows, cols, vals = [], [], []
    for s in range(n):
        d = V * nadj(s) - mu * popc(s)
        rows.append(s); cols.append(s); vals.append(d)
        for j in range(L - 1):
            for ops, amp in (([('d', j + 1), ('a', j)], w), ([('d', j), ('a', j + 1)], np.conj(w)),
                             ([('a', j), ('a', j + 1)], -Dl), ([('d', j + 1), ('d', j)], -Dl)):
                r = apply_seq(s, ops)
                if r is not None:
                    rows.append(r[1]); cols.append(s); vals.append(amp * r[0])
    return sp.csr_matrix((vals, (rows, cols)), shape=(n, n), dtype=complex)


if __name__ == "__main__":
    import sys
    # Schrieffer-Wolff check against the full chain
    L = 12
    for V in (20.0, 40.0, 80.0):
        Dl = 1.0
        t2 = -Dl ** 2 / (2 * V)
        w = 0.7j * abs(t2) + 0.3 * abs(t2)
        mu = -3 * Dl ** 2 / V + 0.5 * abs(t2)   # place the band near zero
        Hf = full_open(L, Dl, V, w, mu)
        ef = np.sort(np.linalg.eigvalsh(Hf.toarray()))[:30]
        b, He = eff_rules(L, Dl, V, w, mu)
        ee = np.sort(np.linalg.eigvalsh(He.toarray()))[:30]
        b2, Ha = eff_analytic(L, Dl, V, w, mu)
        ea = np.sort(np.linalg.eigvalsh(Ha.toarray()))[:30]
        print(f"V={V}: |t2|={abs(t2):.4g}  max|E_full - E_rules| (30 lowest) = {np.max(np.abs(ef - ee)):.3e}"
              f"  rel to |t2|: {np.max(np.abs(ef - ee)) / abs(t2):.3e};  max|E_rules - E_analytic| = {np.max(np.abs(ee - ea)):.2e}")
