"""Entropy curves s(T) of the benchmark chains (Ising, XX exactly in the thermodynamic limit;
XXZ, Potts, Takhtajan-Babujian by sector-resolved exact diagonalization on rings), and the
conformal window T*_eps. Output: data/bench_entropy.npz and a printed table.

Usage: python3 bench_thermo.py [Lmax_xxz] [Lmax_potts] [Lmax_tb]
"""
import sys
import numpy as np
import scipy.sparse as sp
from models_ed import models
from thermo import entropy_from_levels, ff_entropy

Ts = np.round(np.concatenate([np.arange(0.05, 1.0, 0.025), np.arange(1.0, 4.01, 0.1)]), 4)


def basis_with_charge(d, L, q, Q, mod):
    idx = np.arange(d ** L)
    digits = np.array([(idx // d ** (L - 1 - i)) % d for i in range(L)]).T
    tot = q[digits].sum(1)
    keep = (tot % mod == Q % mod) if mod else (tot == Q)
    return digits[keep]


def sector_H(h, d, L, digits):
    """Sparse H on the basis 'digits' (rows = configurations), periodic ring."""
    pw = d ** np.arange(L - 1, -1, -1)
    codes = digits @ pw
    order = np.argsort(codes); cs = codes[order]
    D = len(digits)
    rows, cols, vals = [], [], []
    hr = h.real if np.allclose(h.imag, 0) else h
    for i in range(L):
        j = (i + 1) % L
        ab = digits[:, i] * d + digits[:, j]
        for out in range(d * d):
            amp = hr[out, ab]
            nz = np.abs(amp) > 1e-14
            if not nz.any():
                continue
            nd = digits[nz].copy()
            nd[:, i] = out // d; nd[:, j] = out % d
            nc = nd @ pw
            pos = np.searchsorted(cs, nc)
            pos = np.minimum(pos, D - 1)
            ok = cs[pos] == nc
            rows.append(order[pos[ok]]); cols.append(np.nonzero(nz)[0][ok]); vals.append(amp[nz][ok])
    if not vals:
        return sp.csr_matrix((D, D))
    return sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(D, D))


def ed_entropy(h, d, L, q, mod):
    E = []
    Qs = range(mod) if mod else range(int(q.min()) * L, int(q.max()) * L + 1)
    for Q in Qs:
        dig = basis_with_charge(d, L, q, Q, mod)
        if len(dig) == 0:
            continue
        H = sector_H(h, d, L, dig).toarray()
        assert np.allclose(H, H.conj().T)
        E.append(np.linalg.eigvalsh(H))
    E = np.concatenate(E)
    assert len(E) == d ** L
    return entropy_from_levels(E, np.ones_like(E), Ts) / L, E.min() / L


def potts_tau_basis():
    """Potts bond in the basis where tau is diagonal (real symmetric), charge n mod 3."""
    h, d, c, v = models()["Potts3"]
    F = np.array([[np.exp(2j * np.pi * a * b / 3) for b in range(3)] for a in range(3)]) / np.sqrt(3)
    tau = np.roll(np.eye(3), 1, axis=0)
    ev, U = np.linalg.eig(tau)  # diagonalize tau
    U = np.linalg.qr(U)[0]
    UU = np.kron(U, U)
    hb = UU.conj().T @ h @ UU
    # charge of each tau eigenvector: eigenvalue omega^n
    n = np.round(np.angle(np.diag(U.conj().T @ tau @ U)) / (2 * np.pi / 3)).astype(int) % 3
    return hb, n


def main():
    Lx = int(sys.argv[1]) if len(sys.argv) > 1 else 14
    Lp = int(sys.argv[2]) if len(sys.argv) > 2 else 9
    Lt = int(sys.argv[3]) if len(sys.argv) > 3 else 9
    M = models()
    res = {"T": Ts}
    k = (np.arange(200000) + 0.5) * 2 * np.pi / 200000 - np.pi
    res["Ising_inf"] = ff_entropy(4 * np.abs(np.sin(k / 2)), Ts)  # one Bogoliubov mode per k
    res["XX_inf"] = ff_entropy(2 * np.cos(k), Ts)
    h, d, c, v = M["XXZ_D0.5"]
    for L in (Lx - 2, Lx):
        res[f"XXZ_L{L}"], _ = ed_entropy(h, 2, L, np.array([0, 1]), 0)
    hb, n = potts_tau_basis()
    for L in (Lp - 1, Lp):
        res[f"Potts_L{L}"], _ = ed_entropy(hb, 3, L, n, 3)
    h, d, c, v = M["TB_spin1"]
    for L in (Lt - 1, Lt):
        res[f"TB_L{L}"], _ = ed_entropy(h, 3, L, np.array([1, 0, -1]), 0)
    np.savez("data/bench_entropy.npz", **res)
    for key in res:
        if key != "T":
            print(key, " ".join(f"{x:.4f}" for x in res[key][::8]))


if __name__ == "__main__":
    main()
