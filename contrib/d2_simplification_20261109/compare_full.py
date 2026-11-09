"""Compare rational simplifications of C183's full V=16 qubit bond.

Uses the project's momentum-sector assembler (C183/est.py) with a new
literal bond constructor and exact-diagonalization comparison.  The
small-ring data are diagnostics, not estimates of a thermodynamic c.
"""

import json
from pathlib import Path
import sys
import time

import numpy as np
from scipy.sparse.linalg import eigsh
from scipy.optimize import minimize_scalar

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "codex/work_C183"))
import est

I = np.eye(2)
X = np.array([[0, 1], [1, 0]], complex)
Y = np.array([[0, -1j], [1j, 0]], complex)
Z = np.diag([1, -1])


def bond(c):
    return (c["Jx"] * np.kron(X, X) + c["Jy"] * np.kron(Y, Y)
            + c["Jz"] * np.kron(Z, Z)
            + c["Dz"] * (np.kron(X, Y) - np.kron(Y, X))
            + c["Dy"] * (np.kron(X, Z) - np.kron(Z, X))
            + c["hz"] / 2 * (np.kron(Z, I) + np.kron(I, Z)))


def choices():
    base = dict(Jx=.4979983, Jy=-.5020017, Jz=4., Dz=.0152122, Dy=.13693064)
    binary = dict(Jx=255/512, Jy=-257/512, Jz=4., Dz=1/64, Dy=9/64)
    no_real = dict(base, Jx=.5, Jy=-.5)
    w_comp = -(2 * base["Dy"]) ** 2 / 16
    compensation_only = dict(base, Jx=(1+w_comp)/2, Jy=(-1+w_comp)/2)
    return {"C183": base, "binary": binary, "zero_real_hop": no_real,
            "compensation_only": compensation_only}


def lower(h, L, m, count=2):
    table = est.table(L)
    mat = table.sector_matrix(h, m)
    if mat.shape[0] <= 400:
        values, vectors = np.linalg.eigh(mat.toarray())
        values, vectors = values[:count], vectors[:, :count]
    else:
        values, vectors = eigsh(mat, k=count, which="SA", tol=2e-11)
        order = np.argsort(values)
        values, vectors = values[order], vectors[:, order]
    particles = np.array([int(s).bit_count() for s in table.sector(m)["R"]])
    density = (np.abs(vectors) ** 2 * particles[:, None]).sum(axis=0) / L
    return values, density


def harmonic_basis(k):
    return np.array([np.ones_like(k), np.sin(k), np.cos(k), np.sin(2*k),
                     np.cos(2*k), np.sin(3*k), np.cos(3*k), np.sin(4*k), np.cos(4*k)])


def far_band(c, L):
    h = bond(dict(c, hz=-8.09375))
    eps = np.empty(L)
    vacuum_energy = None
    vacuum_density = None
    for m in range(L):
        e, n = lower(h, L, m)
        eps[m] = e[0] if m else e[1]
        if not m:
            vacuum_energy, vacuum_density = e[0], n[0]
    eps -= vacuum_energy
    k = np.arange(L)*2*np.pi/L
    coef = np.linalg.lstsq(harmonic_basis(k).T, eps, rcond=None)[0]
    fun = lambda x: float(coef @ harmonic_basis(np.array(x)))
    minima = [minimize_scalar(fun, bounds=window, method="bounded")
              for window in ((-1., .5), (2.4, 4.2))]
    t2 = -coef[4]/2
    mismatch = float(np.angle(np.exp(2j*(minima[0].x+minima[1].x))))
    return dict(L=L, hz=-8.09375, energies=eps.tolist(), vacuum_density=float(vacuum_density),
                fit_coefficients=coef.tolist(), fit_residual=float(np.max(np.abs(harmonic_basis(k).T@coef-eps))),
                t2=float(t2), valleys=[dict(k=float(v.x), energy=float(v.fun)) for v in minima],
                signed_split_in_t2=float((minima[1].fun-minima[0].fun)/t2),
                quartet_momentum_mismatch=mismatch)


def density_point(c, L, hz):
    h = bond(dict(c, hz=hz))
    records = []
    for m in range(L):
        e, n = lower(h, L, m)
        records.append(dict(momentum_index=m, energy=float(e[0]), density=float(n[0])))
    winner = min(records, key=lambda x: x["energy"])
    return dict(L=L, hz=hz, ground=winner, sectors=records)


def run():
    start = time.monotonic()
    output = {"scope": "full nearest-neighbour qubit rings; no c estimate", "candidates": []}
    for name, c in choices().items():
        result = dict(name=name, couplings=c, far_bands=[], density_ladder=[])
        for L in (12, 14):
            row = far_band(c, L)
            result["far_bands"].append(row)
            print(name, L, "far-band split/t2", row["signed_split_in_t2"], "quartet mismatch", row["quartet_momentum_mismatch"], flush=True)
        for hz in (-7.950, -7.948, -7.9463, -7.944, -7.942):
            row = density_point(c, 12, hz)
            result["density_ladder"].append(row)
            print(name, "hz", hz, "L12 density", row["ground"]["density"], flush=True)
        output["candidates"].append(result)
    output["seconds"] = time.monotonic()-start
    return output


if __name__ == "__main__":
    result = run()
    target = Path(__file__).with_name("full_comparison.json")
    target.write_text(json.dumps(result, indent=2)+"\n")
    print(target)
