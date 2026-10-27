#!/usr/bin/env python3
"""Closed-book small-port swap gadgets. Requires only numpy and scipy.

Run: python3 gadget.py
The default run checks (d,N)=(2,3),(2,4),(3,4), both signs analytically
at second order, and the antiferromagnetic sign by exact diagonalization.
No input files are read. Results are printed and saved to gadget_results.json.
"""

import os

# Set before numpy is imported, for predictable small-matrix performance.
for _name in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS",
              "VECLIB_MAXIMUM_THREADS"):
    os.environ[_name] = "1"

import argparse
import json
import time
import numpy as np
import scipy.linalg as la
from scipy.sparse.linalg import LinearOperator, eigsh


def traceless_basis(n):
    """Hermitian basis with tr(T_mu T_nu)=delta_mu,nu."""
    out = []
    for j in range(n):
        for k in range(j + 1, n):
            x = np.zeros((n, n), complex)
            x[j, k] = x[k, j] = 1 / np.sqrt(2)
            out.append(x)
            y = np.zeros((n, n), complex)
            y[j, k], y[k, j] = -1j / np.sqrt(2), 1j / np.sqrt(2)
            out.append(y)
    for j in range(1, n):
        z = np.zeros((n, n), complex)
        z[np.arange(j), np.arange(j)] = 1 / np.sqrt(j * (j + 1))
        z[j, j] = -j / np.sqrt(j * (j + 1))
        out.append(z)
    return np.array(out)


def complete_isometry(w):
    """A full unitary whose first columns are the supplied isometry."""
    q, _ = la.qr(w, mode="full")
    q[:, :w.shape[1]] = w
    return q


def sqrt_psd(m):
    vals, vecs = la.eigh(m)
    assert vals.min() > 0
    return (vecs * np.sqrt(vals)) @ vecs.conj().T


def make_effects(d, n, sign=-1):
    """sign=-1 produces AF exchange; sign=+1 produces ferromagnetic."""
    r, m = d * d - 1, n * n - 1
    capacity = r * r - 1
    q = (m + capacity - 1) // capacity
    ell = q + 1
    energies = 3.0 ** np.arange(ell)
    kernel = 1 / (energies[:, None] + energies[None, :])
    # Helmert columns form an orthonormal basis for the sum-zero subspace.
    z = np.zeros((ell, q))
    for j in range(1, ell):
        z[:j, j - 1] = 1 / np.sqrt(j * (j + 1))
        z[j, j - 1] = -j / np.sqrt(j * (j + 1))
    vals, vecs = la.eigh(z.T @ kernel @ z)
    assert vals.min() > 0
    x = z @ ((vecs / np.sqrt(vals)) @ vecs.T)
    syndrome_basis, logical_basis = traceless_basis(r), traceless_basis(n)
    da = np.zeros((ell, r * n, r * n), complex)
    db = np.zeros_like(da)
    for mu in range(m):
        p, j = mu % capacity, mu // capacity
        da += x[:, j, None, None] * np.kron(syndrome_basis[p], logical_basis[mu])
        db += x[:, j, None, None] * np.kron(syndrome_basis[p].T, logical_basis[mu])
    w = np.ones(ell) / ell
    maxnorm = max(np.max(np.abs(la.eigvalsh(v))) for v in [*da, *db])
    t = 0.75 / (ell * maxnorm)
    ma = w[:, None, None] * np.eye(r * n) + t * da
    mb = w[:, None, None] * np.eye(r * n) + sign * t * db
    diagnostics = {
        "effect_min_eigenvalue": float(min(la.eigvalsh(v).min() for v in [*ma, *mb])),
        "effect_sum_error": float(max(la.norm(ma.sum(axis=0) - np.eye(r*n)),
                                       la.norm(mb.sum(axis=0) - np.eye(r*n)))),
        "whitening_error": float(la.norm(x.T @ kernel @ x - np.eye(q))),
        "sum_x_error": float(la.norm(x.sum(axis=0))),
    }
    a2 = -r * (w @ kernel @ w) + sign * t*t / n
    b2 = -sign * t*t
    return ma, mb, energies, float(a2), float(b2), float(t), diagnostics


def realize_block(d, n, effects, energies):
    """Explicit physical H on C^d tensor C^K, and its specified ground basis."""
    ell, rn, _ = effects.shape
    r = d*d - 1
    dim = d * ((n + ell * rn + d - 1) // d)
    aux = dim // d
    assert aux >= d*n
    g = np.zeros((dim, n), complex)
    for a in range(d):
        for i in range(n):
            g[a * aux + a * n + i, i] = 1 / np.sqrt(d)
    ports = np.sqrt(d) * traceless_basis(d)
    s = np.concatenate([np.kron(f, np.eye(aux)) @ g for f in ports], axis=1)
    gs = np.concatenate([g, s], axis=1)
    assert la.norm(gs.conj().T @ gs - np.eye(n + rn)) < 1e-11
    physical_frame = complete_isometry(gs)
    qdim = dim - n
    w = np.zeros((qdim, rn), complex)
    for e in range(ell):
        w[e*rn:(e+1)*rn] = sqrt_psd(effects[e])
    u = complete_isometry(w)
    ev = np.full(qdim, energies[-1])
    ev[:ell*rn] = np.repeat(energies, rn)
    # Columns are the excited eigenvectors in physical coordinates.
    excited_vectors = physical_frame[:, n:] @ u.conj().T
    hlocal = (excited_vectors * ev) @ excited_vectors.conj().T
    hlocal = (hlocal + hlocal.conj().T) / 2
    transitions = [excited_vectors.conj().T @ np.kron(f, np.eye(aux)) @ g
                   for f in ports]
    assert la.norm(hlocal @ g) < 1e-11
    return {"H": hlocal, "G": g, "E": ev, "transitions": transitions,
            "ports": ports, "dim": dim, "aux": aux,
            "unitarity_error": float(la.norm(u.conj().T @ u - np.eye(qdim)))}


def swap(n):
    return np.eye(n*n).reshape(n, n, n, n).transpose(1, 0, 2, 3).reshape(n*n, n*n)


def second_order(a, b):
    coupling = sum(np.kron(x, y) for x, y in zip(a["transitions"], b["transitions"]))
    denominator = (a["E"][:, None] + b["E"][None, :]).ravel()
    return -coupling.conj().T @ (coupling / denominator[:, None])


def exact_low_energy(d, n, a, b, lam, method="iterative"):
    """Diagonalize the full H_A+H_B+lambda*h, not a truncated effective H."""
    dim, aux = a["dim"], a["aux"]
    size = dim * dim
    ha, hb = a["H"], b["H"]

    def mv(v):
        mat = v.reshape(dim, dim)
        # sum_a F_a tensor F_a = d*P_d - I on the two ports.
        portswap = mat.reshape(d, aux, d, aux).transpose(2, 1, 0, 3).reshape(dim, dim)
        return (ha @ mat + mat @ hb.T + lam * (d * portswap - mat)).ravel()

    op = LinearOperator((size, size), matvec=mv, dtype=np.complex128)
    count = n*n + 2
    if method == "dense":
        full = np.column_stack([mv(v) for v in np.eye(size)])
        values, vectors = la.eigh(full, subset_by_index=(0, count-1), driver="evr")
    else:
        rng = np.random.default_rng(2000 + 100*d + n)
        v0 = rng.normal(size=size) + 1j*rng.normal(size=size)
        values, vectors = eigsh(op, k=count, which="SA", v0=v0,
                               ncv=min(size, max(100, 5*count)), tol=3e-14, maxiter=10000)
        order = np.argsort(values)
        values, vectors = values[order], vectors[:, order]
    # Complex ARPACK can lose mutual orthogonality in these tight clusters.
    # Rayleigh-Ritz in the returned subspace removes that numerical ambiguity.
    raw_gram_error = la.norm(vectors.conj().T @ vectors - np.eye(count))
    basis, _ = la.qr(vectors, mode="economic")
    h_basis = np.column_stack([mv(basis[:, j]) for j in range(count)])
    rayleigh = basis.conj().T @ h_basis
    values, rotation = la.eigh((rayleigh + rayleigh.conj().T) / 2)
    vectors = basis @ rotation
    gram_error = la.norm(vectors.conj().T @ vectors - np.eye(count))
    residual = max(la.norm(mv(vectors[:, j]) - values[j]*vectors[:, j]) for j in range(count))
    isolation_gap = values[n*n] - values[n*n-1]
    values, vectors = values[:n*n], vectors[:, :n*n]
    ground = np.kron(a["G"], b["G"])
    overlap = ground.conj().T @ vectors
    left, singular, right = la.svd(overlap, full_matrices=False)
    polar = left @ right
    heff = (polar * values) @ polar.conj().T
    heff = (heff + heff.conj().T) / 2
    p = swap(n)
    f = p - np.eye(n*n)/n
    bfit = float(np.trace(f @ heff).real / (n*n - 1))
    afit = float(np.trace(heff).real / (n*n) - bfit/n)
    fitres = float(la.norm(heff - afit*np.eye(n*n) - bfit*p, "fro"))
    # AF: the lower n(n-1)/2 states should approach the antisymmetric space.
    na = n*(n-1)//2
    swap_expectations = np.diag(polar.conj().T @ p @ polar).real
    pa_exact = polar[:, :na] @ polar[:, :na].conj().T
    pa = (np.eye(n*n) - p)/2
    anti, sym = values[:na], values[na:]
    return {
        "lambda": lam, "a_fit": afit, "b_fit": bfit,
        "b_over_lambda2": bfit / lam**2,
        "fit_residual_fro": fitres, "fit_residual_over_lambda3": fitres / lam**3,
        "anti_width": float(np.ptp(anti)), "sym_width": float(np.ptp(sym)),
        "mean_splitting": float(sym.mean() - anti.mean()),
        "anti_swap_mean": float(swap_expectations[:na].mean()),
        "sym_swap_mean": float(swap_expectations[na:].mean()),
        "anti_projector_error_fro": float(la.norm(pa_exact - pa, "fro")),
        "overlap_min_singular": float(singular.min()),
        "isolation_gap": float(isolation_gap), "eigen_residual": float(residual),
        "raw_eigenvector_gram_error": float(raw_gram_error),
        "eigenvector_gram_error": float(gram_error),
        "levels": [float(v) for v in values],
    }, heff


def run_case(d, n, lambdas, method):
    started = time.monotonic()
    print(f"\nCASE d={d}, N={n}", flush=True)
    ma, mb, en, a2, b2, t, diagnostics = make_effects(d, n)
    a, b = realize_block(d, n, ma, en), realize_block(d, n, mb, en)
    target = a2*np.eye(n*n) + b2*swap(n)
    h2 = second_order(a, b)
    h2error = float(la.norm(h2 - target, "fro"))
    h1 = sum(np.kron(a["G"].conj().T @ np.kron(f, np.eye(a["aux"])) @ a["G"],
                     b["G"].conj().T @ np.kron(f, np.eye(b["aux"])) @ b["G"])
             for f in a["ports"])
    diagnostics.update({"h1_norm": float(la.norm(h1)), "h2_error_fro": h2error,
                        "unitarity_error_A": a["unitarity_error"],
                        "unitarity_error_B": b["unitarity_error"]})
    assert h2error < 1e-10
    # The opposite sign has the same positivity bound and physical construction.
    maf, mbf, _, a2f, b2f, _, _ = make_effects(d, n, sign=1)
    bf = realize_block(d, n, mbf, en)
    ferro_error = float(la.norm(second_order(a, bf) - a2f*np.eye(n*n) - b2f*swap(n)))
    assert ferro_error < 1e-10
    diagnostics["ferro_h2_error_fro"] = ferro_error
    record = {"d": d, "N": n, "local_dimension": a["dim"], "aux_dimension": a["aux"],
              "energies": en.tolist(), "t": t, "a2": a2, "b2": b2,
              "diagnostics": diagnostics, "exact": []}
    print(f"D={a['dim']}, K={a['aux']}, energies={en}, a2={a2:.12g}, b2={b2:.12g}", flush=True)
    print(f"H1 norm={la.norm(h1):.3e}; H2 error={h2error:.3e}; ferro error={ferro_error:.3e}", flush=True)
    for lam in lambdas:
        row, heff = exact_low_energy(d, n, a, b, lam, method)
        row["second_order_error_fro"] = float(la.norm(heff - lam**2*target, "fro"))
        row["relative_b_error"] = float(abs(row["b_over_lambda2"]/b2 - 1))
        record["exact"].append(row)
        print(f"lambda={lam:.5g}: b/lambda^2={row['b_over_lambda2']:.10g}; "
              f"fitres={row['fit_residual_fro']:.3e}; "
              f"widths=({row['anti_width']:.3e},{row['sym_width']:.3e}); "
              f"Panti error={row['anti_projector_error_fro']:.3e}; "
              f"eigres={row['eigen_residual']:.2e}", flush=True)
        assert row["isolation_gap"] > .5
        assert row["eigen_residual"] < 1e-10
    if len(record["exact"]) >= 2:
        for key in ("mean_splitting", "fit_residual_fro", "second_order_error_fro",
                    "anti_width", "sym_width", "anti_projector_error_fro"):
            rows = record["exact"]
            record[key + "_log_slope"] = float(np.polyfit(
                np.log([v["lambda"] for v in rows]), np.log([v[key] for v in rows]), 1)[0])
        print("slopes: splitting={:.5f}, residual={:.5f}, H2 error={:.5f}".format(
            record["mean_splitting_log_slope"], record["fit_residual_fro_log_slope"],
            record["second_order_error_fro_log_slope"]), flush=True)
    record["runtime_seconds"] = time.monotonic() - started
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cases", default="2,3;2,4;3,4")
    parser.add_argument("--lambdas", default="0.02,0.01,0.005,0.0025,0.00125")
    parser.add_argument("--method", choices=("iterative", "dense"), default="iterative")
    parser.add_argument("--output", default="gadget_results.json")
    args = parser.parse_args()
    started = time.monotonic()
    cases = [tuple(map(int, v.split(","))) for v in args.cases.split(";")]
    lambdas = [float(v) for v in args.lambdas.split(",")]
    records = [run_case(d, n, lambdas, args.method) for d, n in cases]
    result = {"cases": records, "total_runtime_seconds": time.monotonic()-started,
              "convention": "H(lambda)=H_A+H_B+lambda*h; canonical/des Cloizeaux H_eff"}
    with open(args.output, "w") as f:
        json.dump(result, f, indent=2)
    print(f"\nTotal runtime {result['total_runtime_seconds']:.2f} s; wrote {args.output}", flush=True)


if __name__ == "__main__":
    main()
