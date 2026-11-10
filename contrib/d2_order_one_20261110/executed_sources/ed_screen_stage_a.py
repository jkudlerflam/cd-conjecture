#!/usr/bin/env python3
"""Independent full-spin ED for the order-one-coupling search (2026-11-10).

No particle-number or blockade truncation is made. Site j is integer bit L-1-j,
so vectors reshape in the usual tensor-product order. Pauli normalization is
H=sum[Jx XX+Jy YY+Dz(XY-YX)+Jz ZZ]+hz sum Z. The optional distance=2 is
ONLY an estimator control (two interleaved XX chains), never a target model.

JSON config: {"points":[{"name":"...", "couplings":{...}}], "sizes":[10,12],
"sectors":[0,1], "periodic":true, "k":2}. Point values override global values
for sizes, periodic, distance, sectors and k. Each completed point/size is
atomically appended to JSONL. Existing (name,L) rows are skipped on resume.
"""
import os
for _key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_key] = "1"

import argparse
import hashlib
import json
import time
from pathlib import Path
import numpy as np
import scipy
from scipy import sparse
from scipy.sparse.linalg import eigsh


def coupling_metrics(p):
    """Literal Pauli coefficients plus rotation-invariant exchange diagnostics."""
    x, y, z, d, h = (float(p.get(k, 0)) for k in ("Jx", "Jy", "Jz", "Dz", "hz"))
    # XY and YX are distinct Pauli strings with equal coefficient magnitudes.
    values = [abs(v) for v in (x, y, z, d, -d, h) if v != 0]
    exchange = [abs(v) for v in (x, y, z, d, -d) if v != 0]
    J = np.array([[x, d, 0], [-d, y, 0], [0, 0, z]])
    singular = np.linalg.svd(J, compute_uv=False)
    nzsingular = singular[singular > 1e-13 * max(1, singular.max())]
    sym = np.linalg.norm((J + J.T) / 2)
    asym = np.linalg.norm((J - J.T) / 2)
    return dict(pauli_string_count=len(values), grouped_term_count=sum(v != 0 for v in (x, y, z, d, h)),
                raw_coefficient_spread=max(values) / min(values) if values else None,
                exchange_coefficient_spread=max(exchange) / min(exchange) if exchange else None,
                exchange_singular_values=singular.tolist(),
                nonzero_singular_spread=float(nzsingular.max() / nzsingular.min()) if len(nzsingular) else None,
                exchange_symmetric_frobenius=float(sym), exchange_antisymmetric_frobenius=float(asym),
                antisymmetric_to_symmetric=float(asym / sym) if sym else None,
                field_norm=abs(h), field_to_exchange_frobenius=float(abs(h) / np.linalg.norm(J)) if np.linalg.norm(J) else None,
                coefficient_convention="Pauli strings; full on-site hz counted once; absent terms excluded")


def sector_basis(L, parity):
    full = np.arange(1 << L, dtype=np.int64)
    par = np.zeros(len(full), dtype=np.int8)
    for j in range(L):
        par ^= ((full >> j) & 1).astype(np.int8)
    return full[par == parity]


def build_sector(L, couplings, parity=0, periodic=True, distance=1):
    """Assemble CSR independently by Pauli actions on bit strings."""
    if L < 3 or distance < 1 or 2 * distance >= L:
        raise ValueError("Need L>2*distance to count every undirected bond once")
    if parity not in (0, 1):
        raise ValueError("parity must be 0 or 1")
    p = {k: float(couplings.get(k, 0)) for k in ("Jx", "Jy", "Jz", "Dz", "hz")}
    unknown = set(couplings) - set(p)
    if unknown:
        raise ValueError(f"Unsupported couplings {unknown}")
    basis = sector_basis(L, parity)
    dim = len(basis)
    inverse = np.full(1 << L, -1, dtype=np.int64)
    inverse[basis] = np.arange(dim)
    z = np.array([1 - 2 * ((basis >> (L - 1 - j)) & 1) for j in range(L)], dtype=np.int8)
    diag = p["hz"] * z.sum(axis=0).astype(float)
    rows, cols, data = [], [], []
    source = np.arange(dim)
    bonds = [(j, (j + distance) % L) for j in range(L if periodic else L - distance)]
    for j, k in bonds:
        diag += p["Jz"] * z[j] * z[k]
        amp = p["Jx"] - p["Jy"] * z[j] * z[k] + 1j * p["Dz"] * (z[k] - z[j])
        dest = inverse[basis ^ (1 << (L - 1 - j)) ^ (1 << (L - 1 - k))]
        keep = amp != 0
        rows.append(dest[keep]); cols.append(source[keep]); data.append(amp[keep])
    rows.append(source); cols.append(source); data.append(diag)
    H = sparse.coo_matrix((np.concatenate(data), (np.concatenate(rows), np.concatenate(cols))), shape=(dim, dim)).tocsr()
    if not p["Dz"]:
        H = H.real
    H.eliminate_zeros()
    return H, basis


def entropy_profile(psi, L):
    vals = []
    for cut in range(1, L):
        s = np.linalg.svd(psi.reshape(1 << cut, -1), compute_uv=False)
        prob = s * s
        prob = prob[prob > 1e-15]
        vals.append(float(-np.dot(prob, np.log(prob))))
    return np.array(vals)


def state_metrics(psi_sector, basis, L, periodic=True):
    psi = np.zeros(1 << L, dtype=psi_sector.dtype)
    psi[basis] = psi_sector
    prob = np.abs(psi_sector) ** 2
    bits = np.array([((basis >> (L - 1 - j)) & 1) for j in range(L)], dtype=float)
    n = bits @ prob
    nn = (bits * prob) @ bits.T
    connected_z = 4 * (nn - np.outer(n, n))
    qs = 2 * np.pi * np.arange(L) / L
    fourier = np.exp(-1j * np.outer(qs, np.arange(L)))
    structure = np.einsum('qi,ij,qj->q', fourier, connected_z, fourier.conj()).real / L
    S = entropy_profile(psi, L)
    cuts = np.arange(1, L)
    central = (cuts >= L / 4) & (cuts <= 3 * L / 4)
    xx = np.log((L / np.pi) * np.sin(np.pi * cuts[central] / L))
    fit = np.polyfit(xx, S[central], 1)
    fit_error = S[central] - np.polyval(fit, xx)
    N = bits.sum(axis=0)
    # Translation moves the last site to the first (orientation immaterial for modulus).
    shifted = (basis >> 1) | ((basis & 1) << (L - 1))
    translation = np.vdot(psi[shifted], psi[basis])
    return dict(density_mean=float(n.mean()), density_profile=n.tolist(),
                density_profile_std=float(n.std()), density_profile_range=float(np.ptp(n)),
                N_mean=float(prob @ N), N_variance=float(prob @ (N * N) - (prob @ N) ** 2),
                entropy_profile=S.tolist(), entropy_half=float(S[L // 2 - 1]),
                entropy_central_mean=float(np.mean(S[central])),
                c_profile=float((3 if periodic else 6) * fit[0]),
                entropy_fit_rms=float(np.sqrt(np.mean(fit_error ** 2))),
                structure_q=qs.tolist(), connected_Szz=structure.tolist(),
                connected_Z_correlation=connected_z[0].tolist(),
                translation_expectation=[float(translation.real), float(translation.imag)],
                translation_modulus=float(abs(translation)))


def solve_sector(L, couplings, parity, periodic=True, distance=1, k=2, tol=1e-11, checkpoint_dir=None, name="point"):
    started = time.time()
    H, basis = build_sector(L, couplings, parity, periodic, distance)
    rng = np.random.default_rng(20261110 + 997 * L + parity)
    v0 = rng.normal(size=H.shape[0])
    if np.iscomplexobj(H.data):
        v0 = v0 + 1j * rng.normal(size=H.shape[0])
    evals, vecs = eigsh(H, k=k, which="SA", v0=v0, tol=tol, maxiter=100000)
    order = np.argsort(evals)
    evals, vecs = evals[order], vecs[:, order]
    residuals = [float(np.linalg.norm(H @ vecs[:, a] - evals[a] * vecs[:, a])) for a in range(k)]
    hermiticity = H - H.getH()
    result = dict(parity=parity, dimension=H.shape[0], energies=evals.tolist(), residuals=residuals,
                  hermiticity_error=float(np.max(np.abs(hermiticity.data))) if hermiticity.nnz else 0.0,
                  same_parity_gap=float(evals[1] - evals[0]) if k > 1 else None,
                  ground_degenerate_warning=bool(k > 1 and evals[1] - evals[0] < 1e-8),
                  **state_metrics(vecs[:, 0], basis, L, periodic))
    if checkpoint_dir and L >= 16:
        directory = Path(checkpoint_dir); directory.mkdir(parents=True, exist_ok=True)
        path = directory / f"{name}_L{L}_P{parity}.npz"
        temp = path.with_name(path.name + ".tmp.npz")
        np.savez_compressed(temp, basis=basis, eigenvectors=vecs, eigenvalues=evals,
                            metadata=json.dumps(dict(L=L, parity=parity, couplings=couplings, periodic=periodic, distance=distance)))
        temp.replace(path)
        result["checkpoint"] = str(path)
    result["seconds"] = time.time() - started
    return result


def solve_point(L, couplings, periodic=True, sectors=(0, 1), k=2, checkpoint_dir=None, distance=1, name="point", tol=1e-11):
    results = [solve_sector(L, couplings, p, periodic, distance, k, tol, checkpoint_dir, name) for p in sectors]
    selected = min(results, key=lambda r: r["energies"][0])
    levels = sorted((e, r["parity"], i) for r in results for i, e in enumerate(r["energies"]))
    return dict(name=name, L=L, couplings=couplings, periodic=periodic, distance=distance,
                target_geometry=distance == 1, metrics=coupling_metrics(couplings), sectors=results,
                selected_parity=selected["parity"], ground_energy=selected["energies"][0],
                ground_gap=levels[1][0] - levels[0][0], lowest_levels=levels,
                entropy_half=selected["entropy_half"], c_profile=selected["c_profile"],
                density_mean=selected["density_mean"],
                parity_splitting=(results[1]["energies"][0] - results[0]["energies"][0]) if tuple(sectors) == (0, 1) else None)


def dense_reference(L, couplings, periodic=True, distance=1):
    """Independent Kronecker assembly used only for tiny-L validation."""
    I = np.eye(2)
    X = np.array([[0, 1], [1, 0]])
    Y = np.array([[0, -1j], [1j, 0]])
    Z = np.diag([1, -1])
    def op(items):
        out = np.array([[1.]])
        for j in range(L):
            out = np.kron(out, items.get(j, I))
        return out
    H = np.zeros((1 << L, 1 << L), dtype=complex)
    for j in range(L if periodic else L - distance):
        k = (j + distance) % L
        for a, b, v in ((X, X, couplings.get("Jx", 0)), (Y, Y, couplings.get("Jy", 0)),
                        (Z, Z, couplings.get("Jz", 0)), (X, Y, couplings.get("Dz", 0)),
                        (Y, X, -couplings.get("Dz", 0))):
            H += v * op({j: a, k: b})
    for j in range(L):
        H += couplings.get("hz", 0) * op({j: Z})
    return H


def validate():
    """Assembly check is independent of physical screen selection."""
    p = dict(Jx=.5, Jy=-.5, Jz=.75, Dz=.125, hz=-1.25)
    checks = []
    for periodic in (False, True):
        dense = dense_reference(6, p, periodic)
        for parity in (0, 1):
            H, basis = build_sector(6, p, parity, periodic)
            error = float(np.max(np.abs(H.toarray() - dense[np.ix_(basis, basis)])))
            assert error < 1e-13, error
            checks.append(dict(L=6, periodic=periodic, parity=parity, max_matrix_error=error))
    return dict(assembly_checks=checks, numpy=np.__version__, scipy=scipy.__version__,
                source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--validate", action="store_true")
    parser.add_argument("--checkpoint-dir", default="/private/tmp/d2-order-one-checkpoints")
    args = parser.parse_args()
    if args.validate:
        print(json.dumps(validate(), indent=2), flush=True)
        return
    if not args.config or not args.output:
        parser.error("--config and --output are required except with --validate")
    config = json.loads(args.config.read_text())
    done = set()
    if args.output.exists():
        for line in args.output.read_text().splitlines():
            r = json.loads(line)
            done.add((r["name"], r["L"]))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    provenance = dict(config_sha256=hashlib.sha256(args.config.read_bytes()).hexdigest(),
                      source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                      numpy=np.__version__, scipy=scipy.__version__)
    for point in config["points"]:
        for L in point.get("sizes", config["sizes"]):
            if (point["name"], L) in done:
                continue
            options = {key: point.get(key, config.get(key, default)) for key, default in
                       (("periodic", True), ("sectors", [0, 1]), ("distance", 1), ("k", 2), ("tol", 1e-11))}
            result = solve_point(L, point["couplings"], name=point["name"], checkpoint_dir=args.checkpoint_dir, **options)
            result["provenance"] = provenance
            # One line per completed point. flush + fsync preserves completed work.
            with args.output.open("a") as handle:
                handle.write(json.dumps(result) + "\n"); handle.flush(); os.fsync(handle.fileno())
            print(json.dumps({key: result[key] for key in ("name", "L", "ground_energy", "ground_gap", "density_mean", "entropy_half", "c_profile", "selected_parity")}), flush=True)


if __name__ == "__main__":
    main()
