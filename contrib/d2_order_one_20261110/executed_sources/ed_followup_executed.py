#!/usr/bin/env python3
"""Registered larger-ring ED with momentum resolution of ground multiplets.

Keep raw Lanczos entropies. Resolve translation on the span of all lattice
translates of returned ground vectors; this recovers momentum partners which
single-vector Lanczos can miss. Discard orbit singular values below 1e-6 of
the largest, check projected energy residuals and translation closure, and
report every retained momentum branch. A repeated momentum can retain an
unresolved internal degeneracy, which is explicitly labeled. Select the
lowest momentum index within the ground manifold, never the largest entropy.
"""
import ed_screen as ed  # also limits BLAS threads before NumPy import
import argparse
import hashlib
import json
import os
import time
from pathlib import Path
import numpy as np
from scipy.linalg import schur
from scipy.sparse.linalg import eigsh


def translation_permutation(basis, L):
    inverse = np.full(1 << L, -1, dtype=np.int64)
    inverse[basis] = np.arange(len(basis))
    dest = inverse[(basis >> 1) | ((basis & 1) << (L - 1))]
    assert np.all(dest >= 0)
    return dest


def momentum_ground_basis(H, basis, evals, vecs, residuals, L):
    tolerance = max(1e-8, 20 * max(residuals))
    selected = np.flatnonzero(evals - evals[0] <= tolerance)
    dest = translation_permutation(basis, L)
    def translate(V):
        out = np.empty_like(V)
        out[dest] = V
        return out
    raw_ground = vecs[:, selected]
    raw_Q, _ = np.linalg.qr(raw_ground)
    translated = translate(raw_Q)
    raw_closure = float(np.linalg.norm(translated - raw_Q @ (raw_Q.conj().T @ translated)))
    orbit = []
    v = raw_ground
    for _ in range(L):
        orbit.append(v)
        v = translate(v)
    orbit = np.concatenate(orbit, axis=1)
    gram = orbit.conj().T @ orbit
    eig, U = np.linalg.eigh((gram + gram.conj().T) / 2)
    keep = eig > max(eig[-1] * 1e-12, 1e-13)
    Q = orbit @ (U[:, keep] / np.sqrt(eig[keep]))
    Q, _ = np.linalg.qr(Q)
    del orbit
    TQ = translate(Q)
    small_T = Q.conj().T @ TQ
    closure = float(np.linalg.norm(TQ - Q @ small_T))
    triangular, rotation = schur(small_T, output='complex')
    states = Q @ rotation
    branches = []
    for a in range(states.shape[1]):
        psi = states[:, a]
        psi /= np.linalg.norm(psi)
        tv = np.vdot(psi, translate(psi))
        energy = float(np.vdot(psi, H @ psi).real)
        residual = float(np.linalg.norm(H @ psi - energy * psi))
        phase = float(np.angle(tv))
        momentum = int(round(phase * L / (2 * np.pi))) % L
        branches.append(dict(momentum_index=momentum, momentum_phase=phase,
                             translation_modulus=float(abs(tv)), energy=energy, residual=residual,
                             psi=psi))
    branches.sort(key=lambda b: (b['momentum_index'], b['energy']))
    momenta = [b['momentum_index'] for b in branches]
    return branches, dict(energy_degeneracy_tolerance=tolerance,
                          returned_ground_vectors=len(selected), translation_orbit_rank=len(branches),
                          raw_translation_closure_error=raw_closure,
                          completed_translation_closure_error=closure,
                          translation_schur_offdiagonal=float(np.linalg.norm(triangular - np.diag(np.diag(triangular)))),
                          ground_multiplet_reaches_returned_spectrum_end=bool(len(selected) == len(evals)),
                          repeated_momentum_unresolved=len(set(momenta)) != len(momenta),
                          orbit_singular_relative_cutoff=1e-6)


def solve_sector(L, couplings, parity, periodic=True, distance=1, k=4, tol=1e-11,
                 checkpoint_dir=None, name='point'):
    if not periodic:
        raise ValueError('Momentum follow-up requires PBC')
    if k < 4:
        raise ValueError('Registered follow-up needs at least four eigenvectors')
    start = time.time()
    H, basis = ed.build_sector(L, couplings, parity, periodic, distance)
    rng = np.random.default_rng(20261110 + 997 * L + parity)
    v0 = rng.normal(size=H.shape[0])
    if np.iscomplexobj(H.data):
        v0 = v0 + 1j * rng.normal(size=H.shape[0])
    evals, vecs = eigsh(H, k=k, which='SA', v0=v0, tol=tol, maxiter=100000, ncv=max(4 * k + 1, 25))
    order = np.argsort(evals); evals = evals[order]; vecs = vecs[:, order]
    residuals = [float(np.linalg.norm(H @ vecs[:, i] - evals[i] * vecs[:, i])) for i in range(k)]
    raw = ed.state_metrics(vecs[:, 0], basis, L, periodic)
    branches, diagnostics = momentum_ground_basis(H, basis, evals, vecs, residuals, L)
    checkpoint_states = np.column_stack([b['psi'] for b in branches])
    for branch in branches:
        psi = branch.pop('psi')
        branch.update(ed.state_metrics(psi, basis, L, periodic))
    # state_metrics uses the conjugate translation convention; momentum_phase
    # from the actual T action above remains the stated momentum convention.
    chosen = branches[0]
    good = (diagnostics['completed_translation_closure_error'] < 1e-7
            and all(b['residual'] < 1e-7 and abs(b['translation_modulus'] - 1) < 1e-7 for b in branches)
            and not diagnostics['repeated_momentum_unresolved'])
    result = dict(parity=parity, dimension=H.shape[0], energies=evals.tolist(), residuals=residuals,
                  same_parity_gap=float(evals[1] - evals[0]), raw=raw,
                  momentum_branches=branches, selected_momentum_index=chosen['momentum_index'],
                  momentum_resolution_pass=good, momentum_diagnostics=diagnostics,
                  **{key: value for key, value in chosen.items() if key not in ('energy', 'residual')})
    if checkpoint_dir and L >= 16:
        directory = Path(checkpoint_dir); directory.mkdir(parents=True, exist_ok=True)
        path = directory / f'followup_{name}_L{L}_P{parity}.npz'
        tmp = path.with_name(path.name + '.tmp.npz')
        np.savez_compressed(tmp, basis=basis, eigenvectors=vecs, eigenvalues=evals,
                            momentum_ground_vectors=checkpoint_states,
                            metadata=json.dumps(dict(L=L, parity=parity, couplings=couplings, periodic=periodic, distance=distance)))
        tmp.replace(path); result['checkpoint'] = str(path)
    result['seconds'] = time.time() - start
    return result


def solve_point(L, couplings, periodic=True, sectors=(0, 1), k=4, checkpoint_dir=None,
                distance=1, name='point', tol=1e-11):
    sectors = [solve_sector(L, couplings, p, periodic, distance, k, tol, checkpoint_dir, name) for p in sectors]
    ground = min(sectors, key=lambda r: r['energies'][0])
    levels = sorted((e, r['parity'], i) for r in sectors for i, e in enumerate(r['energies']))
    return dict(name=name, L=L, couplings=couplings, periodic=periodic, distance=distance, k=k, tol=tol,
                requested_sectors=[r['parity'] for r in sectors], target_geometry=distance == 1,
                metrics=ed.coupling_metrics(couplings), sectors=sectors,
                selected_parity=ground['parity'], ground_energy=ground['energies'][0],
                ground_gap=levels[1][0] - levels[0][0], lowest_levels=levels,
                entropy_half=ground['entropy_half'], raw_entropy_half=ground['raw']['entropy_half'],
                c_profile=ground['c_profile'], raw_c_profile=ground['raw']['c_profile'],
                density_mean=ground['density_mean'], selected_momentum_index=ground['selected_momentum_index'],
                momentum_resolution_pass=ground['momentum_resolution_pass'],
                parity_splitting=sectors[1]['energies'][0] - sectors[0]['energies'][0])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--checkpoint-dir', default='/private/tmp/d2-order-one-checkpoints')
    args = parser.parse_args()
    config = json.loads(args.config.read_text())
    done = {}
    for record in ed.read_jsonl_resumable(args.output):
        key = (record['name'], record['L'])
        if key in done:
            raise ValueError(f'Duplicate completed point {key}')
        done[key] = ed.record_signature(record)
    provenance = dict(config_sha256=hashlib.sha256(args.config.read_bytes()).hexdigest(),
                      source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                      ed_screen_sha256=hashlib.sha256(Path(ed.__file__).read_bytes()).hexdigest(),
                      numpy=np.__version__, scipy=ed.scipy.__version__)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    for point in config['points']:
        for L in point.get('sizes', config['sizes']):
            options = {key: point.get(key, config.get(key, default)) for key, default in
                       (('periodic', True), ('sectors', [0, 1]), ('distance', 1), ('k', 4), ('tol', 1e-11))}
            key = (point['name'], L)
            if key in done:
                if done[key] != ed.request_signature(point, L, options):
                    raise ValueError(f'Resume configuration mismatch for {key}')
                continue
            result = solve_point(L, point['couplings'], name=point['name'], checkpoint_dir=args.checkpoint_dir, **options)
            result['provenance'] = provenance
            with args.output.open('a') as handle:
                handle.write(json.dumps(result) + '\n'); handle.flush(); os.fsync(handle.fileno())
            print(json.dumps({key: result[key] for key in ('name', 'L', 'ground_energy', 'ground_gap',
                  'density_mean', 'entropy_half', 'raw_entropy_half', 'c_profile', 'selected_parity',
                  'selected_momentum_index', 'momentum_resolution_pass')}), flush=True)


if __name__ == '__main__':
    main()
