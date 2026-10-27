"""Closed-book port-idle FK experiment. Run with python3 experiment.py.
All BLAS thread counts are set before importing numpy/scipy.
"""
import os
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
            'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS', 'BLIS_NUM_THREADS'):
    os.environ[key] = '1'
import json
import numpy as np
from scipy.linalg import eigh, qr

SEED = 20261027
rng = np.random.default_rng(SEED)
I2 = np.eye(2, dtype=complex)
pauli = [np.array([[0, 1], [1, 0]], complex),
         np.array([[0, -1j], [1j, 0]], complex),
         np.diag([1, -1]).astype(complex)]
omega = np.array([1, 0, 0, 1], complex) / np.sqrt(2)
syndrome = np.column_stack([np.kron(f, I2) @ omega for f in pauli])
Qom = np.eye(4) - np.outer(omega, omega.conj())
T, L, Dw, N = 5, 6, 4, 2
D = 4 * Dw
# Workspace ordering: logical qubit, initialized ancilla qubit.
Ew = np.kron(I2, np.array([[1], [0]], complex))
Ein = np.kron(np.eye(4), Ew)
Ginput = np.kron(omega[:, None], Ew)
Gtarget = np.kron(omega[:, None], I2)
Pbad = np.kron(np.eye(4), np.kron(I2, np.diag([0, 1])))
Ffull = [np.kron(np.kron(f, I2), np.eye(Dw)) for f in pauli]
Ftarget = [np.kron(np.kron(f, I2), I2) for f in pauli]


def opnorm(a):
    return float(np.linalg.norm(a, 2))


def haar(d):
    a = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
    q, r = qr(a)
    return q @ np.diag(np.diag(r) / np.abs(np.diag(r)))


def positive(d, floor):
    u = haar(d)
    return (u * np.linspace(floor, 1.0, d)) @ u.conj().T


Ugates = [haar(Dw) for _ in range(T)]
# Each h_t acts on p,s and workspace's logical qubit c, never on the ancilla.
hidle = []
for t in range(L):
    k = syndrome @ positive(3, 0.25) @ syndrome.conj().T
    y = positive(2, 0.2)
    hidle.append(np.kron(k, np.kron(y, I2)))


def construct(active=False):
    """The active model has one Z gate on p and conjugated later penalties.
    Its effective target and exact ground code agree with the idle target.
    """
    Ws = [np.eye(D, dtype=complex)]
    gates = []
    hs = []
    zport = np.kron(np.kron(pauli[2], I2), np.eye(Dw))
    for t, u in enumerate(Ugates, start=1):
        gate = np.kron(np.eye(4), u)
        if active and t == 3:
            gate = zport @ gate
        gates.append(gate)
        Ws.append(gate @ Ws[-1])
    for t, h in enumerate(hidle):
        hs.append(zport @ h @ zport if active and t >= 3 else h)
    H0 = np.zeros((L * D, L * D), complex)
    for t, u in enumerate(gates, start=1):
        st, sp = slice(t * D, (t + 1) * D), slice((t - 1) * D, t * D)
        H0[st, st] += np.eye(D) / 2
        H0[sp, sp] += np.eye(D) / 2
        H0[st, sp] -= u / 2
        H0[sp, st] -= u.conj().T / 2
    H0[:D, :D] += Pbad
    V = np.zeros_like(H0)
    for t, h in enumerate(hs):
        V[t * D:(t + 1) * D, t * D:(t + 1) * D] = L * h
    eta = np.vstack([w @ Ein for w in Ws]) / np.sqrt(L)
    G = np.vstack([w @ Ginput for w in Ws]) / np.sqrt(L)
    HA = sum(Ein.conj().T @ w.conj().T @ h @ w @ Ein
             for w, h in zip(Ws, hs))
    Fs = [np.kron(np.eye(L), f) for f in Ffull]
    source = np.concatenate([f @ G for f in Fs], axis=1)
    return H0, V, eta, G, HA, Fs, source, Ws, hs


def boundary(e, u, source, z):
    x = u.conj().T @ source
    return x.conj().T @ (x / (e + z)[:, None])


def exchange(e, u, Fs, G):
    # All excited states are retained; physical ground states are omitted.
    keep = e > 1e-7
    ee, ue = e[keep], u[:, keep]
    x = np.stack([ue.conj().T @ f @ G for f in Fs], axis=1)
    v = np.einsum('eai,faj->efij', x, x).reshape(len(ee)**2, N*N)
    denom = (ee[:, None] + ee[None, :]).reshape(-1)
    return -v.conj().T @ (v / denom[:, None])


zs = [0.05, 0.2, 1.0, 3.0]
Js = [10.0**k for k in range(1, 7)]
models = {name: construct(active) for name, active in [('idle', False), ('active', True)]}
HA = models['idle'][4]
ea, ua = eigh(HA)
sourceA = np.concatenate([f @ Gtarget for f in Ftarget], axis=1)
RA = {z: boundary(ea, ua, sourceA, z) for z in zs}
KA = exchange(ea, ua, Ftarget, Gtarget)
Delta = float(ea[ea > 1e-8][0])
# Independent target check: diagonalize the 64-dimensional two-cell Hamiltonian.
pairH = np.kron(HA, np.eye(len(HA))) + np.kron(np.eye(len(HA)), HA)
pairv = sum(np.kron(f, f) for f in Ftarget)
pairsource = pairv @ np.kron(Gtarget, Gtarget)
ep, up = eigh(pairH)
pairkeep = ep > 1e-8
pairx = up[:, pairkeep].conj().T @ pairsource
Kdirect = -pairx.conj().T @ (pairx / ep[pairkeep, None])
assert opnorm(Kdirect-KA) < 1e-10
ktr = KA.reshape(N, N, N, N)
local_left = np.einsum('ijkj->ik', ktr) / N
local_right = np.einsum('ijil->jl', ktr) / N
connected_KA = (KA - np.kron(local_left, np.eye(N))
                - np.kron(np.eye(N), local_right)
                + np.trace(KA) / N**2 * np.eye(N*N))
# With three idle and three post-Z times, X,Y have zero history projection.
lowActive = sourceA.copy()
lowActive[:, :2*N] = 0
Rlim = {z: boundary(ea, ua, lowActive, z) for z in zs}
activeFs = [0 * Ftarget[0], 0 * Ftarget[1], Ftarget[2]]
Klim = exchange(ea, ua, activeFs, Gtarget)
summary = {
    'seed': SEED, 'T': T, 'workspace_qubits': 2, 'full_dimension': L*D,
    'target_dimension': len(HA), 'logical_dimension': N,
    'z_grid': zs, 'target_gap': Delta, 'target_norm': opnorm(HA),
    'penalty_norm_M': opnorm(models['idle'][1]),
    'target_K2_eigenvalues': np.linalg.eigvalsh(KA).tolist(),
    'target_K2_connected_norm': opnorm(connected_KA),
    'target_K2_direct_pair_check': opnorm(Kdirect-KA),
    'example_lambda': 1e-3,
    'active_limit_R_error': max(opnorm(Rlim[z]-RA[z]) for z in zs),
    'active_limit_K2_error': opnorm(Klim-KA), 'models': {}
}
lines = [f"seed={SEED}, T={T}, workspace=2 qubits, dimension={L*D}, N={N}",
         f"Delta={Delta:.12g}, B={opnorm(HA):.12g}, M={opnorm(models['idle'][1]):.12g}",
         'K2 target eigenvalues: ' + str(np.linalg.eigvalsh(KA)),
         f'K2 connected interaction norm: {opnorm(connected_KA):.12g}',
         f'K2 independent two-cell target check: {opnorm(Kdirect-KA):.6e}',
         'Errors use the spectral norm; R is the 3N by 3N block operator.',
         'K2 includes every excited state. Actual second-order energy = lambda^2 K2.']
for name, (H0, V, eta, G, target, Fs, source, Ws, hs) in models.items():
    diagnostics = {
        'target_match_error': opnorm(target-HA),
        'target_ground_residual': opnorm(HA@Gtarget),
        'term_ground_residual': max(opnorm(h@w@Ginput) for h, w in zip(hs, Ws)),
        'history_isometry_error': opnorm(eta.conj().T@eta-np.eye(8)),
        'history_ground_orthonormality_error': opnorm(G.conj().T@G-np.eye(N)),
        'history_source_leakage': opnorm(source-eta@(eta.conj().T@source)),
        'xi': opnorm(V@eta-eta@target),
    }
    lines += ['', name, json.dumps(diagnostics, sort_keys=True),
              'Jc       nground  splitting       compression      R_error          K2_error         gap              R_to_limit       K2_to_limit']
    rows = []
    for J in Js:
        H = J*H0+V
        e, u = eigh(H)
        ground = e[np.abs(e) < 1e-7]
        compression = max(opnorm(G.conj().T@f@G) for f in Fs)
        rerror = max(opnorm(boundary(e, u, source, z)-RA[z]) for z in zs)
        ksim = exchange(e, u, Fs, G)
        rlimit = RA if name == 'idle' else Rlim
        klimit = KA if name == 'idle' else Klim
        row = {'Jc': J, 'ground_dimension': len(ground),
               'ground_splitting': float(np.ptp(ground)) if len(ground) else None,
               'max_abs_ground_energy': float(np.max(np.abs(ground))) if len(ground) else None,
               'ground_residual': opnorm(H@G),
               'port_compression_error': compression,
               'R_error': rerror, 'K2_error': opnorm(ksim-KA),
               'gap': float(e[e > 1e-7][0]),
               'R_to_limit': max(opnorm(boundary(e, u, source, z)-rlimit[z]) for z in zs),
               'K2_to_limit': opnorm(ksim-klimit)}
        if name == 'idle':
            gamma = J/(2*L**2)
            d, d2 = gamma-opnorm(V), gamma-2*opnorm(V)
            xi2 = diagnostics['xi']**2
            if d > 0 and xi2/d < Delta:
                delta = xi2/d
                row['uniform_R_bound'] = delta/(Delta*(Delta-delta))
                assert row['R_error'] <= row['uniform_R_bound'] + 1e-8
            if d2 > 0 and 2*xi2/d2 < 2*Delta:
                delta2 = 2*xi2/d2
                row['K2_bound'] = 3*delta2/(2*Delta*(2*Delta-delta2))
                assert row['K2_error'] <= row['K2_bound'] + 1e-8
        rows.append(row)
        lines.append(f"{J:8.1e} {len(ground):3d} {row['ground_splitting']:15.6e} {compression:15.6e} {rerror:15.6e} {row['K2_error']:15.6e} {row['gap']:15.6e} {row['R_to_limit']:15.6e} {row['K2_to_limit']:15.6e}")
    tail = rows[-3:]
    slopes = {key: float(np.polyfit(np.log([r['Jc'] for r in tail]),
                                   np.log([r[key] for r in tail]), 1)[0])
              for key in ['R_error', 'K2_error', 'R_to_limit', 'K2_to_limit']}
    summary['models'][name] = {'diagnostics': diagnostics, 'rows': rows, 'tail_loglog_slopes': slopes}
    lines.append('Last-three-point log-log slopes: ' + json.dumps(slopes, sort_keys=True))
lines += ['', f"Active limiting R error: {summary['active_limit_R_error']:.12g}",
          f"Active limiting K2 error: {summary['active_limit_K2_error']:.12g}"]
with open('numerics.json', 'w') as f:
    json.dump(summary, f, indent=2)
with open('numerics.txt', 'w') as f:
    f.write('\n'.join(lines)+'\n')
print('\n'.join(lines))
