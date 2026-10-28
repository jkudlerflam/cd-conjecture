"""Independent implementation of C45 PROOF Sec. 3-4 (local SW and threshold).

Written from the formulas of C45 PROOF (8)-(16), not from C45's checks.py.
Cells of dimension 4: code P (2 states, energy 0) and Q (2 states, energies
1 and 9, so a multi-scale complement with Delta = 1). Open chains of N = 2..5.

Tests
 (a) P D P = lam^2 (PFP - P V Q A^-1 Q V P), and P V^d P = 0, eq. (10).
 (b) ||R|| / lam^3 converges as lam -> 0 and grows linearly in N (a
     volume-independent per-cell remainder), eq. (13).
 (c) The explicit bound ||R|| <= (N-1) * 2 * C_SW lam^3 in the regime (12),
     with C45's constants (11) evaluated in the a=0 interaction norm.
 (d) The ground-relative threshold D >= E_P + (Delta - 2J) sum_j Q_j, (15).
 (e) Whether P V Q A^-1 Q V P is a sum of two-cell bonds: true for a
     port-like V (P x P -> Q x Q only), false for a generic V with PVP = 0.
"""
import json
import numpy as np

rng = np.random.default_rng(4445)
d = 4
Pc = np.diag([1., 1., 0., 0.])
Ac = np.diag([0., 0., 1., 9.])
DELTA = 1.0


def emb(op, start, nsites, N):
    return np.kron(np.kron(np.eye(d ** start), op), np.eye(d ** (N - start - nsites)))


def herm(n):
    m = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    return (m + m.conj().T) / 2


def local_ops(n):
    P = np.ones((1, 1))
    A = np.zeros((1, 1))
    for _ in range(n):
        A = np.kron(A, np.eye(d)) + np.kron(np.eye(P.shape[0]), Ac)
        P = np.kron(P, Pc)
    Q = np.eye(d ** n) - P
    e = np.diag(A).real
    Rinv = np.diag(np.where(np.diag(Q).real > 0.5, 1 / np.where(e > 0, e, 1), 0.))
    return P, Q, Rinv


def split(terms):
    dg, od, gen = {}, {}, {}
    for (s, n), B in terms.items():
        P, Q, R = local_ops(n)
        dg[(s, n)] = P @ B @ P + Q @ B @ Q
        od[(s, n)] = P @ B @ Q + Q @ B @ P
        gen[(s, n)] = R @ B @ P - P @ B @ R
    return dg, od, gen


def comm(T1, T2):
    out = {}
    for (s1, n1), B in T1.items():
        for (s2, n2), C in T2.items():
            lo, hi = min(s1, s2), max(s1 + n1, s2 + n2)
            if max(s1, s2) >= min(s1 + n1, s2 + n2):
                continue
            Bz = np.kron(np.kron(np.eye(d ** (s1 - lo)), B), np.eye(d ** (hi - s1 - n1)))
            Cz = np.kron(np.kron(np.eye(d ** (s2 - lo)), C), np.eye(d ** (hi - s2 - n2)))
            key = (lo, hi - lo)
            out[key] = out.get(key, 0) + Bz @ Cz - Cz @ Bz
    return out


def add(*pairs):
    out = {}
    for w, T in pairs:
        for k, v in T.items():
            out[k] = out.get(k, 0) + w * v
    return out


def full(T, N):
    return sum(emb(B, s, n, N) for (s, n), B in T.items())


def norm0(T, N):
    """sup_j sum_{X contains j} ||T_X||"""
    acc = np.zeros(N)
    for (s, n), B in T.items():
        acc[s:s + n] += np.linalg.norm(B, 2)
    return acc.max()


def make_bonds(port_like):
    PP = np.kron(Pc, Pc)
    QQ = np.kron(np.eye(d) - Pc, np.eye(d) - Pc)
    if port_like:   # from the code, only both-cells-excited states are reached
        X = PP @ herm(16) @ QQ
        rest = np.eye(16) - PP
        V = X + X.conj().T + rest @ herm(16) @ rest
    else:           # generic, with PVP = 0 only
        V = herm(d * d)
        V -= PP @ V @ PP
    F = herm(d * d) * 0.7
    f1 = herm(d) * 0.3
    return V / np.linalg.norm(V, 2), F, f1


def run(N, V, F, f1, lams):
    Vt = {(j, 2): V for j in range(N - 1)}
    Ft = add((1, {(j, 2): F for j in range(N - 1)}), (1, {(j, 1): f1 for j in range(N)}))
    Vd, Vo, S1 = split(Vt)
    G2 = add((1, Ft), (1, comm(S1, Vd)), (0.5, comm(S1, Vo)))
    G2d, G2o, S2 = split(G2)
    A = sum(emb(Ac, j, 1, N) for j in range(N))
    P, Q, Rg = local_ops(N)
    Vf, Ff, S1f, S2f, Vdf, G2df = (full(T, N) for T in (Vt, Ft, S1, S2, Vd, G2d))
    res = dict(N=N, dim=d ** N)
    res['PVdP'] = float(np.linalg.norm(P @ Vdf @ P, 2))
    target = P @ (Ff - Vf @ Rg @ Vf) @ P
    res['identity10_error'] = float(np.linalg.norm(P @ G2df @ P - target, 2))
    v0 = max(1., norm0(Vt, N))
    f0 = norm0(Ft, N)
    G = f0 + 48 * v0 ** 2 / DELTA
    C_SW = 2 ** 16 * (v0 * G / DELTA + v0 ** 3 / DELTA ** 2)
    res.update(v0=v0, f0=f0, C_SW=C_SW, lam_max_12=min(v0 / G, DELTA / (192 * v0)))
    runs = []
    nQ = sum(emb(np.eye(d) - Pc, j, 1, N) for j in range(N))
    code = np.flatnonzero(np.diag(P).real > 0.5)
    for lam in lams:
        H = A + lam * Vf + lam ** 2 * Ff
        S = lam * S1f + lam ** 2 * S2f
        w, U = np.linalg.eigh(-1j * S)           # S anti-Hermitian
        eS = (U * np.exp(1j * w)) @ U.conj().T
        Hp = eS @ H @ eS.conj().T
        Dm = A + lam * Vdf + lam ** 2 * G2df
        R = Hp - Dm
        Wt = add((lam, Vd), (lam ** 2, G2d))
        J = norm0(Wt, N)
        EP = np.linalg.eigvalsh(Dm[np.ix_(code, code)])[0]
        thr = np.linalg.eigvalsh(Dm - EP * np.eye(d ** N) - (DELTA - 2 * J) * nQ)[0]
        runs.append(dict(lam=lam, R_over_lam3=float(np.linalg.norm(R, 2) / lam ** 3),
                         PRQ_over_lam3=float(np.linalg.norm(P @ R @ Q, 2) / lam ** 3),
                         bound_ratio=float(np.linalg.norm(R, 2) / (2 * (N - 1) * C_SW * lam ** 3)),
                         threshold_min_eig=float(thr), J=float(J)))
    res['runs'] = runs
    return res


def bond_sum_test(V):
    """Three-cell connected part of P V Q A^-1 Q V P."""
    N = 3
    P, Q, Rg = local_ops(N)
    Vf = emb(V, 0, 2, N) + emb(V, 1, 2, N)
    K3 = P @ Vf @ Rg @ Vf @ P
    P2, Q2, R2 = local_ops(2)
    Kb = P2 @ V @ R2 @ V @ P2
    Kbonds = np.kron(Kb, Pc) + np.kron(Pc, Kb)
    return float(np.linalg.norm(K3 - Kbonds, 2) / np.linalg.norm(K3, 2))


out = {}
for label, port in (('generic_V', False), ('port_like_V', True)):
    V, F, f1 = make_bonds(port)
    out[label] = dict(three_cell_part_of_second_order=bond_sum_test(V), chains=[])
    for N in (2, 3, 4, 5):
        lams = (0.04, 0.02, 0.01, 0.002)
        out[label]['chains'].append(run(N, V, F, f1, lams))
json.dump(out, open('ref_localsw.json', 'w'), indent=1)
for label in out:
    print(label, 'three-cell part', out[label]['three_cell_part_of_second_order'])
    for c in out[label]['chains']:
        print(' N=%d id10=%.1e PVdP=%.1e lam12=%.2e' % (c['N'], c['identity10_error'], c['PVdP'], c['lam_max_12']),
              ' '.join('[%.3g: R/l3=%.4f PRQ/l3=%.4f bnd=%.1e thr=%.1e]' % (r['lam'], r['R_over_lam3'], r['PRQ_over_lam3'],
                                                                         r['bound_ratio'], r['threshold_min_eig']) for r in c['runs']))
