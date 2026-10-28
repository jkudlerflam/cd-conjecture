"""Independent rebuild of the C44 vacancy-charge sector (PROOFS Sec. 3-4).

Written from the prose of C44 PROOFS only, no C44/C32 code imported.
A cell of N = k + p sites holds k occupied labels in fixed order and p holes.
Diagonal: U per occupied-hole adjacent pair inside the cell.
Off-diagonal: -w for every adjacent occupied/hole swap.

Checks: the two charge gaps quoted by the C44 cell ED (n=1, rho=0.2), the
explicit lower bound (5), the asymptotic coefficient (6), the adjacency bound
1-d_(p,j) <= zeta^2 and the endpoint weight r_p >= (1-zeta^2)/2 of (7).
Small gaps are resolved by 80-digit inertia bisection (stdlib Decimal).
"""
import itertools, json, math
from decimal import Decimal, getcontext
import numpy as np

getcontext().prec = 90


def build(k, p, U, w):
    N = k + p
    confs = [frozenset(c) for c in itertools.combinations(range(N), p)]
    idx = {c: i for i, c in enumerate(confs)}
    m = len(confs)
    H = [[0 for _ in range(m)] for _ in range(m)]
    for c, i in idx.items():
        occ = [s not in c for s in range(N)]
        H[i][i] += U * sum(occ[s] != occ[s + 1] for s in range(N - 1))
        for s in range(N - 1):
            if occ[s] != occ[s + 1]:
                hole = s if not occ[s] else s + 1
                other = s + 1 if hole == s else s
                j = idx[(c - {hole}) | {other}]
                H[i][j] += -w
    return confs, H


def count_below(H, E):
    """Number of eigenvalues of symmetric H below E (Sylvester inertia, LDL^T)."""
    m = len(H)
    A = [[H[i][j] - (E if i == j else 0) for j in range(m)] for i in range(m)]
    neg = 0
    for i in range(m):
        piv = A[i][i]
        if piv == 0:
            piv = Decimal(10) ** -80
        if piv < 0:
            neg += 1
        for r in range(i + 1, m):
            if A[r][i] != 0:
                f = A[r][i] / piv
                for c in range(i + 1, m):
                    A[r][c] -= f * A[i][c]
    return neg


def level(H, which, lo, hi, iters=280):
    lo, hi = Decimal(lo), Decimal(hi)
    for _ in range(iters):
        mid = (lo + hi) / 2
        if count_below(H, mid) > which:
            hi = mid
        else:
            lo = mid
    return (lo + hi) / 2


def gap(k, p, rho, U=1):
    U, w = Decimal(U), Decimal(rho) * Decimal(U)
    _, H = build(k, p, U, w)
    top = 2 * U * (k + p) + 2 * w * (k + p)
    e0 = level(H, 0, -top, top)
    e1 = level(H, 1, -top, top)
    return e0, e1 - e0


def dense_checks(k, p, rho, U=1.0):
    w = rho * U
    confs, H = build(k, p, U, w)
    N = k + p
    # The splitting can be below double precision, so diagonalize in the
    # reflection-even sector (the Perron vector is reflection even).
    idx = {c: i for i, c in enumerate(confs)}
    m = len(confs)
    Rm = np.zeros((m, m))
    for c, i in idx.items():
        Rm[idx[frozenset(N - 1 - s for s in c)], i] = 1
    Pe = (np.eye(m) + Rm) / 2
    ev, eb = np.linalg.eigh(Pe)
    B = eb[:, ev > 0.5]
    Hf = np.array(H, float)
    vals, vecs = np.linalg.eigh(B.T @ Hf @ B)
    chi = B @ vecs[:, 0]
    chi *= np.sign(chi.sum())
    # adjacency of occupied labels j, j+1 (positions of occupied sites in order)
    d = []
    for j in range(k - 1):
        prob = 0.0
        for c, a in zip(confs, chi):
            occ = [s for s in range(N) if s not in c]
            prob += a * a * (occ[j + 1] - occ[j] == 1)
        d.append(prob)
    r_left = sum(a * a for c, a in zip(confs, chi) if 0 not in c)
    r_right = sum(a * a for c, a in zip(confs, chi) if N - 1 not in c)
    zeta = rho * (k + 1) / (1 - rho * (k + 1))
    return dict(min_chi=float(chi.min()), positive=bool(chi.min() > 0),
                e_p=float(vals[0]), max_one_minus_d=float(1 - min(d)),
                zeta2=zeta ** 2, r_left=float(r_left), r_right=float(r_right),
                r_bound=(1 - zeta ** 2) / 2)


out = {}
# 1. reproduce C44's cell-ED next levels (n=1: k=3, rho=.2, U=1)
for p, quoted in ((1, 0.01475806347), (2, 0.00003363095421)):
    e0, g = gap(3, p, '0.2')
    out[f'C44_cellED_n1_p{p}'] = dict(ours=float(g), quoted=quoted,
                                      rel_diff=float(abs(g - Decimal(quoted)) / Decimal(quoted)))
# 2. bound (5), asymptotic (6), and (7) for n = 1..4
rows = []
for n in range(1, 5):
    k = n + 2
    for p in (1, 2):
        m = math.comb(k + p, p)
        for rho in ('0.1', '0.03', '0.01', '0.001'):
            e0, g = gap(k, p, rho)
            r = Decimal(rho)
            a = r / (k + p - 1)
            bound5 = r * a ** (2 * k * p) / (m * k * p)
            row = dict(n=n, p=p, rho=float(r), gap=float(g),
                       bound5=float(bound5), bound5_holds=bool(g >= bound5),
                       coeff=float(g / r ** (k * p)), e0=float(e0))
            if float(r) * (k + 1) < 0.5:
                row.update(dense_checks(k, p, float(r)))
                row['d_bound_holds'] = row['max_one_minus_d'] <= row['zeta2'] + 1e-14
                row['r_bound_holds'] = min(row['r_left'], row['r_right']) >= row['r_bound'] - 1e-14
            rows.append(row)
out['rows'] = rows
out['all_bound5'] = all(r['bound5_holds'] for r in rows)
out['all_d_bounds'] = all(r.get('d_bound_holds', True) for r in rows)
out['all_r_bounds'] = all(r.get('r_bound_holds', True) for r in rows)
out['asymptotic_coeff_rho_1e-3'] = {f'n{r["n"]}_p{r["p"]}': r['coeff'] for r in rows if r['rho'] == 0.001}
json.dump(out, open('ref_charge.json', 'w'), indent=1)
print(json.dumps({k: v for k, v in out.items() if k != 'rows'}, indent=1))
