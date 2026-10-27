#!/usr/bin/env python3
"""Ground-state Koo-Saleur test for four critical spin chains.

The implementation is closed-book and matrix-free.  Spin-flip models use the
full computational basis, while U(1)-symmetric models use Sz=0.  Local energy
densities are also represented as Pauli strings, which is used both for h_q
and for the energy current

    J = i sum_j sum_{r>0} r [h_j,h_{j+r}].

This is the coefficient of q in [H,h_q] for h_q=sum_j exp(i q j)h_j.
"""
import argparse
import csv
import time
from itertools import combinations

import numpy as np
from numba import njit
from scipy.sparse.linalg import LinearOperator, eigsh


@njit(cache=True)
def poppar(x):
    p = 0
    while x:
        p ^= 1
        x &= x - 1
    return p


@njit(cache=True)
def apply_paulis(states, inv, xs, zs, cs, v):
    out = np.zeros(v.size, dtype=np.complex128)
    for a in range(states.size):
        s = states[a]
        va = v[a]
        for t in range(xs.size):
            bstate = s ^ xs[t]
            b = inv[bstate]
            if b >= 0:
                sign = 1.0 if poppar(s & zs[t]) == 0 else -1.0
                out[b] += cs[t] * sign * va
    return out


@njit(cache=True)
def ising_mv(v, n, lam):
    d = v.size
    out = np.empty_like(v)
    mask = (1 << n) - 1
    for s in range(d):
        zsum = 0.0
        for j in range(n):
            zsum += 1.0 if ((s >> j) & 1) == 0 else -1.0
        val = -zsum * v[s]
        for j in range(n):
            k = (j + 1) % n
            zl = 1.0 if ((s >> ((j - 1) % n)) & 1) == 0 else -1.0
            zr = 1.0 if ((s >> ((j + 2) % n)) & 1) == 0 else -1.0
            dest = s ^ (1 << j) ^ (1 << k)
            val += (-1.0 + lam * (zl + zr)) * v[dest]
        out[s] = val
    return out


@njit(cache=True)
def heis_mv(v, states, inv, n, j2):
    out = np.zeros_like(v)
    for a in range(states.size):
        s = states[a]
        diag = 0.0
        for j in range(n):
            for dist in range(1, 3):
                coupling = 1.0 if dist == 1 else j2
                if coupling == 0.0:
                    continue
                k = (j + dist) % n
                bj = (s >> j) & 1
                bk = (s >> k) & 1
                diag += coupling * (0.25 if bj == bk else -0.25)
                if bj != bk:
                    dest = s ^ (1 << j) ^ (1 << k)
                    out[inv[dest]] += 0.5 * coupling * v[a]
        out[a] += diag * v[a]
    return out


def sector(n, u1):
    if not u1:
        states = np.arange(1 << n, dtype=np.int64)
        inv = np.arange(1 << n, dtype=np.int32)
    else:
        ss = []
        for comb in combinations(range(n), n // 2):
            x = 0
            for j in comb:
                x |= 1 << j
            ss.append(x)
        states = np.array(sorted(ss), dtype=np.int64)
        inv = np.full(1 << n, -1, dtype=np.int32)
        inv[states] = np.arange(states.size, dtype=np.int32)
    return states, inv


def add_term(dic, x, z, c, tol=1e-14):
    key = (int(x), int(z))
    dic[key] = dic.get(key, 0.0j) + c
    if abs(dic[key]) < tol:
        del dic[key]


def spin_dot(dic, a, b, c):
    # S.S=(XX+YY+ZZ)/4, with Y=iXZ in the canonical X^x Z^z form.
    add_term(dic, (1 << a) | (1 << b), 0, c / 4)
    add_term(dic, (1 << a) | (1 << b), (1 << a) | (1 << b), -c / 4)
    add_term(dic, 0, (1 << a) | (1 << b), c / 4)


def local_terms(kind, n, param, j):
    d = {}
    if kind == "ising":
        jm, jp = (j - 1) % n, (j + 1) % n
        add_term(d, (1 << jm) | (1 << j), 0, -0.5)
        add_term(d, (1 << j) | (1 << jp), 0, -0.5)
        add_term(d, 0, 1 << j, -1.0)
        if param:
            add_term(d, (1 << jm) | (1 << j), 1 << jp, param)
            add_term(d, (1 << j) | (1 << jp), 1 << jm, param)
    else:
        spin_dot(d, j, (j + 1) % n, 1.0)
        if param:
            spin_dot(d, j, (j + 2) % n, param)
    return d


def multiply(a, b):
    (xa, za), ca = a
    (xb, zb), cb = b
    sign = -1.0 if ((za & xb).bit_count() & 1) else 1.0
    return (xa ^ xb, za ^ zb), ca * cb * sign


def commutator_into(out, A, B, pref):
    for aa in A.items():
        for bb in B.items():
            k1, c1 = multiply(aa, bb)
            k2, c2 = multiply(bb, aa)
            add_term(out, k1[0], k1[1], pref * (c1 - c2))


def make_hq_and_current(kind, n, param, q):
    locals_ = [local_terms(kind, n, param, j) for j in range(n)]
    hq = {}
    for j, hj in enumerate(locals_):
        w = np.exp(1j * q * j)
        for (x, z), c in hj.items():
            add_term(hq, x, z, w * c)
    cur = {}
    # The densities have diameter at most two, so only r=1,2 contribute.
    for r in (1, 2):
        for j in range(n):
            commutator_into(cur, locals_[j], locals_[(j + r) % n], 1j * r)
    return hq, cur


def arrays(dic):
    items = list(dic.items())
    return (np.array([k[0] for k, c in items], dtype=np.int64),
            np.array([k[1] for k, c in items], dtype=np.int64),
            np.array([c for k, c in items], dtype=np.complex128))


def translate(v, states, inv, n):
    mask = (1 << n) - 1
    deststates = ((states << 1) & mask) | (states >> (n - 1))
    out = np.empty_like(v)
    out[inv[deststates]] = v
    return out


def lanczos_measure(matvec, start, e0, steps=55):
    """Ritz energies and spectral weights for the measure generated by start."""
    norm = np.linalg.norm(start)
    q = start / norm
    qprev = np.zeros_like(q)
    alpha, beta = [], []
    Q = []
    bprev = 0.0
    for k in range(steps):
        Q.append(q.copy())
        z = matvec(q) - (bprev * qprev if k else 0)
        a = np.vdot(q, z).real
        z -= a * q
        # Full reorthogonalization is cheap relative to H at these sizes and
        # prevents ghost stress-tensor levels.
        for qq in Q:
            z -= np.vdot(qq, z) * qq
        b = np.linalg.norm(z)
        alpha.append(a)
        if k == steps - 1 or b < 1e-12:
            break
        beta.append(b)
        qprev, q, bprev = q, z / b, b
    T = np.diag(alpha)
    if beta:
        bb = np.array(beta[:len(alpha)-1])
        T += np.diag(bb, 1) + np.diag(bb, -1)
    ev, U = np.linalg.eigh(T)
    weights = norm * norm * np.abs(U[0, :]) ** 2
    gaps = ev - e0
    good = np.where(gaps > 1e-7)[0]
    good = good[:min(10, good.size)]
    pick = good[np.argmax(weights[good])]
    return gaps[pick], weights[pick], ev, weights


def run_one(kind, param, n, lanczos_steps):
    u1 = kind == "heis"
    states, inv = sector(n, u1)
    if kind == "ising":
        def mv_real(v): return ising_mv(v, n, param)
    else:
        def mv_real(v): return heis_mv(v, states, inv, n, param)
    dim = states.size
    H = LinearOperator((dim, dim), matvec=mv_real, dtype=np.float64)
    vals, vecs = eigsh(H, k=2, which="SA", tol=2e-10, maxiter=1200,
                       ncv=min(28, dim), v0=np.ones(dim))
    order = np.argsort(vals)
    vals, vecs = vals[order], vecs[:, order]
    e0, gap0 = vals[0], vals[1] - vals[0]
    g = vecs[:, 0]
    tg = translate(g, states, inv, n)
    tovl = np.vdot(g, tg)
    tres = np.linalg.norm(tg - tovl * g)

    # h_-2 has exp(-i 4 pi j/N), as in the question.
    hminus, cur = make_hq_and_current(kind, n, param, -4 * np.pi / n)
    hx, hz, hc = arrays(hminus)
    jx, jz, jc = arrays(cur)
    v = apply_paulis(states, inv, hx, hz, hc, g.astype(np.complex128))
    m0 = np.vdot(v, v).real
    Hv = mv_real(v)
    m1 = np.vdot(v, Hv - e0 * v).real
    jv = apply_paulis(states, inv, jx, jz, jc, g.astype(np.complex128))
    jmean = np.vdot(g, jv)
    jvar = np.vdot(jv, jv).real - abs(jmean) ** 2
    c2 = 8 * m0 ** 3 / m1 ** 2

    # H is real-linear on complex vectors as well.  Lanczos remains in the
    # momentum sector of v and directly resolves its spectral measure.
    dT, wT, _, _ = lanczos_measure(mv_real, v, e0, lanczos_steps)
    cmv = 8 * wT / dT ** 2
    return dict(N=n, dim=dim, E0=e0, unique_gap=gap0,
                trans_abs=abs(tovl), trans_res=tres, m0N2=m0*n*n,
                m1N3=m1*n**3, JvarN=jvar/n, c2=c2, cMV=cmv,
                Tgap=dT, Tweight=wT)


def fit_report(rows):
    N = np.array([r["N"] for r in rows], float)
    c2 = np.array([r["c2"] for r in rows])
    cmv = np.array([r["cMV"] for r in rows])
    # General intercept fits and the theoretically motivated N*c2 fit.
    deg = min(2, len(N)-1)
    pp = np.polyfit(1/N, c2, deg)
    p2 = np.pad(pp, (3-len(pp), 0))
    pmv = np.polyfit(1/N**2, cmv, min(1, len(N)-1))
    if len(pmv) == 1:
        pmv = np.array([0.0, pmv[0]])
    tail = N >= min(N[-1], max(12, N[-1]-6))
    if np.count_nonzero(tail) < 2:
        tail = np.arange(len(N)) >= max(0, len(N)-2)
    amp = np.polyfit(1/N[tail], c2[tail], 1)
    X0 = np.column_stack((1/N, 1/N**2))
    zero_a1, zero_a2 = np.linalg.lstsq(X0, c2, rcond=None)[0]
    power = np.polyfit(np.log(N[tail]), np.log(c2[tail]), 1)[0]
    return dict(c2_intercept=p2[2], c2_a1=p2[1], c2_a2=p2[0],
                c2_tail_intercept=amp[1], c2_tail_a1=amp[0],
                c2_zero_a1=zero_a1, c2_zero_a2=zero_a2,
                c2_tail_power=power,
                cmv_intercept=pmv[1], cmv_a2=pmv[0])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-n", type=int, default=20)
    ap.add_argument("--min-n", type=int, default=8)
    ap.add_argument("--lanczos", type=int, default=55)
    ap.add_argument("--models", default="tfim,of02,of035,heis,j1j2")
    args = ap.parse_args()
    models = {
        "tfim": ("ising", 0.0), "of02": ("ising", 0.2),
        "of035": ("ising", 0.35), "heis": ("heis", 0.0),
        "j1j2": ("heis", 0.2),
    }
    chosen = [x for x in args.models.split(",") if x]
    allrows = {}
    t0 = time.time()
    for name in chosen:
        kind, param = models[name]
        rows = []
        for n in range(args.min_n, args.max_n + 1, 2):
            t = time.time()
            r = run_one(kind, param, n, args.lanczos)
            r["model"] = name
            rows.append(r)
            print(name, n, " ".join(f"{k}={r[k]:.9g}" for k in
                  ("c2", "m0N2", "m1N3", "JvarN", "cMV")),
                  f"gap={r['unique_gap']:.3g} transres={r['trans_res']:.2g}",
                  f"sec={time.time()-t:.1f}", flush=True)
        allrows[name] = rows
        print("FIT", name, fit_report(rows), flush=True)
    fields = ["model", "N", "dim", "E0", "unique_gap", "trans_abs",
              "trans_res", "m0N2", "m1N3", "JvarN", "c2", "cMV",
              "Tgap", "Tweight"]
    with open("ks_results.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for rows in allrows.values():
            w.writerows(rows)
    with open("ks_fits.txt", "w") as f:
        for name, rows in allrows.items():
            f.write(name + " " + repr(fit_report(rows)) + "\n")
    print(f"TOTAL_SECONDS {time.time()-t0:.1f}")


if __name__ == "__main__":
    main()
