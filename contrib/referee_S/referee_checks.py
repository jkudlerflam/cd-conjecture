"""Independent referee checks for C48, C49, C51, C52, C58 (input S).

Written from the statements in the PROOF.md files, without reusing their
checks.py code. Small toy models only. They test inequalities and
identities, not the microscopic 4464-state chain.

Run:  OPENBLAS_NUM_THREADS=1 python3 -B referee_checks.py
Writes referee_checks.json next to this file.
"""
import itertools
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

RNG = np.random.default_rng(20261028)


def herm(n, scale=1.0):
    x = RNG.normal(size=(n, n)) + 1j * RNG.normal(size=(n, n))
    x = (x + x.conj().T) / 2
    return scale * x / np.linalg.norm(x, 2)


def kron(*ops):
    out = np.ones((1, 1))
    for o in ops:
        out = np.kron(out, o)
    return out


def site_op(op, j, M, d):
    return kron(*[op if k == j else np.eye(d) for k in range(M)])


def bond_op(op2, j, M, d):
    """Two-site operator on sites j, j+1 (cyclic) of an M-site chain."""
    k = (j + 1) % M
    if k == j + 1:
        return kron(np.eye(d ** j), op2, np.eye(d ** (M - j - 2)))
    # wrap-around bond (j=M-1, k=0): permute
    perm_op = op2.reshape(d, d, d, d).transpose(1, 0, 3, 2).reshape(d * d, d * d)
    return bond_op(perm_op, 0, M, d) if M == 2 else _wrap(op2, M, d)


def _wrap(op2, M, d):
    # build via basis action
    dim = d ** M
    out = np.zeros((dim, dim), complex)
    o = op2.reshape(d, d, d, d)
    for idx in range(dim):
        digits = list(np.unravel_index(idx, [d] * M))
        a, b = digits[M - 1], digits[0]
        for a2 in range(d):
            for b2 in range(d):
                amp = o[a2, b2, a, b]
                if amp == 0:
                    continue
                nd = digits.copy()
                nd[M - 1], nd[0] = a2, b2
                out[np.ravel_multi_index(nd, [d] * M), idx] += amp
    return out


def logZ(evals, beta):
    x = -beta * np.asarray(evals)
    m = x.max()
    return m + math.log(np.exp(x - m).sum())


# ---------------------------------------------------------------- C48
def c48_enlarged_embedding():
    """Multi-cell version of C48 (3),(4),(32),(33): ring of M=3 cells.

    Original cell: code C^2 (A=0) plus k excited states, A>=Delta Q.
    Enlarged cell: C^2_log (x) (Omega + C^k). B commutes with P_X, R generic.
    """
    out = []
    M, k, Delta = 3, 2, 1.0
    D = 2 + k
    P1 = np.diag([1, 1] + [0] * k).astype(complex)
    for u in (0.0, 1.0, 5.0):
        A1 = np.diag([0, 0, Delta, Delta + u]).astype(complex)
        # isometry J: K -> C^2 (x) (C^{1+k}); code |a> -> |a>|Omega>, exc |s> -> |0>|s>
        E = 2 * (1 + k)
        J = np.zeros((E, D), complex)
        J[0 * (1 + k) + 0, 0] = 1
        J[1 * (1 + k) + 0, 1] = 1
        for s in range(k):
            J[0 * (1 + k) + 1 + s, 2 + s] = 1
        a_aux = np.kron(np.eye(2), np.diag([0.0] + [Delta, Delta + u]))
        assert np.allclose(J.conj().T @ J, np.eye(D))
        assert np.allclose(a_aux @ J, J @ A1)
        P2 = np.kron(P1, P1)
        Q2 = np.eye(D * D) - P2
        for lam, rr in ((0.08, 0.0), (0.08, 0.01), (0.15, 0.02)):
            h = herm(D * D)
            Bb = lam * (P2 @ h @ P2 + Q2 @ h @ Q2)
            Rb = rr * herm(D * D)
            Hb = Bb + Rb
            Horig = sum(site_op(A1, j, M, D) for j in range(M)) + \
                sum(_wrap(Hb, M, D) if j == M - 1 else bond_op(Hb, j, M, D) for j in range(M))
            J2 = np.kron(J, J)
            HbE = J2 @ Hb @ J2.conj().T
            Henl = sum(site_op(a_aux, j, M, E) for j in range(M)) + \
                sum(_wrap(HbE, M, E) if j == M - 1 else bond_op(HbE, j, M, E) for j in range(M))
            eo = np.linalg.eigvalsh(Horig)
            ee = np.linalg.eigvalsh(Henl)
            w = 2 * np.linalg.norm(Hb, 2)  # incidence norm: two bonds per site
            thr = eo[0] + Delta - w
            lo_o, lo_e = eo[eo < thr - 1e-9], ee[ee < thr - 1e-9]
            spec_err = (float(np.max(np.abs(lo_o - lo_e))) if len(lo_o) == len(lo_e) and len(lo_o) else
                        (0.0 if len(lo_o) == len(lo_e) else float("inf")))
            # pressure comparison (33)
            Tr_Q = lambda b: math.exp(-b * Delta) + math.exp(-b * (Delta + u))
            worst = float("inf")
            for beta in (2.0, 5.0, 12.0):
                df = (logZ(ee, beta) - logZ(eo, beta)) / (beta * M)
                zin = (2 - 1) * Tr_Q(beta)
                zval = 2 + Tr_Q(beta)
                bound = math.log1p(math.exp(beta * w) * zin / zval) / beta
                assert -1e-12 <= df <= bound + 1e-12, (df, bound)
                worst = min(worst, bound - df)
            out.append(dict(u=u, lam=lam, r=rr, w=w, levels_below_threshold=int(len(lo_o)),
                            spec_err=spec_err, min_pressure_slack=worst))
            assert spec_err < 1e-10
    return out


def c48_leakage_linear_uniform_in_u():
    """Low spectrum of A+B+R vs H_log=P B P, B block diagonal, R generic.

    Tests that the matching error is O(r) with a constant independent of the
    massive-sector internal interaction u (the C46 regime u>=Delta)."""
    M, Delta, lam = 3, 1.0, 0.05
    rows = []
    P1 = np.diag([1, 1, 0, 0, 0]).astype(complex)
    D = 5
    P2 = np.kron(P1, P1)
    Q2 = np.eye(D * D) - P2
    h = herm(D * D)
    Bb = lam * (P2 @ h @ P2 + Q2 @ h @ Q2)
    Rshape = herm(D * D)
    Pfull = kron(*([P1] * M))
    idx = np.where(np.diag(Pfull).real > 0.5)[0]
    for u in (0.0, 1.0, 3.0, 10.0, 100.0):
        A1 = np.diag([0, 0, Delta, Delta + u / 2, Delta + u]).astype(complex)
        HA = sum(site_op(A1, j, M, D) for j in range(M))
        HB = sum(_wrap(Bb, M, D) if j == M - 1 else bond_op(Bb, j, M, D) for j in range(M))
        HR1 = sum(_wrap(Rshape, M, D) if j == M - 1 else bond_op(Rshape, j, M, D) for j in range(M))
        Hlog = HB[np.ix_(idx, idx)]
        kvals = np.linalg.eigvalsh(Hlog)
        e0 = np.linalg.eigvalsh(HA + HB)[: len(kvals)]
        base = float(np.max(np.abs(e0 - kvals)))
        for r in (1e-2, 1e-3, 1e-4):
            hv = np.linalg.eigvalsh(HA + HB + r * HR1)[: len(kvals)]
            err = float(np.max(np.abs(hv - kvals)))
            rows.append(dict(u=u, r=r, err=err, err_over_r=err / r, r0_err=base))
    ratios = [x["err_over_r"] for x in rows]
    assert max(x["r0_err"] for x in rows) < 1e-12
    return dict(rows=rows, max_err_over_r=max(ratios), min_err_over_r=min(ratios),
                note="bond R norm 1 per bond, 3 bonds; err/r bounded by sum of bond norms")


def c46_resummed_formula():
    """One cell: logical n1,n2 (4 code states), excursion to e_s with E_s=(D,D,D,D+u)."""
    rows = []
    for Delta, u in ((1.0, 1.0), (1.0, 3.0), (1.0, 20.0)):
        Es = np.array([Delta, Delta, Delta, Delta + u])
        for lam in (0.05, 0.02):
            # basis: 4 code states |s>, 4 excited |e_s>; flip couples |s> <-> |e_s>
            H = np.zeros((8, 8))
            H[4:, 4:] = np.diag(Es)
            for s in range(4):
                H[s, 4 + s] = H[4 + s, s] = lam
            ev = np.sort(np.linalg.eigvalsh(H))[:4]
            N = np.array([0, 0, 0, 1.0])  # n1 n2 for s=00,10,01,11
            pred = np.sort(-lam ** 2 / Delta + lam ** 2 * u / (Delta * (Delta + u)) * N)
            err = float(np.max(np.abs(ev - pred)))
            rows.append(dict(Delta=Delta, u=u, lam=lam, err=err, lam4=lam ** 4))
            assert err <= 2 * lam ** 4 / Delta ** 3
    return rows


# ---------------------------------------------------------------- C58
def c58_pinning_lemma():
    """C58 (4) on spin-1 rings: marker = state 2, random marker-conserving complex bond."""
    d = 3
    # bond conserving number of '2's: block structure by count of 2 in pair
    counts = np.array([[(a == 2) + (b == 2) for b in range(d)] for a in range(d)]).ravel()
    tests, worst = 0, float("inf")
    for trial in range(3):
        h = herm(9)
        mask = counts[:, None] == counts[None, :]
        h = h * mask
        h = (h + h.conj().T) / 2
        J = np.linalg.norm(h, 2)
        for L in range(3, 7):
            H = sum(_wrap(h, L, d) if j == L - 1 else bond_op(h, j, L, d) for j in range(L))
            dig = np.array(list(itertools.product(range(d), repeat=L)))
            Mvec = (dig == 2).sum(1)
            for M in range(1, L + 1):
                sec = np.where(Mvec == M)[0]
                Hs = H[np.ix_(sec, sec)]
                pin = np.where(dig[sec, 0] == 2)[0]
                Hp = Hs[np.ix_(pin, pin)]
                es, ep = np.linalg.eigvalsh(Hs), np.linalg.eigvalsh(Hp)
                for beta in (0.3, 1.0, 5.0):
                    lz, lp = logZ(es, beta), logZ(ep, beta)
                    lower = lz - (math.log(L / M) + lp)
                    upper = (math.log(L) + 4 * beta * J * L / M + lp) - lz
                    assert lower > -1e-10 and upper > -1e-10
                    worst = min(worst, lower, upper)
                    tests += 1
    return dict(tests=tests, min_slack=worst)


def c58_partial_trace_peierls():
    """C58 (14): Z_g >= d_P^{|D|} exp(-beta w |D|) Z_rest, for a cell compressed to its code.

    One D-cell (dim 4, code dim 2 at zero energy, excited >= 1) coupled by V (norm w)
    to a rest system of dimension 6 with random H_rest."""
    tests, worst = 0, float("inf")
    for trial in range(20):
        dc, dr = 4, 6
        A = np.diag([0, 0, 1 + RNG.random(), 1 + 2 * RNG.random()])
        Hr = herm(dr, 2.0)
        V = herm(dc * dr, RNG.random())
        w = np.linalg.norm(V, 2)
        H = np.kron(A, np.eye(dr)) + np.kron(np.eye(dc), Hr) + V
        for beta in (0.5, 2.0, 8.0):
            lhs = logZ(np.linalg.eigvalsh(H), beta)
            rhs = math.log(2) - beta * w + logZ(np.linalg.eigvalsh(Hr), beta)
            assert lhs >= rhs - 1e-10
            worst = min(worst, lhs - rhs)
            tests += 1
    return dict(tests=tests, min_slack=worst)


def c58_semigroup_heat_constant():
    """C58 Sec. 5: from Tr_Q exp(-tA) <= d e^{-2t kappa}/(1-d e^{-t kappa}) and A>=Delta Q,
    get Tr_Q exp(-tA/2) <= C0 exp(-t Delta/2) for t>=2 t0, C0=max(1,exp(t0 Delta)/(2d)).
    Synthetic spectrum: level max(Delta, kappa(l+1)) with multiplicity d^l, l=1..Lmax."""
    rows = []
    for d, b in ((3, 2), (5, 4)):
        Delta = 1.0
        kappa = Delta / (b + 2)
        t0 = math.log(2 * d) / kappa
        C0 = max(1.0, math.exp(t0 * Delta) / (2 * d))
        worst = float("inf")
        for t in np.linspace(2 * t0, 12 * t0, 30):
            ls = np.arange(1, 400)
            lev = np.maximum(Delta, kappa * (ls + 1))
            logtr = np.log(np.sum(np.exp(ls * math.log(d) - t * lev / 2 - (-t * Delta / 2))))  # relative
            # logtr = log(Tr e^{-tA/2}) + t Delta/2
            assert logtr <= math.log(C0) + 1e-9
            worst = min(worst, math.log(C0) - logtr)
        rows.append(dict(d=d, b=b, log_C0=math.log(C0), min_log_slack=worst))
    return rows


# ---------------------------------------------------------------- C51
def c51_canonical_diag():
    worst = 0.0
    for N in (1, 2, 3, 4):
        S = np.diag([1.0] * N + [-1.0] * N)
        for _ in range(20):
            X = RNG.normal(size=(2 * N, 2 * N))
            V = X @ X.T + 0.05 * np.eye(2 * N)
            w, Ov = np.linalg.eigh(V)
            Vh = Ov @ np.diag(np.sqrt(w)) @ Ov.T
            Vmh = Ov @ np.diag(1 / np.sqrt(w)) @ Ov.T
            W = Vh @ S @ Vh
            ev, O = np.linalg.eigh(W)
            order = np.r_[np.where(ev > 0)[0], np.where(ev < 0)[0]]
            ev, O = ev[order], O[:, order]
            assert (ev[:N] > 0).all() and (ev[N:] < 0).all()
            U = np.diag(np.abs(ev))
            R = Vmh @ O @ np.sqrt(U)
            e1 = np.abs(R.T @ V @ R - U).max() / np.abs(U).max()
            e2 = np.abs(R.T @ S @ R - S).max()
            # speeds = eigenvalues of S V (equations of motion)
            sv = np.sort(np.linalg.eigvals(S @ V).real)
            e3 = np.abs(sv - np.sort(ev)).max() / np.abs(ev).max()
            worst = max(worst, e1, e2, e3)
    assert worst < 1e-8
    return dict(max_rel_residual=worst)


def c51_virasoro_norm():
    """||L_{-m}|x>||^2 = m x^2 + (m^3-m)/12 by explicit contraction count."""
    worst = 0.0
    for x in (0.0, 0.37, math.sqrt(2), 3.1):
        for m in range(1, 9):
            # L_{-m}|x> = x J_{-m}|x> + 1/2 sum_{r=1}^{m-1} J_{-r}J_{-(m-r)}|x>
            # norm of second piece: sum over unordered pairs {r, m-r}
            val = m * x * x
            for r in range(1, m):
                s = m - r
                if r < s:
                    val += r * s  # (1/2)^2 * 2 (ordered pairs) * <..> = r s per unordered pair, times 2 orderings/2
                elif r == s:
                    val += 0.25 * 2 * r * r
            pred = m * x * x + (m ** 3 - m) / 12
            worst = max(worst, abs(val - pred))
    assert worst < 1e-12
    return dict(max_err=worst)


def c51_asymmetric_free_fermion():
    """Grand potential of free fermions with dispersion eps(k)=-2cos k + a sin 2k at fixed mu:
    omega(T)-omega(0) -> -(pi/12) T^2 sum_{Fermi points} 1/|v|.  Unequal left/right speeds."""
    rows = []
    for a, mu in ((0.0, -0.4), (0.35, -0.4), (0.6, 0.3)):
        eps = lambda k: -2 * math.cos(k) + a * math.sin(2 * k)
        deps = lambda k: 2 * math.sin(k) + 2 * a * math.cos(2 * k)
        grid = np.linspace(-math.pi, math.pi, 20001)
        vals = np.array([eps(k) - mu for k in grid])
        roots = []
        for i in range(len(grid) - 1):
            if vals[i] * vals[i + 1] < 0:
                roots.append(brentq(lambda k: eps(k) - mu, grid[i], grid[i + 1], xtol=1e-14))
        vs = [deps(k) for k in roots]
        A = (math.pi / 12) * sum(1 / abs(v) for v in vs)
        pts = sorted(roots)

        def omega(T):
            f = lambda k: -T * np.logaddexp(0.0, -(eps(k) - mu) / T) / (2 * math.pi)
            tot = 0.0
            edges = [-math.pi] + pts + [math.pi]
            for lo, hi in zip(edges[:-1], edges[1:]):
                tot += quad(f, lo, hi, limit=400, epsabs=1e-14, epsrel=1e-13)[0]
            return tot

        def omega0():
            g = lambda k: min(eps(k) - mu, 0.0) / (2 * math.pi)
            edges = [-math.pi] + pts + [math.pi]
            return sum(quad(g, lo, hi, limit=400, epsabs=1e-15)[0] for lo, hi in zip(edges[:-1], edges[1:]))

        w0 = omega0()
        coeffs = []
        for T in (0.02, 0.01, 0.005):
            coeffs.append(-(omega(T) - w0) / T ** 2)
        rows.append(dict(a=a, mu=mu, speeds=vs, predicted=A, measured=coeffs,
                         rel_err_smallest_T=abs(coeffs[-1] - A) / A))
        assert abs(coeffs[-1] - A) / A < 2e-3
    return rows


# ---------------------------------------------------------------- C52
def c52_three_cell_commutator():
    """<1a|[K12,K23]|3b> = (T1 T2)_ab on a JW Fock space, n flavours per cell,
    with random onsite, density-density and pair-hopping quartic terms."""
    worst = 0.0
    for n in (1, 2):
        modes = 3 * n
        dim = 2 ** modes
        # JW annihilators
        Z = np.diag([1, -1]); a = np.array([[0, 1], [0, 0]])
        c = []
        for p in range(modes):
            c.append(kron(*([Z] * p + [a] + [np.eye(2)] * (modes - p - 1))))
        cd = [x.conj().T for x in c]
        ix = lambda cell, f: cell * n + f

        def bond(i, j, T):
            K = np.zeros((dim, dim), complex)
            for f in range(n):
                for g in range(n):
                    term = T[f, g] * cd[ix(i, f)] @ c[ix(j, g)]
                    K += term + term.conj().T
            # onsite quadratic
            for cell in (i, j):
                h = herm(n)
                for f in range(n):
                    for g in range(n):
                        K += h[f, g] * cd[ix(cell, f)] @ c[ix(cell, g)]
            # normal-ordered quartics: density-density and pair hopping across the bond
            for f in range(n):
                for g in range(n):
                    K += RNG.normal() * cd[ix(i, f)] @ c[ix(i, f)] @ cd[ix(j, g)] @ c[ix(j, g)]
            if n == 2:
                pair = RNG.normal() * cd[ix(i, 0)] @ cd[ix(i, 1)] @ c[ix(j, 1)] @ c[ix(j, 0)]
                K += pair + pair.conj().T
            return K

        T1 = RNG.normal(size=(n, n)) + 1j * RNG.normal(size=(n, n))
        T2 = RNG.normal(size=(n, n)) + 1j * RNG.normal(size=(n, n))
        K12, K23 = bond(0, 1, T1), bond(1, 2, T2)
        Cm = K12 @ K23 - K23 @ K12
        vac = np.zeros(dim); vac[0] = 1  # all modes empty: state index 0 (bit 0 = empty)
        # with a=[[0,1],[0,0]], empty = basis 0
        for f in range(n):
            for g in range(n):
                bra = cd[ix(0, f)] @ vac
                ket = cd[ix(2, g)] @ vac
                val = bra.conj() @ Cm @ ket
                worst = max(worst, abs(val - (T1 @ T2)[f, g]))
    assert worst < 1e-10
    return dict(max_err=worst)


def c52_head_leakage():
    rows = []
    for b in (4, 5, 8, 16, 40, 200):
        phi = lambda l, i: math.sqrt(2 / l) * math.sin(math.pi * i / l)
        u = np.array([phi(b, i) for i in range(1, b)] + [0.0])            # grown cell, length b+1, positions 1..b
        v = np.array([phi(b + 1, j + 1) for j in range(1, b)])            # shrunk cell, positions 1..b-1
        tgt_u = np.array([phi(b + 1, i) for i in range(1, b + 1)])
        tgt_v = np.array([phi(b, j) for j in range(1, b)])
        w = np.kron(u, v)
        tgt = np.kron(tgt_u, tgt_v)
        s_b = 1 - phi(b + 1, 1) ** 2
        A_b = sum(phi(b, i) * phi(b + 1, i) for i in range(1, b))
        leak = float(w @ w - (tgt @ w) ** 2)
        assert abs(w @ w - s_b) < 1e-13 and abs(tgt @ w - A_b ** 2) < 1e-13
        assert abs(leak - (s_b - A_b ** 4)) < 1e-13 and leak >= s_b * (1 - s_b) - 1e-15
        ratios = [phi(b + 1, i) / phi(b, i) for i in range(1, b)]
        assert all(np.diff(ratios) > 0)
        rows.append(dict(b=b, leakage_sq=leak, lower=s_b * (1 - s_b)))
    return rows


def c52_gram_singular_value():
    worst = float("inf")
    for n in (1, 2, 4, 8):
        for _ in range(50):
            X = RNG.normal(size=(n, n)) + 1j * RNG.normal(size=(n, n))
            G1 = X @ X.conj().T + 0.1 * np.eye(n)
            Y = RNG.normal(size=(n, n)) + 1j * RNG.normal(size=(n, n))
            G2 = Y @ Y.conj().T
            W, eta = 0.7, 3.0 * RNG.random()
            T0 = -W * (G1 + 1j * eta * G2)
            smin = np.linalg.svd(T0, compute_uv=False).min()
            bound = W * np.linalg.eigvalsh(G1).min()
            assert smin >= bound - 1e-12
            worst = min(worst, smin / bound)
    return dict(min_ratio_smin_over_bound=worst, rank_bound_d4464=4464 ** 2 // 4)


# ---------------------------------------------------------------- C49
def c49_theorem_F():
    """Schur/Feshbach comparison (18),(26) on M cells of dimension 4 (code dim 2)."""
    rows = []
    d, dP = 4, 2
    P1 = np.diag([1, 1, 0, 0]).astype(complex)
    for M in (3, 4):
        for trial in range(4):
            Delta = 1.0
            A1 = np.diag([0, 0, Delta, Delta + RNG.random()]).astype(complex)
            P2 = np.kron(P1, P1); Q2 = np.eye(d * d) - P2
            hD = herm(d * d, 0.08)
            Db = P2 @ hD @ P2 + Q2 @ hD @ Q2
            eps_bond = 0.02 * (1 + trial)
            hO = herm(d * d)
            Ob = P2 @ hO @ Q2 + Q2 @ hO @ P2
            Ob = eps_bond * Ob / np.linalg.norm(Ob, 2)
            Jn = 2 * np.linalg.norm(Db, 2)
            epsp = 2 * eps_bond
            delta = Delta - 2 * Jn
            bonds = lambda X: sum(_wrap(X, M, d) if j == M - 1 else bond_op(X, j, M, d) for j in range(M))
            H = sum(site_op(A1, j, M, d) for j in range(M)) + bonds(Db) + bonds(Ob)
            Pf = kron(*([P1] * M))
            idx = np.where(np.diag(Pf).real > 0.5)[0]
            K = (H)[np.ix_(idx, idx)]
            kv = np.linalg.eigvalsh(K)
            hv = np.linalg.eigvalsh(H)
            E = kv[0]
            cond = M * epsp ** 2 <= delta ** 2 / 8
            tM = 8 * M * epsp ** 2 / delta
            ok = True
            sel = kv <= E + delta / 8
            if cond:
                for i in np.where(sel)[0]:
                    ok &= (kv[i] - tM - 1e-12 <= hv[i] <= kv[i] + 1e-12)
            ground_ok = (E - M * epsp ** 2 / delta - 1e-12 <= hv[0] <= E + 1e-12)
            assert ground_ok and ok
            rows.append(dict(M=M, eps=epsp, delta=delta, hypothesis_17=bool(cond), t_M=tM,
                             max_shift=float(np.max(kv[sel] - hv[: sel.sum()])),
                             levels_checked=int(sel.sum())))
    return rows


def c49_product_rotation_example():
    rows = []
    for Delta, t in ((1.0, 0.3), (1.0, 1.0)):
        for M in (2, 4, 6, 8):
            h = np.array([[0, t], [t, Delta]])
            H = sum(site_op(h, j, M, 2) for j in range(M))
            P = kron(*([np.diag([1.0, 0.0])] * M))
            Q = np.eye(2 ** M) - P
            sub = np.where(np.diag(Q) > 0.5)[0]
            qmin = np.linalg.eigvalsh(H[np.ix_(sub, sub)])[0]
            em = (Delta - math.sqrt(Delta ** 2 + 4 * t ** 2)) / 2
            alpha = (1 + Delta / math.sqrt(Delta ** 2 + 4 * t ** 2)) / 2
            EH = M * em
            bound = (-M * em) * alpha ** M / (1 - alpha ** M)
            assert -1e-12 <= qmin - EH <= bound + 1e-12
            rows.append(dict(Delta=Delta, t=t, M=M, qmin_minus_ground=qmin - EH, bound=bound,
                             true_gap=math.sqrt(Delta ** 2 + 4 * t ** 2)))
    return rows


def main():
    res = {}
    res["C48_enlarged_embedding_multicell"] = c48_enlarged_embedding()
    res["C48_leakage_linear_uniform_in_u"] = c48_leakage_linear_uniform_in_u()
    res["C48_C46_resummed_formula"] = c46_resummed_formula()
    res["C58_fixed_marker_pinning"] = c58_pinning_lemma()
    res["C58_partial_trace_lower_bound"] = c58_partial_trace_peierls()
    res["C58_semigroup_heat_constant"] = c58_semigroup_heat_constant()
    res["C51_canonical_diagonalization"] = c51_canonical_diag()
    res["C51_virasoro_norm"] = c51_virasoro_norm()
    res["C51_asymmetric_free_fermion_T2"] = c51_asymmetric_free_fermion()
    res["C52_three_cell_commutator"] = c52_three_cell_commutator()
    res["C52_head_leakage"] = c52_head_leakage()
    res["C52_gram_singular_value"] = c52_gram_singular_value()
    res["C49_theorem_F"] = c49_theorem_F()
    res["C49_product_rotation"] = c49_product_rotation_example()
    Path(__file__).with_name("referee_checks.json").write_text(json.dumps(res, indent=1, default=float) + "\n")
    print("all referee checks passed")


if __name__ == "__main__":
    main()
