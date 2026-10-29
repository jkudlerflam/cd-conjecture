"""Independent referee checks for the collaborator lemmas imported by C58.

numpy and scipy only. Nothing is imported from the collaborator or Astra
check scripts. Run in this folder with

    OPENBLAS_NUM_THREADS=1 python3 -B referee_imported_checks.py

It writes referee_imported_checks.json. These are small toy instances that
test the identities and inequalities the proofs rely on. They are not
simulations of the d0 = 4464 chain.
"""
from __future__ import annotations

import itertools
import json
import math
from pathlib import Path

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh

RNG = np.random.default_rng(20261029)


def opnorm_sparse(m):
    """Operator norm of a Hermitian sparse matrix."""
    m = sp.csr_matrix(m)
    if m.nnz == 0:
        return 0.0
    vals = eigsh(m, k=1, which="LM", return_eigenvectors=False, tol=1e-12)
    return float(abs(vals[0]))


def maxabs(m):
    if sp.issparse(m):
        m = sp.csr_matrix(m)
        return float(abs(m.data).max()) if m.nnz else 0.0
    return float(np.max(np.abs(m))) if m.size else 0.0


def logsumexp(x):
    x = np.asarray(x, dtype=float)
    top = np.max(x)
    return float(top + np.log(np.exp(x - top).sum()))


# ---------------------------------------------------------------------------
# 1. SEPARATE_JW: separate payload and length strings with graded auxiliaries
# ---------------------------------------------------------------------------

def separate_jw_checks():
    ncell = 3
    a = sp.csr_matrix(np.array([[0., 1.], [0., 0.]]))   # lowering, a|1>=|0>
    z = sp.csr_matrix(np.diag([1., -1.]))
    i2 = sp.identity(2, format="csr")
    # Auxiliary space: vacuum Omega plus four excited states with integer
    # (payload, length) charges, including a negative length charge.
    aux_charges = [(0, 0), (0, 0), (1, 0), (0, 1), (-1, 1)]
    daux = len(aux_charges)
    iaux = sp.identity(daux, format="csr")
    pi_loc = sp.diags([(-1.) ** qp for qp, _ in aux_charges]).tocsr()
    lam_loc = sp.diags([(-1.) ** ql for _, ql in aux_charges]).tocsr()
    dcell = 2 * 2 * daux
    icell = sp.identity(dcell, format="csr")
    dim = dcell ** ncell

    def cellop(payload=None, length=None, aux=None):
        return sp.kron(sp.kron(payload if payload is not None else i2,
                               length if length is not None else i2),
                       aux if aux is not None else iaux, format="csr")

    def at(op, j):
        ops = [icell] * ncell
        ops[j] = op
        out = ops[0]
        for o in ops[1:]:
            out = sp.kron(out, o, format="csr")
        return out

    c = [at(cellop(payload=a), j) for j in range(ncell)]
    b = [at(cellop(length=a), j) for j in range(ncell)]
    p = [at(cellop(payload=z), j) for j in range(ncell)]
    q = [at(cellop(length=z), j) for j in range(ncell)]
    pi = [at(cellop(aux=pi_loc), j) for j in range(ncell)]
    lam = [at(cellop(aux=lam_loc), j) for j in range(ncell)]
    one = sp.identity(dim, format="csr")
    P = [p[j] @ pi[j] for j in range(ncell)]
    Q = [q[j] @ lam[j] for j in range(ncell)]

    def prod(ops):
        out = one
        for o in ops:
            out = out @ o
        return out

    S = [prod(P[:j]) for j in range(ncell)]
    T = [prod(Q[:j]) for j in range(ncell)]
    K = prod(P)
    psi = [S[j] @ c[j] for j in range(ncell)]
    chi = [K @ T[j] @ b[j] for j in range(ncell)]
    fields = psi + chi

    car = 0.0
    for i, f in enumerate(fields):
        for k, g in enumerate(fields):
            car = max(car, maxabs(f @ g + g @ f))
            anti = f @ g.T + g.T @ f
            if i == k:
                anti = anti - one
            car = max(car, maxabs(anti))

    def homogeneous(o):
        comps = {}
        for r in (0, 1):
            for s in (0, 1):
                acc = np.zeros_like(o)
                for u in (0, 1):
                    for v in (0, 1):
                        g = (np.linalg.matrix_power(pi_loc.toarray(), u)
                             @ np.linalg.matrix_power(lam_loc.toarray(), v))
                        acc = acc + (-1) ** (r * u + s * v) * g @ o @ g
                comps[(r, s)] = acc / 4
        return comps

    def xi(o_rs, j, r, s):
        out = at(cellop(aux=sp.csr_matrix(o_rs)), j)
        pre = one
        if r:
            pre = pre @ S[j] @ p[j]
        if s:
            pre = pre @ K @ pi[j] @ T[j] @ q[j]
        return pre @ out

    raw = [RNG.normal(size=(daux, daux)) + 1j * RNG.normal(size=(daux, daux))
           for _ in range(2)]
    comps = [homogeneous(o) for o in raw]
    graded, homomorph, isometry = 0.0, 0.0, 0.0
    images = {}
    for j in range(ncell):
        for idx, cp in enumerate(comps):
            for (r, s), o_rs in cp.items():
                images[(j, idx, r, s)] = xi(o_rs, j, r, s)
    for (j, idx, r, s), img in images.items():
        sign = (-1) ** (r + s)
        for f in fields:
            graded = max(graded, maxabs(f @ img - sign * img @ f))
        for (j2, idx2, u, v), img2 in images.items():
            if j2 != j:
                sgn = (-1) ** ((r + s) * (u + v))
                graded = max(graded, maxabs(img @ img2 - sgn * img2 @ img))
            elif idx2 != idx:
                target = xi(comps[idx][(r, s)] @ comps[idx2][(u, v)], j,
                            r ^ u, s ^ v)
                homomorph = max(homomorph, maxabs(img @ img2 - target))
    # Isometry of Xi on a full (inhomogeneous) auxiliary operator.
    for j in range(ncell):
        herm = raw[0] + raw[0].conj().T
        cp = homogeneous(herm)
        img = sum(xi(o, j, r, s) for (r, s), o in cp.items())
        isometry = max(isometry, abs(opnorm_sparse(img) - np.max(np.abs(np.linalg.eigvalsh(herm)))))

    # Local charges for neutrality.
    pay_cell = np.zeros(dcell)
    len_cell = np.zeros(dcell)
    for np_ in (0, 1):
        for nl in (0, 1):
            for k, (qp, ql) in enumerate(aux_charges):
                idx = np_ * 2 * daux + nl * daux + k
                pay_cell[idx] = np_ + qp
                len_cell[idx] = nl + ql

    def neutral_random(ncells_x):
        d = dcell ** ncells_x
        m = RNG.normal(size=(d, d)) + 1j * RNG.normal(size=(d, d))
        m = m + m.conj().T
        tp = np.zeros(d)
        tl = np.zeros(d)
        for word in itertools.product(range(dcell), repeat=ncells_x):
            idx = 0
            for w in word:
                idx = idx * dcell + w
            tp[idx] = sum(pay_cell[w] for w in word)
            tl[idx] = sum(len_cell[w] for w in word)
        mask = (tp[:, None] == tp[None, :]) & (tl[:, None] == tl[None, :])
        return m * mask

    def embed2(op2, first):
        if first == 0:
            return sp.kron(sp.csr_matrix(op2), icell, format="csr")
        return sp.kron(icell, sp.csr_matrix(op2), format="csr")

    outside_gens = lambda k: ([psi[k], chi[k]] +
                              [images[key] for key in images if key[0] == k])
    locality = 0.0
    for first, outside in ((0, 2), (1, 0)):
        o = embed2(neutral_random(2), first)
        for g in outside_gens(outside):
            locality = max(locality, maxabs(o @ g - g @ o))
    # Charge-exchange vertex: logical payload lowering times auxiliary raising.
    raise_aux = np.zeros((daux, daux))
    raise_aux[2, 0] = 1.0                 # Omega -> payload-charge-one state
    vertex = at(cellop(payload=a, aux=sp.csr_matrix(raise_aux)), 1)
    vertex = vertex + vertex.T
    for g in outside_gens(2) + outside_gens(0):
        locality = max(locality, maxabs(vertex @ g - g @ vertex))
    # Negative control: even under combined parity, odd under each parity.
    bad = at(cellop(payload=a.T, length=a), 1)
    bad = bad + bad.T
    nonlocal_control = maxabs(bad @ psi[2] - psi[2] @ bad)

    # Reference extension E_JW and the K_D cancellation.
    nmodes = 2 * ncell
    dl = 2 ** nmodes
    ai = np.array([[0., 1.], [0., 0.]])
    zi = np.diag([1., -1.])

    def lop(op, mode):
        ops = [np.eye(2)] * nmodes
        ops[mode] = op
        out = ops[0]
        for o in ops[1:]:
            out = np.kron(out, o)
        return out

    lc = [lop(ai, 2 * j) for j in range(ncell)]
    lb = [lop(ai, 2 * j + 1) for j in range(ncell)]
    lp = [lop(zi, 2 * j) for j in range(ncell)]
    lq = [lop(zi, 2 * j + 1) for j in range(ncell)]
    lone = np.eye(dl)
    K0 = lone.copy()
    for j in range(ncell):
        K0 = K0 @ lp[j]
    psi0 = []
    chi0 = []
    for j in range(ncell):
        sj = lone.copy()
        tj = lone.copy()
        for k in range(j):
            sj = sj @ lp[k]
            tj = tj @ lq[k]
        psi0.append(sj @ lc[j])
        chi0.append(K0 @ tj @ lb[j])
    # Mode ordering: (psi_0, chi_0, psi_1, chi_1, ...).
    f0 = []
    ffull = []
    for j in range(ncell):
        f0 += [psi0[j], chi0[j]]
        ffull += [psi[j], chi[j]]
    maj0 = []
    majf = []
    for f, g in zip(f0, ffull):
        maj0 += [f + f.T, 1j * (f - f.T)]
        majf += [g + g.T, 1j * (g - g.T)]
    monos0 = []
    subsets = []
    for mask in range(2 ** (2 * nmodes)):
        sub = [i for i in range(2 * nmodes) if (mask >> i) & 1]
        m = lone.astype(complex)
        for i in sub:
            m = m @ maj0[i]
        monos0.append(m)
        subsets.append(sub)

    def ejw(o_log):
        # coef_S = Tr(Gamma_S^dag O)/dl
        coef = np.array([np.trace(monos0[i].conj().T @ o_log) / dl
                         for i in range(len(monos0))])
        keep = np.flatnonzero(np.abs(coef) > 1e-13)
        recon = sum(coef[i] * monos0[i] for i in keep)
        out = sp.csr_matrix((dim, dim), dtype=complex)
        support_modes = set()
        for i in keep:
            m = one.astype(complex)
            for k in subsets[i]:
                m = m @ majf[k]
                support_modes.add(k // 4)
            out = out + coef[i] * m
        return out, float(np.max(np.abs(recon - o_log))), support_modes

    # Physical cell: 4 code states |pay,len> then 4 excited states.
    phys_charges = [(0, 0), (0, 1), (1, 0), (1, 1)] + aux_charges[1:]
    dphys = len(phys_charges)
    pcode = np.diag([1.] * 4 + [0.] * 4)
    J = np.zeros((dcell, dphys))
    for np_ in (0, 1):
        for nl in (0, 1):
            J[np_ * 2 * daux + nl * daux + 0, np_ * 2 + nl] = 1.0
    for k in range(1, daux):
        J[0 * 2 * daux + 0 * daux + k, 3 + k] = 1.0
    out = {"CAR_error": car, "graded_relation_error": graded,
           "Xi_homomorphism_error": homomorph, "Xi_isometry_error": isometry,
           "neutral_locality_error": locality,
           "nonneutral_parity_even_commutator": nonlocal_control}
    kd_rows = []
    for first in (0, 1):
        d2 = dphys ** 2
        m = RNG.normal(size=(d2, d2)) + 1j * RNG.normal(size=(d2, d2))
        m = m + m.conj().T
        px = np.kron(pcode, pcode)
        m = px @ m @ px + (np.eye(d2) - px) @ m @ (np.eye(d2) - px)
        tp = np.array([phys_charges[i][0] + phys_charges[k][0]
                       for i in range(dphys) for k in range(dphys)])
        tl = np.array([phys_charges[i][1] + phys_charges[k][1]
                       for i in range(dphys) for k in range(dphys)])
        m = m * ((tp[:, None] == tp[None, :]) & (tl[:, None] == tl[None, :]))
        bnorm = float(np.max(np.abs(np.linalg.eigvalsh(m))))
        cidx = [i * dphys + k for i in range(4) for k in range(4)]
        o2 = m[np.ix_(cidx, cidx)]           # P_X B P_X on code tensor code
        onorm = float(np.max(np.abs(np.linalg.eigvalsh(o2))))
        o_log = np.kron(o2, np.eye(4)) if first == 0 else np.kron(np.eye(4), o2)
        e_op, recon_err, modes = ejw(o_log)
        jj = np.kron(J, J)
        jbj = embed2(jj @ m @ jj.T, first)
        kd = jbj - e_op
        # Projector onto auxiliary vacuum at both cells of X.
        vac = np.zeros((daux, daux))
        vac[0, 0] = 1.
        pv = [at(cellop(aux=sp.csr_matrix(vac)), j) for j in range(ncell)]
        pvx = pv[first] @ pv[first + 1]
        annih = max(maxabs(kd @ pvx), maxabs(pvx @ kd))
        kdnorm = opnorm_sparse(kd)
        enorm = opnorm_sparse(e_op)
        # Sectorwise unitary equivalence: E_JW in each auxiliary basis
        # configuration has the spectrum of O_X tensor identity.
        ref = np.sort(np.linalg.eigvalsh(o_log))
        sector_err = 0.0
        e_dense_rows = e_op.tocsr()
        for auxcfg in itertools.product(range(daux), repeat=ncell):
            idx = []
            for logword in itertools.product((0, 1), repeat=nmodes):
                ix = 0
                for j in range(ncell):
                    ix = ix * dcell + (logword[2 * j] * 2 * daux
                                       + logword[2 * j + 1] * daux + auxcfg[j])
                idx.append(ix)
            blk = e_dense_rows[idx][:, idx].toarray()
            sector_err = max(sector_err, float(np.max(np.abs(
                np.sort(np.linalg.eigvalsh(blk)) - ref))))
        # E_JW, J B J^dag and K_D are local on X.
        outside = 2 if first == 0 else 0
        loc = 0.0
        for g in outside_gens(outside):
            loc = max(loc, maxabs(e_op @ g - g @ e_op), maxabs(jbj @ g - g @ jbj))
        kd_rows.append({"X": [first, first + 1], "B_norm": bnorm,
                        "O_X_norm": onorm, "E_JW_norm": enorm,
                        "K_D_norm": kdnorm, "K_D_over_2B": kdnorm / (2 * bnorm),
                        "two_sided_vacuum_annihilation": annih,
                        "sector_spectrum_error": sector_err,
                        "Majorana_reconstruction_error": recon_err,
                        "coefficient_modes_cells": sorted(modes),
                        "locality_error": loc})
        assert annih < 1e-10 and kdnorm <= 2 * bnorm + 1e-9
        assert abs(enorm - onorm) < 1e-9 and sector_err < 1e-9 and loc < 1e-10
    out["K_D_rows"] = kd_rows

    # Equation (9): CAR hoppings in tensor variables, and the naive extension.
    hop_psi = psi[0].T @ psi[1]
    hop_chi = chi[0].T @ chi[1]
    eq9 = max(maxabs(hop_psi - c[0].T @ p[0] @ pi[0] @ c[1]),
              maxabs(hop_chi - b[0].T @ lam[0] @ b[1]))
    naive = c[0].T @ p[0] @ c[1] + b[0].T @ b[1]
    naive_gap = opnorm_sparse((hop_psi + hop_chi) - naive + ((hop_psi + hop_chi) - naive).T)
    out["eq9_identity_error"] = eq9
    out["naive_identity_extension_deviation_norm"] = naive_gap

    # Equation (13): normal-ordered coefficient bound on r = 2 cells.
    m4 = 4
    d4 = 2 ** m4
    loc_ops = []
    for mode in range(m4):
        ops = [np.eye(2)] * m4
        for k in range(mode):
            ops[k] = zi
        ops[mode] = ai
        out_op = ops[0]
        for o in ops[1:]:
            out_op = np.kron(out_op, o)
        loc_ops.append(out_op)
    monos = []
    degs = []
    for A in range(2 ** m4):
        for B in range(2 ** m4):
            mm = np.eye(d4)
            for k in range(m4):
                if (A >> k) & 1:
                    mm = mm @ loc_ops[k].T
            for k in reversed(range(m4)):
                if (B >> k) & 1:
                    mm = mm @ loc_ops[k]
            monos.append(mm.reshape(-1))
            degs.append(bin(A).count("1") + bin(B).count("1"))
    Mat = np.array(monos).T
    degs = np.array(degs)
    ratios = []
    for _ in range(10):
        o = RNG.normal(size=(d4, d4)) + 1j * RNG.normal(size=(d4, d4))
        o = o + o.conj().T
        coef = np.linalg.solve(Mat, o.reshape(-1))
        on = float(np.max(np.abs(np.linalg.eigvalsh(o))))
        for H in (1.0, 2.0, 3.0):
            lhs = float(np.sum(np.abs(coef) * H ** degs))
            ratios.append(lhs / ((8 * H * H) ** m4 * on))
    out["eq13_max_ratio_to_bound"] = max(ratios)
    assert max(ratios) <= 1.0
    assert car < 1e-12 and graded < 1e-12 and homomorph < 1e-12
    assert locality < 1e-12 and nonlocal_control > 0.5 and eq9 < 1e-12
    return out


# ---------------------------------------------------------------------------
# 2. LENGTH_FUGACITY Section 2
# ---------------------------------------------------------------------------

def length_fugacity_classical():
    """Interacting renewal model: pressure zero equals fixed-L growth rate."""
    rows = []
    lmax = 40
    ell = np.arange(1, lmax + 1)
    for beta, coupling in ((1.0, 0.4), (3.0, 0.8), (6.0, 1.5)):
        g = np.where(ell <= 2, 1, 2) * 1.0
        energy = 0.7 * np.abs(ell - 4) + 0.2 * np.sin(ell)
        v = coupling * np.cos(0.7 * ell[:, None] - 1.3 * ell[None, :] + 0.4)

        def tmat(s):
            w = g * np.exp(-beta * (energy + s * ell))
            return w[:, None] * np.exp(-beta * v)

        def pcell(s):
            lam = np.max(np.abs(np.linalg.eigvals(tmat(s))))
            return math.log(lam) / beta

        lo, hi = -0.5, 0.5
        while pcell(lo) <= 0:
            lo -= 0.5
        while pcell(hi) >= 0:
            hi += 0.5
        assert pcell(lo) > 0 > pcell(hi)
        for _ in range(200):
            mid = (lo + hi) / 2
            if pcell(mid) > 0:
                lo = mid
            else:
                hi = mid
        sstar = (lo + hi) / 2
        # Monotonicity (8) with ell_min = 1.
        mono = min(pcell(s) - ell[0] * t - pcell(s + t)
                   for s in np.linspace(sstar - 1, sstar + 1, 9) for t in (0.05, 0.3))
        # Fixed-L open sum over all cell numbers, rescaled by exp(-beta s* L).
        L = 4000
        w = g * np.exp(-beta * (energy + sstar * ell))
        boltz = np.exp(-beta * v)       # boltz[l_prev, l_next]
        f = np.zeros((L + 1, lmax))
        logscale = np.zeros(L + 1)
        for tot in range(1, L + 1):
            vec = np.zeros(lmax)
            if tot <= lmax:
                vec[tot - 1] += w[tot - 1]
            for k in range(lmax):
                prev = tot - ell[k]
                if prev >= 1:
                    vec[k] += w[k] * (f[prev] @ boltz[:, k]) * math.exp(logscale[prev] - logscale[tot - 1])
            norm = vec.sum()
            f[tot] = vec / norm
            logscale[tot] = logscale[tot - 1] + math.log(norm)
        # log S_L e^{-beta s* L} = logscale[L]; its growth per unit length.
        growth = (logscale[L] - logscale[L // 2]) / (beta * (L - L // 2))
        # Control: S_L grows like exp(beta s* L), so rescaling by a shifted
        # tilt leaves a visible residual growth rate.
        shifted = 0.0
        ws = g * np.exp(-beta * (energy + (sstar + 0.01) * ell))
        logs = [0.0] * (L + 1)
        fvec = np.zeros((L + 1, lmax))
        for tot in range(1, L // 4 + 1):
            vec = np.zeros(lmax)
            if tot <= lmax:
                vec[tot - 1] += ws[tot - 1]
            for k in range(lmax):
                prev = tot - ell[k]
                if prev >= 1:
                    vec[k] += ws[k] * (fvec[prev] @ boltz[:, k]) * math.exp(logs[prev] - logs[tot - 1])
            nrm = vec.sum()
            fvec[tot] = vec / nrm
            logs[tot] = logs[tot - 1] + math.log(nrm)
        shifted = (logs[L // 4] - logs[L // 8]) / (beta * (L // 4 - L // 8))
        assert abs(shifted + 0.01) < 1e-6
        rows.append({"beta": beta, "coupling": coupling, "s_star": sstar,
                     "fixed_L_growth_minus_s_star": growth,
                     "control_growth_with_s_star_plus_0.01": shifted,
                     "monotonicity_slack_min": mono})
        assert abs(growth) < 5e-3 and mono >= -1e-12
    return rows


def length_fugacity_quantum():
    """Finite quantum checks of (7), (8) and (9b) with length-transferring bonds."""
    lengths = [1, 2, 2, 3, 3]          # cell basis: one length-1 state, two of each other
    dc = len(lengths)
    A = np.zeros((dc, dc))
    for l in set(lengths):
        idx = [i for i in range(dc) if lengths[i] == l]
        blk = RNG.normal(size=(len(idx), len(idx)))
        A[np.ix_(idx, idx)] = (blk + blk.T) / 4 + np.eye(len(idx)) * (0.5 * abs(l - 2) + 0.3)
    pair_len = np.array([lengths[i] + lengths[k] for i in range(dc) for k in range(dc)])
    Bm = RNG.normal(size=(dc * dc, dc * dc))
    Bm = (Bm + Bm.T) / 2
    Bm = Bm * (pair_len[:, None] == pair_len[None, :])
    Bm *= 0.6 / np.max(np.abs(np.linalg.eigvalsh(Bm)))
    ccut = 0.6
    # Length transfer is present.
    transfer = 0.0
    for i in range(dc * dc):
        for k in range(dc * dc):
            if abs(Bm[i, k]) > 1e-12 and lengths[i // dc] != lengths[k // dc]:
                transfer = max(transfer, abs(Bm[i, k]))
    assert transfer > 1e-3

    memo = {}

    def words_with_length(M, L):
        if (M, L) in memo:
            return memo[(M, L)]
        if M == 0:
            res = [()] if L == 0 else []
        elif L < M or L > 3 * M:
            res = []
        else:
            res = [w + (x,) for x in range(dc) for w in words_with_length(M - 1, L - lengths[x])]
        memo[(M, L)] = res
        return res

    def sector_h(M, L):
        states = words_with_length(M, L)
        if not states:
            return None
        index = {w: i for i, w in enumerate(states)}
        h = np.zeros((len(states), len(states)))
        for col, w in enumerate(states):
            for j in range(M):
                for x2 in range(dc):
                    if A[x2, w[j]] != 0:
                        w2 = w[:j] + (x2,) + w[j + 1:]
                        h[index[w2], col] += A[x2, w[j]]
            for j in range(M - 1):
                pin = w[j] * dc + w[j + 1]
                for pout in np.flatnonzero(Bm[:, pin]):
                    w2 = w[:j] + (pout // dc, pout % dc) + w[j + 2:]
                    h[index[w2], col] += Bm[pout, pin]
        return h

    def logz_sector(M, L, beta):
        h = sector_h(M, L)
        if h is None:
            return -np.inf
        return logsumexp(-beta * np.linalg.eigvalsh(h))

    results = {"superadditivity_9b": [], "cut_7": []}
    Lmax = 11
    for beta in (0.5, 2.0, 5.0):
        logS = {}
        for L in range(1, Lmax + 1):
            vals = [logz_sector(M, L, beta) for M in range(1, L + 1)]
            logS[L] = logsumexp(vals)
        worst = np.inf
        for L1 in range(1, Lmax):
            for L2 in range(1, Lmax - L1 + 1):
                worst = min(worst, logS[L1 + L2] - (logS[L1] + logS[L2] - beta * ccut))
        results["superadditivity_9b"].append({"beta": beta, "min_slack": float(worst)})
        assert worst >= -1e-10
        for s in (-0.3, 0.0, 0.4):
            def logzM(M):
                vals = [logz_sector(M, L, beta) - beta * s * L for L in range(M, 3 * M + 1)]
                return logsumexp(vals)
            worst7 = np.inf
            for M, N in ((1, 1), (1, 2), (2, 2), (1, 3), (2, 3), (3, 3), (2, 4)):
                dev = abs(logzM(M + N) - logzM(M) - logzM(N))
                worst7 = min(worst7, beta * ccut - dev)
            results["cut_7"].append({"beta": beta, "s": s, "min_slack": float(worst7)})
            assert worst7 >= -1e-10
    return results


def length_fugacity_9d():
    """Lemma 2.2: zero of p0 + A T^2 + r_T moves by A(s0) T^2 / ell0."""
    ell0, s0 = 2.5, 0.1
    p0 = lambda s: -ell0 * (s - s0) - 0.8 * (s - s0) ** 2 + 0.3 * abs(s - s0) ** 1.5
    Afun = lambda s: 0.9 + 2.0 * s
    rows = []
    for T in (0.1, 0.03, 0.01, 0.003):
        r = lambda s: T ** 3 * math.sin(s / T)
        f = lambda s: p0(s) + Afun(s) * T * T + r(s)
        lo, hi = s0 - 0.5, s0 + 0.5
        for _ in range(200):
            mid = (lo + hi) / 2
            if f(mid) > 0:
                lo = mid
            else:
                hi = mid
        sT = (lo + hi) / 2
        pred = s0 + Afun(s0) * T * T / ell0
        rows.append({"T": T, "relative_error_over_T2": (sT - pred) / T ** 2})
    for row in rows:
        assert abs(row["relative_error_over_T2"]) <= 5 * row["T"]
    return rows


# ---------------------------------------------------------------------------
# 3. TILTED_MATCHING as used in C58 Section 5
# ---------------------------------------------------------------------------

def tilted_matching_checks():
    out = {}
    rows = []
    for b in (4, 5, 9, 30):
        Delta = 1.0
        kappa = Delta / (b + 2)
        sa = min(kappa / 4, Delta / (4 * b))
        ells = np.arange(1, 50 * b)
        A = np.maximum(Delta, kappa * (ells + 1))   # saturates C58 (18)
        qv = ells - b
        bound_ok = np.all(np.abs(qv) <= b + A / kappa + 1e-12)
        ratio = min(np.min((A + sgn * sa * qv) / A) for sgn in (1, -1))
        Ca = b / Delta + 1 / kappa
        crit = min(np.min(A[qv > 0] / qv[qv > 0]), np.min(A[qv < 0] / (-qv[qv < 0])))
        rows.append({"b": b, "s_a": sa, "TILTED_MATCHING_s_star": 1 / (2 * Ca),
                     "min_ReA_s_over_A": float(ratio),
                     "accretivity_loss_tilt": float(crit),
                     "abs_q_bound_holds": bool(bound_ok)})
        assert bound_ok and ratio >= 0.5 - 1e-12 and sa <= 1 / (2 * Ca) + 1e-15
    out["strip"] = rows

    # Complex-tilt Schatten-Hoelder trace bound with noncommuting insertions.
    d = 40
    b = 5
    Delta = 1.0
    kappa = Delta / (b + 2)
    sa = min(kappa / 4, Delta / (4 * b))
    ells = RNG.integers(1, 3 * b, size=d)
    ells[0] = b                              # vacuum-like state
    A = np.maximum(Delta, kappa * (ells + 1)) + RNG.uniform(0, 1, size=d)
    A[0] = 0.0
    qv = (ells - b).astype(float)
    qv[0] = 0.0
    worst = np.inf
    worst_half = np.inf
    for trial in range(200):
        s = complex(RNG.uniform(-sa, sa), RNG.uniform(-5, 5))
        beta = RNG.uniform(1, 20)
        n = RNG.integers(1, 5)
        cuts = np.sort(RNG.uniform(0, beta, size=n - 1))
        times = np.diff(np.concatenate([[0], cuts, [beta]]))
        ops = [RNG.normal(size=(d, d)) + 1j * RNG.normal(size=(d, d)) for _ in range(n)]
        lhs = np.eye(d, dtype=complex)
        for t, o in zip(times, ops):
            lhs = lhs @ np.diag(np.exp(-t * (A + s * qv))) @ o
        lhs = abs(np.trace(lhs))
        prodnorm = np.prod([np.linalg.norm(o, 2) for o in ops])
        rhs = np.sum(np.exp(-beta * (A + s.real * qv))) * prodnorm
        rhs_half = np.sum(np.exp(-beta * A / 2)) * prodnorm
        worst = min(worst, rhs - lhs)
        worst_half = min(worst_half, rhs_half - rhs)
    out["complex_hoelder_min_slack"] = float(worst)
    out["half_heat_trace_domination_min_slack"] = float(worst_half)
    assert worst >= -1e-9 and worst_half >= -1e-9

    # Uniform resolvent bound on the strip and Cauchy derivative bounds for a
    # Feshbach-type massive influence map W(s) = P K Q (a_s + QKQ)^{-1} Q K P.
    nlog = 2
    dq = d - 1
    Kc = RNG.normal(size=(dq, nlog)) + 1j * RNG.normal(size=(dq, nlog))
    Kc *= 0.1 / np.linalg.norm(Kc, 2)
    Kq = RNG.normal(size=(dq, dq)) + 1j * RNG.normal(size=(dq, dq))
    Kq = (Kq + Kq.conj().T)
    Kq *= 0.1 / np.linalg.norm(Kq, 2)
    aq = A[1:]
    qq = qv[1:]

    def W(s):
        R = np.linalg.inv(np.diag(aq + s * qq) + Kq)
        return Kc.conj().T @ R @ Kc, R

    uniform = 0.1 ** 2 / (Delta / 2 - 0.1)
    sup = 0.0
    for _ in range(400):
        s = complex(RNG.uniform(-sa, sa), RNG.uniform(-50, 50))
        sup = max(sup, np.linalg.norm(W(s)[0], 2))
    assert sup <= uniform
    deriv_rows = []
    for s0 in (-sa / 2, 0.0, sa / 2):
        _, R = W(s0)
        Qd = np.diag(qq)
        term = R.copy()
        for m in range(1, 6):
            term = term @ Qd @ R
            dm = math.factorial(m) * np.linalg.norm(Kc.conj().T @ term @ Kc, 2)
            bound = math.factorial(m) * (2 / sa) ** m * uniform
            deriv_rows.append({"s0": s0, "m": m, "ratio_to_cauchy_bound": dm / bound})
            assert dm <= bound
    out["strip_sup_over_uniform_bound"] = sup / uniform
    out["derivatives"] = deriv_rows
    return out


# ---------------------------------------------------------------------------
# 4. C44 pin structure as used in C58 Section 3
# ---------------------------------------------------------------------------

ARROWS = {"#": (1, 0), "L": (0, 0), "X": (0, 1), "R": (1, 1)}


def cells_of(word):
    n = len(word)
    marks = [i for i, ch in enumerate(word) if ch == "#"]
    return marks


def cell_flags(word, marks):
    """Mismatch count per ordered cell, including its outgoing bond."""
    n = len(word)
    out = []
    for k, start in enumerate(marks):
        nxt = marks[(k + 1) % len(marks)]
        size = (nxt - start) % n or n
        z = 0
        for off in range(size):
            i = (start + off) % n
            j = (i + 1) % n
            if ARROWS[word[i]][1] != ARROWS[word[j]][0]:
                z += 1
        out.append(z)
    return out


def pin_structure_checks():
    rewrites = 0
    pin_failures = 0
    legal_shape_failures = 0
    words = 0
    for n in range(2, 10):
        for tup in itertools.product("#LXR", repeat=n):
            if "#" not in tup:
                continue
            word = list(tup)
            words += 1
            marks = cells_of(word)
            flags = cell_flags(word, marks)
            m = len(marks)
            for k in range(m):
                start = marks[k]
                nxt = marks[(k + 1) % m]
                last = (nxt - 1) % n
                prev_last = (start - 1) % n
                # Outgoing pin of a good cell, incoming pin after a good predecessor.
                if flags[k] == 0 and word[last] not in "XR":
                    pin_failures += 1
                if flags[(k - 1) % m] == 0 and word[prev_last] not in "XR":
                    pin_failures += 1
                if flags[k] == 0:
                    size = (nxt - start) % n or n
                    cell = "".join(word[(start + o) % n] for o in range(size))
                    body = cell[1:]
                    ok = (cell[0] == "#" and body.count("X") == 1
                          and set(body.split("X")[0]) <= {"L"}
                          and set(body.split("X")[1]) <= {"R"})
                    if not ok:
                        legal_shape_failures += 1
            # All literal selector rewrites, applied at every bond.
            for i in range(n):
                j = (i + 1) % n
                pair = word[i] + word[j]
                new = None
                if pair == "LX":
                    new = "XR"
                elif pair == "XR":
                    new = "LX"
                elif pair == "R#" and n > 1:
                    new = "#L"
                elif pair == "#L":
                    new = "R#"
                if new is None:
                    continue
                w2 = word[:]
                w2[i], w2[j] = new[0], new[1]
                if pair in ("R#", "#L"):
                    # Ordered marker k moves by one site; keep its label.
                    k = marks.index(j if pair == "R#" else i)
                    m2 = marks[:]
                    m2[k] = i if pair == "R#" else j
                    if len(set(m2)) != len(m2):
                        continue
                    f2 = cell_flags(w2, m2)
                else:
                    f2 = cell_flags(w2, marks)
                rewrites += 1
                assert f2 == flags, (word, w2, flags, f2)
    # Fragment arithmetic and the kappa bound of C58 (18).
    frag_min = np.inf
    kappa_min = np.inf
    for _ in range(20000):
        U = RNG.uniform(0.1, 3)
        t = RNG.uniform(0.1, 1e3)
        Ksel = 4 * t + 9 * U
        for zc in range(1, 8):
            frag_min = min(frag_min, (Ksel * zc - (2 * t + 4 * U) * (zc + 1)) / U - 1)
        b = int(RNG.integers(4, 60))
        gamma = RNG.uniform(0.01, 5)
        e = RNG.uniform(0.01, 5)
        Delta = RNG.uniform(0, 1) * min(gamma, U, e)
        kap = Delta / (b + 2)
        for l in range(1, 20 * b):
            need = kap * (l + 1)
            legal = Delta if l in (b, b + 1) else (gamma * (1 + abs(l - b)) if l >= 2 else np.inf)
            illegal = U + e * l
            kappa_min = min(kappa_min, legal - need, illegal - need)
    assert pin_failures == 0 and legal_shape_failures == 0
    assert frag_min >= -1e-9 and kappa_min >= -1e-9
    return {"ring_words": words, "rewrites_checked": rewrites,
            "pin_enable_failures": pin_failures,
            "legal_cell_shape_failures": legal_shape_failures,
            "fragment_bound_min_slack_over_U": float(frag_min),
            "kappa_bound_min_slack": float(kappa_min)}


def main():
    result = {
        "separate_jw": separate_jw_checks(),
        "length_fugacity_classical": length_fugacity_classical(),
        "length_fugacity_quantum": length_fugacity_quantum(),
        "length_fugacity_lemma_2_2": length_fugacity_9d(),
        "tilted_matching": tilted_matching_checks(),
        "c44_pin_structure": pin_structure_checks(),
        "status": "all assertions passed",
    }
    path = Path(__file__).with_name("referee_imported_checks.json")
    path.write_text(json.dumps(result, indent=2, default=float) + "\n")
    print(json.dumps({k: (v if k == "status" else "ok") for k, v in result.items()}, indent=2))


if __name__ == "__main__":
    main()
