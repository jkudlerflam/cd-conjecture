#!/usr/bin/env python3
# Copied unchanged (apart from this line) from ../counterexample_search/idmrg.py on 2026-10-27 for contrib/d2_two_fluid.
"""Two-site infinite DMRG (McCulloch growth with wavefunction prediction) for a
translation-invariant nearest-neighbour bond h on C^d. numpy/scipy only.

The converged state has a two-site unit cell, so a spontaneously dimerized
ground state shows up as unequal energies on the two inequivalent bonds.
Outputs per bond dimension chi: energy per site, both bond energies,
dimerization, both entanglement entropies, correlation length from the
two-site transfer matrix.
"""
import os
for _k in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_k, '8')
import json, time, argparse
import numpy as np
from scipy import linalg
from scipy.sparse.linalg import LinearOperator, eigsh, eigs


def mpo_from_bond(h, d, tol=1e-12):
    """W[a,b,s',s] with h = sum_i A_i (x) B_i ; left state 0, right state w-1."""
    M = h.reshape(d, d, d, d).transpose(0, 2, 1, 3).reshape(d * d, d * d)   # (s1' s1),(s2' s2)
    U, S, Vh = linalg.svd(M)
    r = int(np.sum(S > tol * S[0]))
    A = [(U[:, i] * np.sqrt(S[i])).reshape(d, d) for i in range(r)]
    B = [(Vh[i, :] * np.sqrt(S[i])).reshape(d, d) for i in range(r)]
    w = r + 2
    W = np.zeros((w, w, d, d), dtype=complex)
    W[0, 0] = np.eye(d)
    W[w - 1, w - 1] = np.eye(d)
    for i in range(r):
        W[0, i + 1] = A[i]
        W[i + 1, w - 1] = B[i]
    return W


def heff_matvec(Lenv, W, R, x, shape):
    th = x.reshape(shape)                       # (a, s1, s2, b)
    t = np.tensordot(Lenv, th, axes=([2], [0]))           # (a, w, s1, s2, b)
    t = np.tensordot(t, W, axes=([1, 2], [0, 3]))         # (a, s2, b, w', s1')
    t = np.tensordot(t, W, axes=([3, 1], [0, 3]))         # (a, b, s1', w'', s2')
    t = np.tensordot(t, R, axes=([3, 1], [1, 2]))         # (a, s1', s2', b)
    return t.reshape(-1)


def grow_L(Lenv, A, W):
    t = np.tensordot(Lenv, A, axes=([2], [0]))            # (a, w, s, b')
    t = np.tensordot(t, W, axes=([1, 2], [0, 3]))         # (a, b', w', s')
    t = np.tensordot(A.conj(), t, axes=([0, 1], [0, 3]))  # (b, b', w')
    return t.transpose(0, 2, 1)


def grow_R(R, B, W):
    t = np.tensordot(B, R, axes=([2], [2]))               # (a', s, b, w)
    t = np.tensordot(t, W, axes=([1, 3], [3, 1]))         # (a', b, w, s')   W[w,w'',s',s] contract s and w''
    t = np.tensordot(B.conj(), t, axes=([1, 2], [3, 1]))  # (a, a', w)
    return t.transpose(0, 2, 1)


def bond_energy(theta, h, d):
    a, _, _, b = theta.shape
    t = theta.transpose(1, 2, 0, 3).reshape(d * d, a * b)
    return float(np.real(np.vdot(t, h @ t)) / np.real(np.vdot(t, t)))


def entropy(lam):
    p = lam ** 2
    p = p[p > 1e-16]
    p = p / p.sum()
    return float(-np.sum(p * np.log(p)))


def run(h, d, chis, steps_per_chi=200, tol_e=1e-9, seed=0, log=None, lanczos_tol=1e-10):
    W = mpo_from_bond(h, d)
    w = W.shape[0]
    Lenv = np.zeros((1, w, 1), complex); Lenv[0, 0, 0] = 1
    R = np.zeros((1, w, 1), complex); R[0, w - 1, 0] = 1
    rng = np.random.default_rng(seed)
    lam_prev = np.ones(1)
    A = B = None
    lam = np.ones(1)
    Eprev = 0.0
    results = []
    nstep = 0
    for chi in chis:
        e_hist = []
        for it in range(steps_per_chi):
            nstep += 1
            a_dim, b_dim = Lenv.shape[0], R.shape[0]
            shape = (a_dim, d, d, b_dim)
            if A is not None and A.shape[2] == len(lam) and len(lam_prev) == B.shape[2] and B.shape[2] == A.shape[0]:
                # prediction: (Lambda B) Lambda_prev^{-1} (A Lambda)
                LB = lam[:, None, None] * B
                AL = A * lam[None, None, :]
                inv = np.where(lam_prev > 1e-7, 1 / np.maximum(lam_prev, 1e-300), 0.0)
                guess = np.tensordot(LB * inv[None, None, :], AL, axes=([2], [0]))
                if guess.shape != shape:
                    guess = rng.normal(size=shape) + 0j
            else:
                guess = rng.normal(size=shape) + 0j
            n = int(np.prod(shape))
            op = LinearOperator((n, n), matvec=lambda x: heff_matvec(Lenv, W, R, x, shape), dtype=complex)
            if n <= 64:
                Hd = np.array([op.matvec(col) for col in np.eye(n)]).T
                ev, vv = linalg.eigh((Hd + Hd.conj().T) / 2)
                E, th = ev[0], vv[:, 0]
            else:
                g = guess.reshape(-1)
                if not np.all(np.isfinite(g)) or np.linalg.norm(g) < 1e-12:
                    g = rng.normal(size=n) + 0j
                ev, vv = eigsh(op, k=2, which='SA', v0=g, tol=lanczos_tol, ncv=min(n - 1, 30))
                i0 = int(np.argmin(ev))
                E, th = ev[i0], vv[:, i0]
            th = th.reshape(a_dim * d, d * b_dim)
            U, S, Vh = linalg.svd(th, full_matrices=False)
            keep = min(chi, int(np.sum(S > 1e-14)))
            trunc = float(np.sum(S[keep:] ** 2))
            U, S, Vh = U[:, :keep], S[:keep], Vh[:keep]
            S = S / np.linalg.norm(S)
            A_new = U.reshape(a_dim, d, keep)
            B_new = Vh.reshape(keep, d, b_dim)
            Lenv = grow_L(Lenv, A_new, W)
            R = grow_R(R, B_new, W)
            lam_prev, lam = lam, S
            A, B = A_new, B_new
            e_site = (E - Eprev) / 2
            Eprev = E
            e_hist.append(e_site)
            # period-2 convergence (a dimerized state alternates strong and weak centre bonds)
            if it > 20 and abs(e_hist[-1] - e_hist[-3]) < tol_e and abs(e_hist[-2] - e_hist[-4]) < tol_e:
                break
        # observables on the converged two-site cell: Lp Ga L Gb Lp
        res = measure(A, B, lam, lam_prev, h, d)
        e2 = float(np.real(e_hist[-1] + e_hist[-2]) / 2)
        res.update(chi=chi, steps=it + 1, e_site=e2, trunc=trunc, nstep=nstep,
                   e_bonds=(res['eAB'] + res['eBA']) / 2, consistent=bool(abs(e2 - (res['eAB'] + res['eBA']) / 2) < 1e-3))
        results.append(res)
        if log:
            print(json.dumps(res), flush=True)
    return results


def measure(A, B, lam, lam_prev, h, d):
    """Two-site unit cell in Vidal form: Ga = lam_prev^{-1} A, Gb = B lam_prev^{-1}."""
    k = len(lam_prev)
    good = lam_prev > 1e-10
    inv = np.where(good, 1 / np.maximum(lam_prev, 1e-300), 0.0)
    Ga = inv[:, None, None] * A
    Gb = B * inv[None, None, :]
    # bond AB: lam_prev Ga lam Gb lam_prev
    thAB = np.tensordot(lam_prev[:, None, None] * Ga * lam[None, None, :], Gb * lam_prev[None, None, :], axes=([2], [0]))
    # bond BA: lam Gb lam_prev Ga lam
    thBA = np.tensordot(lam[:, None, None] * Gb * lam_prev[None, None, :], Ga * lam[None, None, :], axes=([2], [0]))
    eAB = bond_energy(thAB, h, d)
    eBA = bond_energy(thBA, h, d)
    # transfer matrix of the cell, right-canonical form: Gb lam_prev Ga lam  -> use X = Ga lam Gb lam_prev
    X = np.tensordot(Ga * lam[None, None, :], Gb * lam_prev[None, None, :], axes=([2], [0]))  # (a, s1, s2, b)
    chi = X.shape[0]
    Xm = X.reshape(chi, d * d, chi)

    def tm(v):
        v = v.reshape(chi, chi)
        t = np.tensordot(Xm, v, axes=([2], [0]))           # (a, s, b')
        t = np.tensordot(t, Xm.conj(), axes=([1, 2], [1, 2]))  # (a, a')
        return t.reshape(-1)
    try:
        nev = min(4, chi * chi - 2)
        if chi * chi <= 64:
            T = np.array([tm(c) for c in np.eye(chi * chi)]).T
            vals = np.sort(np.abs(linalg.eigvals(T)))[::-1]
        else:
            op = LinearOperator((chi * chi, chi * chi), matvec=tm, dtype=complex)
            vals = np.sort(np.abs(eigs(op, k=nev, which='LM', tol=1e-8, return_eigenvectors=False)))[::-1]
        xi = float(-2 / np.log(vals[1] / vals[0])) if vals[1] > 1e-14 else 0.0
        tmvals = (vals / vals[0]).tolist()
    except Exception as exc:
        xi, tmvals = float('nan'), [str(exc)]
    return dict(eAB=eAB, eBA=eBA, dimer=eAB - eBA, S_AB=entropy(lam), S_BA=entropy(lam_prev), xi=xi, tm=tmvals[:4])


if __name__ == '__main__':
    import models
    p = argparse.ArgumentParser()
    p.add_argument('--model', required=True)
    p.add_argument('--theta', type=float, default=0.0)
    p.add_argument('--coeff', type=float, nargs='+')
    p.add_argument('--chis', type=int, nargs='+', default=[16, 32, 64, 96, 128])
    p.add_argument('--steps', type=int, default=300)
    p.add_argument('--out', required=True)
    p.add_argument('--label', default='')
    a = p.parse_args()
    if a.model == 'sp4':
        m = models.sp4(); d = 4
        h = np.cos(a.theta) * m['hs'][0] + np.sin(a.theta) * m['hs'][1]
    elif a.model == 'heis':
        m = models.spin_s(0.5); d = 2; h = m['X']
    elif a.model == 'uls':
        from ced import swap
        d = int(a.coeff[0]); h = swap(d)
    elif a.model == 'orbital':
        m = models.spin_orbital(); d = 4
        h = sum(c * x for c, x in zip(a.coeff, m['hs']))
    elif a.model == 'spin32':
        m = models.spin_s(1.5); d = 4
        h = sum(c * x for c, x in zip(a.coeff, m['P']))
    elif a.model == 'spin1':
        m = models.spin_s(1); d = 3
        h = sum(c * x for c, x in zip(a.coeff, m['P']))
    else:
        raise ValueError(a.model)
    h = (h + h.conj().T) / 2
    t = time.time()
    res = run(h, d, a.chis, steps_per_chi=a.steps, log=True)
    with open(a.out, 'a') as f:
        for r in res:
            r.update(model=a.model, theta=a.theta, coeff=a.coeff, label=a.label)
            f.write(json.dumps(r) + '\n')
    print('done %.0fs' % (time.time() - t))
