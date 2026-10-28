#!/usr/bin/env python3
"""The full route-A chain: every cell is ONE site of dimension D = 48 carrying the exact C8 cell
Hamiltonian H0, and neighbouring cells are coupled by lam * sum_a FR_a(j) FL_a(j+1).
    H = sum_j H0_j + lam sum_j sum_a FR_a(j) FL_a(j+1)        (no projection, no perturbation theory)

The logical exchange (~ 3.3e-3 lam^2) is 1e-4 to 1e-6 of the cell gap, so plain Lanczos inside
DMRG/VUMPS stalls.  Every local eigenproblem is solved here by Davidson with the preconditioner
(H0_local + delta)^{-1} on the physical legs (H0 is diagonal in the basis used), which removes the
cell-excitation scale.  The shift delta is of the order of the logical exchange, so that the logical
sector is rescaled to O(1) as well and Davidson converges as for an ordinary d = 3 problem.

Modes
  fes   VUMPS on the infinite chain, n-site unit cell, over a chi list (finite-entanglement scaling).
  obc   two-site finite DMRG on an open chain of L cells, for comparison with the logical chain.
usage: bigsite.py LAM fes --n 2 --chis 9 16 24 32 48 64 --out data/big_fes_LAM.json
       bigsite.py LAM obc --L 12 --chis 27 54 81 --out data/big_obc_LAM_L12.json
Written by Claude (Opus 5.5) for contrib/routeA_numerics, 2026-10-27.
"""
import os, sys, json, time, argparse
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'd2_followup_220_1282'))
import numpy as np
from scipy import linalg
import vumps as V
import fdmrg as F
from cell import load_c8, hermitian_basis

DELTA = 0.3
EVALS = None          # cell energies in the working basis, set in main()


def davidson(apply, x0, prec, tol=1e-10, maxit=3000, mmax=40, keep=4):
    """Lowest eigenpair by preconditioned Davidson with thick restart (incremental Rayleigh matrix)."""
    shp = x0.shape
    n = x0.size
    x = x0.ravel().astype(complex)
    if np.linalg.norm(x) < 1e-300:
        x = np.random.default_rng(0).normal(size=n) + 0j
    Vm = np.zeros((mmax + 1, n), complex)
    Am = np.zeros((mmax + 1, n), complex)
    Hs = np.zeros((mmax + 1, mmax + 1), complex)
    Vm[0] = x / np.linalg.norm(x)
    Am[0] = apply(Vm[0].reshape(shp)).ravel()
    Hs[0, 0] = np.vdot(Vm[0], Am[0])
    m = 1
    theta, rn = None, np.inf
    for it in range(maxit):
        w, s = linalg.eigh(0.5 * (Hs[:m, :m] + Hs[:m, :m].conj().T))
        theta = w[0]
        u = s[:, 0] @ Vm[:m]
        r = s[:, 0] @ Am[:m] - theta * u
        rn = np.linalg.norm(r)
        if rn < tol:
            break
        if m >= mmax:                                  # thick restart on the lowest Ritz vectors
            S = s[:, :keep]
            Vm[:keep] = S.T @ Vm[:m]; Am[:keep] = S.T @ Am[:m]
            Hs[:keep, :keep] = np.diag(w[:keep])
            m = keep
        # Olsen correction: t = M^-1 r - eps M^-1 u with u^dagger t = 0 (avoids stagnation when M ~ H - theta)
        Mr = prec(r.reshape(shp), theta).ravel()
        Mu = prec(u.reshape(shp), theta).ravel()
        t = Mr - (np.vdot(u, Mr) / np.vdot(u, Mu)) * Mu
        for _ in range(2):
            t -= Vm[:m].T @ (Vm[:m].conj() @ t)
        nt = np.linalg.norm(t)
        if nt < 1e-14:
            break
        Vm[m] = t / nt
        Am[m] = apply(Vm[m].reshape(shp)).ravel()
        Hs[:m + 1, m] = Vm[:m + 1].conj() @ Am[m]
        Hs[m, :m] = Hs[:m, m].conj()
        m += 1
    return theta, (u / np.linalg.norm(u)).reshape(shp), it, rn


PAIR = None          # eigen-decomposition of the two-cell Hamiltonian, set in main()


def prec_pair(r, theta):
    """(h_pair - eps_0 + DELTA)^{-1} on the two physical legs of a two-site tensor (a, s1, s2, b),
    h_pair = H0 (x) 1 + 1 (x) H0 + lam V: exact for everything inside the pair."""
    eps, U = PAIR
    a, d1, d2, b = r.shape
    x = np.tensordot(U.conj().T, r.transpose(1, 2, 0, 3).reshape(d1 * d2, a * b), axes=1)
    x = x / (eps - eps[0] + DELTA)[:, None]
    return (U @ x).reshape(d1, d2, a, b).transpose(2, 0, 1, 3)


def prec_phys(axes):
    """(H0 on the physical legs `axes` + DELTA)^{-1}."""
    def p(r, theta):
        den = np.full(r.shape, DELTA)
        for ax in axes:
            sh = [1] * r.ndim; sh[ax] = len(EVALS)
            den = den + EVALS.reshape(sh)
        return r / den
    return p


STATS = dict(calls=0, iters=0, maxres=0.0)
MAXIT = 200          # local solves need not be exact in early sweeps; residuals are recorded


def lowest_patched(op_apply, x0, tol):
    n = x0.size
    if x0.ndim == 3 and x0.shape[1] == len(EVALS):
        e, v, it, rn = davidson(op_apply, x0, prec_phys([1]), tol=max(tol, 1e-12), maxit=MAXIT)
        STATS['calls'] += 1; STATS['iters'] += it; STATS['last_res'] = float(rn)
        return e, v
    if n <= 4096:
        H = np.zeros((n, n), complex)
        for i in range(n):
            ee = np.zeros(n, complex); ee[i] = 1
            H[:, i] = op_apply(ee.reshape(x0.shape)).ravel()
        w, U = linalg.eigh((H + H.conj().T) / 2)
        return w[0], U[:, 0].reshape(x0.shape)
    e, v, it, rn = davidson(op_apply, x0, lambda r, t: r, tol=max(tol, 1e-12))
    return e, v


def two_site_patched(Lenv, W1, W2, Renv, th, tol):
    f = lambda v: V.hac2(Lenv, W1, W2, Renv, v)
    e, v, it, rn = davidson(f, th, prec_pair, tol=max(tol, 1e-11), maxit=MAXIT)
    STATS['calls'] += 1; STATS['iters'] += it; STATS['maxres'] = max(STATS['maxres'], float(rn))
    return e, v


def big_bond(cell, lam):
    D = cell['D']
    H0 = cell['H0']
    I = np.eye(D)
    h = 0.5 * (np.kron(H0, I) + np.kron(I, H0))
    for a in range(3):
        h = h + lam * np.kron(cell['FR'][a], cell['FL'][a])
    return h


def logical_projectors(cell):
    """Generators of su(3) on the code space (the first N basis states of the cell)."""
    N, D = cell['N'], cell['D']
    out = []
    for T in hermitian_basis(N)[1:]:
        O = np.zeros((D, D), complex); O[:N, :N] = T; out.append(O)
    Pc = np.zeros((D, D), complex); Pc[:N, :N] = np.eye(N)
    return out, Pc


def main():
    global EVALS, PAIR
    ap = argparse.ArgumentParser()
    ap.add_argument('lam', type=float)
    ap.add_argument('mode', choices=['fes', 'obc'])
    ap.add_argument('--n', type=int, default=2)
    ap.add_argument('--chis', type=int, nargs='+', default=[9, 16, 24, 32])
    ap.add_argument('--L', type=int, default=12)
    ap.add_argument('--sweeps', type=int, default=10)
    ap.add_argument('--maxit', type=int, default=600)
    ap.add_argument('--tol', type=float, default=1e-9)
    ap.add_argument('--vmode', default='vumps')
    ap.add_argument('--delta', type=float, default=None, help='preconditioner shift (default 10 x 3.3e-3 lam^2)')
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    global DELTA
    DELTA = a.delta if a.delta is not None else 10 * 3.3e-3 * a.lam ** 2
    cell = load_c8()
    EVALS = cell['e'].copy()
    D = cell['D']
    h = big_bond(cell, a.lam)
    W = V.mpo_from_term(h, 2, d=D)
    hp = np.kron(cell['H0'], np.eye(D)) + np.kron(np.eye(D), cell['H0'])
    for x in range(3):
        hp = hp + a.lam * np.kron(cell['FR'][x], cell['FL'][x])
    PAIR = linalg.eigh(hp)
    V.lowest = lowest_patched
    F.two_site_min = two_site_patched
    gens, Pc = logical_projectors(cell)
    data = dict(lam=a.lam, mode=a.mode, results=[])
    t0 = time.time()
    if a.mode == 'fes':
        def save(psi, rec):
            rec['time_s'] = time.time() - t0
            rec['ebond'] = []
            Hb = h.reshape(D, D, D, D)
            for k in range(psi.n):
                AC = np.tensordot(psi.AL[k], psi.C[(k + 1) % psi.n], ([2], [0]))
                th = np.tensordot(AC, psi.AR[(k + 1) % psi.n], ([2], [0]))
                rec['ebond'].append(float(np.vdot(th, np.einsum('stuv,auvb->astb', Hb, th)).real / np.vdot(th, th).real))
            rec['code_weight'] = V.local_expect(psi, [Pc]).tolist()
            rec['loc'] = V.local_expect(psi, gens).tolist()
            rec['davidson'] = dict(STATS)
            data['results'].append(rec)
            json.dump(data, open(a.out + '.tmp', 'w')); os.replace(a.out + '.tmp', a.out)
            print(f"big lam={a.lam} chi={rec['chi']} e={rec['e']:.14f} err={rec['err']:.1e} it={rec['iters']} "
                  f"S={np.round(rec['S'], 5).tolist()} xi={rec['xi']:.3f} t={rec['time_s']:.0f}s "
                  f"eb={[round(x, 10) for x in rec['ebond']]} dav={STATS}", flush=True)
        # start from a random logical uMPS embedded in the code space (zero weight on cell excitations)
        # plus a small random admixture of cell excitations (an exact code-only uMPS is a stationary
        # point of the single-site update, because every port compresses to zero on the code)
        p3 = V.UMPS.random(a.n, cell['N'], a.chis[0], 0)
        rng = np.random.default_rng(5)
        pad = lambda A: np.concatenate([A, 0.05 * (rng.normal(size=(A.shape[0], D - A.shape[1], A.shape[2]))
                                                   + 1j * rng.normal(size=(A.shape[0], D - A.shape[1], A.shape[2])))], axis=1)
        AL = [np.linalg.qr(pad(x).reshape(-1, x.shape[2]))[0].reshape(x.shape[0], D, x.shape[2]) for x in p3.AL]
        AR = [np.linalg.qr(pad(x).reshape(x.shape[0], -1).T)[0].T.reshape(x.shape[0], D, x.shape[2]).conj() for x in p3.AR]
        psi0 = V.UMPS(AL, AR, p3.C)
        V.run([W] * a.n, a.n, a.chis, d=D, tol=a.tol, maxit=a.maxit, mode=a.vmode, log=lambda r: None,
              on_chi=save, psi=psi0)
    else:
        Ws = [W] * a.L
        dm = F.FDMRG(W, a.L, d=D, seed=1, chi0=4, edge=0.5 * cell['H0'])
        # restrict the random initial MPS to the code space, then restore right-canonical form
        M = dm.M
        for i in range(a.L):
            M[i][:, cell['N']:, :] = 0
        for i in range(a.L - 1, 0, -1):
            q, r = np.linalg.qr(M[i].reshape(M[i].shape[0], -1).T)
            M[i] = q.T.reshape(-1, D, M[i].shape[2])
            M[i - 1] = np.tensordot(M[i - 1], r.T, ([2], [0]))
        M[0] /= np.linalg.norm(M[0])
        for i in range(a.L - 1, 0, -1):
            dm.Re[i] = V.rt(dm.Re[i + 1], M[i], dm.Ws[i])
        for chi in a.chis:
            Eprev = None
            for sw in range(a.sweeps):
                E, tr = dm.sweep(chi, tol=1e-10)
                print(f"  sweep {sw} chi={chi} E={E:.14f} trunc={tr:.1e} dav={STATS} t={time.time()-t0:.0f}s", flush=True)
                STATS['maxres'] = 0.0
                if Eprev is not None and abs(E - Eprev) < 1e-13 * a.L:
                    break
                Eprev = E
            rec = dict(chi=chi, E=E, trunc=tr, S=dm.S.tolist(), sweeps=sw + 1, time_s=time.time() - t0,
                       code_weight=dm.local_ops([Pc])[0].tolist())
            data['results'].append(rec)
            json.dump(data, open(a.out + '.tmp', 'w')); os.replace(a.out + '.tmp', a.out)
            print(f"big lam={a.lam} L={a.L} chi={chi} E={E:.14f} Smid={dm.S[a.L//2-1]:.6f}", flush=True)


if __name__ == '__main__':
    main()
