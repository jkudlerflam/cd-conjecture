#!/usr/bin/env python3
"""Logical (effective) chains and their calibration models, with three estimators of c.

Models (MODEL argument):
  heis            spin-1/2 Heisenberg, h = S.S                          (c = 1, SU(2)_1)
  uls             SU(3) Uimin-Lai-Sutherland, h = P_3                   (c = 2, SU(3)_1)
  blbq:T          spin-1 cos(T pi) S.S + sin(T pi) (S.S)^2              (T = theta/pi)
  eff:LAM         route-A logical chain from data/heff_c8*.npz at coupling LAM, two-cell bond h2 only
  eff3:LAM        the same plus the connected three-cell term h3 (next-nearest-neighbour, O(lam^4))
  amp:LAM:K       (P - 1/3) + K [eff3:LAM - (P - 1/3)]: the measured SU(3) breaking amplified K times
  file:PATH       a bond (r = 2) or three-site term (r = 3) saved as 'h' in an npz (with 'r', 'd')
Effective Hamiltonians are divided by their swap coefficient b (so the logical exchange is ~1).

Modes:
  fes   VUMPS (damped, n-site cell) over a chi list: e, S on every bond, xi, bond energies.
  obc   two-site finite DMRG on an open chain of length L: entanglement profile S(l), bond energies.
  ed    periodic chains, Lanczos with momentum labels: E0(L), low levels, momenta.

usage:  effchain.py MODEL fes --n 6 --chis 16 24 32 48 64 96 128 --out data/fes_MODEL.json
        effchain.py MODEL obc --L 96 --chis 64 128 256 --out ...
        effchain.py MODEL ed  --Ls 6 9 12 --nev 24 --out ...
Reuses contrib/d2_followup_220_1282/{vumps.py, fdmrg.py}.  Written by Claude (Opus 5.5), 2026-10-27.
"""
import os, sys, json, time, argparse
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'd2_followup_220_1282'))
import numpy as np
from scipy import linalg
from scipy.sparse.linalg import LinearOperator, eigsh
import vumps as V
from cell import swap, hermitian_basis

S12 = [np.array([[0, 1], [1, 0]]) / 2, np.array([[0, -1j], [1j, 0]]) / 2, np.diag([0.5, -0.5])]
S1 = [np.array([[0, 1, 0], [1, 0, 1], [0, 1, 0]]) / np.sqrt(2),
      np.array([[0, -1j, 0], [1j, 0, -1j], [0, 1j, 0]]) / np.sqrt(2), np.diag([1.0, 0, -1])]


def heff_file(lam):
    """h2 (and h3 if available) at coupling lam, searching data/heff_c8*.npz."""
    import glob
    h2 = h3 = None
    for f in sorted(glob.glob(os.path.join(HERE, 'data', 'heff_c8*.npz'))):
        z = np.load(f)
        for k in z.files:
            if k == f'h2_{lam}' and h2 is None:
                h2 = z[k]
            if k == f'h3_{lam}':
                h3 = z[k]
    return h2, h3


def model(name):
    """Returns h (d^r x d^r), r, d, c_expected, 2-site bond used for bond energies, note."""
    if name == 'heis':
        h = sum(np.kron(s, s) for s in S12).real.astype(complex)
        return h, 2, 2, 1.0, h, 'Heisenberg S=1/2'
    if name == 'uls':
        h = swap(3).astype(complex)
        return h, 2, 3, 2.0, h, 'SU(3) ULS'
    if name.startswith('blbq:'):
        t = float(name.split(':')[1]) * np.pi
        ss = sum(np.kron(s, s) for s in S1)
        h = np.cos(t) * ss + np.sin(t) * ss @ ss
        return h, 2, 3, None, h, f'BLBQ theta={name.split(":")[1]} pi'
    if name.startswith('eff:') or name.startswith('eff3:'):
        lam = float(name.split(':')[1])
        h2, h3 = heff_file(lam)
        P = swap(3)
        Pc = P - np.eye(9) / 3
        b = np.trace(Pc @ h2).real / np.trace(Pc @ Pc).real
        h2n = (h2 - np.trace(h2) / 9 * np.eye(9)) / b
        if name.startswith('eff3:'):
            assert h3 is not None, 'no h3 for this lambda'
            h3n = (h3 - np.trace(h3) / 27 * np.eye(27)) / b
            h = np.kron(h2n, np.eye(3)) + h3n
            return h, 3, 3, None, h2n, f'route-A logical chain h2+h3, lam={lam}, b={b:.6e}'
        return h2n, 2, 3, None, h2n, f'route-A logical chain h2 only, lam={lam}, b={b:.6e}'
    if name.startswith('amp:'):
        # amplified breaking: h = (P - 1/3) (x) 1 + kappa [h_eff3(lam) - (P - 1/3) (x) 1]
        _, lam, kap = name.split(':')
        h, r, d, c, hb, note = model(f'eff3:{lam}')
        base = np.kron(swap(3) - np.eye(9) / 3, np.eye(3))
        kap = float(kap)
        hb2 = (swap(3) - np.eye(9) / 3) + kap * (hb - (swap(3) - np.eye(9) / 3))
        return base + kap * (h - base), 3, 3, None, hb2, f'amplified route-A breaking, lam={lam}, kappa={kap}'
    if name.startswith('file:'):
        z = np.load(name[5:])
        h = z['h']; r = int(z['r']); d = int(z['d'])
        hb = z['hb'] if 'hb' in z.files else (h if r == 2 else None)
        return h, r, d, None, hb, f'file {name[5:]}'
    raise ValueError(name)


def gens(d):
    if d == 2:
        return [2 * s for s in S12]
    return hermitian_basis(d)[1:]


# ---------------------------------------------------------------- VUMPS measurements
def bond_energies_umps(psi, hb):
    d = psi.AL[0].shape[1]
    H = hb.reshape(d, d, d, d)
    out = []
    for k in range(psi.n):
        AC = np.tensordot(psi.AL[k], psi.C[(k + 1) % psi.n], ([2], [0]))
        th = np.tensordot(AC, psi.AR[(k + 1) % psi.n], ([2], [0]))           # a s1 s2 b
        nrm = np.vdot(th, th).real
        Ht = np.einsum('stuv,auvb->astb', H, th)
        out.append(float(np.vdot(th, Ht).real / nrm))
    return out


def run_fes(a, h, r, d, hb):
    W = V.mpo_from_term(h, r, d=d)
    data = dict(model=a.model, mode='fes', n=a.n, results=[])
    t0 = time.time()

    def save(psi, rec):
        rec['time_s'] = time.time() - t0
        if hb is not None:
            rec['ebond'] = bond_energies_umps(psi, hb)
        rec['loc'] = V.local_expect(psi, gens(d)).tolist()
        if a.blocks:
            rec['Sblock'] = V.block_entropies(psi, 0, a.blocks, every=max(1, a.blocks // 24))
        data['results'].append(rec)
        json.dump(data, open(a.out + '.tmp', 'w'))
        os.replace(a.out + '.tmp', a.out)
        print(f"{a.model} chi={rec['chi']} e={rec['e']:.12f} err={rec['err']:.1e} it={rec['iters']} "
              f"S={np.round(rec['S'], 5).tolist()} xi={rec['xi']:.3f} t={rec['time_s']:.0f}s "
              f"eb={np.round(rec.get('ebond', []), 6).tolist()}", flush=True)
    V.run([W] * a.n, a.n, a.chis, d=d, tol=a.tol, maxit=a.maxit, seed=a.seed, mode=a.vmode,
          log=lambda rec: None, on_chi=save, stol=a.stol)


# ---------------------------------------------------------------- finite DMRG
def run_obc(a, h, r, d, hb):
    from fdmrg import FDMRG
    W = V.mpo_from_term(h, r, d=d)
    dm = FDMRG(W, a.L, d=d, seed=a.seed)
    data = dict(model=a.model, mode='obc', L=a.L, results=[])
    t0 = time.time()
    for chi in a.chis:
        Eprev = None
        for sw in range(a.sweeps):
            E, tr = dm.sweep(chi, tol=a.tol)
            if Eprev is not None and abs(E - Eprev) < a.etol * a.L and sw >= 2:
                break
            Eprev = E
        rec = dict(chi=chi, E=E, sweeps=sw + 1, trunc=tr, S=dm.S.tolist(), time_s=time.time() - t0)
        if hb is not None and r == 2:
            X = dm.left_norm_envs()
            Hb = hb.reshape(d, d, d, d)
            eb = []
            for j in range(a.L - 1):
                th = np.tensordot(dm.M[j], dm.M[j + 1], ([2], [0]))
                t = np.tensordot(X[j], th, ([1], [0]))
                eb.append(float(np.tensordot(th.conj(), np.einsum('stuv,auvb->astb', Hb, t), 4).real))
            rec['ebond'] = eb
        data['results'].append(rec)
        json.dump(data, open(a.out + '.tmp', 'w'))
        os.replace(a.out + '.tmp', a.out)
        print(f"{a.model} L={a.L} chi={chi} E={E:.12f} sweeps={sw+1} trunc={tr:.1e} Smid={dm.S[a.L//2-1]:.5f} "
              f"t={rec['time_s']:.0f}s", flush=True)


# ---------------------------------------------------------------- periodic ED
def ring_op(h, r, d, L):
    dim = d ** L
    hr = h.reshape(d ** r, d ** r)

    def mv(x):
        x = x.reshape([d] * L)
        out = np.zeros_like(x)
        for j in range(L):
            ax = [(j + i) % L for i in range(r)]
            y = np.moveaxis(x, ax, list(range(r))).reshape(d ** r, -1)
            y = (hr @ y).reshape([d] * L)
            out += np.moveaxis(y, list(range(r)), ax)
        return out.ravel()
    return LinearOperator((dim, dim), matvec=mv, dtype=complex)


def run_ed(a, h, r, d, hb):
    data = dict(model=a.model, mode='ed', results=[])
    for L in a.Ls:
        t0 = time.time()
        op = ring_op(h, r, d, L)
        k = min(a.nev, d ** L - 2)
        w, U = eigsh(op, k=k, which='SA', tol=1e-12, ncv=max(2 * k + 10, 40), maxiter=100000)
        idx = np.argsort(w); w, U = w[idx], U[:, idx]
        # momenta: diagonalize translation inside clusters of degenerate levels
        Tm = lambda v: np.moveaxis(v.reshape([d] * L), 0, L - 1).ravel()
        TU = np.array([Tm(U[:, i]) for i in range(k)]).T
        Tmat = U.conj().T @ TU
        mom = np.full(k, np.nan)
        i = 0
        while i < k:
            j = i + 1
            while j < k and abs(w[j] - w[j - 1]) < 1e-7 * max(1, abs(w[i])):
                j += 1
            if j < k:                       # only clusters known to be complete
                ev = linalg.eigvals(Tmat[i:j, i:j])
                mom[i:j] = np.sort(np.angle(ev) * L / (2 * np.pi))
            i = j
        rec = dict(L=L, E=w.tolist(), k=mom.tolist(), time_s=time.time() - t0)
        data['results'].append(rec)
        json.dump(data, open(a.out + '.tmp', 'w'))
        os.replace(a.out + '.tmp', a.out)
        print(f"{a.model} L={L} E0/L={w[0]/L:.12f} gaps*L/2pi=" +
              ' '.join(f'{(e - w[0]) * L / (2 * np.pi):.4f}[k={q:+.1f}]' for e, q in zip(w[1:10], mom[1:10])) +
              f" t={time.time()-t0:.0f}s", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('model')
    ap.add_argument('mode', choices=['fes', 'obc', 'ed'])
    ap.add_argument('--n', type=int, default=6)
    ap.add_argument('--chis', type=int, nargs='+', default=[16, 24, 32, 48, 64])
    ap.add_argument('--vmode', default='damped')
    ap.add_argument('--maxit', type=int, default=2000)
    ap.add_argument('--tol', type=float, default=1e-9)
    ap.add_argument('--stol', type=float, default=1e-7)
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--blocks', type=int, default=0)
    ap.add_argument('--L', type=int, default=48)
    ap.add_argument('--sweeps', type=int, default=12)
    ap.add_argument('--etol', type=float, default=1e-10)
    ap.add_argument('--Ls', type=int, nargs='+', default=[6, 9])
    ap.add_argument('--nev', type=int, default=24)
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    h, r, d, c, hb, note = model(a.model)
    if a.out is None:
        a.out = os.path.join(HERE, 'data', f"{a.mode}_{a.model.replace(':', '_').replace('/', '_')}.json")
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    print(note, flush=True)
    dict(fes=run_fes, obc=run_obc, ed=run_ed)[a.mode](a, h, r, d, hb)


if __name__ == '__main__':
    main()
