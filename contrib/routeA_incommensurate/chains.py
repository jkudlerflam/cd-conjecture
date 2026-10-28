#!/usr/bin/env python3
"""Logical (effective) chains and calibration chains, with VUMPS finite-entanglement scaling,
open-chain DMRG and ring exact diagonalization.

MODEL:
  xx            one complex fermion (= XX chain at filling 0.962/pi), c = 1
  xx2           two decoupled XX spin chains (spin (x) spin per site, no string between them),
                fillings and phases of INC, c = 2
  free:T        the target logical chain T in {INC, SYM, INC3} (JW across cells), c = n
  eff:T:LAM     h2(LAM) from two-cell ED of the full cell chain (all orders, data/heff_T.npz)
  eff3:T:LAM    the same plus the three-cell term h3(LAM)
All effective bonds are divided by J0 LAM^2 so the logical hopping is ~1.

MODE:
  fes   VUMPS (damped), n-site unit cell, chi list.  Records S on every bond, xi (sites), energy,
        bond energies and flavour densities per site (translation breaking), transfer spectrum.
  obc   two-site DMRG on an open chain of L sites: S(l), bond energies, densities.
  ed    rings of L sites with fermionic (JW) boundary terms in each parity sector: lowest levels and
        momenta (translation eigenvalues), for velocities and z.

usage:  chains.py MODEL fes --n 2 --chis 16 24 32 48 64 96 128 --out data/fes_MODEL.json
        chains.py MODEL obc --L 64 --chis 64 128 --out ...
        chains.py MODEL ed --Ls 6 8 10 --nev 12 --out ...
Reuses contrib/d2_followup_220_1282/{vumps.py, fdmrg.py} (read only).
Written by Claude (Opus 5.5) for contrib/routeA_incommensurate, 2026-10-27.
"""
import os, sys, json, time, argparse
for _k in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ.setdefault(_k, '4')
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'd2_followup_220_1282'))
import numpy as np
from scipy import linalg
from scipy.sparse.linalg import LinearOperator, eigsh
import vumps as V
import cell as CE


def model(name):
    """Returns h (d^r x d^r), r, d, n_flavours, number operators per site, note."""
    if name == 'xx':
        par = dict(t=[1.0], th=[0.0], kF=[0.962])
        h = CE.target_bond(par)
        a, P = CE.fermion_ops(1)
        return h, 2, 2, 1, [a[0].conj().T @ a[0]], 'XX chain, filling 0.962/pi'
    if name == 'xx2':
        par = CE.INC
        h1 = CE.target_bond(dict(t=[1.0], th=[par['th'][0]], kF=[par['kF'][0]]))
        h2 = CE.target_bond(dict(t=[1.0], th=[par['th'][1]], kF=[par['kF'][1]]))
        I2 = np.eye(2)
        # site = spin1 (x) spin2; bond ordering (s1, s2)_j (s1, s2)_{j+1}
        h = np.zeros((16, 16), complex)
        h1f = h1.reshape(2, 2, 2, 2); h2f = h2.reshape(2, 2, 2, 2)
        H1 = np.einsum('abcd,ef,gh->aebgcfdh', h1f, I2, I2).reshape(16, 16)
        H2 = np.einsum('ef,gh,abcd->eagbfchd', I2, I2, h2f).reshape(16, 16)
        h = H1 + H2
        n = np.diag([0, 1.0])
        return h, 2, 4, 2, [np.kron(n, I2), np.kron(I2, n)], 'two decoupled XX chains (INC fillings)'
    parts = name.split(':')
    par = getattr(CE, parts[1])
    nf = len(par['t'])
    a, P = CE.fermion_ops(nf)
    nops = [x.conj().T @ x for x in a]
    N = 2 ** nf
    if parts[0] == 'free':
        return CE.target_bond(par), 2, N, nf, nops, f'free target {parts[1]}'
    lam = float(parts[2])
    z = np.load(os.path.join(HERE, 'data', f'heff_{parts[1]}.npz'))
    J0 = float(z['J0'])
    h2 = z[f'h2_{lam}']
    h2 = (h2 - np.trace(h2) / N ** 2 * np.eye(N ** 2)) / (J0 * lam ** 2)
    if parts[0] == 'eff':
        return h2, 2, N, nf, nops, f'two-cell effective bond at lam={lam}'
    h3 = z[f'h3_{lam}'] / (J0 * lam ** 2)
    return np.kron(h2, np.eye(N)) + h3, 3, N, nf, nops, f'two- and three-cell effective terms at lam={lam}'


def bond_energies_umps(psi, h, d):
    out = []
    Hb = h.reshape(d, d, d, d)
    for k in range(psi.n):
        AC = np.tensordot(psi.AL[k], psi.C[(k + 1) % psi.n], ([2], [0]))
        th = np.tensordot(AC, psi.AR[(k + 1) % psi.n], ([2], [0]))
        out.append(float(np.vdot(th, np.einsum('stuv,auvb->astb', Hb, th)).real / np.vdot(th, th).real))
    return out


# ----------------------------------------------------------------------------- ring ED
def ring_ed(h, d, L, nev, pdiag):
    """Ring of L sites with the same bond h on every bond, including (L-1, 0).  This is the ring
    of cells itself; in fermion language the boundary condition is twisted by -Q, Q the global
    parity, which is conserved.  Each Q sector separately; momenta from the cyclic shift."""
    dim = d ** L
    hb = h.reshape(d, d, d, d)
    par = np.ones([d] * L)
    for j in range(L):
        sh = [1] * L; sh[j] = d
        par = par * pdiag.reshape(sh)
    par = par.ravel()

    def mv(x):
        x = x.reshape([d] * L + [-1])
        y = np.zeros_like(x)
        for j in range(L):
            j1 = (j + 1) % L
            t = np.tensordot(hb, x, ([2, 3], [j, j1]))
            y += np.moveaxis(t, [0, 1], [j, j1])
        return y.reshape(dim, -1)

    def shift(x):                     # (T x)(s_0..s_{L-1}) = x(s_1..s_{L-1}, s_0): site j -> j+1
        x = x.reshape([d] * L + [-1])
        return np.moveaxis(x, 0, L - 1).reshape(dim, -1)
    res = []
    for Q in (1, -1):
        sec = np.where(par == Q)[0]
        dsec = len(sec)

        def mvs(v):
            full = np.zeros((dim, v.shape[1] if v.ndim > 1 else 1), complex)
            full[sec] = v.reshape(dsec, -1)
            return mv(full)[sec]
        op = LinearOperator((dsec, dsec), matvec=lambda v: mvs(v)[:, 0], matmat=mvs, dtype=complex)
        if dsec <= 3000:
            Hd = op.matmat(np.eye(dsec, dtype=complex))
            w, U = linalg.eigh(0.5 * (Hd + Hd.conj().T))
            w, U = w[:nev], U[:, :nev]
        else:
            w, U = eigsh(op, k=nev, which='SA', tol=1e-11, ncv=max(3 * nev, 40))
            o = np.argsort(w); w, U = w[o], U[:, o]
        full = np.zeros((dim, U.shape[1]), complex); full[sec] = U
        TU = shift(full)[sec]
        ks = np.zeros(len(w))
        i = 0
        while i < len(w):                      # momenta inside clusters of (near-)degenerate levels
            j = i
            while j + 1 < len(w) and abs(w[j + 1] - w[i]) < 1e-8 * max(1, abs(w[i])):
                j += 1
            Mt = U[:, i:j + 1].conj().T @ TU[:, i:j + 1]
            ks[i:j + 1] = np.sort(np.angle(linalg.eigvals(Mt)))
            i = j + 1
        res.append(dict(Q=Q, E=w.tolist(), k=ks.tolist()))
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('model')
    ap.add_argument('mode', choices=['fes', 'obc', 'ed'])
    ap.add_argument('--n', type=int, default=2)
    ap.add_argument('--chis', type=int, nargs='+', default=[16, 24, 32, 48, 64, 96, 128])
    ap.add_argument('--L', type=int, default=64)
    ap.add_argument('--Ls', type=int, nargs='+', default=[6, 8, 10])
    ap.add_argument('--nev', type=int, default=12)
    ap.add_argument('--sweeps', type=int, default=8)
    ap.add_argument('--vmode', default='damped')
    ap.add_argument('--maxit', type=int, default=600)
    ap.add_argument('--tol', type=float, default=1e-9)
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    h, r, d, nf, nops, note = model(a.model)
    data = dict(model=a.model, mode=a.mode, note=note, args=vars(a), results=[])
    t0 = time.time()
    if a.mode == 'fes':
        W = V.mpo_from_term(h, r, d=d)
        h2site = h if r == 2 else None

        def save(psi, rec):
            rec['time_s'] = time.time() - t0
            if h2site is not None:
                rec['ebond'] = bond_energies_umps(psi, h2site, d)
            rec['dens'] = V.local_expect(psi, nops).tolist()
            data['results'].append(rec)
            json.dump(data, open(a.out + '.tmp', 'w')); os.replace(a.out + '.tmp', a.out)
            print(f"{a.model} chi={rec['chi']} e={rec['e']:.12f} err={rec['err']:.1e} it={rec['iters']} "
                  f"S={np.round(rec['S'], 5).tolist()} xi={rec['xi']:.3f} dens={np.round(rec['dens'], 5).tolist()} "
                  f"t={rec['time_s']:.0f}s", flush=True)
        V.run([W] * a.n, a.n, a.chis, d=d, tol=a.tol, maxit=a.maxit, mode=a.vmode, log=lambda r_: None,
              on_chi=save, seed=a.seed)
    elif a.mode == 'obc':
        import fdmrg as F
        W = V.mpo_from_term(h, r, d=d)
        dm = F.FDMRG(W, a.L, d=d, seed=1, chi0=8)
        # open-chain boundary: the bond's one-site parts o_L (on its left site) and o_R (on its right
        # site) reach every bulk site once each, but site 0 misses o_R and site L-1 misses o_L.
        # Add them, so the open chain has no boundary potential.
        h2b = h if r == 2 else None
        if h2b is not None:
            hb = h2b.reshape(d, d, d, d)
            oL = np.einsum('aibi->ab', hb) / d
            oR = np.einsum('iaib->ab', hb) / d
            oL = oL - np.trace(oL) / d * np.eye(d); oR = oR - np.trace(oR) / d * np.eye(d)
            D0 = W.shape[0]
            W0 = W.copy(); W0[0, D0 - 1] = W0[0, D0 - 1] + oR
            WL = W.copy(); WL[0, D0 - 1] = WL[0, D0 - 1] + oL
            dm.Ws[0] = W0; dm.Ws[a.L - 1] = WL
            for i in range(a.L - 1, 0, -1):
                dm.Re[i] = F.rt(dm.Re[i + 1], dm.M[i], dm.Ws[i])
        for chi in a.chis:
            Eprev = None
            for sw in range(a.sweeps):
                E, tr = dm.sweep(chi, tol=1e-10)
                if Eprev is not None and abs(E - Eprev) < 1e-10 * a.L:
                    break
                Eprev = E
            rec = dict(chi=chi, E=E, trunc=tr, S=dm.S.tolist(), sweeps=sw + 1, time_s=time.time() - t0,
                       dens=dm.local_ops(nops).tolist())
            data['results'].append(rec)
            json.dump(data, open(a.out + '.tmp', 'w')); os.replace(a.out + '.tmp', a.out)
            print(f"{a.model} L={a.L} chi={chi} E={E:.12f} trunc={tr:.1e} Smid={dm.S[a.L // 2 - 1]:.5f} t={rec['time_s']:.0f}s", flush=True)
    else:
        pdiag = np.real(np.diag(CE.fermion_ops(nf)[1])) if a.model != 'xx2' else np.array([1., -1., -1., 1.])
        for L in a.Ls:
            res = ring_ed(h, d, L, a.nev, pdiag)
            data['results'].append(dict(L=L, sectors=res, time_s=time.time() - t0))
            json.dump(data, open(a.out + '.tmp', 'w')); os.replace(a.out + '.tmp', a.out)
            print(f"{a.model} L={L} " + ' '.join(f"Q={s['Q']}: {np.round(s['E'][:4], 6).tolist()}" for s in res), flush=True)


if __name__ == '__main__':
    main()
