#!/usr/bin/env python3
"""C166 task 2 (Claude Opus 5.5 agent, 2026-11-08): damped VUMPS (vumps.py, unchanged) in an n-site cell for a
named bond (c166lib), growing chi through a list, with a checkpoint so that a job can resume.
Output JSON has the fesx.py layout (args.n, results[...] with chi, e, S, xi, lam, err, iters), so
fesfit_c158.analyse reads it, plus per chi:
  mag[a][k]      <sigma^a> on site k of the cell
  chir[a][k]     <(sigma_k x sigma_{k+1})^a>  (vector chirality, Theta-even, parity-odd)
  bond_e[k]      <h_{k,k+1}> including f . sigma_k (one copy of the field per site)
  tr_radius      spectral radius per cell of the mixed transfer matrix <psi|Theta psi>, Theta = i sigma^y K on
                 every site; 1 iff the uMPS is time-reversal invariant (up to phase).
usage: vumps_cells.py NAME --n N --chis 8 16 ... --seed S --out f.json [--npz f.npz] [--maxit 600]"""
import os, sys, json, time, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import vumps as V
import fesx as F
import c166lib as X

SY = V.PAULI[1]


def tr_radius(psi):
    from scipy.sparse.linalg import LinearOperator, eigs
    n = psi.n
    A = psi.AL
    U = 1j * SY
    B = [np.einsum('st,atb->asb', U, A[k].conj()) for k in range(n)]
    ca = A[0].shape[0]

    def mv(v):
        X_ = v.reshape(ca, ca)
        for k in range(n):
            t = np.tensordot(X_, A[k], ([1], [0]))                # b, s, a'
            X_ = np.tensordot(B[k].conj(), t, ([0, 1], [0, 1]))   # b', a'
        return X_.ravel()
    N = ca * ca
    if N <= 300:
        M = np.array([mv(e) for e in np.eye(N, dtype=complex)]).T
        return float(np.max(np.abs(np.linalg.eigvals(M))))
    w = eigs(LinearOperator((N, N), matvec=mv, dtype=complex), k=1, which='LM', tol=1e-10, maxiter=3000,
             return_eigenvectors=False)
    return float(np.abs(w[0]))


def two_site(psi, k, O):
    """<O> for a 4x4 operator O on sites (k, k+1) of the cell."""
    n = psi.n
    k1 = (k + 1) % n
    A2 = np.tensordot(psi.AL[k], psi.AC[k1], ([2], [0]))     # a s t b
    O4 = O.reshape(2, 2, 2, 2)
    t = np.einsum('uvst,astb->auvb', O4, A2)
    return complex(np.vdot(A2, t) / np.vdot(A2, A2))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('name')
    ap.add_argument('--n', type=int, default=2)
    ap.add_argument('--chis', type=int, nargs='+', default=[8, 16, 32])
    ap.add_argument('--maxit', type=int, default=600)
    ap.add_argument('--tol', type=float, default=1e-9)
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--nev', type=int, default=12)
    ap.add_argument('--out', required=True)
    ap.add_argument('--npz', default=None)
    a = ap.parse_args()
    J, f = X.JB(a.name)
    Js = [J] * a.n; Bs = [f] * a.n
    W = F.mpo_cell(Js, Bs)
    hb = sum(J[p, q] * np.kron(V.PAULI[p], V.PAULI[q]) for p in range(3) for q in range(3)) \
        + sum(f[p] * np.kron(V.PAULI[p], np.eye(2)) for p in range(3))
    chirops = []
    for c in range(3):
        p, q = (c + 1) % 3, (c + 2) % 3
        chirops.append(np.kron(V.PAULI[p], V.PAULI[q]) - np.kron(V.PAULI[q], V.PAULI[p]))
    data = dict(args=vars(a), Js=[np.asarray(J).tolist()] * a.n, Bs=[np.asarray(f).tolist()] * a.n, exact_e=None,
                results=[], v=X.load_vec(a.name).tolist())
    t0 = time.time()
    psi0 = None
    chis = list(a.chis)
    if a.npz and os.path.exists(a.npz):
        z = np.load(a.npz)
        nn = int(z['n'])
        psi0 = V.UMPS([z[f'AL{k}'] for k in range(nn)], [z[f'AR{k}'] for k in range(nn)], [z[f'C{k}'] for k in range(nn)])
        if os.path.exists(a.out):
            data['results'] = json.load(open(a.out))['results']
        done = {r['chi'] for r in data['results']}
        chis = [c for c in chis if c not in done and c > psi0.AL[0].shape[0]]
        if not chis:
            print('nothing to do', flush=True); return
        print('restart from', a.npz, 'chi', psi0.AL[0].shape[0], '->', chis, flush=True)
        _, _, _, envs = V.vumps_step(psi0, W, 1e-6, 1e-8, {}, tau=1e-9, krylov=4)
        psi0 = V.expand(psi0, W, envs, chis[0])

    def on_chi(psi, rec):
        rec['time_s'] = time.time() - t0
        lam = V.transfer_spectrum(psi, nev=a.nev)
        rec['lam'] = [[float(z.real), float(z.imag)] for z in lam]
        rec['schmidt0'] = [float(s) for s in V.schmidt(psi, 0)[:40]]
        rec['mag'] = V.local_expect(psi, [V.PAULI[0], V.PAULI[1], V.PAULI[2]]).tolist()
        rec['chir'] = [[two_site(psi, k, O).real for k in range(psi.n)] for O in chirops]
        rec['bond_e'] = [two_site(psi, k, hb).real for k in range(psi.n)]
        try:
            rec['tr_radius'] = tr_radius(psi)
        except Exception as ex:
            rec['tr_radius'] = str(ex)
        data['results'].append(rec)
        json.dump(data, open(a.out + '.tmp', 'w'), default=float); os.replace(a.out + '.tmp', a.out)
        m = np.array(rec['mag'])
        print(f"{a.name} n={a.n} seed={a.seed} chi={rec['chi']} e={rec['e']:.12f} err={rec['err']:.1e} it={rec['iters']} "
              f"S={np.round(rec['S'], 4).tolist()} xi={rec['xi']:.2f} |m|max={np.max(np.abs(m)):.4f} "
              f"TR={rec['tr_radius'] if isinstance(rec['tr_radius'], str) else round(rec['tr_radius'], 6)} "
              f"arg1/pi={np.angle(complex(*rec['lam'][1])) / np.pi:.4f} t={rec['time_s']:.0f}s", flush=True)
        if a.npz:
            dd = dict(n=psi.n)
            for k in range(psi.n):
                dd[f'AL{k}'] = psi.AL[k]; dd[f'AR{k}'] = psi.AR[k]; dd[f'C{k}'] = psi.C[k]
            np.savez(a.npz + '.tmp.npz', **dd); os.replace(a.npz + '.tmp.npz', a.npz)

    V.run(W, a.n, chis, d=2, tol=a.tol, maxit=a.maxit, seed=a.seed, mode='damped', on_chi=on_chi, psi=psi0,
          log=lambda r: None)
    print('done', a.out, f'{time.time() - t0:.0f}s', flush=True)


if __name__ == '__main__':
    main()
