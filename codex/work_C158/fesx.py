#!/usr/bin/env python3
"""C158 (Claude Opus agent, 2026-11-07): finite-entanglement scaling with the damped VUMPS of vumps.py
(copied unchanged from codex/work_C153) for qubit chains near the XZ point.

Bond convention (C143/C144): h_{j,j+1} = sum_ab J[a,b] s^a_j s^b_{j+1} + sum_a B[a] s^a_j (one copy of the
field per site), Pauli matrices.  A cell of n sites may carry site-dependent bonds (Js[k] for bond (k,k+1))
and fields (Bs[k] on site k).

Models (--model):
  tfi            -XX - Z                              (critical Ising, c = 1/2)
  xx             XX + YY                              (c = 1)
  heis           XX + YY + ZZ                         (c = 1, marginal log corrections)
  xzzx           XZ + ZX                              (two decoupled critical Majoranas, c = 1)
  xzalt   --eps  XZ + eps (-1)^j ZX                   (sign alternation is a gauge: B critical, c = 1)
  xzdim   --eps  XZ + eps ZX on even bonds only        (A critical, B in a commuting stabilizer state: c = 1/2 exactly)
  xzeps   --eps  XZ + eps ZX                          (two decoupled Majoranas, velocities 1 : eps, c = 1)
  bis     --hp   XZ + B (x + z)/sqrt2, B = -2 hp      (Kitaev chain in a [110] field, hp in units of K)
  dir     --dirs 'yx=0.3,zy=0.1,...'                  (XZ + sum_t c_t * term t; terms xx..zz, Bx, By, Bz)
  gen     --J 9 numbers --B 3 numbers                 (general bond)
usage: fesx.py MODEL [params] --n 2 --chis 8 12 16 ... --out file.json [--npz state.npz]
"""
import os, sys, json, time, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import vumps as V

P4 = [np.eye(2, dtype=complex)] + [V.PAULI[a] for a in range(3)]
NM = 'xyz'


def term_JB(name):
    J = np.zeros((3, 3)); B = np.zeros(3)
    if name[0] == 'B':
        B[NM.index(name[1])] = 1.0
    else:
        J[NM.index(name[0]), NM.index(name[1])] = 1.0
    return J, B


def mpo_cell(Js, Bs):
    """Upper-triangular MPO (D = 5) for an n-site cell. W[k] at site k opens bond (k,k+1) and closes bond (k-1,k)."""
    n = len(Js)
    D = 5
    Ws = []
    for k in range(n):
        W = np.zeros((D, D, 2, 2), complex)
        W[0, 0] = P4[0]; W[D - 1, D - 1] = P4[0]
        Jprev = np.asarray(Js[(k - 1) % n], float)
        for a in range(3):
            W[0, 1 + a] = P4[1 + a]
            W[1 + a, D - 1] = sum(Jprev[a, b] * P4[1 + b] for b in range(3))
        W[0, D - 1] = sum(Bs[k][a] * P4[1 + a] for a in range(3))
        Ws.append(W)
    return Ws


HAD = np.array([[1, 1], [1, -1]], complex) / np.sqrt(2)


def sym_radius(psi, U=HAD, shift=0):
    """Spectral radius per cell of the mixed transfer matrix between |psi> and (prod_j U_j) R |psi>, R the spatial
    reflection (cell order reversed, optionally rolled by `shift` sites).  Equal to 1 iff the uMPS is invariant
    (up to a phase) under this symmetry; < 1 means the state breaks it (or is a different translate)."""
    from scipy.sparse.linalg import LinearOperator, eigs
    n = psi.n
    A = psi.AL
    Bt = [np.einsum('st,atb->bsa', U, A[n - 1 - k]) for k in range(n)]
    Bt = Bt[shift:] + Bt[:shift]
    ca = A[0].shape[0]; cb = Bt[0].shape[0]

    def mv(v):
        X = v.reshape(cb, ca)
        for k in range(n):
            t = np.tensordot(X, A[k], ([1], [0]))                 # b, s, a'
            X = np.tensordot(Bt[k].conj(), t, ([0, 1], [0, 1]))   # b', a'
        return X.ravel()
    N = ca * cb
    if N <= 300:
        M = np.array([mv(e) for e in np.eye(N, dtype=complex)]).T
        return float(np.max(np.abs(np.linalg.eigvals(M))))
    w = eigs(LinearOperator((N, N), matvec=mv, dtype=complex), k=1, which='LM', tol=1e-10, maxiter=3000,
             return_eigenvectors=False)
    return float(np.abs(w[0]))


def model(a):
    n = a.n
    Z3 = np.zeros((3, 3))
    XZ = Z3.copy(); XZ[0, 2] = 1
    ZX = Z3.copy(); ZX[2, 0] = 1
    Bs = [np.zeros(3) for _ in range(n)]
    if a.model == 'tfi':
        J = Z3.copy(); J[0, 0] = -1; Js = [J] * n; Bs = [np.array([0, 0, -1.])] * n
        exact = -4 / np.pi
    elif a.model == 'xx':
        Js = [np.diag([1., 1, 0])] * n; exact = -4 / np.pi
    elif a.model == 'heis':
        Js = [np.diag([1., 1, 1])] * n; exact = 1 - 4 * np.log(2)
    elif a.model == 'xzzx':
        Js = [XZ + ZX] * n; exact = -4 / np.pi
    elif a.model == 'xzalt':
        assert n % 2 == 0
        Js = [XZ + a.eps * (-1) ** k * ZX for k in range(n)]
        exact = -(2 / np.pi) * (1 + a.eps)
    elif a.model == 'xzdim':
        assert n % 2 == 0
        Js = [XZ + a.eps * (k % 2 == 0) * ZX for k in range(n)]
        exact = -(2 / np.pi) - a.eps / 2
    elif a.model == 'xzeps':
        Js = [XZ + a.eps * ZX] * n; exact = -(2 / np.pi) * (1 + a.eps)
    elif a.model == 'bis':
        Js = [XZ] * n; Bs = [-2 * a.hp * np.array([1, 0, 1]) / np.sqrt(2)] * n; exact = None
    elif a.model == 'dir':
        J = XZ.copy(); B = np.zeros(3)
        for it in a.dirs.split(','):
            nm, c = it.split('=')
            dJ, dB = term_JB(nm); J = J + float(c) * dJ; B = B + float(c) * dB
        Js = [J] * n; Bs = [B] * n; exact = None
    elif a.model == 'gen':
        Js = [np.array(a.J, float).reshape(3, 3)] * n; Bs = [np.array(a.B, float)] * n; exact = None
    else:
        raise SystemExit('unknown model')
    return Js, Bs, exact


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('model')
    ap.add_argument('--eps', type=float, default=0.0)
    ap.add_argument('--hp', type=float, default=0.0)
    ap.add_argument('--dirs', default='')
    ap.add_argument('--J', type=float, nargs=9)
    ap.add_argument('--B', type=float, nargs=3, default=[0, 0, 0])
    ap.add_argument('--n', type=int, default=2)
    ap.add_argument('--chis', type=int, nargs='+', default=[8, 12, 16, 24, 32])
    ap.add_argument('--maxit', type=int, default=600)
    ap.add_argument('--tol', type=float, default=1e-9)
    ap.add_argument('--mode', default='damped')
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--nev', type=int, default=12)
    ap.add_argument('--out', required=True)
    ap.add_argument('--npz', default=None)
    a = ap.parse_args()
    Js, Bs, exact = model(a)
    W = mpo_cell(Js, Bs)
    data = dict(args=vars(a), Js=[np.asarray(J).tolist() for J in Js], Bs=[np.asarray(B).tolist() for B in Bs],
                exact_e=exact, results=[])
    t0 = time.time()
    psi0 = None
    if a.npz and os.path.exists(a.npz):
        z = np.load(a.npz)
        nn = int(z['n'])
        psi0 = V.UMPS([z[f'AL{k}'] for k in range(nn)], [z[f'AR{k}'] for k in range(nn)], [z[f'C{k}'] for k in range(nn)])
        print('restarting from', a.npz, 'chi', psi0.AL[0].shape[0], flush=True)
        if os.path.exists(a.out):
            data['results'] = json.load(open(a.out))['results']
        done = {r['chi'] for r in data['results']}
        a.chis = [c for c in a.chis if c not in done and c > psi0.AL[0].shape[0]]
        if not a.chis:
            print('nothing to do', flush=True); return
        _, _, _, envs = V.vumps_step(psi0, W, 1e-6, 1e-8, {}, tau=1e-9, krylov=4)
        psi0 = V.expand(psi0, W, envs, a.chis[0])

    def on_chi(psi, rec):
        rec['time_s'] = time.time() - t0
        lam = V.transfer_spectrum(psi, nev=a.nev)
        rec['lam'] = [[float(z.real), float(z.imag)] for z in lam]
        rec['schmidt0'] = [float(s) for s in V.schmidt(psi, 0)[:40]]
        loc = V.local_expect(psi, [V.PAULI[0], V.PAULI[1], V.PAULI[2]])
        rec['mag'] = loc.tolist()
        try:
            rec['sym_HR'] = [sym_radius(psi, HAD, sh) for sh in range(min(psi.n, 2))]
        except Exception as ex:
            rec['sym_HR'] = [str(ex)]
        data['results'].append(rec)
        json.dump(data, open(a.out + '.tmp', 'w'), default=float); os.replace(a.out + '.tmp', a.out)
        print(f"chi={rec['chi']} e={rec['e']:.12f} err={rec['err']:.1e} it={rec['iters']} S={np.round(rec['S'],5).tolist()} "
              f"xi={rec['xi']:.2f} symHR={np.round(rec['sym_HR'],6).tolist() if not isinstance(rec['sym_HR'][0], str) else rec['sym_HR']} arg/pi={np.angle(complex(*rec['lam'][1]))/np.pi:.4f} t={rec['time_s']:.0f}s", flush=True)
        if a.npz:
            dd = dict(n=psi.n)
            for k in range(psi.n):
                dd[f'AL{k}'] = psi.AL[k]; dd[f'AR{k}'] = psi.AR[k]; dd[f'C{k}'] = psi.C[k]
            np.savez(a.npz + '.tmp.npz', **dd); os.replace(a.npz + '.tmp.npz', a.npz)

    V.run(W, a.n, a.chis, d=2, tol=a.tol, maxit=a.maxit, seed=a.seed, mode=a.mode, on_chi=on_chi, psi=psi0,
          log=lambda r: None)
    print('done', a.out, f'{time.time() - t0:.0f}s', flush=True)


if __name__ == '__main__':
    main()
