#!/usr/bin/env python3
"""Finite-entanglement scaling (VUMPS, ../d2_followup_220_1282/vumps.py) for the d = 11 vacancy chain
and its calibrations.

usage: run_fes.py NAME --model {vac,xx,twoxx} [--t 1 --J 0.25 --mu_v 0.6 --lam 0.15] [--r 0.01]
                  --n 1 --chis 16 24 32 ... [--out data/fes/NAME.json]
  vac   : the full d = 11 chain of vacbond.py (random initial state; nothing imposed)
  xx    : XX chain, c = 1
  twoxx : two decoupled XX chains with velocities 1 and r on a d = 4 site (c = 2, two velocities)
Records per chi: energy per site, gradient error, entanglement, the 8 leading transfer-matrix eigenvalues
(xi = -n / ln|l1/l0|), occupied density, B density, adjacency <n_i n_i+1>, and logical charge per site.
Claude (Opus 5.5), 2026-10-28.
"""
import os, sys, json, time, argparse
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'd2_followup_220_1282'))
sys.path.insert(0, HERE)
import numpy as np
import vumps as V
import vacbond as VB


def model(a):
    if a.model == 'vac':
        h = VB.vacancy_lift(VB.demo_bond(a.lam), t=a.t, J=a.J, mu_v=a.mu_v)
        d = 11
        n = np.diag(np.r_[0, np.ones(10)])
        nB = np.diag(np.r_[0, 0, 0, np.ones(8)])
        q = np.diag(np.r_[0, VB.Q2_LETTER / 2.0])
        return h, d, [n, nB, q], np.kron(n, n)
    sp = np.array([[0., 1.], [0., 0.]]); sm = sp.T
    hxx = -(np.kron(sp, sm) + np.kron(sm, sp))
    if a.model == 'xx':
        return hxx, 2, [np.diag([1., -1.])], None
    if a.model == 'twoxx':
        I2 = np.eye(2)
        # site = (a, b); bond on (a1 b1)(a2 b2): XX on a1 a2 plus r XX on b1 b2
        A = np.kron(np.kron(sp, I2), np.kron(sm, I2)) + np.kron(np.kron(sm, I2), np.kron(sp, I2))
        B = np.kron(np.kron(I2, sp), np.kron(I2, sm)) + np.kron(np.kron(I2, sm), np.kron(I2, sp))
        h = -(A + a.r * B)
        z = np.diag([1., -1.])
        return h, 4, [np.kron(z, I2), np.kron(I2, z)], None
    if a.model == 'lv':
        # logical vacancy model LV_n of logical_vac.py (n = len(taus) flavours)
        import logical_vac as LV
        nf = len(a.taus)
        h = LV.model_bond(nf, a.t, a.J, a.mu_v, a.taus, a.eps)
        d = 1 + 2 ** nf
        n = np.diag(np.r_[0, np.ones(d - 1)])
        fl = [np.diag(np.r_[0, [(b >> (nf - 1 - f)) & 1 for b in range(d - 1)]]).astype(float) for f in range(nf)]
        return h, d, [n] + fl, np.kron(n, n)
    if a.model == 'kxx':
        # k decoupled XX chains with velocities a.vels on a d = 2^k site (c = k)
        k = len(a.vels)
        I = np.eye(2)

        def emb(o, i):
            m = np.ones((1, 1))
            for j in range(k):
                m = np.kron(m, o if j == i else I)
            return m
        h = 0
        for i, v in enumerate(a.vels):
            A = np.kron(emb(sp, i), emb(sm, i))
            h = h - v * (A + A.T)
        z = np.diag([1., -1.])
        return h, 2 ** k, [emb(z, i) for i in range(k)], None
    raise ValueError(a.model)


def canonical_umps(A, tol=1e-12, maxit=5000):
    """Mixed canonical form (A_L, A_R, C) of the one-site uMPS A (chi, d, chi)."""
    chi, d, _ = A.shape
    A = A / np.sqrt(np.abs(V.linalg.eigvals(np.einsum('asb,csd->acbd', A, A.conj()).reshape(chi * chi, chi * chi))).max())
    Lm = np.eye(chi)
    for _ in range(maxit):
        AL, R = np.linalg.qr(np.tensordot(Lm, A, ([1], [0])).reshape(chi * d, chi))
        ph = np.sign(np.diag(R)); ph[ph == 0] = 1
        AL = AL * ph[None, :]; R = ph[:, None] * R
        R = R / np.linalg.norm(R)
        if np.linalg.norm(R - Lm) < tol:
            Lm = R; break
        Lm = R
    AL = AL.reshape(chi, d, chi)
    Rm = np.eye(chi)
    for _ in range(maxit):
        Q, Rr = np.linalg.qr(np.tensordot(A, Rm, ([2], [0])).reshape(chi, d * chi).T)
        ph = np.sign(np.diag(Rr)); ph[ph == 0] = 1
        Q = Q * ph[None, :]; Rr = ph[:, None] * Rr
        Lq = Rr.T / np.linalg.norm(Rr)
        if np.linalg.norm(Lq - Rm) < tol:
            Rm = Lq; break
        Rm = Lq
    AR = Q.T.reshape(chi, d, chi)
    C = Lm @ Rm
    U, S, Vh = np.linalg.svd(C)
    AL = np.einsum('ab,bsc,cd->asd', U.conj().T, AL, U)
    AR = np.einsum('ab,bsc,cd->asd', Vh, AR, Vh.conj().T)
    C = np.diag(S / np.linalg.norm(S))
    return V.UMPS([AL.astype(complex)], [AR.astype(complex)], [C.astype(complex)])


def pattern_umps(rho, aux=2, seed=0):
    """chi = 3 aux uMPS supported on the period-3 type pattern P B P of the occupied sites: the bond carries
    the pattern phase s (next occupied site has type pattern[s]); a vacancy keeps s, a letter advances it."""
    rng = np.random.default_rng(seed)
    pat = ['P', 'B', 'P']
    chi = 3 * aux
    A = np.zeros((chi, 11, chi))
    for s in range(3):
        sl = slice(s * aux, (s + 1) * aux)
        A[sl, 0, sl] = np.sqrt(1 - rho) * (np.eye(aux) + 0.1 * rng.normal(size=(aux, aux)))
        s1 = (s + 1) % 3
        sl1 = slice(s1 * aux, (s1 + 1) * aux)
        for a in VB.LETTERS[pat[s]]:
            A[sl, 1 + a, sl1] = np.sqrt(rho / len(VB.LETTERS[pat[s]])) * (np.eye(aux) + 0.3 * rng.normal(size=(aux, aux)))
    return canonical_umps(A)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('name')
    ap.add_argument('--model', default='vac')
    for k, v in dict(t=1.0, J=0.25, mu_v=0.6, lam=0.15, r=0.01).items():
        ap.add_argument('--' + k, type=float, default=v)
    ap.add_argument('--taus', type=float, nargs='+', default=[1.0], help='lv: flavour hoppings')
    ap.add_argument('--eps', type=float, nargs='+', default=[0.0], help='lv: flavour fields')
    ap.add_argument('--vels', type=float, nargs='+', default=[1.0, 0.5], help='kxx: chain velocities')
    ap.add_argument('--chis', type=int, nargs='+', default=[16, 24, 32])
    ap.add_argument('--n', type=int, default=1)
    ap.add_argument('--mode', default='damped')
    ap.add_argument('--maxit', type=int, default=400)
    ap.add_argument('--tol', type=float, default=1e-9)
    ap.add_argument('--stol', type=float, default=1e-8)
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--out', default=None)
    ap.add_argument('--state_dir', default=None)
    ap.add_argument('--init', default='random', help="random, or 'pattern' (vac only): start in the P B P type sector")
    ap.add_argument('--rho0', type=float, default=0.6)
    a = ap.parse_args()
    h, d, ops, nn = model(a)
    W = V.mpo_from_term(h, 2, d=d)
    out = a.out or os.path.join(HERE, 'data', 'fes', a.name + '.json')
    os.makedirs(os.path.dirname(out), exist_ok=True)
    data = dict(name=a.name, args=vars(a), results=[])
    t0 = time.time()

    def on_chi(psi, rec):
        rec['time_s'] = time.time() - t0
        rec['local'] = V.local_expect(psi, ops).tolist()
        if nn is not None:
            vals = []
            for k in range(psi.n):
                AC = np.tensordot(psi.AL[k], psi.C[(k + 1) % psi.n], ([2], [0]))
                th = np.tensordot(AC, psi.AR[(k + 1) % psi.n], ([2], [0])).reshape(AC.shape[0], d * d, -1)
                nrm = np.vdot(th, th).real
                vals.append(float(np.einsum('asb,s,asb->', th.conj(), np.diag(nn), th).real / nrm))
            rec['adjacent_occupied'] = vals
        try:
            lam = V.transfer_spectrum(psi, nev=10)
            rec['lam10'] = [[float(z.real), float(z.imag)] for z in lam]
        except Exception as e:  # noqa: BLE001
            rec['lam10_error'] = str(e)
        data['results'].append(rec)
        json.dump(data, open(out + '.tmp', 'w'))
        os.replace(out + '.tmp', out)
        print(f"{a.name} chi={rec['chi']} e={rec['e']:.12f} err={rec['err']:.1e} it={rec['iters']} "
              f"S={np.round(rec['S'], 5).tolist()} xi={rec['xi']:.3f} local={np.round(rec['local'], 5).tolist()} "
              f"t={rec['time_s']:.0f}s", flush=True)
        if a.state_dir:
            os.makedirs(a.state_dir, exist_ok=True)
            np.savez_compressed(os.path.join(a.state_dir, f"{a.name}_chi{rec['chi']}.npz"),
                                **{f'AL{k}': psi.AL[k] for k in range(psi.n)},
                                **{f'AR{k}': psi.AR[k] for k in range(psi.n)},
                                **{f'C{k}': psi.C[k] for k in range(psi.n)})

    psi0 = None
    if a.init == 'pattern':
        assert a.model == 'vac' and a.n == 1
        psi0 = pattern_umps(a.rho0, aux=max(1, a.chis[0] // 3), seed=a.seed)
        chis = [3 * max(1, a.chis[0] // 3)] + [c for c in a.chis if c > 3 * max(1, a.chis[0] // 3)]
    else:
        chis = a.chis
    V.run([W] * a.n, a.n, chis, d=d, tol=a.tol, maxit=a.maxit, seed=a.seed, mode=a.mode, psi=psi0,
          log=lambda rec: None, on_chi=on_chi, stol=a.stol)


if __name__ == '__main__':
    main()
