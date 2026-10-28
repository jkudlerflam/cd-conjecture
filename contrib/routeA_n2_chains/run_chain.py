#!/usr/bin/env python3
"""Numerics on the n = 2 logical chains (d = 4 per cell) and their calibration chains.

CHAIN = CELL:x:LEVEL[:nu=..|:mu=..]
  CELL   c28 | c29 | spatial | compile  (bonds.py), or xx (one flavour, c = 1)
  x      coupling (g or tau)
  LEVEL  full   the complete order-lambda^2 bond (all generated terms)
         free   its quadratic part, normal ordered with respect to the empty state
         lead   the source's leading-order free chain (c28 leading_bond, c29 'leading')
         hf     the self-consistent Hartree-Fock quadratic chain of the full bond (exact c = 2 when it
                has four Fermi points), at the same chemical potential as the full chain
  nu=    retune the uniform field so that the Hartree-Fock (full, hf) or bare (free) density is nu
  mu=    add -mu (n_0 + n_1) per cell
Every bond is divided by s = max |T_ab| of its own quadratic part, so the hopping is O(1).

MODE
  fes    VUMPS (damped, one-cell unit cell) for a list of chi: S, xi, energy, densities; at the last chi
         the single-particle correlations <c^+_{0a} c_{rb}> and their Fourier transform n(k).
  fes2   the same with a two-cell unit cell at one chi (translation breaking).
  obc    two-site DMRG on an open chain of L cells: S(l), energy, densities, bond energies; for quadratic
         chains also the exact entropies from the correlation matrix.
usage: run_chain.py CHAIN fes --chis 16 24 32 ... --out data/X.json
       run_chain.py CHAIN obc --L 128 --chis 64 128 256 --out data/Y.json
Written by Claude (Opus 5.5) for contrib/routeA_n2_chains, 2026-10-28.
"""
import os, sys
for _k in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ.setdefault(_k, '4')
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import json, time, argparse
import numpy as np
from scipy import linalg
import bonds as B
import bands as BD
import vumps as V

NOPS = [B.jw(2)[a].conj().T @ B.jw(2)[a] for a in (0, 1)]
C2 = B.jw(2)                       # local annihilators on one cell
P2 = np.diag([1.0, -1.0, -1.0, 1.0]).astype(complex)   # cell parity


def quad_op(h1):
    """Quadratic 16 x 16 operator sum h1_ab c_a^+ c_b on the four bond modes."""
    return sum(h1[a, b] * B.CS[a].conj().T @ B.CS[b] for a in range(4) for b in range(4))


def build(spec):
    """Returns K (16x16 bond, normalized), meta."""
    parts = spec.split(':')
    cell, x, level = parts[0], float(parts[1]), parts[2]
    opts = dict(p.split('=') for p in parts[3:])
    meta = dict(spec=spec, cell=cell, x=x, level=level)
    if cell == 'xx':
        # one complex fermion with hopping -1, filling set by mu (default k_F = 0.962), flavour 1 empty
        kF = float(opts.get('kF', 0.962))
        h1 = np.zeros((4, 4), complex)
        h1[0, 2] = h1[2, 0] = -1.0
        mu0 = -2 * np.cos(kF)
        h1[0, 0] = h1[2, 2] = -0.5 * mu0
        h1[1, 1] = h1[3, 3] = 5.0
        K = B.traceless(quad_op(h1))
        meta.update(c_exact=1, kF=kF)
        return K, meta
    if level == 'lead':
        if cell == 'c28':
            z = np.load(os.path.join(HERE, 'src', 'c28', 'logical_chain_n2.npz'))
            K = B.traceless(z['leading_bond'])
        else:
            K, _ = B.c29(x, 'leading')
    else:
        K, info = B.CELLS[cell](x)
        meta.update({k: v for k, v in info.items() if k not in ('level', 'cell')})
    c, h, U, Hq, H4, r = B.normal_order(K)
    T, M = BD.blocks(h)
    s = float(np.abs(T).max())
    K = K / s
    meta['scale'] = s
    c, h, U, Hq, H4, r = B.normal_order(K)
    meta['nonquadratic_norm'] = float(np.linalg.norm(K - Hq, 2))
    meta['beyond_quartic_norm'] = r
    dmu = 0.0
    Fh = None
    if 'mu' in opts:
        dmu = float(opts['mu'])
    if level in ('full', 'hf'):
        nu = float(opts['nu']) if 'nu' in opts else None
        if nu is not None:
            Tf, Mf, dmu, it, Fh = BD.hartree_fock(K, nu=nu)
        Kmu = K + B.split_field(B.site_field([-dmu, -dmu]))
        Tf, Mf, _, it, Fh = BD.hartree_fock(Kmu, nu=None, mu=0.0)
        meta['hf_iters'] = it
        pts, _ = BD.fermi_points(Tf, Mf, 0.0)
        meta['hf_points'] = [(p['k'], p['band'], p['v']) for p in pts]
        meta['hf_density'] = BD.density(Tf, Mf, 0.0)
        if len(pts) >= 2:
            d2, procs = BD.margins(pts)
            meta['hf_centres'] = BD.centres(pts)
            meta['hf_delta2'] = d2
            meta['hf_mismatches'] = [(float(m), lab) for m, lab in procs[:8]]
    elif level in ('free', 'lead'):
        if 'nu' in opts:
            dmu = BD.mu_for_density(T / s, M / s, float(opts['nu']))
        K = B.traceless(Hq)
    meta['dmu'] = float(dmu)
    K = K + B.split_field(B.site_field([-dmu, -dmu]))
    if level == 'hf':
        K = B.traceless(quad_op(Fh))       # Fh already includes the -dmu field
    K = B.traceless(0.5 * (K + K.conj().T))
    c, h, U, Hq, H4, r = B.normal_order(K)
    T, M = BD.blocks(h)
    pts, _ = BD.fermi_points(T, M, 0.0)
    meta['bare_points'] = [(p['k'], p['band'], p['v']) for p in pts]
    meta['bare_density'] = BD.density(T, M, 0.0)
    meta['quadratic'] = bool(np.linalg.norm(K - Hq, 2) < 1e-12)
    return K, meta


# ----------------------------------------------------------------- uMPS measurements
def string_corr(psi, rmax):
    """rho_ab(r) = <c^+_{0,a} c_{r,b}> on a one-site uMPS, r = 0..rmax (JW strings included)."""
    AL, AR = psi.AL[0], psi.AR[0]
    C = psi.C[0]
    AC = np.tensordot(AL, C, ([2], [0]))
    nrm = np.vdot(AC, AC).real
    out = np.zeros((rmax + 1, 2, 2), complex)
    for a in range(2):
        Od = C2[a].conj().T
        for b in range(2):
            out[0, a, b] = np.vdot(AC, np.tensordot(Od @ C2[b], AC, ([1], [1])).transpose(1, 0, 2)) / nrm
        O1 = Od @ P2
        t = np.tensordot(O1, AC, ([1], [1])).transpose(1, 0, 2)          # a, s', b
        Y = np.tensordot(AC.conj(), t, ([0, 1], [0, 1])) / nrm           # (b', b)
        for r in range(1, rmax + 1):
            tb = np.tensordot(Y, AR, ([1], [0]))                          # b', s, c
            for b in range(2):
                t2 = np.tensordot(tb, C2[b], ([1], [1]))                  # b', c, s'
                out[r, a, b] = np.tensordot(AR.conj(), t2, ([0, 1, 2], [0, 2, 1]))
            t2 = np.tensordot(tb, P2, ([1], [1]))
            Y = np.tensordot(AR.conj(), t2, ([0, 1], [0, 2]))
    return out


def nk_from_corr(rho, nk=2048, window=None):
    """n(k) = sum_r rho(r) e^{-ikr} + h.c. part; eigenvalues of the 2x2 matrix n_ab(k), and the
    points where the total occupation jumps (largest |dn/dk|), with a Hann window of length rmax."""
    rmax = rho.shape[0] - 1
    w = np.hanning(2 * rmax + 3)[rmax + 1:-1] if window is None else window
    ks = np.linspace(-np.pi, np.pi, nk, endpoint=False)
    N = np.zeros((nk, 2, 2), complex)
    for i, k in enumerate(ks):
        M = rho[0].copy()
        for r in range(1, rmax + 1):
            M = M + w[r - 1] * (rho[r] * np.exp(-1j * k * r) + rho[r].conj().T * np.exp(1j * k * r))
        N[i] = M.T                                     # n_ab(k) = <c_ka^+ c_kb> (transpose convention)
    tot = np.real(np.einsum('kaa->k', N))
    der = np.gradient(tot, ks)
    # local extrema of |dn/dk| above a threshold
    a = np.abs(der)
    peaks = [i for i in range(nk) if a[i] >= a[i - 1] and a[i] >= a[(i + 1) % nk] and a[i] > 0.25 * a.max()]
    return ks, tot, [(float(ks[i]), float(np.sign(-der[i]))) for i in peaks]


def bond_energy_cell(psi, K, k):
    """<K> on the bond (k, k+1) of a uMPS with an n-site cell."""
    n = psi.n
    A1 = psi.AL[k]
    A2 = np.tensordot(psi.AL[(k + 1) % n], psi.C[(k + 2) % n], ([2], [0]))
    th = np.tensordot(A1, A2, ([2], [0]))                       # a s1 s2 b (left canonical left, C on right)
    d = A1.shape[1]
    H = K.reshape(d, d, d, d)
    Ht = np.einsum('stuv,auvb->astb', H, th)
    return float(np.vdot(th, Ht).real / np.vdot(th, th).real)


def run_fes(K, meta, a):
    W = V.mpo_from_term(K, 2, d=4)
    meta['mpo_dim'] = W.shape[0]
    data = dict(meta=meta, mode=a.mode, results=[])
    t0 = time.time()
    n = 2 if a.mode == 'fes2' else 1

    def on_chi(psi, rec):
        rec['time_s'] = time.time() - t0
        rec['dens'] = V.local_expect(psi, NOPS).tolist()
        rec['ebond'] = [bond_energy_cell(psi, K, k) for k in range(psi.n)]
        if n == 1 and rec['chi'] == a.chis[-1]:
            xi = rec['xi']
            rmax = int(min(max(8 * xi, 100), a.rmax))
            rho = string_corr(psi, rmax)
            ks, tot, jumps = nk_from_corr(rho)
            rec['corr_rmax'] = rmax
            rec['rho_r'] = [[float(z.real), float(z.imag)] for z in rho[:, :, :].reshape(-1)]
            rec['nk_jumps'] = jumps
            rec['nk_total'] = tot[::8].tolist()
        data['results'].append(rec)
        json.dump(data, open(a.out + '.tmp', 'w'))
        os.replace(a.out + '.tmp', a.out)
        print(f"[{meta['spec']}] chi={rec['chi']} e={rec['e']:.12f} S={np.round(rec['S'], 5).tolist()} xi={rec['xi']:.2f} "
              f"err={rec['err']:.1e} it={rec['iters']} dens={np.round(np.array(rec['dens']).sum(0), 5).tolist()} "
              f"t={rec['time_s']:.0f}s", flush=True)

    # chi by chi, with a restart from a fresh random state when the converged uMPS is not injective
    # (|lambda_1| -> 1, a cat of two states) or the gradient norm stays large.
    psi = None
    xi_prev = 1.0
    for chi in a.chis:
        good = None
        for attempt in range(3):
            if psi is None or attempt > 0:
                start = V.UMPS.random(n, 4, chi, a.seed + 17 * attempt + chi)
            else:
                cache = {}
                _, _, _, envs = V.vumps_step(psi, [W] * n, 1e-6, 1e-8, cache, tau=1e-9, krylov=4)
                start = V.expand(psi, [W] * n, envs, chi)
            cand, res = V.run([W] * n, n, [chi], d=4, tol=a.tol, maxit=a.maxit, seed=a.seed, log=lambda r: None,
                              psi=start, mode='damped', stol=a.stol)
            rec = res[0]
            rec['attempt'] = attempt
            ok = np.isfinite(rec['xi']) and rec['xi'] < 50 * max(xi_prev, 10.0) and rec['err'] < 1e-3
            if good is None or (ok and not good[2]) or (ok == good[2] and rec['e'] < good[1]['e']):
                good = (cand, rec, ok)
            if ok:
                break
            print(f"[{meta['spec']}] chi={chi} attempt {attempt} rejected: xi={rec['xi']:.3g} err={rec['err']:.1e}", flush=True)
        psi, rec, ok = good
        rec['accepted'] = bool(ok)
        xi_prev = rec['xi'] if ok else xi_prev
        on_chi(psi, rec)
    return data


def exact_open(K, L):
    """Exact entropies of the open chain sum_{j<L-1} K(j,j+1) for quadratic K."""
    c, h, U, Hq, H4, r = B.normal_order(K)
    n = 2 * L
    H = np.zeros((n, n), complex)
    for j in range(L - 1):
        H[2 * j:2 * j + 4, 2 * j:2 * j + 4] += h
    w, Uv = np.linalg.eigh(H)
    occ = Uv[:, w < 0]
    Cm = occ.conj() @ occ.T
    S = []
    for l in range(1, L):
        e = np.linalg.eigvalsh(Cm[:2 * l, :2 * l]).real
        e = e[(e > 1e-14) & (e < 1 - 1e-14)]
        S.append(float(-(e * np.log(e) + (1 - e) * np.log(1 - e)).sum()))
    E0 = float(w[w < 0].sum()) + (L - 1) * float(c)
    return dict(S=S, E=E0, N=int((w < 0).sum()), gap_sp=float(np.min(np.abs(w))))


def run_obc(K, meta, a):
    from fdmrg import FDMRG
    W = V.mpo_from_term(K, 2, d=4)
    meta['mpo_dim'] = W.shape[0]
    data = dict(meta=meta, mode='obc', L=a.L, results=[])
    if meta.get('quadratic'):
        data['exact'] = exact_open(K, a.L)
    dm = FDMRG(W, a.L, d=4, seed=a.seed)
    t0 = time.time()
    for chi in a.chis:
        Eprev = None
        for sw in range(a.sweeps * 3):
            E, tr = dm.sweep(chi)
            print(f"[{meta['spec']}] L={a.L} chi={chi} sweep={sw} E={E:.12f} trunc={tr:.1e} Smid={dm.S[a.L // 2 - 1]:.5f} "
                  f"t={time.time() - t0:.0f}s", flush=True)
            if sw + 1 >= a.sweeps and Eprev is not None and abs(E - Eprev) < 1e-10 * a.L:
                break
            Eprev = E
        dens = dm.local_ops(NOPS)
        X = dm.left_norm_envs()
        eb = []
        for j in range(a.L - 1):
            th = np.tensordot(dm.M[j], dm.M[j + 1], ([2], [0]))
            t = np.tensordot(X[j], th, ([1], [0]))
            Ht = np.einsum('stuv,auvb->astb', K.reshape(4, 4, 4, 4), t)
            eb.append(float(np.tensordot(th.conj(), Ht, ([0, 1, 2, 3], [0, 1, 2, 3])).real))
        rec = dict(chi=chi, E=E, sweeps=sw + 1, trunc=tr, S=dm.S.tolist(), dens=dens.tolist(), ebond=eb,
                   N=float(dens.sum()), time_s=time.time() - t0)
        data['results'].append(rec)
        json.dump(data, open(a.out + '.tmp', 'w'))
        os.replace(a.out + '.tmp', a.out)
    return data


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('chain')
    ap.add_argument('mode', choices=['fes', 'fes2', 'obc', 'meta'])
    ap.add_argument('--chis', type=int, nargs='+', default=[16, 24, 32, 48, 64, 96, 128])
    ap.add_argument('--L', type=int, default=128)
    ap.add_argument('--sweeps', type=int, default=6)
    ap.add_argument('--tol', type=float, default=1e-8)
    ap.add_argument('--stol', type=float, default=2e-5)
    ap.add_argument('--maxit', type=int, default=250)
    ap.add_argument('--rmax', type=int, default=3000)
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    K, meta = build(a.chain)
    if a.mode == 'meta':
        print(json.dumps(meta, default=str, indent=1)); return
    a.out = a.out or os.path.join(HERE, 'data', a.chain.replace(':', '_') + f'_{a.mode}.json')
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    if a.mode in ('fes', 'fes2'):
        run_fes(K, meta, a)
    else:
        run_obc(K, meta, a)


if __name__ == '__main__':
    main()
