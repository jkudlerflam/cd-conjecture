#!/usr/bin/env python3
"""Finite-entanglement scaling with VUMPS for the encoded chain and its controls.

usage:
  run_fes.py ising   --ncells 2 --chis 8 16 32 64
  run_fes.py logical --q 2 --qp 1 --n 1 --ncells 4 --chis ...          (first-order chain, t = 1)
  run_fes.py encoded --q 2 --qp 1 --n 1 --g 0.5 --lam 0.1 --ncells 2 --chis ... [--mode damped]
  run_fes.py bigsite --q 2 ...                                         (literal D x D site)
Output: JSON with one record per chi (energy per cell, entropies on every bond, transfer-matrix
spectrum, xi in cells, bond/cell energies of each cell in the unit cell, leakage).

Written by Claude (Opus 5.5) for contrib/encoded_chain_numerics, 2026-10-27.
"""
import os, sys, json, time, argparse
for _k in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ.setdefault(_k, '4')
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np
import cells as C
import umps_md as U


def build(a):
    info = {}
    if a.model == 'ising':
        d, t = C.ising_terms()
        dims, terms = [], []
        for j in range(a.ncells):
            dims += d
            terms += [((s + j) % a.ncells, m, H) for (s, m, H) in t]
        W = C.mpo_from_terms(dims, terms)
        return dims, W, [1] * len(dims), 1, info, None
    c = C.cell(a.q, a.qp, a.n, 1.0, a.g)
    if a.model == 'logical':
        dims, terms = C.logical_chain_terms(c, 1.0, a.ncells)
        return dims, C.mpo_from_terms(dims, terms), [1] * len(dims), 1, info, c
    if a.model == 'core':
        import effective as EF
        so = EF.second_order(c)
        hb = np.load(os.path.join(HERE, 'data', 'core_bonds.npz'))[f'g{a.g}_l{a.lam}']
        dims = [c['N']] * a.ncells
        terms = [(j, 2, hb) for j in range(a.ncells)] + [(j, 3, a.lam ** 2 * so['T3']) for j in range(a.ncells)]
        info.update(alpha=float(so['alpha']), t1=float(a.lam * so['alpha'] ** 2))
        return dims, C.mpo_from_terms(dims, terms), [1] * len(dims), 1, info, c
    if a.model == 'effective':
        import effective as EF
        so = EF.second_order(c)
        so['T2'] = so['T2'] * a.kappa; so['T3'] = so['T3'] * a.kappa
        dims, terms = EF.eff_terms(c, a.lam, a.ncells, so)
        info.update(alpha=float(so['alpha']), t1=float(a.lam * so['alpha'] ** 2), kappa=a.kappa)
        return dims, C.mpo_from_terms(dims, terms), [1] * len(dims), 1, info, c
    o = C.cell_spectrum_and_alpha(c)
    info.update(alpha=float(o['alphaL'][0]), cell_gap=float(o['gap']), cell_E0=float(o['E'][0]),
                t1=float(a.lam * o['alphaL'][0] ** 2))
    if a.model == 'bigsite':
        dims, W = C.bigsite_mpo(c, a.lam, a.ncells)
        return dims, W, [1] * len(dims), 1, info, c
    if a.model == 'sparse':
        import sparse_mpo as SPM
        SPM.install()
        dims, W = SPM.bigsite_spmpo(c, a.lam, a.ncells, C.bigsite_H, C.bigsite_ports)
        return dims, W, [1] * len(dims), 1, info, c
    if a.model == 'trunc':
        import sparse_mpo as SPM
        SPM.install()
        Em, FLm, FRm, Qm, gcut = C.truncated_cell(c, a.M)
        info.update(M=a.M, gap_at_cut=float(gcut))
        c['trunc'] = (Em, FLm, FRm, Qm)
        dims, W = SPM.bigsite_spmpo(c, a.lam, a.ncells, lambda cc: np.diag(Em), lambda cc: (FLm, FRm))
        return dims, W, [1] * len(dims), 1, info, c
    dims, W = C.chain_mpo(c, a.lam, a.ncells)
    m = a.mult if a.mult > 0 else a.q          # q*chi on inner bonds = exact big-site MPS
    mults = [1, m, m] * a.ncells
    return dims, W, mults, 3, info, c


def bond_expect(psi, j, FR, FL):
    """sum_a <F^R_a(j) F^L_a(j+1)> for big sites, without forming the two-site tensor."""
    n = psi.n
    AC, AR = psi.AC[j], psi.AR[(j + 1) % n]
    nrm = np.real(np.vdot(AC, AC))
    tot = 0.0
    for fr, fl in zip(FR, FL):
        t = np.tensordot(fr, AC, ([1], [1])).transpose(1, 0, 2)          # a s' c
        X = np.tensordot(AC.conj(), t, ([0, 1], [0, 1]))                 # c', c
        t2 = np.tensordot(fl, AR, ([1], [1])).transpose(1, 0, 2)         # c s' b
        t2 = np.tensordot(X, t2, ([1], [0]))                             # c' s' b
        tot += np.real(np.tensordot(AR.conj(), t2, ([0, 1, 2], [0, 1, 2])))
    return float(tot / nrm)


def charge_ops_q2():
    Z = np.diag([1.0, -1.0]); I2 = np.eye(2)
    k = lambda *o: __import__('functools').reduce(np.kron, o)
    return Z, I2, k


def measure(psi, a, c, sub):
    out = {}
    if a.model == 'trunc':
        Em, FLm, FRm, Qm = c['trunc']
        out['e_bond'] = [a.lam * bond_expect(psi, j, FRm, FLm) for j in range(a.ncells)]
        out['charge'] = [U.one_site_expect(psi, j, Qm) for j in range(a.ncells)]
        out['p_code'] = [U.one_site_expect(psi, j, np.diag((np.arange(len(Em)) < c['N']).astype(float))) for j in range(a.ncells)]
        return out
    if a.model in ('bigsite', 'sparse'):
        H = C.bigsite_H(c)
        FL, FR = C.bigsite_ports(c)
        q, N = c['q'], c['N']
        leakop = np.kron(np.eye(q * q) - c['PomL'], np.eye(N * q * q))
        out['e_bond'] = [a.lam * bond_expect(psi, j, FR, FL) for j in range(a.ncells)]
        out['e_cell_op'] = [U.one_site_expect(psi, j, H) for j in range(a.ncells)]
        out['leakL'] = [U.one_site_expect(psi, j, leakop) for j in range(a.ncells)]
        Q = C.u1_charge(c)[0]
        out['charge'] = [U.one_site_expect(psi, j, Q) for j in range(a.ncells)]
        return out
    if a.model in ('effective', 'core'):
        gam = c['gam']
        Ql = 0.5j * sum(gam[2 * k] @ gam[2 * k + 1] for k in range(len(gam) // 2))
        out['charge'] = [U.one_site_expect(psi, j, Ql) for j in range(a.ncells)]
        return out
    if a.model != 'encoded':
        return out
    q, N = c['q'], c['N']
    bond = sum(np.kron(f, f) for f in c['F'])
    leakop = np.kron(np.eye(q * q) - c['PomL'], np.eye(N * q))     # on (pL, mid): pair L not in |Om>
    eb, eL, eR, leak = [], [], [], []
    for j in range(a.ncells):
        eb.append(a.lam * U.two_site_expect(psi, 3 * j + 2, bond))
        eL.append(U.two_site_expect(psi, 3 * j, c['HLs']))
        eR.append(U.two_site_expect(psi, 3 * j + 1, c['HRs']))
        leak.append(U.two_site_expect(psi, 3 * j, leakop))
    out.update(e_bond=eb, e_HL=eL, e_HR=eR, leakL=leak)
    if q == 2:   # U(1) charge per cell, Q = (Z_pL - Z_sL)/2 + Z_log/2 + (Z_pR - Z_sR)/2 (exact symmetry)
        Z = np.diag([1.0, -1.0]); I2 = np.eye(2)
        Qmid = 0.5 * (-np.kron(np.kron(Z, I2), I2) + np.kron(np.kron(I2, Z), I2) - np.kron(np.kron(I2, I2), Z))
        out['charge'] = [U.one_site_expect(psi, 3 * j, Z / 2) + U.one_site_expect(psi, 3 * j + 1, Qmid)
                         + U.one_site_expect(psi, 3 * j + 2, Z / 2) for j in range(a.ncells)]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('model', choices=['ising', 'logical', 'encoded', 'bigsite', 'sparse', 'effective', 'trunc', 'core'])
    ap.add_argument('--q', type=int, default=2)
    ap.add_argument('--qp', type=int, default=1)
    ap.add_argument('--n', type=int, default=1)
    ap.add_argument('--g', type=float, default=0.5)
    ap.add_argument('--lam', type=float, default=0.1)
    ap.add_argument('--ncells', type=int, default=2)
    ap.add_argument('--mult', type=float, default=0, help='chi multiplier on the two inner bonds (default q)')
    ap.add_argument('--chis', type=int, nargs='+', default=[8, 12, 16, 24, 32, 48, 64])
    ap.add_argument('--mode', default='vumps', choices=['vumps', 'damped'])
    ap.add_argument('--maxit', type=int, default=400)
    ap.add_argument('--tol', type=float, default=1e-8)
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--out', required=True)
    ap.add_argument('--save_state', action='store_true')
    ap.add_argument('--lanczos', action='store_true', help='restarted Lanczos instead of ARPACK')
    ap.add_argument('--kappa', type=float, default=1.0, help='effective model: scale of the lam^2 terms')
    ap.add_argument('--stol', type=float, default=1e-7, help='stop when S changes less than this over 20 iterations')
    ap.add_argument('--M', type=int, default=68, help='trunc: number of cell eigenstates kept')
    ap.add_argument('--verbose', type=int, default=0, help='print progress every k iterations (stderr)')
    ap.add_argument('--load_state', default=None, help='restart from a saved *_state.npz (first chi = its chi)')
    ap.add_argument('--eig_factor', type=float, default=1e-2, help='eigensolver tol = eig_factor * gradient')
    a = ap.parse_args()
    dims, W, mults, sub, info, c = build(a)
    data = dict(args=vars(a), dims=dims, mpo_dims=[w.shape[0] for w in W], info=info, results=[])
    t0 = time.time()
    print(json.dumps(dict(dims=dims, mpo=data['mpo_dims'], info=info)), flush=True)

    def on_chi(psi, rec):
        rec.update(measure(psi, a, c, sub))
        rec['xi_cells'] = rec['xi'] * a.ncells
        rec['e_cell'] = rec['e'] * len(dims) / a.ncells
        rec['time_s'] = time.time() - t0
        data['results'].append(rec)
        json.dump(data, open(a.out + '.tmp', 'w'))
        os.replace(a.out + '.tmp', a.out)
        if a.save_state:
            np.savez(a.out.replace('.json', '') + '_state.npz',
                     **{f'AL{k}': x for k, x in enumerate(psi.AL)}, **{f'AR{k}': x for k, x in enumerate(psi.AR)},
                     **{f'C{k}': x for k, x in enumerate(psi.C)})

    def log(rec):
        S0 = [rec['S'][sub * j] for j in range(len(dims) // sub)]
        extra = ''
        if 'e_bond' in rec:
            extra = f" eb={np.round(rec['e_bond'], 7).tolist()}"
            for key in ('leakL', 'p_code'):
                if key in rec:
                    extra += f" {key}={np.round(rec[key], 5).tolist()}"
            if 'charge' in rec:
                extra += f" Q={np.round(rec['charge'], 5).tolist()}"
        print(f"chi={rec['chi']} {rec['chis_actual']} e_cell={rec['e_cell']:.12f} err={rec['err']:.1e} it={rec['iters']} "
              f"S0={np.round(S0, 6).tolist()} xi_cells={rec['xi_cells']:.3f}{extra} t={rec['time_s']:.0f}s", flush=True)

    expand = None
    if a.model in ('sparse', 'trunc'):
        import sparse_mpo as SPM
        expand = SPM.expand_sp
        if a.lanczos:
            SPM.install_lanczos()
    psi0 = None
    if a.load_state:
        z = np.load(a.load_state)
        nn = len(dims)
        import vumps as V
        psi0 = V.UMPS([z[f'AL{k}'] for k in range(nn)], [z[f'AR{k}'] for k in range(nn)], [z[f'C{k}'] for k in range(nn)])
        if os.path.exists(a.out):
            data['results'] = json.load(open(a.out))['results']
    U.run_md(W, dims, a.chis, mults=mults, tol=a.tol, maxit=a.maxit, seed=a.seed, log=log, on_chi=on_chi,
             mode=a.mode, expand=expand, eig_factor=a.eig_factor, psi=psi0, verbose=a.verbose, stol=a.stol)


if __name__ == '__main__':
    main()
