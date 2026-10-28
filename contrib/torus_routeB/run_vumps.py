#!/usr/bin/env python3
"""VUMPS finite-entanglement scaling for the task-1 chains (local dimension 16) and for the
task-2 effective logical chains (torus_cell.py, local dimension 16, two- and three-site terms).

usage:
  run_vumps.py model --case ii --s 1 --g 0.3 --seed 0 --ncells 2 --chis 16 24 32 ... --out f.json
  run_vumps.py eff   --tag tuned --lam 0.2 --ncells 2 --chis ... --out f.json   (needs data/eff_<tag>.npz)
Reuses ../encoded_chain_numerics/umps_md.py (multi-site VUMPS) and cells.mpo_from_terms.
Per chi it records the energy per site, entropies on every bond, xi (sites) from the transfer
matrix, and flavour densities <n_a>.

Written by Claude (Opus 5.5) for contrib/torus_routeB, 2026-10-27.
"""
import os, sys, json, time, argparse
for _k in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ.setdefault(_k, '4')
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'encoded_chain_numerics'))
import numpy as np
import cells as C          # noqa: E402  (encoded_chain_numerics)
import umps_md as U        # noqa: E402
import fermions as F       # noqa: E402
import tcases as M         # noqa: E402


def load_eff(tag, here):
    """T1, T2, T3 (rebuilt from its factors), Q2 of data/eff_<tag>.npz (routeB_torus.py)."""
    z = np.load(os.path.join(here, 'data', f'eff_{tag}.npz'))
    r = z['M'].shape[0]
    T3 = sum(np.kron(np.kron(z['CRw'][A], z['M'][A, B]), z['CLw'][B]) for A in range(r) for B in range(r))
    T3 = (T3 + T3.conj().T) / 2
    return z['T1'], z['T2'], T3, z['Q2'], z['alpha2']


def build(a):
    if a.mode == 'model':
        k, h, inter = M.two_site_term(a.case, a.s, a.g, a.seed, no=bool(a.no))
        dims = [16] * a.ncells
        terms = [(j, 2, h) for j in range(a.ncells)]
        return dims, C.mpo_from_terms(dims, terms), dict(k=list(map(float, k)))
    T1, T2, T3, Q2, al2 = load_eff(a.tag, HERE)
    T2, T3 = a.kappa * T2, a.kappa * T3                  # kappa: amplification of the lam^2 terms
    mu0 = a.mu0 * 4 * a.lam * float(np.mean(al2))       # same convention as routeB_torus (mu0_rel)
    dims = [16] * a.ncells
    sc = 1.0 / (a.lam * float(np.mean(al2)))             # energy unit: the mean first-order hopping
    if a.order >= 2:
        h2 = sc * (a.lam * T1 + a.lam ** 2 * T2 + mu0 * Q2)
        terms = [(j, 2, h2) for j in range(a.ncells)] + [(j, 3, sc * a.lam ** 2 * T3) for j in range(a.ncells)]
    else:
        h2 = sc * (a.lam * T1 + mu0 * Q2)
        terms = [(j, 2, h2) for j in range(a.ncells)]
    return dims, C.mpo_from_terms(dims, terms), dict(alpha2=al2.tolist(), mu0=mu0, energy_unit=1 / sc)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('mode', choices=['model', 'eff'])
    ap.add_argument('--case', default='ii')
    ap.add_argument('--s', type=float, default=1.0)
    ap.add_argument('--g', type=float, default=0.3)
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--no', type=int, default=0, help='model: 1 = Hartree-subtracted densities')
    ap.add_argument('--tag', default='tuned')
    ap.add_argument('--lam', type=float, default=0.2)
    ap.add_argument('--mu0', type=float, default=0.0, help='eff: mu0_rel, total-charge term mu0_rel * 4 lam mean(t_A) * Q')
    ap.add_argument('--order', type=int, default=2)
    ap.add_argument('--kappa', type=float, default=1.0)
    ap.add_argument('--ncells', type=int, default=2)
    ap.add_argument('--chis', type=int, nargs='+', default=[16, 24, 32, 48, 64])
    ap.add_argument('--mode_v', default='vumps', choices=['vumps', 'damped'])
    ap.add_argument('--maxit', type=int, default=400)
    ap.add_argument('--tol', type=float, default=1e-8)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    dims, W, info = build(a)
    data = dict(args=vars(a), dims=dims, mpo=[w.shape[0] for w in W], info=info, results=[])
    print(json.dumps(dict(mpo=data['mpo'], info=info)), flush=True)
    c, P = F.site_ops(4)
    nops = [x.conj().T @ x for x in c]
    t0 = time.time()

    def on_chi(psi, rec):
        rec['n'] = [[U.one_site_expect(psi, j, n) for n in nops] for j in range(a.ncells)]
        rec['xi_sites'] = rec['xi'] * a.ncells
        rec['e_site'] = rec['e']
        rec['time_s'] = time.time() - t0
        data['results'].append(rec)
        json.dump(data, open(a.out + '.tmp', 'w'))
        os.replace(a.out + '.tmp', a.out)

    def log(rec):
        print(f"chi={rec['chi']} {rec['chis_actual']} e={rec['e']:.12f} err={rec['err']:.1e} it={rec['iters']} "
              f"S={np.round(rec['S'], 5).tolist()} xi={rec['xi_sites']:.3f} "
              f"n={np.round(np.mean(rec['n'], 0), 4).tolist()} t={rec['time_s']:.0f}s", flush=True)
    U.run_md(W, dims, a.chis, tol=a.tol, maxit=a.maxit, log=log, on_chi=on_chi, mode=a.mode_v)


if __name__ == '__main__':
    main()
