#!/usr/bin/env python3
"""Task 2 driver: torus route-B cells, their second-order chains, symmetry checks, current
projection and one-loop flows.

Configurations (Delta = 1, q = 4, q+ = 2, n = 4, cell dimension 4096):
  untuned  g = 0.5 for every flavour, bond weights 1, mu0 = 0    (the red team's (4,2) cell)
  uniform_mu  the untuned cell plus a common chemical potential mu0 (equal Fermi points)
  tuned    flavour-dependent dressing g_a (torus-symmetric, U(2) x U(2) broken, code still exactly
           degenerate), so every flavour has its own hopping, plus a common chemical potential mu0
           on the total charge, so that every flavour has its own Fermi points.
For each: exact 2^n degeneracy, compressions, [T_s, H_cell] and [T_s, bond]; T1, T2, T3 and
their logical torus symmetry; Majorana bilinears -> Fermi points of each flavour; quartic part
-> so(8) current couplings (route_b_second_order.current_couplings, c_{j,A} -> psi_{A,0} +
(-1)^j psi_{A,pi}); decomposition of the couplings by flavour charge (q at the k ~ 0 point, q' at
the k ~ pi point), with torus charge W(q + q') and momentum mismatch q.k_0 + q'.k_pi; one-loop
flows dG/dl = (1/4 pi v) f f G G with the components whose mismatch exceeds the running cutoff
e^{-l} removed.  kappa multiplies the lam^2 terms (kappa = 1 is the physical chain; kappa = 10 is a
stress test).  Writes data/eff_<tag>.npz (T1, T2, factors of T3, Q2) for run_vumps / run_obc.

usage: python3 routeB_torus.py [--cache DIR] [--tags ...] [--lams ...] [--kappas ...]
       (eigendecompositions cached in DIR; about 25 minutes for the three configurations)
Written by Claude (Opus 5.5) for contrib/torus_routeB, 2026-10-27.
"""
import os, sys, json, argparse, itertools, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np
from scipy import linalg
import torus_cell as TC
import route_b_second_order as RB
import cells as C

CONFIGS = {
    'untuned': dict(g=[0.5, 0.5, 0.5, 0.5], w=[1.0, 1.0, 1.0, 1.0], mu0_rel=0.0),
    'uniform_mu': dict(g=[0.5, 0.5, 0.5, 0.5], w=[1.0, 1.0, 1.0, 1.0], mu0_rel=0.3),
    'tuned': dict(g=[0.35, 0.5, 0.65, 0.8], w=[1.0, 1.0, 1.0, 1.0], mu0_rel=0.3),
}


def current_basis_charges(c):
    """ad-charges of the so(8) currents J^{AB} = i psi_A psi_B (A<B, itertools order) under the
    flavour numbers n_a = (1 + i gamma_2a gamma_2a+1)/2, computed in the logical Clifford module.
    Returns U (columns: charge eigencurrents in the J basis), Uinv, q (n_eig x nf)."""
    g = c['gam']
    m = len(g)
    pairs = list(itertools.combinations(range(m), 2))
    J = [1j * g[A] @ g[B] for (A, B) in pairs]
    nf = m // 2
    nops = [(np.eye(c['N']) + 1j * g[2 * a] @ g[2 * a + 1]) / 2 for a in range(nf)]
    nrm = [np.trace(x.conj().T @ x).real for x in J]
    K = []
    for a in range(nf):
        M = np.zeros((len(J), len(J)), complex)
        for Aidx, JA in enumerate(J):
            comm = nops[a] @ JA - JA @ nops[a]
            for Bidx, JB in enumerate(J):
                M[Bidx, Aidx] = np.trace(JB.conj().T @ comm) / nrm[Bidx]
        K.append(M)
    rng = np.random.default_rng(7)
    coef = rng.normal(size=nf) + 3.0
    w, U = np.linalg.eig(sum(cf * Km for cf, Km in zip(coef, K)))
    Uinv = np.linalg.inv(U)
    q = np.array([[(Uinv @ Km @ U)[mu, mu] for Km in K] for mu in range(len(J))])
    assert np.abs(q.imag).max() < 1e-9 and np.abs(q.real - np.round(q.real)).max() < 1e-9
    return U, Uinv, np.round(q.real).astype(int), pairs


def fermi_points(quad, m, nk=20001):
    """For each flavour a, the two zero crossings (k1 near 0, k2 near pi) of the particle branch
    of the flavour block of h(k), and the velocities there."""
    ks = np.linspace(-np.pi, np.pi, nk, endpoint=False)
    H = TC.flavour_dispersion(quad, m, ks)
    nf = m // 2
    out = []
    for a in range(nf):
        blk = H[:, 2 * a:2 * a + 2, 2 * a:2 * a + 2]
        # a complex fermion f = (gamma_2a + i gamma_2a+1)/2: eps_f(k) is the eigenvalue of h(k) along (1, -i)/sqrt2 or (1, i)
        u = np.array([1, 1j]) / np.sqrt(2)
        eps = np.real(np.einsum('i,kij,j->k', u.conj(), blk, u))
        z = []
        for i in range(nk):
            j = (i + 1) % nk
            if eps[i] == 0 or eps[i] * eps[j] < 0:
                kk = ks[i] + (ks[1] - ks[0]) * eps[i] / (eps[i] - eps[j])
                vv = (eps[j] - eps[i]) / (ks[1] - ks[0])
                z.append((kk, vv))
        out.append(z)
    return out


def flow(Gm, f, v, Uinv, U, P, lmax=40.0, Lam0=1.0, h=0.002, Gstop=None):
    """One-loop flow in the J basis; component (mu, nu) of Uinv G Uinv^T is removed once
    |P[mu, nu]| > Lam0 e^{-l}.  Runaway: max|G| > max(Gstop, 5 max|G0|)."""
    kap = 1.0 / (4 * np.pi * v)
    nA = Gm.shape[0]
    fr = f.reshape(nA * nA, nA)

    def B(G):
        X = np.tensordot(f, G, ([1], [1]))
        t = np.tensordot(G, X, ([1], [0]))
        return kap * (fr.T @ t.transpose(0, 2, 1).reshape(nA * nA, nA))

    def proj(G, l):
        Gt = Uinv @ G @ Uinv.T
        Gt = Gt * (np.abs(P) <= Lam0 * np.exp(-l) + 1e-12)
        return np.real(U @ Gt @ U.T)
    G = proj(Gm, 0.0)
    g0 = np.abs(G).max()
    gstop = max(Gstop if Gstop is not None else 0.0, 5 * g0)
    l = 0.0
    while l < lmax:
        if np.abs(G).max() > gstop:
            return 'runaway', l, G
        dl = min(h / (kap * max(np.abs(G).max(), 1e-12)), 0.05)
        k1 = B(G); k2 = B(G + 0.5 * dl * k1)
        G = proj(G + dl * k2, l + dl)
        l += dl
    return 'bounded', l, G


def analyse(tag, cfg, lams, cache=None, log=print, kappas=(1.0,)):
    t0 = time.time()
    c = TC.build_cell(4, 2, 4, 1.0, cfg['g'], cfg['g'])
    w = np.asarray(cfg['w'], float)
    sign, terr, berr, _ = TC.check_cell_torus(c, w)
    log(f'\n=== {tag}: g = {cfg["g"]}, bond weights w = {cfg["w"]} ===')
    log(f'torus: logical charge sign {sign:+d}; max_s |[T_s, H_cell]| = {terr:.1e}, |[T_s, bond]| = {berr:.1e}')
    fn = os.path.join(cache, f'eig_{tag}.npz') if cache else None
    if fn and os.path.exists(fn):
        z = np.load(fn); E, V, G = z['E'], z['V'], z['G']
    else:
        E, V, G = TC.code_and_frame(c)
        if fn:
            np.savez(fn, E=E, V=V, G=G)
    so = TC.second_order(c, w, E, V, G)
    N, r = c['N'], c['r']
    lev = np.unique(np.round(E - E[0], 8))
    log(f'cell: code splitting {so["split"]:.1e}, gap {so["gap"]:.4f}, lowest levels {lev[:5]}, '
        f'compression residual {so["resid"]:.1e} ({time.time()-t0:.0f}s)')
    aL, aR = so['aL'], so['aR']
    tA = w[np.arange(r) // 2] * aL * aR
    log('compressions per flavour: aL = ' + str(np.round(aL[::2], 5)) + ', aR = ' + str(np.round(aR[::2], 5))
        + ', first-order Majorana hopping / lam = w aL aR = ' + str(np.round(tA[::2], 5)))
    eT1 = TC.check_logical_torus(c, so['T1'], 2)
    eT2 = TC.check_logical_torus(c, so['T2'], 2)
    eT3 = TC.check_logical_torus(c, so['T3'], 3)
    log(f'logical torus: relative |[T, T1]| = {eT1:.1e}, |[T, T2]| = {eT2:.1e}, |[T, T3]| = {eT3:.1e}')
    nl = TC.logical_numbers(c)
    Q = sum(nl)
    Q2 = np.kron(Q, np.eye(N))
    Ts = TC.bigsite_torus(c, sign)
    X = G.conj().T @ (Ts[0] + Ts[1]) @ G
    cst = np.trace(X + Q).real / N
    log(f'physical upper-species charge on the code = -Q + {cst:.3f}, residual {np.linalg.norm(X + Q - cst * np.eye(N)):.1e} '
        f'(so mu0 Q is a physical, torus-symmetric term, added exactly after SW)')
    # T3 factors (weighted)
    Mfac = np.zeros((r, r, N, N), complex)
    for A in range(r):
        for Bb in range(r):
            Mfac[A, Bb] = -(so['YRL'][Bb, A] + so['YLR'][A, Bb])
    CRw = np.array([so['wA'][A] * so['CR'][A] for A in range(r)])
    CLw = np.array([so['wA'][A] * so['CL'][A] for A in range(r)])
    T3re = sum(np.kron(np.kron(CRw[A], Mfac[A, Bb]), CLw[Bb]) for A in range(r) for Bb in range(r))
    T3re = (T3re + T3re.conj().T) / 2
    log(f'T3 from factors: |T3 - rebuilt| = {np.linalg.norm(T3re - so["T3"]):.1e}')
    np.savez(os.path.join(HERE, 'data', f'eff_{tag}.npz'), T1=so['T1'], T2=so['T2'], CRw=CRw, CLw=CLw, M=Mfac,
             Q2=Q2, alpha2=tA, g=np.asarray(cfg['g']), w=w)
    # Majorana expansion
    view = TC.RBView(c, so, 1.0)
    quart = RB.unit_cell_terms(view, 4)
    m = 2 * c['n']
    Gm, pairs = RB.current_couplings(quart, m)
    f = RB.so_structure(m)
    U, Uinv, q, pairs2 = current_basis_charges(c)
    Wt = TC.torus_charge_matrix(2, 4)               # (flavour, species)
    Gt = Uinv @ Gm @ Uinv.T
    tor = np.array([[np.abs((q[mu] + q[nu]) @ Wt).sum() for nu in range(len(q))] for mu in range(len(q))])
    viol = np.abs(Gt[tor > 0]).max() if (tor > 0).any() else 0.0
    log(f'quartic monomials {len(quart)}; |G| (per lam^2) = {np.linalg.norm(Gm):.4f}; largest torus-charged '
        f'component {viol:.1e} (allowed components {int((tor == 0).sum())} of {tor.size})')
    # classify allowed components
    kinds = {}
    for mu in range(len(q)):
        for nu in range(len(q)):
            if tor[mu, nu] or abs(Gt[mu, nu]) < 1e-10:
                continue
            a, b = tuple(q[mu]), tuple(q[nu])
            if not any(a) and not any(b):
                kd = 'Cartan (forward)'
            elif sum(q[mu]) == 0 and np.all(q[mu] == -q[nu]):
                kd = 'exchange e_a - e_b | e_b - e_a'
            elif sum(q[mu]) != 0 and np.all(q[mu] == -q[nu]):
                kd = 'umklapp e_a + e_b | -(e_a + e_b)'
            else:
                kd = 'rectangle'
            kinds[kd] = kinds.get(kd, 0.0) + abs(Gt[mu, nu]) ** 2
    log('  allowed couplings by kind (Frobenius, per lam^2): '
        + ', '.join(f'{k}: {np.sqrt(v):.4f}' for k, v in kinds.items()))
    res = dict(tag=tag, lams={})
    for lam in lams:
      for kap in kappas:
        mu0 = cfg['mu0_rel'] * 4 * lam * float(np.mean(tA))
        quad = TC.quadratic_terms(c, so, lam, mu0, kap)
        fp = fermi_points(quad, m)
        k0 = np.array([min(z, key=lambda x: abs(x[0]))[0] for z in fp])
        kpi = np.array([min(z, key=lambda x: abs(abs(x[0]) - np.pi))[0] for z in fp])
        vs = np.array([abs(min(z, key=lambda x: abs(x[0]))[1]) for z in fp])
        P = np.array([[q[mu] @ k0 + q[nu] @ kpi for nu in range(len(q))] for mu in range(len(q))])
        P = np.abs((P + np.pi) % (2 * np.pi) - np.pi)
        allowed = (tor == 0) & (np.abs(Gt) > 1e-10)
        nonres = allowed & (P > 1e-6)
        pmin = P[nonres].min() if nonres.any() else np.nan
        vbar = float(np.mean(vs))
        Gabs = kap * lam ** 2 * Gm
        fa = flow(Gabs, f, vbar, Uinv, U, np.where(tor == 0, 0.0, 10.0), lmax=80)          # nothing detuned
        fb = flow(Gabs, f, vbar, Uinv, U, np.where(tor == 0, P, 10.0), lmax=80)            # detuned
        ncart = (np.abs(q).sum(1)[:, None] + np.abs(q).sum(1)[None, :]) > 0
        log(f'lam = {lam}, kappa = {kap}: mu0 = {mu0:.5f}; Fermi points k0/pi = {np.round(k0 / np.pi, 4)}, '
            f'kpi/pi = {np.round(kpi / np.pi, 4)}, velocities {np.round(vs, 4)}; smallest nonzero mismatch '
            f'of an allowed coupling {pmin:.4f}; max|G|/(2 pi v) = {np.abs(Gabs).max() / (2 * np.pi * vbar):.4f}')
        log(f'   one-loop, all allowed couplings resonant: {fa[0]} (l = {fa[1]:.2f}); '
            f'with the actual detuning: {fb[0]} (l = {fb[1]:.2f}), final max non-Cartan coupling '
            f'{np.abs((Uinv @ fb[2] @ Uinv.T)[ncart]).max():.2e} (initial {np.abs((Uinv @ Gabs @ Uinv.T)[ncart & (tor == 0)]).max():.2e})')
        res['lams'][f'{lam}_k{kap}'] = dict(mu0=mu0, k0=k0.tolist(), kpi=kpi.tolist(), v=vs.tolist(), pmin=float(pmin),
                                            flow_all=[fa[0], fa[1]], flow_detuned=[fb[0], fb[1]])
    return res


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--cache', default=None)
    ap.add_argument('--tags', nargs='+', default=['untuned', 'tuned'])
    ap.add_argument('--lams', type=float, nargs='+', default=[0.1, 0.2, 0.4])
    ap.add_argument('--kappas', type=float, nargs='+', default=[1.0, 10.0])
    a = ap.parse_args()
    np.set_printoptions(linewidth=160)
    out = {}
    for tag in a.tags:
        out[tag] = analyse(tag, CONFIGS[tag], a.lams, a.cache, log=lambda s: print(s, flush=True), kappas=a.kappas)
    json.dump(out, open(os.path.join(HERE, 'data', 'routeB_torus_summary.json'), 'w'), indent=1)
