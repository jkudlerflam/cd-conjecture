#!/usr/bin/env python3
"""VUMPS finite-entanglement scaling for the two-flavour demonstrator and its calibrations.

models
  micro : the full d = 13 bond (cell2.bond), n-site unit cell, random start (types not imposed), damped steps
  cell3 : the same Hamiltonian restricted to its conserved period-3 type sector, unit cell (P_L, B, P_R) with
          site dimensions (9, 4, 9); exact, since no term changes a letter's type
  cellbig: the same sector as one 324-dimensional site per cell (P_L B P_R); energies per cell
  eff   : the logical chain on d = 4 (n1, n2) per cell, with the exact two-cell bond at --lam (--order 2 for
          the second-order bond)
  xx1   : one XX chain (c = 1 control), d = 2
  xx2   : two decoupled XX chains with velocities 1 and --ratio (c = 2 calibration), d = 4
Records per chi: energy, gradient error, entropies on every bond, transfer-matrix spectrum, xi (in unit cells),
and <n1>, <n2> on the logical site.  Saves the MPS (npz) when --state_dir is given.
Claude (Opus 5.5), 2026-10-28.
"""
import os, sys, json, time, argparse
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np
import vumps2 as V
import cell2 as C

PAR = dict(DN=1.0, DC=1.0, kap=(0.6, 0.6), sg=-1.0, DCd=(0.0, -0.2), Jp=1.0, Jz=1.0)


def gauge_phases(lam):
    """phi_a = arg(t_a) + pi from the exact two-cell bond: in this frame both logical hoppings are real and
    negative (a ferromagnetic sign, which a one-site uMPS represents without a sublattice cat)."""
    import effective2 as E
    He, _, _ = E.exact_eff(dict(PAR), lam)
    c, _ = E.fit(He)
    return tuple(float(np.angle(c[f'hop{a}'] + 1j * c[f'hopi{a}']) + np.pi) for a in (1, 2))


def gauge_bond16(H, ph):
    from scipy.linalg import expm
    n1 = np.kron(np.diag([0, 1.]), np.eye(2)); n2 = np.kron(np.eye(2), np.diag([0, 1.]))
    Dm = expm(1j * np.kron(np.eye(4), ph[0] * n1 + ph[1] * n2))
    return Dm @ H @ Dm.conj().T


def eff_bond(lam, order):
    import effective2 as E
    par = dict(PAR)
    if order == 2:
        _, H2, _, _ = E.second_order(par)
        return lam ** 2 * H2
    He, _, _ = E.exact_eff(par, lam)
    return He


def model(a):
    nop = [np.kron(C.NOP, np.eye(2)), np.kron(np.eye(2), C.NOP)]
    if a.model in ('micro', 'cell3'):
        g = gauge_phases(a.lam) if a.gauge else (0.0, 0.0)
        print('gauge phases', g, flush=True)
        h = C.bond(lam=a.lam, mu=a.mu, K=a.K, **PAR, h2=0.0, gauge=g)
        W = V.mpo_from_term(h, 2, d=13)
        ops13 = []
        for o in nop:
            m = np.zeros((13, 13)); m[9:, 9:] = o; ops13.append(m)
        nB = np.diag([0.] * 9 + [1.] * 4)
        if a.model == 'micro':
            return [W] * a.n, 13, [nB] + ops13
        WP, WB = W[:, :, :9, :9].copy(), W[:, :, 9:, 9:].copy()
        return [WP, WB, WP], [9, 4, 9], None
    if a.model == 'cellbig':
        H0, ports, code, qc = C.cell_ops(PAR['DN'], PAR['DC'], PAR['kap'], PAR['sg'], PAR['DCd'], 0.0)
        terms = C.bond_terms(PAR['Jp'], PAR['Jz'])
        D = 2 + len(terms)
        Wc = np.zeros((D, D, 324, 324), complex)
        Wc[0, 0] = np.eye(324); Wc[D - 1, D - 1] = np.eye(324)
        Wc[0, D - 1] = H0 - a.mu * np.eye(324)
        for i, (cf, r, l) in enumerate(terms):
            Wc[0, 1 + i] = a.lam * cf * ports['R'][r]
            Wc[1 + i, D - 1] = ports['L'][l]
        lo = []
        for o in nop:
            lo.append(np.kron(np.kron(np.eye(9), o), np.eye(9)))
        return [Wc], 324, lo
    if a.model == 'eff':
        h = eff_bond(a.lam, a.order)
        if a.gauge:
            h = gauge_bond16(h, gauge_phases(a.lam))
        h = h - np.trace(h) / 16 * np.eye(16)
        if not a.raw:
            h = h / np.abs(h).max()      # overall scale is irrelevant for c; keeps solver tolerances meaningful
        return [V.mpo_from_term(h, 2, d=4)] * a.n, 4, nop
    if a.model == 'flav':
        # one flavour of the effective chain: the exact two-cell bond with the other flavour traced out
        # (inter-flavour terms are O(lam^4)); used to build warm starts with the right fillings
        H = eff_bond(a.lam, a.order).reshape(2, 2, 2, 2, 2, 2, 2, 2)      # (n1 n2 n1' n2' | ...)
        if a.flavour == 1:
            h = np.einsum('abcdebgd->aceg', H) / 4
        else:
            h = np.einsum('abcdafch->bdfh', H) / 4
        h = h.reshape(4, 4); h = h - np.trace(h) / 4 * np.eye(4); h = h / np.abs(h).max()
        return [V.mpo_from_term(h, 2, d=2)] * a.n, 2, [C.NOP]
    if a.model == 'xx1':
        hx = -(np.kron(C.SP, C.SM) + np.kron(C.SM, C.SP))
        return [V.mpo_from_term(hx, 2, d=2)] * a.n, 2, [C.NOP]
    if a.model == 'xx2':
        hx = -(np.kron(C.SP, C.SM) + np.kron(C.SM, C.SP))
        h4 = np.kron(hx, np.eye(4)) + a.ratio * np.kron(np.eye(4), hx)      # (a1 a2 b1 b2)
        h = h4.reshape([2] * 8).transpose(0, 2, 1, 3, 4, 6, 5, 7).reshape(16, 16)
        return [V.mpo_from_term(h, 2, d=4)] * a.n, 4, nop
    raise ValueError(a.model)


def product_start(a, n):
    """chi = 1 uMPS: Bell states on the P letters and a U(1)-breaking superposition on each logical mode,
    so that the fillings can relax continuously (a random start can stall in an empty logical sector)."""
    th = [np.arcsin(np.sqrt(f)) for f in a.fill]
    q = [np.array([np.cos(t), np.sin(t)]) for t in th]
    logical = np.kron(q[0], q[1]).astype(complex)
    if a.model in ('eff', 'xx2'):
        vecs = [logical] * n
    elif a.model == 'cell3':
        vecs = [C.bell('L'), logical, C.bell('R')]
    elif a.model == 'micro':
        def emb(v, off):
            w = np.zeros(13, complex); w[off:off + len(v)] = v; return w
        pat = [emb(C.bell('L'), 0), emb(logical, 9), emb(C.bell('R'), 0)]
        vecs = [pat[k % 3] for k in range(n)]
    A = [v.reshape(1, -1, 1) / np.linalg.norm(v) for v in vecs]
    return V.UMPS(A, [x.copy() for x in A], [np.ones((1, 1), complex) for _ in range(n)])


def pair_start(files, model):
    """Product of two single-flavour uMPS (bond dimension chi1 chi2), on B alone ('eff') or embedded in the
    microscopic cell with Bell states on the P letters ('cell3')."""
    z1, z2 = np.load(files[0]), np.load(files[1])
    def comb(X, Y):
        T = np.einsum('asb,ctd->acstbd', X, Y)
        return T.reshape(X.shape[0] * Y.shape[0], 4, X.shape[2] * Y.shape[2])
    AL, AR = comb(z1['AL0'], z2['AL0']), comb(z1['AR0'], z2['AR0'])
    Cc = np.kron(z1['C0'], z2['C0'])
    if model == 'eff':
        return V.UMPS([AL], [AR], [Cc])
    chi = AL.shape[0]; I = np.eye(chi)
    PL = np.einsum('s,ab->asb', C.bell('L'), I); PR = np.einsum('s,ab->asb', C.bell('R'), I)
    return V.UMPS([PL, AL, PR], [PL.copy(), AR, PR.copy()], [Cc, Cc.copy(), Cc.copy()])


def embed_eff(f):
    """Microscopic (P_L, B, P_R) uMPS from a logical uMPS: Bell states on the P letters (identity on the bond)
    and the logical tensors on B.  This is the zeroth-order code embedding of the logical ground state; the
    microscopic VUMPS then adds the port dressing."""
    z = np.load(f)
    AL, AR, Cc = z['AL0'], z['AR0'], z['C0']
    chi = AL.shape[0]
    I = np.eye(chi)
    PL = np.einsum('s,ab->asb', C.bell('L'), I); PR = np.einsum('s,ab->asb', C.bell('R'), I)
    return V.UMPS([PL, AL, PR], [PL.copy(), AR, PR.copy()], [Cc, Cc.copy(), Cc.copy()])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('name')
    ap.add_argument('--model', default='cell3')
    ap.add_argument('--lam', type=float, default=0.15)
    ap.add_argument('--mu', type=float, default=1.2)
    ap.add_argument('--K', type=float, default=3.0)
    ap.add_argument('--order', type=int, default=0)
    ap.add_argument('--ratio', type=float, default=0.7)
    ap.add_argument('--flavour', type=int, default=1)
    ap.add_argument('--gauge', type=int, default=1, help='work in the frame where both logical hoppings are real and negative')
    ap.add_argument('--raw', type=int, default=0, help='eff: keep the physical scale (needed for warm starts)')
    ap.add_argument('--load_pair', nargs=2, default=None, help='two single-flavour npz files: start from their product')
    ap.add_argument('--chis', type=int, nargs='+', default=[16, 24, 32])
    ap.add_argument('--n', type=int, default=3)
    ap.add_argument('--mode', default='vumps')
    ap.add_argument('--switch_chi', type=int, default=None)
    ap.add_argument('--maxit', type=int, default=400)
    ap.add_argument('--tol', type=float, default=1e-10)
    ap.add_argument('--stol', type=float, default=1e-9)
    ap.add_argument('--tau_max', type=float, default=2000.0)
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--state_dir', default=None)
    ap.add_argument('--load_eff', default=None, help='npz of a logical (eff) uMPS to embed as the start (cell3)')
    ap.add_argument('--init', default='random', help="'random' or 'product' (start at chi = 1; put 1 first in --chis)")
    ap.add_argument('--fill', type=float, nargs=2, default=[0.5, 0.4], help='initial logical fillings for --init product')
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    Ws, d, ops = model(a)
    n = len(Ws)
    out = a.out or os.path.join(HERE, 'data', 'fes', a.name + '.json')
    os.makedirs(os.path.dirname(out), exist_ok=True)
    data = dict(name=a.name, args=vars(a), mpo_D=int(Ws[0].shape[0]), results=[])
    t0 = time.time()

    def logical_ops(psi):
        if a.model == 'cell3':
            # B is site 1 of the cell, dimension 4 = (n1, n2)
            res = np.zeros((2, n))
            for i, o in enumerate([np.kron(C.NOP, np.eye(2)), np.kron(np.eye(2), C.NOP)]):
                AC = np.tensordot(psi.AL[1], psi.C[2 % n], ([2], [0]))
                nrm = np.vdot(AC, AC).real
                res[i, 1] = (np.vdot(AC, np.tensordot(o, AC, ([1], [1])).transpose(1, 0, 2)) / nrm).real
            return res.tolist()
        return V.local_expect(psi, ops).tolist()

    def save(psi, rec):
        rec['time_s'] = time.time() - t0
        rec['local'] = logical_ops(psi)
        rec['xi_cells'] = rec['xi'] / (3 if a.model in ('micro', 'cell3') else 1)
        data['results'].append(rec)
        json.dump(data, open(out + '.tmp', 'w'))
        os.replace(out + '.tmp', out)
        if a.state_dir:
            os.makedirs(a.state_dir, exist_ok=True)
            np.savez_compressed(os.path.join(a.state_dir, f"{a.name}_chi{rec['chi']}.npz"),
                                **{f'AL{k}': psi.AL[k] for k in range(psi.n)}, **{f'AR{k}': psi.AR[k] for k in range(psi.n)},
                                **{f'C{k}': psi.C[k] for k in range(psi.n)})
        print(f"{a.name} chi={rec['chi']} e={rec['e']:.12f} err={rec['err']:.1e} it={rec['iters']} "
              f"S={np.round(rec['S'], 5).tolist()} xi_cells={rec['xi_cells']:.3f} local={np.round(rec['local'], 4).tolist()} "
              f"t={rec['time_s']:.0f}s", flush=True)

    kw = dict(tol=a.tol, maxit=a.maxit, seed=a.seed, log=lambda r: None, stol=a.stol, tau_max=a.tau_max)
    psi0 = product_start(a, n) if a.init == 'product' else None
    if a.load_eff:
        psi0 = embed_eff(a.load_eff)
    if a.load_pair:
        psi0 = pair_start(a.load_pair, a.model)
    if a.switch_chi is None:
        V.run(Ws, n, a.chis, d=d, mode=a.mode, on_chi=save, psi=psi0, **kw)
    else:
        lo = [c for c in a.chis if c <= a.switch_chi]; hi = [c for c in a.chis if c > a.switch_chi]
        psi, _ = V.run(Ws, n, lo, d=d, mode=a.mode, on_chi=save, psi=psi0, **kw)
        V.run(Ws, n, [lo[-1]] + hi, d=d, mode='damped', psi=psi,
              on_chi=lambda p, r: save(p, r) if r['chi'] in hi else None, **kw)


if __name__ == '__main__':
    main()
