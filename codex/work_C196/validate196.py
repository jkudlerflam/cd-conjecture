"""C196 local validation (numpy, small chi).  python3 -E -s validate196.py [part]

part ed:    iDMRG without truncation reproduces open-chain ED (same bulk MPO, P = 0) for L = 2 to 10.
part tm:    Arnoldi transfer-matrix eigenvalues against dense eig.
part ctrl:  XX, two-leg XX, TFIM energy densities against the exact values at small chi.
part cand:  p153c and sym40_n316 at small chi against the C190 open-chain bulk energy and density.
"""
import os
import sys
import json
import time
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import z2dmrg as Z0     # noqa: E402
import z2dmrg_xp as X   # noqa: E402
import idmrg196 as I    # noqa: E402

HZ = dict(p153c=-19.97589125, sym40_n316=-19.976047250613814, rat_n316=-19.975988986408243)
MU_C = float(-2 * np.cos(0.1957 * np.pi))


def dense_H(W, L):
    w = W.shape[0]
    l = np.zeros(w, dtype=complex); l[0] = 1
    r = np.zeros(w, dtype=complex); r[2] = 1
    T = np.einsum('w,wvst->vst', l, W)                      # (v, s', s)
    T = T.reshape(w, 2, 2)
    for _ in range(L - 1):
        T = np.einsum('vab,vust->uasbt', T, W)
        d = T.shape[1] * 2
        T = T.reshape(w, d, d)
    H = np.einsum('v,vab->ab', r, T)
    return H


def part_ed():
    print('# part ed: untruncated iDMRG against ED of the open chain, P = 0 sector')
    cases = [dict(model='p153', hz=HZ['p153c']), dict(model='sym40', hz=HZ['sym40_n316']),
             dict(model='xx2', mu=MU_C), dict(model='tfim', h=1.05)]
    worst = 0.0
    for cfg in cases:
        mpo, W, e_off = I.build_model(cfg)
        st = I.IDMRG(mpo, W, cutoff=0.0, seed=3)
        line = []
        for t in range(5):
            rec = st.step(chi=10 ** 6)
            L = 2 * (t + 1)
            H = dense_H(W, L)
            par = np.array([bin(i).count('1') % 2 for i in range(2 ** L)])
            sel = par == 0
            e0 = np.linalg.eigvalsh(H[np.ix_(sel, sel)])[0]
            d = abs(rec['E'] - e0)
            worst = max(worst, d)
            line.append(f'L={L} dE={d:.1e} res={rec["res"]:.1e}')
        print(cfg['model'], '; '.join(line))
    print(f'worst |E_idmrg - E_ED| = {worst:.2e}')
    return worst


def run_stages(cfg, chis, nmin=20, win=20, tolE=1e-10, tolS=1e-5, cap=4000, verbose=True):
    mpo, W, e_off = I.build_model(cfg)
    st = I.IDMRG(mpo, W, cutoff=1e-14, seed=0)
    t0 = time.time()
    fh = st.init_finite(L0=32, chi0=32)
    if verbose:
        print(f"  finite init L0=32: E={fh[-1]['E']:.10f} res={fh[-1]['res_max']:.1e} t={time.time() - t0:.0f}s", flush=True)
    pts = []
    for chi in chis:
        hist = []
        conv = False
        t0 = time.time()
        for it in range(cap):
            rec = st.step(chi, alpha=I.warmup_alpha(st.nstep))
            E = st.E
            rec['e'] = (E[-1] - E[-2]) / 2 + e_off if len(E) >= 2 else None
            hist.append(rec)
            if len(hist) >= 3 and hist[-3]['e'] is not None:
                rec['de2'] = abs(rec['e'] - hist[-3]['e'])
                rec['dS2'] = abs(rec['S'] - hist[-3]['S'])
            if it + 1 >= nmin and len(hist) >= win + 2 and all('de2' in h and h['alpha'] == 0 for h in hist[-win:]):
                if max(h['de2'] for h in hist[-win:]) <= tolE and max(h['dS2'] for h in hist[-win:]) <= tolS:
                    conv = True
                    break
        xi = st.transfer_spectrum()
        last = hist[-win:]
        pt = dict(chi=chi, it=len(hist), conv=conv, e=hist[-1]['e'], S=(hist[-1]['S'] + hist[-1]['S_prev']) / 2,
                  S_t=hist[-1]['S'], S_prev=hist[-1]['S_prev'],
                  de2=max(h.get('de2', np.inf) for h in last), dS2=max(h.get('dS2', np.inf) for h in last),
                  res=max(h['res'] for h in last), r_pred=hist[-1]['r_pred'], nmv=np.mean([h['nmv'] for h in last]),
                  xi_e=xi['xi_even'], xi_o=xi['xi_odd'], lam0=xi['lam0'], disc=max(h['disc'] for h in last),
                  n=(hist[-1]['n1'] + hist[-1]['n2'] + hist[-2]['n1'] + hist[-2]['n2']) / 4, t=time.time() - t0,
                  D=st.D)
        pts.append(pt)
        if verbose:
            print(f"  chi={chi} D={st.D} it={pt['it']} conv={conv} e={pt['e']:.12f} de2={pt['de2']:.1e} "
                  f"dS2={pt['dS2']:.1e} res={pt['res']:.1e} rpred={pt['r_pred']:.1e} nmv={pt['nmv']:.1f} "
                  f"S={pt['S']:.5f} ({pt['S_t']:.5f},{pt['S_prev']:.5f}) xi_e={pt['xi_e']:.2f} xi_o={pt['xi_o']:.2f} "
                  f"lam0={pt['lam0']:.9f} n={pt['n']:.5f} disc={pt['disc']:.1e} t={pt['t']:.0f}s", flush=True)
    return pts, st


def part_tm():
    print('# part tm: Arnoldi against dense eig of the cell transfer matrix (xx1, chi = 24)')
    cfg = dict(model='xx1', mu=MU_C)
    pts, st = run_stages(cfg, [24], cap=300, verbose=False)
    Ms = st.cell_tensor()
    n0, D = st.n0_prev, st.D_prev
    for r in (0, 1):
        blocks = [(q, q ^ r) for q in (0, 1)]
        idx = []
        for a, b in blocks:
            ia = np.arange(D)[I.bsl(n0, D, a)]
            ib = np.arange(D)[I.bsl(n0, D, b)]
            idx += [(i, j) for i in ia for j in ib]
        n = len(idx)
        T = np.zeros((n, n), dtype=complex)
        for col, (i, j) in enumerate(idx):
            Xf = np.zeros((D, D), dtype=complex); Xf[i, j] = 1
            Y = sum(M.conj().T @ Xf @ M for M in Ms)
            T[:, col] = [Y[a, b] for a, b in idx]
        lam = np.linalg.eigvals(T)
        lam = lam[np.argsort(-np.abs(lam))][:6]

        def mv(v):
            Xf = np.zeros((D, D), dtype=complex)
            for c, (i, j) in enumerate(idx):
                Xf[i, j] = v[c]
            Y = sum(M.conj().T @ Xf @ M for M in Ms)
            return np.array([Y[a, b] for a, b in idx])
        ev, res, nmv, conv = I.arnoldi_top(mv, np.random.default_rng(5).standard_normal(n) + 0j, k=6, m=min(30, n - 1),
                                           tol=1e-12)
        print(f'sector {r}: n={n} dense |lam| {np.round(np.abs(lam), 10).tolist()}')
        print(f'          arnoldi |lam| {np.round(np.abs(ev[:6]), 10).tolist()} nmv={nmv} conv={conv} '
              f'max diff {np.max(np.abs(np.abs(ev[:6]) - np.abs(lam))):.1e}')


def part_ctrl(chis=(16, 32, 64)):
    print('# part ctrl: controls against exact energies')
    for cfg in [dict(model='xx1', mu=MU_C), dict(model='xx2', mu=MU_C), dict(model='tfim', h=1.05)]:
        ex = I.exact_energy(cfg)
        print(cfg, 'exact', ex)
        pts, st = run_stages(cfg, list(chis))
        for p in pts:
            print(f"    chi={p['chi']} e - e_exact = {p['e'] - ex['e']:.3e}")
        if len(pts) >= 2:
            S = np.array([p['S'] for p in pts]); lx = np.log([p['xi_o'] if p['xi_o'] > 0 else 1 for p in pts])
            lxm = np.log([max(p['xi_o'], p['xi_e']) for p in pts])
            print('    local c (6 dS/dln xi):', np.round(6 * np.diff(S) / np.diff(lxm), 3).tolist())


def part_cand(chis=(16, 32, 64)):
    print('# part cand: candidates at small chi')
    ref = {'p153c': (-650.0414853399, -1290.1134157022, -2570.2571184393, 0.2153),
           'sym40_n316': (-650.0460842395, -1290.1231205671, -2570.2770066290, 0.1917)}
    for key, model in (('p153c', 'p153'), ('sym40_n316', 'sym40')):
        cfg = dict(model=model, hz=HZ[key])
        r = ref[key]
        print(key, f'C190 bulk (E128-E64)/64 = {(r[1] - r[0]) / 64:.10f}, (E256-E128)/128 = {(r[2] - r[1]) / 128:.10f},'
              f' n_mid(256) = {r[3]}')
        pts, st = run_stages(cfg, list(chis), cap=3000)


if __name__ == '__main__':
    part = sys.argv[1] if len(sys.argv) > 1 else 'ed'
    if part == 'ed':
        part_ed()
    elif part == 'tm':
        part_tm()
    elif part == 'ctrl':
        part_ctrl(tuple(int(x) for x in sys.argv[2].split(',')) if len(sys.argv) > 2 else (16, 32, 64))
    elif part == 'cand':
        part_cand(tuple(int(x) for x in sys.argv[2].split(',')) if len(sys.argv) > 2 else (16, 32, 64))
