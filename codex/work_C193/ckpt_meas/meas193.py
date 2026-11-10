#!/usr/bin/env python3
"""C193b correlators on stored C190 matrix-product states (Claude Opus 5.5 agent, 2026-11-10).

One left-to-right pass over a right-canonical Z2 MPS (centre at site 0, the state a C190 checkpoint stores),
on numpy or cupy, complex128 throughout.  With E_k the left environment of the identity on bond k (sites
0..k-1 contracted), every right part is the identity because sites >= 1 are right-canonical, so

  <O_i ... O'_j> = close(propagate(E_i, M_i, O_i) ... , M_j, O'_j),

and the Schmidt weights on bond k are the eigenvalues of E_k.  Rows of many correlators are propagated
together as one batched stack.  Measured (operators in the frame of the Pauli Hamiltonian, i.e. rotated back by
the gauge phases, as in z2dmrg.measure):

  Gb_rows   <sm_i sp_j> (sm adds a particle, sp removes one; hard-core boson one-body function) for every
            i < j with i + j <= L - 1 + DSYM (all pairs symmetric about the centre and about centres up to
            DSYM/2 sites to the right of it, and everything left of them), plus every pair inside the stored
            48-site window (for the check against the stored Gb).
  Cnn       <n_i0 n_j> - <n_i0><n_j>, i0 in REFS = L/2-3..L/2, j = i0+1..L-1 (C176's geometry).
  Dpair     <(b_i0 b_i0+d)^dag (b_j b_j+d)> = <sm_i0 sm_i0+d sp_j sp_j+d>, d = 2, 3, i0 in REFS,
            j = i0+d+1..L-1-d (C176's geometry), with the one-point pieces <sp_j sp_j+d> (pp_d) and
            <sm_i0 sm_i0+d> (mm_d) for the connected part.
  S, Sbar1/6/12, dens, energy (through z2dmrg_xp.env_left with the run's MPO), right-canonical error.

Usage
  python3 meas193.py run <run.json> <ckpt.npz> <out.json> [--numpy]
  python3 meas193.py batch <list.txt> <outdir> [--numpy]      (one "<run.json> <ckpt.npz>" per line; skips
                                                               outputs that already have status done)
The module directory lib/ holds unchanged copies of z2dmrg.py, z2dmrg_xp.py and models.py from codex/work_C190.
"""
import os
import sys
import json
import time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, 'lib'))
import z2dmrg as Z            # noqa: E402
import models as Mo           # noqa: E402

DSYM = 9
CHUNK = int(os.environ.get('C193_CHUNK', '48'))
NOP = np.diag([0.0, 1.0])


def get_xp(use_numpy):
    if use_numpy:
        return np
    import cupy
    return cupy


def build(cfg):
    L = cfg['L']
    m = cfg['model']
    if m in ('xx1', 'xx2'):
        mpo, gph = Mo.xx_control(L, cfg['mu'], legs=1 if m == 'xx1' else 2)
        return mpo, gph
    mpo, gph, shift, info = Mo.candidate(m, L, cfg['hz'], gauge=cfg.get('gauge'))
    return mpo, gph


class Stack:
    """Batched row propagators on the current bond: E (n, D, D), with per-row metadata."""
    def __init__(self):
        self.E = []
        self.meta = []


def _opY(Y, O, xp):
    """(O . Y)[n, a, t, b] = sum_s O[t, s] Y[n, a, s, b] for a 2x2 host matrix O."""
    parts = []
    for t in (0, 1):
        acc = None
        for s in (0, 1):
            c = complex(O[t, s])
            if c == 0:
                continue
            term = c * Y[:, :, s, :]
            acc = term if acc is None else acc + term
        if acc is None:
            acc = xp.zeros_like(Y[:, :, 0, :])
        parts.append(acc[:, :, None, :])
    return xp.concatenate(parts, axis=2)


def tonp(a, xp):
    return a if xp is np else xp.asnumpy(a)


def measure_state(M, cfg, xp, gwin_lo=None, gwin=None, log=print):
    """M: list of host arrays (D_l, 2, D_r) of a right-canonical MPS with centre at site 0."""
    t0 = time.time()
    L = len(M)
    mpo, gph = build(cfg)
    sp = [gph[j] * Z.SP for j in range(L)]
    sm = [np.conj(gph[j]) * Z.SM for j in range(L)]
    I2 = np.eye(2)
    Md = [xp.asarray(np.ascontiguousarray(m), dtype=xp.complex128) for m in M]
    REFS = [L // 2 - 3, L // 2 - 2, L // 2 - 1, L // 2]
    DS = (2, 3)
    if gwin is None:
        gwin = min(48, L - 4)
        gwin_lo = L // 2 - gwin // 2
    win = range(gwin_lo, gwin_lo + gwin)

    def jmax_b(i):
        jm = min(L - 1, L - 1 + DSYM - i)
        if i in win:
            jm = max(jm, gwin_lo + gwin - 1)
        return jm

    # right-canonical check (sites >= 1)
    rc_err = 0.0
    for k in range(1, L):
        A = Md[k].reshape(Md[k].shape[0], -1)
        G = A @ A.conj().T
        rc_err = max(rc_err, float(xp.max(xp.abs(G - xp.eye(G.shape[0])))))

    Gb = {}            # (i, j) -> value
    Cnn = {i0: [] for i0 in REFS}
    Dp = {(d, i0): [] for d in DS for i0 in REFS}
    mm = {}            # (d, i0) -> <sm_i0 sm_i0+d>
    pp = {d: {} for d in DS}   # j -> <sp_j sp_j+d>
    dens = np.zeros(L)
    S = np.zeros(L - 1)
    Eid = xp.ones((1, 1), dtype=xp.complex128)
    rows = []          # list of dicts: E (D, D) on the current bond, kind, start, ops
    t_close = 0.0

    def R_pair(j, d):
        """R[a', a] on bond j: sites j..j+d with sp at j and j+d, closed by the right-canonical identity."""
        R = None
        for k in range(j + d, j - 1, -1):
            O = sp[k] if k in (j, j + d) else I2
            A = Md[k]
            OA = xp.asarray(O, dtype=xp.complex128)
            T = xp.tensordot(OA, A, axes=([1], [1]))              # t a b  (O[t,s] A[a,s,b])
            if R is None:
                X = T                                             # right identity
            else:
                X = xp.tensordot(T, R, axes=([2], [1]))           # t a b'
            # R_new[a', a] = sum conj(A[a', t, b']) X[t, a, b']
            R = xp.tensordot(A.conj(), X, axes=([1, 2], [0, 2]))  # a' a
        return R

    for k in range(L):
        A = Md[k]
        Dl, _, Dr = A.shape
        # Schmidt weights on bond k (k >= 1)
        if k >= 1:
            H = (Eid + Eid.conj().T) / 2
            p = tonp(xp.linalg.eigvalsh(H), xp)
            S[k - 1] = Z.entropy(np.clip(p, 0, None))
        # pair closes on bond k
        for d in DS:
            if k + d <= L - 1 and k >= min(REFS) + d + 1:
                tc = time.time()
                R = R_pair(k, d)
                pp[d][k] = complex(xp.sum(Eid * R))
                for r in rows:
                    if r['kind'] == 'pair' and r['d'] == d and k >= r['start'] + d + 1:
                        Dp[(d, r['start'])].append(complex(xp.sum(r['E'] * R)))
                t_close += time.time() - tc
        # new rows starting at site k
        new = []
        if k <= L - 1:
            if jmax_b(k) > k:
                new.append(dict(kind='b', start=k, ops={k: sm[k]}, close=sp, jmax=jmax_b(k)))
            if k in REFS:
                new.append(dict(kind='n', start=k, ops={k: NOP}, close=None, jmax=L - 1))
                for d in DS:
                    if k + d < L:
                        new.append(dict(kind='pair', start=k, d=d, ops={k: sm[k], k + d: sm[k + d]}, close=None,
                                        jmax=L - 1 - d))
        for r in new:
            r['E'] = Eid
        allrows = [dict(kind='id', E=Eid, ops={}, start=-1, jmax=L)] + rows + new
        # batched X = E @ A for all rows, in chunks
        Am = A.reshape(Dl, 2 * Dr)
        Ac = A.conj().reshape(Dl * 2, Dr)
        out_rows = []
        for c0 in range(0, len(allrows), CHUNK):
            chunk = allrows[c0:c0 + CHUNK]
            Es = xp.stack([r['E'] for r in chunk])                       # n Dl Dl
            X = (Es @ Am).reshape(len(chunk), Dl, 2, Dr)                  # n a' s b
            # closes (values at site k), batched per operator
            Acv = A.conj().reshape(-1)
            for kind, O in (('id', NOP), ('b', sp[k]), ('n', NOP)):
                idxs = [i for i, r in enumerate(chunk) if r['kind'] == kind and
                        (kind == 'id' or (k > r['start'] and k <= r['jmax']))]
                if not idxs:
                    continue
                vals = tonp(_opY(X[idxs], O, xp).reshape(len(idxs), -1) @ Acv, xp)
                for t, i in enumerate(idxs):
                    r = chunk[i]
                    if kind == 'id':
                        dens[k] = float(np.real(vals[t]))
                    elif kind == 'b':
                        Gb[(r['start'], k)] = complex(vals[t])
                    else:
                        Cnn[r['start']].append(complex(vals[t]))
            # propagate, grouped by operator
            groups = {}
            for idx, r in enumerate(chunk):
                if k >= r['jmax']:
                    continue                                              # row finished at this site
                O = r['ops'].get(k, I2)
                key = O.tobytes()
                groups.setdefault(key, (O, []))[1].append(idx)
            for key, (O, idxs) in groups.items():
                Xi = X[idxs] if len(idxs) < len(chunk) else X
                Y = _opY(Xi, O, xp).reshape(len(idxs), Dl * 2, Dr)
                En = xp.matmul(Ac.T[None, :, :], Y)                       # n Dr Dr  (b', b)
                for t, idx in enumerate(idxs):
                    r = chunk[idx]
                    r2 = dict(r)
                    r2['E'] = En[t]
                    out_rows.append(r2)
        Eid = None
        rows = []
        for r in out_rows:
            if r['kind'] == 'id':
                Eid = r['E']
            else:
                rows.append(r)
        # one-point pieces of the pair: after site start+d the row holds <sm sm> on its bond
        for r in rows:
            if r['kind'] == 'pair' and k == r['start'] + r['d']:
                mm[(r['d'], r['start'])] = complex(xp.trace(r['E']))
        if k % 32 == 0:
            log(f'   site {k}/{L} rows {len(rows)} t={time.time() - t0:.1f}s')
    norm = float(xp.real(xp.trace(Eid)))
    # energy through the block-sparse MPO environment of z2dmrg_xp (same code as the DMRG)
    import z2dmrg_xp as ZX
    n0 = cfg['_n0']
    W = [xp.asarray(w) for w in mpo.W]
    Ev = xp.zeros((1, mpo.w, 1), dtype=xp.complex128)
    Ev[0, 0, 0] = 1.0
    for j in range(L):
        Ev = ZX.env_left(Ev, Md[j], W[j], n0[j], n0[j + 1], mpo)
    E_mpo = float(xp.real(Ev[0, 2, 0]))
    out = dict(L=L, norm=norm, rc_err=rc_err, E_mpo=E_mpo, dens=dens.tolist(), S=S.tolist(),
               Sbar1=Z.sbar(S, L, 1), Sbar6=Z.sbar(S, L, 6), Sbar12=Z.sbar(S, L, 12), REFS=REFS, DSYM=DSYM,
               gwin=[gwin_lo, gwin])
    # Gb as rows: Gb_rows[i] = [re list, im list] for j = i+1 .. jmax(i)
    gr = {}
    for (i, j), v in Gb.items():
        gr.setdefault(i, {})[j] = v
    out['Gb_rows'] = {str(i): [[gr[i][j].real for j in sorted(gr[i])], [gr[i][j].imag for j in sorted(gr[i])]]
                      for i in sorted(gr)}
    out['Cnn'] = {str(i0): [float(np.real(v)) - dens[i0] * dens[i0 + 1 + t] for t, v in enumerate(Cnn[i0])]
                  for i0 in REFS}
    out['Dpair'] = {f'{d}_{i0}': [[v.real, v.imag] for v in Dp[(d, i0)]] for d in DS for i0 in REFS}
    out['mm'] = {f'{d}_{i0}': [mm[(d, i0)].real, mm[(d, i0)].imag] for (d, i0) in mm}
    out['pp'] = {str(d): {str(j): [v.real, v.imag] for j, v in pp[d].items()} for d in DS}
    out['t_meas'] = time.time() - t0
    out['t_pairclose'] = t_close
    return out


def gb_matrix(out, lo, W):
    """Hermitian one-body matrix on [lo, lo+W) from Gb_rows and dens (diagonal = n_i, as z2dmrg.measure)."""
    G = np.full((W, W), np.nan + 0j)
    for a in range(W):
        i = lo + a
        G[a, a] = out['dens'][i]
        row = out['Gb_rows'].get(str(i))
        if row is None:
            continue
        for t, (re, im) in enumerate(zip(*row)):
            j = i + 1 + t
            if j < lo + W:
                G[a, j - lo] = re + 1j * im
                G[j - lo, a] = re - 1j * im
    return G


def run_one(fjson, fck, fout, xp, log=print):
    d = json.load(open(fjson))
    assert d['status'] == 'done', fjson
    cfg = dict(d['cfg'])
    t0 = time.time()
    z = np.load(fck, allow_pickle=False)
    L = len([k for k in z.files if k.startswith('M')])
    M = [z[f'M{j}'] for j in range(L)]
    n0 = [int(x) for x in z['n0']]
    P = int(z['P'])
    meta = json.loads(str(z['meta'])) if 'meta' in z.files else None
    t_load = time.time() - t0
    assert L == cfg['L'] and P == cfg['P'], (L, P)
    psi = Z.MPS(M, n0, P)
    assert psi.check(1e-12), 'parity structure'
    cfg['_n0'] = n0
    m0 = d['meas']
    lo, W = m0['gwin']
    log(f'{os.path.basename(fjson)}: L={L} P={P} Dmax={max(psi.dims)} load {t_load:.0f}s meta={meta}')
    out = measure_state(M, cfg, xp, gwin_lo=lo, gwin=W, log=log)
    cfg.pop('_n0')
    E_phys = out['E_mpo'] + d['E_off']
    E_st = d.get('E_phys') if d.get('E_phys') is not None else d['hist'][-1]['E_phys']
    Gst = np.array(m0['Gb_re']) + 1j * np.array(m0['Gb_im'])
    Gme = gb_matrix(out, lo, W)
    chk = dict(E_phys=E_phys, E_stored=E_st, dE=E_phys - E_st, E_stored_is_last_sweep=d.get('E_phys') is None,
               Sbar6=out['Sbar6'], Sbar6_stored=d['Sbar6'], dSbar6=out['Sbar6'] - d['Sbar6'],
               dS_max=float(np.max(np.abs(np.array(out['S']) - np.array(d['S'])))),
               dGb_window=float(np.nanmax(np.abs(Gme - Gst))),
               ddens=float(np.max(np.abs(np.array(out['dens']) - np.array(m0['dens'])))),
               ckpt_sweeps=(meta or {}).get('sweeps'), hist_len=len(d['hist']))
    res = dict(status='done', tag=cfg['tag'], src=os.path.abspath(fjson), ckpt=os.path.abspath(fck), cfg=cfg,
               Dmax=max(psi.dims), check=chk, meas=out, t_load=t_load)
    log('   check ' + json.dumps({k: (float(f'{v:.3e}') if isinstance(v, float) else v) for k, v in chk.items()}))
    tmp = fout + '.tmp'
    json.dump(res, open(tmp, 'w'))
    os.replace(tmp, fout)
    return res


def main():
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    xp = get_xp('--numpy' in sys.argv)
    if args[0] == 'run':
        run_one(args[1], args[2], args[3], xp)
    elif args[0] == 'batch':
        outdir = args[2]
        os.makedirs(outdir, exist_ok=True)
        for line in open(args[1]):
            line = line.split('#')[0].strip()
            if not line:
                continue
            fj, fc = line.split()
            tag = os.path.basename(fj)[:-5]
            fout = os.path.join(outdir, tag + '_meas.json')
            if os.path.exists(fout) and json.load(open(fout)).get('status') == 'done':
                print('skip', tag, flush=True)
                continue
            if not os.path.exists(fc):
                print('MISSING checkpoint', tag, fc, flush=True)
                continue
            try:
                run_one(fj, fc, fout, xp, log=lambda s: print(s, flush=True))
                if xp is not np:
                    xp.get_default_memory_pool().free_all_blocks()
            except Exception as e:                     # keep going with the next state
                import traceback
                traceback.print_exc()
                print('FAILED', tag, repr(e), flush=True)


if __name__ == '__main__':
    main()
