"""Wrapped two-cell second-order bond, rebuilt from the rules.

Route 1 (direct): exact eigenbases of the charge x core block and of the
FULL 2^N helix track of each cell, with the source vectors built from the
physical maps; no use of the reached-space or omega formulas of C44 (10).
Route 2 (literal C44 (11)): the one-flip T_N measure and omega weights.
Both are compared with each other and with codex/work_C44/source_bonds.npz.
"""
import sys
import json
import time
import numpy as np
from core46 import Core, Params, SPECIES, port_terms
from wrap import Charge, helix_track, site_op, compensator, plus, WrappedCell


def littleendian_perm(n):
    d = 2 ** n
    return [int(format(i, '0%db' % n)[::-1], 2) if n else 0 for i in range(d)]


def pair_sum(AL, AR, EL, ER, hL, eL, hR, eR, terms, b, d, chunk=512):
    """-b^2 sum_{alpha,beta,l,m} |hL_l|^2|hR_m|^2 Z^dag Z/(E+E+e+e)."""
    wl, wr = np.abs(hL) ** 2, np.abs(hR) ** 2
    W = np.zeros((len(EL), len(ER)))
    for a0 in range(0, len(EL), chunk):
        Ea = EL[a0:a0 + chunk]
        acc = np.zeros((len(Ea), len(ER)))
        for l, wlv in enumerate(wl):
            for m, wrv in enumerate(wr):
                acc += wlv * wrv / (Ea[:, None] + ER[None, :] + eL[l] + eR[m])
        W[a0:a0 + chunk] = acc
    K = np.zeros((d, d, d, d), complex)    # [i, j, i', j']
    for sl, sr, c in terms:
        for sl2, sr2, c2 in terms:
            X = np.einsum('ai,ak->aik', AL[sl].conj(), AL[sl2]).reshape(len(EL), d * d)
            Y = np.einsum('bj,bm->bjm', AR[sr].conj(), AR[sr2]).reshape(len(ER), d * d)
            M = X.T @ (W @ Y)                      # [(i i'), (j j')]
            K += np.conj(c) * c2 * M.reshape(d, d, d, d).transpose(0, 2, 1, 3)
    return -b * b * K.reshape(d * d, d * d)


def run(n, wvac, ratio, Jh=1.0, U=1.0, comp='small', verbose=True, cutoff=1e-13):
    P = Params(n)
    core = Core(P)
    G = core.code()
    k = n + 2
    bsel = k + 1
    phi = 1.0 / bsel
    d = 2 ** n
    ch = {p: Charge(k, p, U, wvac) for p in (1, 2)}
    delta = min(ch[1].gap, ch[2].gap, Jh / ((bsel + 1) * (2 * bsel + 1)))
    s = ratio * delta
    M_L, (u1, u2) = compensator(bsel, phi, ch[1].r, ch[2].r, comp)
    M_R = M_L.conj()
    info = dict(n=n, wvac=wvac, ratio=ratio, s=s, Jh=Jh,
                charge_gaps=[ch[1].gap, ch[2].gap],
                endpoint_weights=[ch[1].r, ch[2].r],
                endpoint_weights_right=[ch[1].rR, ch[2].rR],
                adjacency=[list(ch[1].adj), list(ch[2].adj)],
                compensator_norm=float(np.linalg.norm(M_L, 2)),
                compensator_expectations=[float((u1.conj() @ M_L @ u1).real),
                                          float((u2.conj() @ M_L @ u2).real)])
    # helix tracks: left cell uses track A (pin site 0, source site N-1);
    # right cell uses track B (pin site N-1, source site 0).
    hel = {}
    for p in (1, 2):
        N = k + p
        for side, pin, src, Mop in (('L', 0, N - 1, M_R), ('R', N - 1, 0, M_L)):
            Hh = Jh * helix_track(N, phi, pin)
            e, V = np.linalg.eigh(Hh)
            g = V[:, 0]
            # check ground state is the product helix of C8 (3)
            if side == 'L':
                prod = plus(0)
                for x in range(1, N):
                    prod = np.kron(prod, plus(x * phi))
            else:
                prod = plus(-(N - 1) * phi)
                for x in range(1, N):
                    prod = np.kron(prod, plus((x - N + 1) * phi))
            ov = abs(np.vdot(prod, g))
            vec = site_op(N, src, Mop) @ prod
            amp = V.conj().T @ vec
            keep = np.abs(amp) > cutoff
            hel[(p, side)] = dict(e=e[keep], amp=amp[keep], ground=e[0], gap=e[1] - e[0],
                                  product_overlap=ov, kept=int(keep.sum()),
                                  dropped_max=float(np.abs(amp[~keep]).max(initial=0)))
    info['helix'] = {'%d%s' % key: dict(ground=float(v['ground']), gap=float(v['gap']),
                                        product_overlap=float(v['product_overlap']),
                                        kept=v['kept'], dropped_max=v['dropped_max'])
                     for key, v in hel.items()}
    # charge x core blocks
    cells = {}
    for p in (1, 2):
        for r in (1, 2):
            wc = WrappedCell(core, ch[p], s, r)
            E, Uv = np.linalg.eigh(wc.H)
            cells[(p, r)] = dict(E=E, U=Uv, srcA=wc.source_vectors('A', G),
                                 srcB=wc.source_vectors('B', G), dim=wc.dim)
    info['source_dims'] = [cells[(p, 1)]['dim'] for p in (1, 2)]
    # first-order compression: P V P on the wrapped pair code
    # (every port term maps a vacuum token to an excited one, so check directly)
    pvp = 0.0
    for p in (1, 2):
        chi = ch[p].chi
        for r in (1, 2):
            for sg in SPECIES:
                for end in ('A', 'B'):
                    F = core.source(end, r, sg).toarray()
                    pvp = max(pvp, np.abs(G.conj().T @ F @ G).max())
    info['first_order_single_cell_max'] = float(pvp)
    out = {}
    perm = littleendian_perm(n)
    for p in (1, 2):
        for q in (1, 2):
            K = np.zeros((d * d, d * d), complex)
            for r in (1, 2):
                b, terms = port_terms(P, r)
                cl, cr = cells[(p, r)], cells[(q, r)]
                AL = {sg: cl['U'].conj().T @ cl['srcB'][sg] for sg in SPECIES}
                AR = {sg: cr['U'].conj().T @ cr['srcA'][sg] for sg in SPECIES}
                hl, hr = hel[(p, 'L')], hel[(q, 'R')]
                K += pair_sum(AL, AR, cl['E'], cr['E'], hl['amp'], hl['e'], hr['amp'], hr['e'], terms, b, d)
            sK = s * K
            # to C44's little-endian data ordering
            Pm = np.zeros((d * d, d * d))
            for i in range(d):
                for j in range(d):
                    Pm[perm[i] * d + perm[j], i * d + j] = 1
            out['%d%d' % (p, q)] = Pm @ sK @ Pm.T
    return info, out, cells, hel, ch, core, G, s


if __name__ == '__main__':
    n = int(sys.argv[1])
    t0 = time.time()
    ref = np.load('../../codex/work_C44/source_bonds.npz')
    res = []
    for wvac in (0.16, 0.1):
        for ratio in (0.03, 0.003):
            info, out, *_ = run(n, wvac, ratio)
            cmp = {}
            for key, M in out.items():
                name = 'n%d_w%s_t%s_%s' % (n, wvac, ratio, key)
                if name in ref.files:
                    R = ref[name]
                    cmp[key] = dict(max_abs_diff=float(np.abs(M - R).max()),
                                    scale=float(np.abs(R).max()))
            info['compare_C44'] = cmp
            res.append(info)
            print(json.dumps({k: info[k] for k in ('wvac', 'ratio', 's', 'charge_gaps', 'endpoint_weights', 'compensator_norm', 'compare_C44')}, default=float))
    print('seconds', time.time() - t0)
