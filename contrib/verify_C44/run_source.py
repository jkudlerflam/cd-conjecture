"""Driver for the rebuild (referee item ii): source formula (11), first-order coupling, checks.

Usage: python3 -B run_source.py [n ...]   (default 1 2 3)
Writes source_results_n<ns>.json and source_bonds_rebuilt_n<ns>.npz in this folder.
"""
import json
import sys
import time
import numpy as np
import scipy.linalg as sla
from core46 import Core, Params, SPECIES, port_terms
from wrap import Charge, compensator, plus, helix_track, site_op, WrappedCell
from precise import precise_bonds, kbar_bonds, helix_measure
from source_bond import run as run_plain, pair_sum, littleendian_perm


def tl(M):
    d = M.shape[0]
    return M - np.trace(M) / d * np.eye(d)


def nrm(M):
    return float(np.linalg.norm(M, 2))


def bare(n, P=None):
    """C39 (6), unwrapped, with the code in C44's little-endian ordering."""
    kb = kbar_bonds(n, None, dvals={0: [1.0] * (n + 1)}, P=P)
    return kb['00']


def decompose(B, Kb):
    E = {key: B[key] - Kb for key in B}
    d2 = Kb.shape[0]
    parts = dict(E0=E['11'], EL=E['21'] - E['11'], ER=E['12'] - E['11'],
                 ELR=E['22'] - E['21'] - E['12'] + E['11'])
    out = {}
    for name, M in parts.items():
        out[name + '_scalar'] = float(np.trace(M).real / d2)
        out[name + '_traceless'] = nrm(tl(M))
    out['max_norm_E'] = max(nrm(M) for M in E.values())
    out['payload_traceless_K11'] = nrm(tl(B['11']))
    return out


def first_order_wrapped(n, wvac, U=1.0, Jh=1.0):
    """<code pair| V |code pair> on the explicit wrapped pair codes, all four lengths."""
    P = Params(n)
    core = Core(P)
    G = core.code()
    k = n + 2
    b = k + 1
    phi = 1 / b
    ch = {p: Charge(k, p, U, wvac) for p in (1, 2)}
    M_L, _ = compensator(b, phi, ch[1].r, ch[2].r)
    worst = 0.0
    for p in (1, 2):
        for q in (1, 2):
            # code vectors: chi x G x helix ground (one relevant track per cell)
            NL, NR = k + p, k + q
            gA = plus(0)
            for x in range(1, NL):
                gA = np.kron(gA, plus(x * phi))
            gB = plus(-(NR - 1) * phi)
            for x in range(1, NR):
                gB = np.kron(gB, plus((x - NR + 1) * phi))
            hL = site_op(NL, NL - 1, M_L.conj())
            hR = site_op(NR, 0, M_L)
            for r in (1, 2):
                bb, terms = port_terms(P, r)
                for sl, sr, c in terms:
                    FL = core.source('B', r, sl).toarray()
                    FR = core.source('A', r, sr).toarray()
                    # compression factorizes: charge x core x helix, per cell
                    cl = (ch[p].chi @ (ch[p].OR * ch[p].chi)) * (G.conj().T @ FL @ G) * (gA.conj() @ hL @ gA)
                    cr = (ch[q].chi @ (ch[q].OL * ch[q].chi)) * (G.conj().T @ FR @ G) * (gB.conj() @ hR @ gB)
                    worst = max(worst, float(np.abs(np.kron(cl, cr)).max()))
    return worst


def filter_factorization(n, wvac, ratio=0.03, U=1.0):
    """||H_full (I x g_f) - (I x g_f) H_reduced|| with the occupation filter included."""
    P = Params(n)
    core = Core(P)
    k = n + 2
    phif = 2 * np.pi / (k + 2)
    kap = 1 / ((k + 2) * (2 * k + 3))
    bet = np.sin(np.pi / (k + 2)) ** 2
    eps = kap * bet / ((1 + np.sqrt(kap)) ** 2 + bet)
    Jf = 2 * U / eps
    sv = np.zeros(4, complex)
    sv[1], sv[2] = 1 / np.sqrt(2), -np.exp(-1j * phif) / np.sqrt(2)
    hf = np.outer(sv, sv.conj())
    af = np.eye(2) - np.outer(plus(0), plus(0).conj())
    bv = plus((k - 1) * phif)
    bf = np.eye(2) - np.outer(bv, bv.conj())
    gf = plus(0)
    for j in range(1, k):
        gf = np.kron(gf, plus(j * phif))

    def lab(j, o2, width=1):
        return np.kron(np.kron(np.eye(2 ** j), o2), np.eye(2 ** (k - j - width)))
    worst = 0.0
    for p in (1, 2):
        ch = Charge(k, p, U, wvac)
        s = ratio * ch.gap
        for r in (1, 2):
            wc = WrappedCell(core, ch, s, r)
            m = len(ch.configs)
            dc = len(wc.idx)
            Hf = np.kron(wc.H, np.eye(2 ** k))
            Hf = Hf + Jf * np.kron(np.eye(m * dc), lab(0, af) + lab(k - 1, bf))
            for j in range(k - 1):
                Hf = Hf + Jf * np.kron(np.kron(np.diag(ch.Dj[:, j]), np.eye(dc)), lab(j, hf, 2))
            emb = np.kron(np.eye(m * dc), gf[:, None])
            worst = max(worst, float(np.abs(Hf @ emb - emb @ wc.H).max()))
    return worst, Jf


def bus_superselection(n):
    core = Core(Params(n))
    H = core.H().toarray()
    worst = 0.0
    for bus in (0, 1, 2):
        Pb = np.zeros(core.dim)
        Pb[core.sector(bus)] = 1
        worst = max(worst, float(np.abs(H * Pb[None, :] - Pb[:, None] * H).max()))
    return worst


def sylvester_check(n=1, wvac=0.16, ratio=0.03, U=1.0, Jh=1.0):
    """Independent linear algebra: solve H_L X + X H_R^T = B (Bartels-Stewart)
    with the FULL 2^N helix track, length pair (1,1)."""
    info, out, cells, hel, ch, core, G, s = run_plain(n, wvac, ratio, Jh=Jh)
    P = Params(n)
    k = n + 2
    b = k + 1
    phi = 1 / b
    d = 2 ** n
    M_L, _ = compensator(b, phi, ch[1].r, ch[2].r)
    res = {}
    perm = littleendian_perm(n)
    for p, q in ((1, 1),):
        NL, NR = k + p, k + q
        HhL = Jh * helix_track(NL, phi, 0)
        HhR = Jh * helix_track(NR, phi, NR - 1)
        gA = plus(0)
        for x in range(1, NL):
            gA = np.kron(gA, plus(x * phi))
        gB = plus(-(NR - 1) * phi)
        for x in range(1, NR):
            gB = np.kron(gB, plus((x - NR + 1) * phi))
        hl = site_op(NL, NL - 1, M_L.conj()) @ gA
        hr = site_op(NR, 0, M_L) @ gB
        K = np.zeros((d * d, d * d), complex)
        for r in (1, 2):
            bb, terms = port_terms(P, r)
            wl = WrappedCell(core, ch[p], s, r)
            wr = WrappedCell(core, ch[q], s, r)
            HL = np.kron(wl.H, np.eye(2 ** NL)) + np.kron(np.eye(wl.dim), HhL)
            HR = np.kron(wr.H, np.eye(2 ** NR)) + np.kron(np.eye(wr.dim), HhR)
            sL = wl.source_vectors('B', G)
            sR = wr.source_vectors('A', G)
            vecs = {}
            for i in range(d):
                for j in range(d):
                    Bm = 0
                    for sl, sr, c in terms:
                        Bm = Bm + c * np.outer(np.kron(sL[sl][:, i], hl), np.kron(sR[sr][:, j], hr))
                    vecs[(i, j)] = Bm
            for (i, j), Bm in vecs.items():
                X = sla.solve_sylvester(HL, HR.T, Bm)
                for (i2, j2), Bm2 in vecs.items():
                    K[i2 * d + j2, i * d + j] += -bb * bb * np.vdot(Bm2, X)
        Pm = np.zeros((d * d, d * d))
        for i in range(d):
            for j in range(d):
                Pm[perm[i] * d + perm[j], i * d + j] = 1
        res['%d%d' % (p, q)] = Pm @ (s * K) @ Pm.T
    return res


if __name__ == '__main__':
    ns = [int(a) for a in sys.argv[1:]] or [1, 2, 3]
    ref = np.load('../../codex/work_C44/source_bonds.npz')
    results = dict(cases=[], checks={})
    saved = {}
    t0 = time.time()
    for n in ns:
        Kb = bare(n)
        saved['bare_n%d' % n] = Kb
        if 'bare_n%d' % n in ref.files:
            results['checks']['bare_n%d_vs_C44' % n] = float(np.abs(Kb - ref['bare_n%d' % n]).max())
        results['checks']['bare_n%d_scalar' % n] = float(np.trace(Kb).real / Kb.shape[0])
        results['checks']['bare_n%d_traceless' % n] = nrm(tl(Kb))
        results['checks']['bus_superselection_n%d' % n] = bus_superselection(n)
        for wvac in (0.16, 0.1):
            kb = kbar_bonds(n, wvac)
            for ratio in (0.03, 0.003):
                t1 = time.time()
                B, ex = precise_bonds(n, wvac, ratio)
                case = dict(n=n, wvac=wvac, ratio=ratio, s=ex['s'],
                            charge_gaps=[ex['charge'][p].gap for p in (1, 2)],
                            endpoint_weights=[ex['charge'][p].r for p in (1, 2)],
                            adjacency=[list(ex['charge'][p].adj) for p in (1, 2)],
                            compensator_norm=nrm(ex['M_L']),
                            first_order_wrapped=first_order_wrapped(n, wvac))
                case.update(decompose(B, Kb))
                case['Kbar_minus_bare_max'] = max(nrm(kb[key] - Kb) for key in kb)
                case['wrapped_minus_Kbar_max'] = max(nrm(B[key] - kb[key]) for key in B)
                if n <= 2:
                    _, plain, *_ = run_plain(n, wvac, ratio)
                    case['plain_eigh_minus_precise'] = max(float(np.abs(plain[key] - B[key]).max()) for key in B)
                    names = ['n%d_w%s_t%s_%s' % (n, wvac, ratio, key) for key in B]
                    case['C44_minus_precise'] = max(float(np.abs(ref[nm] - B[nm[-2:]]).max()) for nm in names)
                    R = {nm[-2:]: ref[nm] for nm in names}
                    case['C44_decomposition'] = decompose(R, ref['bare_n%d' % n])
                case['seconds'] = time.time() - t1
                results['cases'].append(case)
                for key in B:
                    saved['n%d_w%s_t%s_%s' % (n, wvac, ratio, key)] = B[key]
                print(json.dumps({k_: case[k_] for k_ in case if k_ not in ('adjacency', 'C44_decomposition')}, default=float))
    if 1 in ns:
        fac, Jf = filter_factorization(1, 0.16)
        results['checks']['filter_factorization_n1'] = fac
        results['checks']['filter_Jf_n1'] = Jf
        syl = sylvester_check()
        Bp, _ = precise_bonds(1, 0.16, 0.03)
        results['checks']['sylvester_full_helix_minus_precise'] = max(float(np.abs(syl[key] - Bp[key]).max()) for key in syl)
    results['seconds'] = time.time() - t0
    print(json.dumps(results['checks'], default=float))
    tag = '_'.join(map(str, ns))
    json.dump(results, open('source_results_n%s.json' % tag, 'w'), indent=1, default=float)
    np.savez('source_bonds_rebuilt_n%s.npz' % tag, **saved)
