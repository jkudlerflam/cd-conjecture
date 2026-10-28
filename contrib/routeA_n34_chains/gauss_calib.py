#!/usr/bin/env python3
"""Exact Gaussian references for the free parts.
(1) c of the infinite free ground state from block entanglement S(l) of l cells (velocity
    independent): slope of S against (1/3) ln l over l = 64..512.
(2) Rings of L cells with the spin-chain twist pi (N - 1): grand-canonical ground state (best N),
    its energy E0(L), block entropies S(l) and the fit c_S(L); Casimir estimates from pairs of
    sizes with the Hartree-Fock velocities, c_Cas(L1, L2) = (12/pi)[E0(L1)/L1 - E0(L2)/L2] /
    [(1/L2^2 - 1/L1^2) sum_p |v_p|].
usage: gauss_calib.py OUT.json SPEC6 [SPEC6 ...]   (true Hartree-Fock chain and proxy free chain)
Written by Claude (Opus 5.5) for contrib/routeA_n34_chains, 2026-10-28."""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import models34 as Mo
import analyze34 as A


def ring_free(T, M, L, const):
    """Grand-canonical ground state of the free ring (spin-chain twist pi (N-1)).
    Returns best N, E0 (with the per-bond constant), S(l) list, and E0 per N."""
    n = len(T)
    res = {}
    for th in (0.0, np.pi):
        ks = (2 * np.pi * np.arange(L) + th) / L
        eps, vec = [], []
        for k in ks:
            w, U = np.linalg.eigh(A.hk(T, M, k))
            for i in range(n):
                eps.append(w[i]); vec.append((k, U[:, i]))
        o = np.argsort(eps)
        res[th] = ([eps[i] for i in o], [vec[i] for i in o])
    best = None
    EN = {}
    for N in range(0, n * L + 1):
        th = 0.0 if (N - 1) % 2 == 0 else np.pi
        eps, vec = res[th]
        E = float(np.sum(eps[:N])) + L * const
        EN[N] = E
        if best is None or E < best[1] - 1e-12:
            best = (N, E)
    N, E = best
    th = 0.0 if (N - 1) % 2 == 0 else np.pi
    eps, vec = res[th]
    # real-space orbitals phi_{(j,a)} = e^{i k j} u_a / sqrt(L)
    Phi = np.zeros((L * n, N), complex)
    for q in range(N):
        k, u = vec[q]
        for j in range(L):
            Phi[j * n:(j + 1) * n, q] = np.exp(1j * k * j) * u / np.sqrt(L)
    C = Phi.conj() @ Phi.T          # <c^+_x c_y>
    S = []
    for l in range(1, L):
        x = np.linalg.eigvalsh(C[:l * n, :l * n])
        x = x[(x > 1e-13) & (x < 1 - 1e-13)]
        S.append(float(-(x * np.log(x) + (1 - x) * np.log(1 - x)).sum()))
    return N, E, S, EN


def fit_cS(S, L):
    l = np.arange(1, L)
    X = np.log((L / np.pi) * np.sin(np.pi * l / L)) / 3
    Am = np.vstack([X, np.ones_like(X)]).T
    c, s0 = np.linalg.lstsq(Am, np.array(S), rcond=None)[0]
    return float(c)


def calib(T, M, n, label, Ls, const=0.0):
    fp = A.fermi(T, M)
    vsum = float(sum(abs(f['v']) for f in fp))
    ells = [32, 64, 128, 256, 512]
    Sl = Mo.gaussian_block_entropy(T, M, ells)
    cs = [3 * (Sl[i + 1] - Sl[i]) / np.log(ells[i + 1] / ells[i]) for i in range(len(ells) - 1)]
    rings = []
    for L in Ls:
        N, E, S, EN = ring_free(T, M, L, const)
        rings.append(dict(L=L, N=N, E0=E, S=S, cS=fit_cS(S, L), EN={str(k): v for k, v in EN.items() if abs(k - N) <= 2}))
    cas = []
    for i in range(len(rings) - 1):
        r1, r2 = rings[i], rings[i + 1]
        L1, L2 = r1['L'], r2['L']
        val = (12 / np.pi) * (r1['E0'] / L1 - r2['E0'] / L2) / ((1 / L2 ** 2 - 1 / L1 ** 2) * vsum)
        cas.append(dict(L1=L1, L2=L2, cCas=float(val)))
    return dict(label=label, n_fermi_points=len(fp), v_sum=vsum, ells=ells, S_block=Sl, c_block_slopes=cs,
                rings=rings, casimir_pairs=cas)


def main():
    out = []
    for spec6 in sys.argv[2:]:
        K, Ws, Th, Mh, info = Mo.base(spec6)
        n = info['n']
        Kp, Tp, Mp, pinfo = Mo.proxy(spec6)
        rec = dict(spec=spec6)
        for lab, (T, M) in (('true_hf_free', (Th, Mh)), ('proxy_free', (Tp, Mp))):
            Kq = Mo.quad_bond(T, M, n)
            const = -float(np.trace(Kq).real) / len(Kq)
            Ls = list(range(4, 9)) + [12, 16, 24, 32, 48, 64, 96, 128]
            r = calib(T, M, n, lab, Ls, const)
            rec[lab] = r
            print(spec6, lab, 'c(block, slopes over l=32..512):', np.round(r['c_block_slopes'], 3),
                  ' c_S(L):', [(x['L'], round(x['cS'], 3)) for x in r['rings']],
                  ' c_Cas:', [(x['L1'], round(x['cCas'], 2)) for x in r['casimir_pairs']], flush=True)
        # matched XX on the proxy
        Kx, bands = Mo.matched_xx(Tp, Mp, n)
        Tx = np.zeros((n, n), complex); Mx = np.zeros((n, n), complex)
        for b in bands:
            if b.get('gapped'):
                Mx[b['band'], b['band']] = -1.0 if b['filled'] else 1.0
                continue
            t = b['v'] / (2 * np.sin(b['kF']))
            Tx[b['band'], b['band']] = -t * np.exp(-1j * b['centre']); Mx[b['band'], b['band']] = 2 * t * np.cos(b['kF'])
        Kq = Mo.quad_bond(Tx, Mx, n)
        r = calib(Tx, Mx, n, 'proxy_xx', list(range(4, 9)), -float(np.trace(Kq).real) / len(Kq))
        rec['proxy_xx'] = r
        print(spec6, 'proxy_xx c_S(L):', [(x['L'], round(x['cS'], 3)) for x in r['rings']], flush=True)
        out.append(rec)
        json.dump(out, open(sys.argv[1], 'w'), indent=1)


if __name__ == '__main__':
    main()
