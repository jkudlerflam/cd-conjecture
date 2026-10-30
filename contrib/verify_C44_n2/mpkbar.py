"""Extended-precision s/delta -> 0 bond (C44 'Kbar') for the n = 2 member.

Same formula as verify_C44/precise.py:kbar_bonds, evaluated in mpmath.
The core matrices, code vectors and charge adjacencies are the double
precision numbers of the rebuild; they are converted exactly and every
subsequent operation (the sums over bonds, the sector eigendecompositions,
the pair sums and the scalar subtraction) is done at mp.dps digits. This
defines the exact second-order bond of the model whose parameters are those
doubles, so the scalar part (about 3e-2) can be removed without cancellation
against payload hoppings of order 1e-15.

python3 -B mpkbar.py RHO G W A ETA V D [DPS]
  -> kbar_n2_rho{RHO}_g{G}_w{W}_a{A}_eta{ETA}_v{V}_D{D}.npz
(q, K at the Params defaults of verify_C44/core46.py)
The file stores, for each length pair pq, the traceless part of K_pq in
units of 1e-12 (complex128, 16x16), the scalar parts as strings, and the
double-precision traceless parts for comparison.
"""
import sys
import time
sys.path.insert(0, '../verify_C44')
import numpy as np
import mpmath as mp
from core46 import Core, Params, SPECIES, port_terms
from wrap import Charge
from source_bond import littleendian_perm
from precise import kbar_bonds


def tomp(A):
    A = np.asarray(A.toarray() if hasattr(A, 'toarray') else A, dtype=complex)
    return mp.matrix([[mp.mpc(float(x.real), float(x.imag)) for x in row] for row in A])


def mp_kbar(n, wvac, P, U=1.0):
    core = Core(P)
    G = tomp(core.code())
    k = n + 2
    d = 2 ** n
    ch = {p: Charge(k, p, U, wvac) for p in (1, 2)}
    dvals = {p: ch[p].adj for p in (1, 2)}
    ends = tomp(core.aL) + tomp(core.bR)
    bonds = [tomp(b) for b in core.bonds]
    src = {(e, r, sg): tomp(core.source(e, r, sg)) for e in 'AB' for r in (1, 2) for sg in SPECIES}
    perm = littleendian_perm(n)
    eig = {}
    for p in (1, 2):
        Hc = ends.copy()
        for j, bj in enumerate(bonds):
            Hc += mp.mpf(float(dvals[p][j])) * bj
        for r in (1, 2):
            idx = core.sector(r)
            Hr = mp.matrix(len(idx), len(idx))
            for a, ia in enumerate(idx):
                for b, ib in enumerate(idx):
                    Hr[a, b] = Hc[ia, ib]
            E, Uv = mp.eighe(Hr)
            UH = Uv.transpose_conj()

            def amp(end, sg):
                F = src[(end, r, sg)] * G
                sub = mp.matrix(len(idx), d)
                for a, ia in enumerate(idx):
                    for c in range(d):
                        sub[a, c] = F[ia, c]
                return UH * sub
            AL = {sg: amp('B', sg) for sg in SPECIES}
            AR = {sg: amp('A', sg) for sg in SPECIES}
            eig[(p, r)] = (E, AL, AR)
    out = {}
    for p in (1, 2):
        for q in (1, 2):
            K = [[mp.mpc(0) for _ in range(d * d)] for _ in range(d * d)]
            for r in (1, 2):
                b2 = mp.mpf(1) if r == 1 else mp.mpf(float(P.eta))
                ph = {1: mp.pi / 2, 2: mp.pi}[r]
                terms = [('0', '0', mp.mpc(1)), ('+', '-', mp.expj(ph)), ('-', '+', mp.expj(-ph))]
                EL, AL, _ = eig[(p, r)]
                ER, _, AR = eig[(q, r)]
                nL, nR = len(EL), len(ER)
                W = mp.matrix(nL, nR)
                for a in range(nL):
                    for b in range(nR):
                        W[a, b] = 1 / (EL[a] + ER[b])
                for sl, sr, c in terms:
                    for sl2, sr2, c2 in terms:
                        X = mp.matrix(d * d, nL)     # [(i i'), a] = conj(AL[sl][a,i]) AL[sl2][a,i']
                        for i in range(d):
                            for i2 in range(d):
                                for a in range(nL):
                                    X[i * d + i2, a] = mp.conj(AL[sl][a, i]) * AL[sl2][a, i2]
                        Y = mp.matrix(nR, d * d)
                        for b in range(nR):
                            for j in range(d):
                                for j2 in range(d):
                                    Y[b, j * d + j2] = mp.conj(AR[sr][b, j]) * AR[sr2][b, j2]
                        M = X * (W * Y)
                        f = -b2 * mp.conj(c) * c2
                        for i in range(d):
                            for i2 in range(d):
                                for j in range(d):
                                    for j2 in range(d):
                                        # K[(i j),(i' j')] in big-endian data order, then little-endian relabel
                                        row = perm[i] * d + perm[j]
                                        col = perm[i2] * d + perm[j2]
                                        K[row][col] += f * M[i * d + i2, j * d + j2]
            out['%d%d' % (p, q)] = K
    return out


if __name__ == '__main__':
    rho, g, w, a, eta, v, D = (float(x) for x in sys.argv[1:8])
    mp.mp.dps = int(sys.argv[8]) if len(sys.argv) > 8 else 34
    n = 2
    P = Params(n, g=g, w=w, a=a, eta=eta, v=v, D=D)
    t0 = time.time()
    K = mp_kbar(n, rho, P)
    Kd = kbar_bonds(n, rho, P=P)
    D = 4 ** n
    save = {}
    for key, Km in K.items():
        Km = [[(Km[i][j] + mp.conj(Km[j][i])) / 2 for j in range(D)] for i in range(D)]
        sc = mp.re(sum(Km[i][i] for i in range(D)) / D)
        tl = np.array([[complex((Km[i][j] - (sc if i == j else 0)) * mp.mpf('1e12')) for j in range(D)] for i in range(D)])
        save['K' + key] = tl
        save['scalar' + key] = mp.nstr(sc, mp.mp.dps)
        kd = (Kd[key] + Kd[key].conj().T) / 2
        save['Kdouble' + key] = (kd - np.trace(kd).real / D * np.eye(D)) * 1e12
        print(key, 'scalar', mp.nstr(sc, 25), 'max|mp - double| traceless (units 1e-12)',
              float(np.abs(tl - save['Kdouble' + key]).max()), flush=True)
    s = {k: mp.mpf(save['scalar' + k]) for k in ('11', '12', '21', '22')}
    V = s['22'] - s['21'] - s['12'] + s['11']
    save['V'] = mp.nstr(V, 20)
    save['lin'] = [mp.nstr(s['21'] - s['11'], 20), mp.nstr(s['12'] - s['11'], 20)]
    save['params'] = np.array([rho, g, w, a, eta, v, D, mp.mp.dps])
    print('V', mp.nstr(V, 15), 'seconds', time.time() - t0)
    np.savez('kbar_n2_rho%s_g%s_w%s_a%s_eta%s_v%s_D%s.npz' % tuple(sys.argv[1:8]), **save)
