#!/usr/bin/env python3
"""Ring spectra of the full microscopic chain in its period-3 sector, resolved by U(1) charge and momentum.

On a ring of L = 3M sites the letters' types are conserved, and in the period-3 type sector the chain is
exactly a ring of M cells (P B P) of dimension 32 with H = sum_j H0_j + lam sum_j V(pR_j, pL_j+1)
(-mu per B is a constant; it is dropped here and restored in the analysis).  Nothing is truncated.

Storage is block sparse in the cell charges t_j = q_j + 5/2 in {0..5} (block sizes 1,5,10,10,5,1),
at fixed total charge Q.  Momentum k = 2 pi m / M is imposed with the real projector
P = (1/M) sum_r cos(k r) T^r (it spans k and -k, which are degenerate by reflection), and the
eigenproblem is P (H - s) P with a shift s so that the projected-out space sits at zero, above
every physical level.

usage: ring_ed.py --M 4 --lam 0.1 --Q 0 --m 0 1 --nev 4 [--model micro|eff --Jxy --Jzz] [--out f.json]
model eff: the logical XXZ ring H = Jxy(s+s- + h.c.) + Jzz ZZ with the same machinery (cell dim 2).
Claude (Opus 5.5), 2026-10-27.
"""
import os, sys, json, time, argparse, itertools
for _k in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_k, '8')
import numpy as np
from scipy.sparse.linalg import eigsh, LinearOperator
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import cell as C


class CellModel:
    """Charge-graded single-cell data: H0 blocks and port operators as (t_in -> t_out) blocks."""

    def __init__(self, H0, qc, ports, bond_terms):
        self.ts = sorted(set(np.round(qc * 2).astype(int)))
        tval = np.round(qc * 2).astype(int)
        self.tmin = min(self.ts)
        self.idx = {t: np.where(tval == t)[0] for t in self.ts}
        self.dim = {t: len(self.idx[t]) for t in self.ts}
        self.H0 = {t: H0[np.ix_(self.idx[t], self.idx[t])] for t in self.ts}
        # bond terms: list of (coef, opR, opL, dtR, dtL) with ops given as full matrices on the cell
        self.terms = []
        for coef, A, B, dA, dB in bond_terms:
            Ab = {t: A[np.ix_(self.idx[t + dA], self.idx[t])] for t in self.ts if t + dA in self.idx}
            Bb = {t: B[np.ix_(self.idx[t + dB], self.idx[t])] for t in self.ts if t + dB in self.idx}
            self.terms.append((coef, Ab, Bb, dA, dB))


def micro_model(lam, D0, D1, kap, Jp, Jz):
    H0, ports, code, qc = C.cell_ops(D0, D1, kap)
    R, L = ports['R'], ports['L']
    # charge steps in units of 1/2 (t = 2q): s+ on a port raises q by 1 -> dt = 2
    terms = [(lam * Jp, R['sp'], L['sm'], 2, -2), (lam * Jp, R['sm'], L['sp'], -2, 2), (lam * Jz, R['z'], L['z'], 0, 0)]
    return CellModel(H0, qc, ports, terms), code, qc


def eff_model(Jxy, Jzz):
    H0 = np.zeros((2, 2))
    qc = np.array([0.5, -0.5])
    sp, sm, z = C.SP, C.SM, C.Z
    terms = [(Jxy, sp, sm, 2, -2), (Jxy, sm, sp, -2, 2), (Jzz, z, z, 0, 0)]
    return CellModel(H0, qc, None, terms), np.eye(2), qc


class Ring:
    def __init__(self, cm, M, Q2):
        """Q2 = 2Q (integer)."""
        self.cm, self.M = cm, M
        keys = [k for k in itertools.product(cm.ts, repeat=M) if sum(k) == Q2]
        self.keys = keys
        self.shapes = {k: tuple(cm.dim[t] for t in k) for k in keys}
        self.offs = {}
        o = 0
        for k in keys:
            self.offs[k] = o
            o += int(np.prod(self.shapes[k]))
        self.N = o

    def unpack(self, v):
        return {k: v[self.offs[k]:self.offs[k] + int(np.prod(self.shapes[k]))].reshape(self.shapes[k]) for k in self.keys}

    def pack(self, blocks):
        v = np.zeros(self.N)
        for k, b in blocks.items():
            v[self.offs[k]:self.offs[k] + b.size] = b.ravel()
        return v

    def apply_H(self, v):
        cm, M = self.cm, self.M
        X = self.unpack(v)
        out = np.zeros(self.N)
        Y = self.unpack(out)                       # views into out
        for k, x in X.items():
            y = Y[k]
            for j in range(M):
                h = cm.H0[k[j]]
                if not h.any():
                    continue
                y += np.moveaxis(np.tensordot(h, x, ([1], [j])), 0, j)
            for j in range(M):
                jn = (j + 1) % M
                for coef, Ab, Bb, dA, dB in cm.terms:
                    if k[j] not in Ab or k[jn] not in Bb:
                        continue
                    nk = list(k); nk[j] += dA; nk[jn] += dB; nk = tuple(nk)
                    if nk not in Y:
                        continue
                    t = np.moveaxis(np.tensordot(Ab[k[j]], x, ([1], [j])), 0, j)
                    t = np.moveaxis(np.tensordot(Bb[k[jn]], t, ([1], [jn])), 0, jn)
                    Y[nk] += coef * t
        return out

    def translate(self, v, r=1):
        """T^r: cell j -> cell j + r."""
        X = self.unpack(v)
        Y = {}
        for k, x in X.items():
            nk = tuple(np.roll(np.array(k), r))
            Y[nk] = np.moveaxis(x, list(range(self.M)), [(i + r) % self.M for i in range(self.M)])
        return self.pack(Y)

    def proj(self, v, m):
        M = self.M
        k = 2 * np.pi * m / M
        out = np.zeros_like(v)
        w = v.copy()
        for r in range(M):
            out += np.cos(k * r) * w
            w = self.translate(w, 1)
        fac = 1.0 if (m % M == 0 or 2 * m == M) else 2.0      # P_k + P_-k for k != 0, pi
        return fac * out / M


def lowest(ring, m, nev, shift, seed=0, tol=1e-10):
    rng = np.random.default_rng(seed)
    v0 = ring.proj(rng.normal(size=ring.N), m)

    def mv(v):                      # v stays in the range of P up to roundoff; one projection per step
        return ring.proj(ring.apply_H(v) - shift * v, m)
    op = LinearOperator((ring.N, ring.N), matvec=mv, dtype=float)
    w = eigsh(op, k=nev, which='SA', v0=v0, tol=tol, ncv=max(2 * nev + 8, 24), maxiter=20000, return_eigenvectors=False)
    return np.sort(w) + shift


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--model', default='micro')
    ap.add_argument('--M', type=int, default=4)
    ap.add_argument('--Q2', type=int, nargs='+', default=None, help='2Q values (default: lowest |Q| and next)')
    ap.add_argument('--m', type=int, nargs='+', default=None, help='momenta 2 pi m / M (default 0..M//2)')
    ap.add_argument('--nev', type=int, default=3)
    for k, v in dict(lam=0.1, kap=0.6, D0=1.0, D1=1.0, Jp=2.0, Jz=-1.0, Jxy=-1.0, Jzz=0.0).items():
        ap.add_argument('--' + k, type=float, default=v)
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    if a.model == 'micro':
        cm, code, qc = micro_model(a.lam, a.D0, a.D1, a.kap, a.Jp, a.Jz)
    else:
        cm, code, qc = eff_model(a.Jxy, a.Jzz)
    M = a.M
    Q2s = a.Q2 if a.Q2 is not None else ([0, 2] if M % 2 == 0 else [1, 3])
    ms = a.m if a.m is not None else list(range(M // 2 + 1))
    res = dict(args=vars(a), levels={})
    t0 = time.time()
    for Q2 in Q2s:
        ring = Ring(cm, M, Q2)
        for m in ms:
            w = lowest(ring, m, a.nev, shift=10.0)
            res['levels'][f'{Q2}_{m}'] = w.tolist()
            print(f'M={M} 2Q={Q2} m={m} dim={ring.N} levels={np.round(w, 10).tolist()} t={time.time() - t0:.0f}s', flush=True)
            if a.out:
                json.dump(res, open(a.out, 'w'), indent=1)


if __name__ == '__main__':
    main()
