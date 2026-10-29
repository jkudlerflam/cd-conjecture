"""Rayleigh-Schroedinger coefficients of the slow-vacancy fluid (t -> 0).

H(t) = H0 + t H1, H0 = the d=12 bond at t=0 (vacancies frozen; material clock,
contact w), H1 = vacancy hopping.  Zero port, K = infinity, bits zero
(exact for all bits by bit blindness).  Everything is done by exact linear
algebra on the frozen-gap segment, so the t^2 and t^4 coefficients are
obtained without taking small differences of energies.

  segment: R_1 E [nL vacancies] cell(p) [nR vacancies] R_0, gaps frozen
     E(nL,nR) = E0 + t^2 E2 + t^4 E4
     X_k = E_k(1,1)-E_k(1,0)-E_k(0,1)+E_k(0,0)   (gap-gap interaction)
     U_k = E_k(2,0)-2E_k(1,0)+E_k(0,0)          (double occupancy)
     B_k = E_k(0,0)-E_k(1,0)                     (binding of a vacancy to the gap)
  two-cell ring, one vacancy: t_eff = t^(p+1) |<g2|H1 (R H1)^p|g1>|
"""
import os, sys, json
import numpy as np
import scipy.sparse as sp
import bond12 as b
from closure import closure, nbad_open

HERE = os.path.dirname(os.path.abspath(__file__))
R0, R1, E, V = ('R', 0), ('R', 1), ('e', 0), b.V
VI = b.IX[V]
d = 12


def operator(states, bonds):
    index = {w: i for i, w in enumerate(states)}
    L = len(states[0]); rows, cols, vals = [], [], []
    for c, w in enumerate(states):
        for i in range(L):
            h = bonds[i]
            if h is None:
                continue
            j = (i + 1) % L
            col = h[:, w[i] * d + w[j]]
            for r in np.flatnonzero(col):
                z = list(w); z[i], z[j] = divmod(int(r), d); z = tuple(z)
                if z in index:
                    rows.append(index[z]); cols.append(c); vals.append(col[r])
    n = len(states)
    return sp.coo_matrix((vals, (rows, cols)), shape=(n, n)).toarray()


def left_count(w, k):
    m = -1
    for i, x in enumerate(w):
        if x != VI:
            m += 1
            if m == k:
                return sum(1 for y in w[:i] if y == VI)


class BlockResolvent:
    """H0 is block diagonal in the vacancy positions (t=0).  Builds
    R = Q (E0 - H0)^(-1) Q blockwise, with Q removing every H0 eigenvector
    degenerate with E0 (tolerance 1e-9)."""

    def __init__(self, states, H0, e0):
        self.blocks = {}
        for i, s in enumerate(states):
            self.blocks.setdefault(tuple(k for k, x in enumerate(s) if x == VI), []).append(i)
        self.data = []
        self.ground = []
        for key, idx in self.blocks.items():
            idx = np.array(idx)
            e, U = np.linalg.eigh(H0[np.ix_(idx, idx)])
            keep = abs(e - e0) > 1e-9
            self.data.append((idx, U[:, keep], 1. / (e0 - e[keep])))
            for k in np.flatnonzero(~keep):
                v = np.zeros(len(states)); v[idx] = U[:, k]; self.ground.append((key, v))
            assert np.all(e > e0 - 1e-9), 'unperturbed state below E0'

    def __call__(self, x):
        y = np.zeros_like(x)
        for idx, U, inv in self.data:
            y[idx] = U @ (inv * (U.T @ x[idx]))
        return y


def sparse_op(states, bonds):
    return sp.csr_matrix(operator(states, bonds))


def seg_coeffs(p, nL, nR, w):
    seed = [R1, E] + [V] * nL + [R0] * p + [E] + [V] * nR + [R0]
    pr = b.params(p, K=0., musel=0., w=w)
    L = len(seed)
    h1 = b.bond(pr, lam=0., t=1.); h0 = b.bond(pr, lam=0., t=0.)
    _, _, _, st = closure(seed, [h1] * (L - 1) + [None], legal_only=True, penfun=nbad_open)
    st = [s for s in st if left_count(s, 2 + p // 2) == nL]
    H0 = operator(st, [h0] * (L - 1) + [None]).real
    H1 = operator(st, [h1] * (L - 1) + [None]).real - H0
    seedi = tuple(b.IX[x] for x in seed)
    key = tuple(k for k in range(L) if seedi[k] == VI)
    blk = [i for i, s in enumerate(st) if tuple(k for k in range(L) if s[k] == VI) == key]
    e = np.linalg.eigvalsh(H0[np.ix_(blk, blk)])
    e0 = e[0]
    assert e[1] - e[0] > 1e-6
    Rm = BlockResolvent(st, H0, e0)
    assert len(Rm.ground) == 1 and Rm.ground[0][0] == key, 'gap configuration is not the unique unperturbed ground state'
    psi = Rm.ground[0][1]
    v = H1 @ psi
    assert abs(psi @ v) < 1e-12
    Rv = Rm(v)
    E2 = v @ Rv
    E4 = Rv @ H1 @ Rm(H1 @ Rv) - E2 * (Rv @ Rv)
    return e0, E2, E4, len(st)


def ring_teff(p, w):
    """Lowest-order amplitude for one vacancy to cross one cell between two
    gap states of a two-cell ring: t_eff = t^(p+1) |<g'|H1 (R H1)^p|g>|."""
    seed = [R0] * p + [E] + [V] + [R0] * p + [E]
    L = len(seed)
    pr = b.params(p, K=0., musel=0., w=w)
    h1 = b.bond(pr, lam=0., t=1.); h0 = b.bond(pr, lam=0., t=0.)
    _, _, _, st = closure(seed, [h1] * L, legal_only=True)
    H0 = sparse_op(st, [h0] * L)
    H1 = (sparse_op(st, [h1] * L) - H0).real
    H0 = H0.real.toarray()
    seedi = tuple(b.IX[x] for x in seed)
    k1 = (p + 1,); k2 = (L - 1,)
    blk1 = [i for i, s in enumerate(st) if tuple(k for k in range(L) if s[k] == VI) == k1]
    e0 = np.linalg.eigvalsh(H0[np.ix_(blk1, blk1)])[0]
    Rm = BlockResolvent(st, H0, e0)
    g = {key: v for key, v in Rm.ground}
    assert k1 in g and k2 in g
    x = H1 @ g[k1]
    for _ in range(p):
        x = H1 @ Rm(x)
    return float(abs(g[k2] @ x)), len(st), len(Rm.ground)


if __name__ == '__main__':
    ps = [int(x) for x in sys.argv[1].split(',')] if len(sys.argv) > 1 else [5, 6, 8, 10, 12]
    ws = [0., -1e-4, -1e-3, -1e-2]
    out = dict(segment=[], tunnelling=[])
    for p in ps:
        for w in ws:
            c = {}
            for cfg in [(0, 0), (1, 0), (0, 1), (1, 1), (2, 0)]:
                c[cfg] = seg_coeffs(p, cfg[0], cfg[1], w)
            def comb(k, f):
                return float(f(lambda cfg: c[cfg][k]))
            X = [comb(k, lambda E_: E_((1, 1)) - E_((1, 0)) - E_((0, 1)) + E_((0, 0))) for k in range(3)]
            U = [comb(k, lambda E_: E_((2, 0)) - 2 * E_((1, 0)) + E_((0, 0))) for k in range(3)]
            BL = [comb(k, lambda E_: E_((0, 0)) - E_((1, 0))) for k in range(3)]
            BR = [comb(k, lambda E_: E_((0, 0)) - E_((0, 1))) for k in range(3)]
            rec = dict(p=p, w=w, X0=X[0], X2=X[1], X4=X[2], U0=U[0], U2=U[1], U4=U[2],
                       BL0=BL[0], BL2=BL[1], BR2=BR[1], dim_11=c[(1, 1)][3])
            out['segment'].append(rec); print(json.dumps(rec), flush=True)
        for w in (0., -1e-3):
            amp, dim, ng = ring_teff(p, w)
            rec = dict(p=p, w=w, t_eff_coefficient=amp, dim=dim, degenerate_gap_states=ng)
            out['tunnelling'].append(rec); print(json.dumps(rec), flush=True)
    json.dump(out, open(os.path.join(HERE, 'perturbative_scales.json'), 'w'), indent=1)
