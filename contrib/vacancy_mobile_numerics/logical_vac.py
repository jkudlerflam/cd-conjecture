#!/usr/bin/env python3
"""Logical vacancy model: the one-vacancy lift of the route-A logical chain at n flavours, and its references.

Model LV_n.  One site is a vacancy (index 0) or a logical cell (index 1 + b, b in [0, 2^n), bit f of b
(flavour 0 most significant) is the occupation of Jordan-Wigner fermion (cell, f)).  The payload bond is the
Gaussian part of the logical chain of ../explicit_candidate/README.md Section 2.4 with a diagonal hopping matrix,

    h_log = - sum_f tau_f (c+_{j,f} c_{j+1,f} + h.c.) + (1/2) sum_f eps_f (n_{j,f} + n_{j+1,f}),

with fermions ordered by cell, then flavour (so the bond carries the intra-cell Z strings).  It is lifted
exactly as in ../counterexample_campaign/doping/bond24.py (vacbond.vacancy_lift):

    h = -t sum_b (|0 b><b 0| + h.c.) + J h_log (occupied pairs only) - (mu_v/2)(n (x) 1 + 1 (x) n).

Vacancies carry no modes, so the vacancy swap needs no fermion sign.  Conserved: occupied number M and every
flavour number N_f.  Expected c = 1 (charge) + number of partially filled flavours.

References (exactly solvable, same number of modes): free-fermion (XX) rings with the same L and particle
numbers (charge: M particles on L sites; flavour f: N_f particles on L sites, the same Friedel wavevector as
the payload in physical coordinates), and k decoupled XX chains with velocity ratios for VUMPS.

usage: logical_vac.py ring --n 1 --L 12 --M 8 --N 4 [--t 1 --J 0.5 --taus 1 --eps 0] --out data/lv/....json
       logical_vac.py check            (literal-lift vs sector ED, small rings)
Claude (Opus 5.5), 2026-10-28.
"""
import os, sys, json, time, argparse, itertools
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh
HERE = os.path.dirname(os.path.abspath(__file__))


# ----------------------------------------------------------------------------- bonds
def jw_ops(nm):
    """Annihilators c_i on nm modes (mode 0 the most significant kron factor), with Z = (-1)^n strings."""
    sm = np.array([[0., 1.], [0., 0.]])        # |0><1|: annihilates
    Z = np.diag([1., -1.])
    I = np.eye(2)
    out = []
    for i in range(nm):
        f = [Z] * i + [sm] + [I] * (nm - i - 1)
        m = f[0]
        for x in f[1:]:
            m = np.kron(m, x)
        out.append(m)
    return out


def payload_bond(taus, eps, vz=None):
    """vz[f]: Ising term vz (n_{j,f} - 1/2)(n_{j+1,f} - 1/2), XXZ anisotropy Delta = vz / (2 tau) (controls)."""
    n = len(taus)
    vz = [0.0] * n if vz is None else vz
    c = jw_ops(2 * n)
    h = np.zeros((4 ** n, 4 ** n))
    I = np.eye(4 ** n)
    for f in range(n):
        h += vz[f] * (c[f].T @ c[f] - I / 2) @ (c[n + f].T @ c[n + f] - I / 2)
        hop = c[f].T @ c[n + f]
        h -= taus[f] * (hop + hop.T)
        h += 0.5 * eps[f] * (c[f].T @ c[f] + c[n + f].T @ c[n + f])
    return h


def lift(hlog, t, J, mu_v):
    q = int(round(np.sqrt(hlog.shape[0])))
    d = q + 1
    out = np.zeros((d * d, d * d))
    occ = np.array([a * d + b for a in range(1, d) for b in range(1, d)])
    out[np.ix_(occ, occ)] = J * hlog
    for a in range(1, d):
        out[a, a * d] = out[a * d, a] = -t
    nn = np.r_[0, np.ones(q)]
    out -= mu_v / 2 * np.diag(np.add.outer(nn, nn).reshape(-1))
    return out


def model_bond(n, t, J, mu_v, taus, eps, vz=None):
    return lift(payload_bond(taus, eps, vz), t, J, mu_v)


# ----------------------------------------------------------------------------- sector ED on rings
def bits(b, n):
    return [(b >> (n - 1 - f)) & 1 for f in range(n)]


def sector_states(L, n, M, N):
    """All configurations (base d = 1 + 2^n, site 0 most significant) with M cells and flavour numbers N."""
    d = 1 + 2 ** n
    pos = np.array(list(itertools.combinations(range(L), M)), dtype=np.int64)          # (P, M)
    # payload strings of length M with N[f] ones in flavour f
    per = []
    for f in range(n):
        cf = np.zeros((0, M), dtype=np.int64)
        rows = []
        for ones in itertools.combinations(range(M), N[f]):
            r = np.zeros(M, dtype=np.int64); r[list(ones)] = 1; rows.append(r)
        per.append(np.array(rows, dtype=np.int64).reshape(-1, M))
    pay = np.zeros((1, M), dtype=np.int64)
    for f in range(n):
        pay = (pay[:, None, :] + (per[f][None, :, :] << (n - 1 - f))).reshape(-1, M)
    pw = d ** (L - 1 - np.arange(L, dtype=np.int64))                                   # place values
    s = (pw[pos][:, None, :] * (1 + pay)[None, :, :]).sum(-1).reshape(-1)
    return np.sort(s), d


def build_H(L, d, states, h, tol=1e-14):
    dim = len(states)
    rows, cols, vals = [], [], []
    diag = np.zeros(dim)
    pw = d ** (L - 1 - np.arange(L, dtype=np.int64))
    digits = ((states[:, None] // pw[None, :]) % d).astype(np.int8)                    # (dim, L)
    for i in range(L):
        j = (i + 1) % L
        p = digits[:, i].astype(np.int64) * d + digits[:, j]
        diag += h[p, p]
        hoff = h.copy(); np.fill_diagonal(hoff, 0)
        for k in np.unique(p):
            sel = np.nonzero(p == k)[0]
            for o in np.nonzero(np.abs(hoff[:, k]) > tol)[0]:
                oi, oj = divmod(o, d)
                new = states[sel] + (oi - digits[sel, i].astype(np.int64)) * pw[i] + (oj - digits[sel, j].astype(np.int64)) * pw[j]
                idx = np.searchsorted(states, new)
                ok = (idx < dim)
                ok[ok] &= states[idx[ok]] == new[ok]
                assert ok.all(), 'leak out of sector'
                rows.append(idx.astype(np.int32)); cols.append(sel.astype(np.int32)); vals.append(np.full(len(sel), hoff[o, k]))
    ar = np.arange(dim, dtype=np.int32)
    rows = np.concatenate(rows + [ar]); cols = np.concatenate(cols + [ar])
    vals = np.concatenate(vals + [diag])
    return sp.csr_matrix((vals, (rows, cols)), shape=(dim, dim)), digits


def interval_entropies(L, d, n, states, digits, psi):
    """S(l) of the block of sites 0..l-1, l = 1..L-1, from the sector vector (blocks by left quantum numbers)."""
    out = []
    occ = (digits > 0).astype(np.int64)
    digits = digits.astype(np.int64)
    fl = [np.where(digits > 0, ((digits - 1) >> (n - 1 - f)) & 1, 0) for f in range(n)]
    for l in range(1, L // 2 + 1):
        base = d ** (L - l)
        left = states // base
        right = states % base
        key = occ[:, :l].sum(1)
        for f in range(n):
            key = key * (l + 1) + fl[f][:, :l].sum(1)
        S = 0.0
        for kv in np.unique(key):
            sel = np.nonzero(key == kv)[0]
            lu, li = np.unique(left[sel], return_inverse=True)
            ru, ri = np.unique(right[sel], return_inverse=True)
            A = sp.csr_matrix((psi[sel], (li, ri)), shape=(len(lu), len(ru)))
            rho = (A @ A.conj().T if A.shape[0] <= A.shape[1] else A.conj().T @ A).toarray()
            w = np.linalg.eigvalsh(rho)
            w = w[w > 1e-15]
            S -= float((w * np.log(w)).sum())
        out.append(S)
    return out + out[:(L - 1) // 2][::-1]          # pure state: S(l) = S(L - l)


def ring(args):
    n = args.n
    taus = args.taus; eps = args.eps
    h = model_bond(n, args.t, args.J, args.mu_v, taus, eps, args.vz)
    t0 = time.time()
    states, d = sector_states(args.L, n, args.M, args.N)
    H, digits = build_H(args.L, d, states, h)
    t1 = time.time()
    k = 4
    w, U = eigsh(H, k=k, which='SA', tol=1e-12, ncv=40, maxiter=100000)
    o = np.argsort(w); w = w[o]; U = U[:, o]
    psi = U[:, 0]
    S = interval_entropies(args.L, d, n, states, digits, psi)
    occ = (digits > 0)
    # payload density per occupied site, charge density
    res = dict(args=vars(args), dim=int(len(states)), levels=w.tolist(), gap=float(w[1] - w[0]), S=S,
               rho=args.M / args.L, adjacency=float(np.sum(np.abs(psi) ** 2 * (occ & np.roll(occ, -1, 1)).sum(1)) / args.L),
               t_build=t1 - t0, t_total=time.time() - t0)
    return res


# ----------------------------------------------------------------------------- free-fermion references
def ff_ring_entropies(L, N):
    """Interval entropies S(l), l = 1..L-1, of the XX (hard-core boson) ring ground state with N particles.
    Fermion boundary condition (-1)^(N-1): closed shell, unique ground state."""
    if N == 0 or N == L:
        return [0.0] * (L - 1)
    phi = 0.0 if N % 2 else 0.5
    k = 2 * np.pi * (np.arange(L) - (L - 1) // 2 + phi) / L
    e = -2 * np.cos(k)
    occ = np.argsort(e, kind='stable')[:N]
    x = np.arange(L)
    C = (np.exp(1j * np.outer(x, k[occ])) @ np.exp(-1j * np.outer(x, k[occ])).T) / L
    out = []
    for l in range(1, L):
        nu = np.linalg.eigvalsh(C[:l, :l]).clip(1e-15, 1 - 1e-15)
        out.append(float(-(nu * np.log(nu) + (1 - nu) * np.log(1 - nu)).sum()))
    return out


def check():
    """Sector ED against the literal d^L Hamiltonian on small rings."""
    out = []
    for (n, L, M, N) in [(1, 5, 3, [1]), (1, 6, 4, [2]), (2, 4, 3, [1, 2]), (2, 5, 3, [2, 1])]:
        taus = [1.0, 0.5][:n]; eps = [0.3, -0.2][:n]
        h = model_bond(n, 1.0, 0.7, 0.4, taus, eps, [0.6, 1.1][:n])
        d = 1 + 2 ** n
        Hf = np.zeros((d ** L, d ** L))
        I = np.eye(d)
        for i in range(L):
            if i < L - 1:
                Hf += np.kron(np.kron(np.eye(d ** i), h), np.eye(d ** (L - i - 2)))
            else:   # wrap bond (L-1, 0): permute
                P = np.arange(d ** L).reshape([d] * L)
                perm = np.roll(np.arange(L), 1)
                hb = np.kron(h, np.eye(d ** (L - 2)))
                T = hb.reshape([d] * (2 * L))
                # hb acts on sites (0,1); relabel to (L-1, 0)
                order = [L - 1, 0] + list(range(1, L - 1))
                inv = np.argsort(order)
                T = T.transpose(list(inv) + [L + x for x in inv])
                Hf += T.reshape(d ** L, d ** L)
        states, _ = sector_states(L, n, M, N)
        Hs, _ = build_H(L, d, states, h)
        sub = Hf[np.ix_(states, states)]
        leak = np.abs(Hf[:, states]).sum() - np.abs(sub).sum()
        out.append(dict(n=n, L=L, M=M, N=N, dim=len(states), max_diff=float(np.abs(sub - Hs.toarray()).max()),
                        leak=float(leak), herm=float(np.abs(h - h.T).max())))
        print(out[-1], flush=True)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('mode', choices=['ring', 'check', 'ref'])
    ap.add_argument('--n', type=int, default=1)
    ap.add_argument('--L', type=int, default=9)
    ap.add_argument('--M', type=int, default=6)
    ap.add_argument('--N', type=int, nargs='+', default=[3])
    ap.add_argument('--t', type=float, default=1.0)
    ap.add_argument('--J', type=float, default=0.5)
    ap.add_argument('--mu_v', type=float, default=0.0)
    ap.add_argument('--taus', type=float, nargs='+', default=[1.0])
    ap.add_argument('--eps', type=float, nargs='+', default=[0.0])
    ap.add_argument('--vz', type=float, nargs='+', default=None)
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    if a.mode == 'check':
        res = check()
    elif a.mode == 'ring':
        res = ring(a)
        print(json.dumps({k: v for k, v in res.items() if k != 'S'}), flush=True)
        print('S =', np.round(res['S'], 5).tolist())
    if a.out:
        os.makedirs(os.path.dirname(a.out), exist_ok=True)
        json.dump(res, open(a.out, 'w'), indent=1)


if __name__ == '__main__':
    main()
