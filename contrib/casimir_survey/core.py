#!/usr/bin/env python3
"""Casimir survey core: the bond norm ||h|| (half spread minimized over h -> h + o(x)1 - 1(x)o),
ring spectra through the symmetry-resolved ED engine of ../counterexample_search/ced.py,
and the two estimators of A = sum_a c_a v_a.

  A_cas : fit E_0(L)/L = e_0 - pi A/(6 L^2) + b/L^4 over commensurate ring lengths.
  A_cv  : c_half (periodic half-ring entropy, two sizes) times v from the lowest level at
          k_0 +- 2 pi/L in the ground-state charge sector (the current state, x = 1).

numpy and scipy only.
"""
import os, sys
for _k in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ.setdefault(_k, '1')
import numpy as np
from scipy import optimize

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(HERE, '..', 'counterexample_search'))
from ced import Sector  # noqa: E402


# ---------------------------------------------------------------- ||h||
def _herm_basis(d):
    B = []
    for a in range(d):
        for b in range(a + 1, d):
            m = np.zeros((d, d), complex); m[a, b] = m[b, a] = 1 / np.sqrt(2); B.append(m)
            m = np.zeros((d, d), complex); m[a, b] = -1j / np.sqrt(2); m[b, a] = 1j / np.sqrt(2); B.append(m)
    for a in range(1, d):
        m = np.zeros((d, d), complex)
        m[np.arange(a), np.arange(a)] = 1; m[a, a] = -a
        B.append(m / np.sqrt(a * (a + 1)))
    return np.array(B)


def spread(h):
    e = np.linalg.eigvalsh((h + h.conj().T) / 2)
    return e[-1] - e[0]


def hnorm(h, d, return_o=False, betas=(4, 16, 64, 256, 1024, 4096)):
    """Half the spread of spec(h + o(x)1 - 1(x)o), minimized over Hermitian o.
    The spread is convex in o; we minimize a log-sum-exp smoothing with an annealed
    inverse temperature (BFGS, analytic gradient) and return the exact spread at the end."""
    h = (np.asarray(h, complex) + np.asarray(h, complex).conj().T) / 2
    B = _herm_basis(d)
    I = np.eye(d)
    KB = np.array([np.kron(b, I) - np.kron(I, b) for b in B])
    s0 = spread(h)
    if s0 < 1e-14:
        return (0.0, np.zeros((d, d))) if return_o else 0.0

    def soft(x, beta):
        K = h + np.tensordot(x, KB, 1)
        e, U = np.linalg.eigh(K)
        bb = beta / s0
        wp = np.exp(bb * (e - e[-1])); wp /= wp.sum()
        wm = np.exp(-bb * (e - e[0])); wm /= wm.sum()
        f = (e[-1] + np.log(np.sum(np.exp(bb * (e - e[-1])))) / bb) - (e[0] - np.log(np.sum(np.exp(-bb * (e - e[0])))) / bb)
        w = wp - wm
        # d e_i / d x_k = <u_i| KB_k |u_i>
        g = np.einsum('kab,ai,bi,i->k', KB, U.conj(), U, w).real
        return f, g

    x = np.zeros(len(B))
    for beta in betas:
        r = optimize.minimize(soft, x, args=(beta,), jac=True, method='BFGS', options=dict(gtol=1e-12, maxiter=2000))
        x = r.x
    best = spread(h + np.tensordot(x, KB, 1))
    if best > s0:
        best, x = s0, np.zeros(len(B))
    o = np.tensordot(x, B, 1)
    return (best / 2, o) if return_o else best / 2


# ---------------------------------------------------------------- ring spectra
def ring(h, d, L, charges=None, targets=(None,), ks=None, nlev=2, entropy=False, tol=1e-10):
    """Lowest nlev levels in each (charge target, momentum k) block of H = sum_i h_{i,i+1}
    on a ring of L sites. ks=None means all k = 0..L-1. Returns a dict with E0, k0, the
    ground-state target, the spectrum, and (entropy=True) the ground-state block entropies."""
    out = {'L': L, 'spec': {}}
    best = (np.inf, None, None, None, None)
    for t in targets:
        sec = Sector(d, L, [h], charges=charges, target=t, cache=False)
        if len(sec.basis) == 0:
            continue
        kk = range(L) if ks is None else [k % L for k in ks]
        for k in kk:
            e, v, res = sec.solve([1.0], k, nlev, tol=tol)
            if len(e) == 0:
                continue
            out['spec'][(None if t is None else tuple(np.atleast_1d(t).tolist()), int(k))] = [float(x) for x in e]
            if e[0] < best[0] - 1e-12:
                best = (float(e[0]), int(k), t, v[:, 0], sec)
    out['E0'], out['k0'], out['target'] = best[0], best[1], best[2]
    if entropy:
        sec = best[4]
        psi = sec.expand(best[1], best[3])
        out['S'] = sec.entropies(psi, cuts=[L // 2])
    return out


def gap_at(r, dk=1, same_target=True):
    """Lowest level at k0 +- dk (units of 2 pi/L) minus E0."""
    L, k0, t = r['L'], r['k0'], r['target']
    tt = None if t is None else tuple(np.atleast_1d(t).tolist())
    vals = []
    for (tg, k), e in r['spec'].items():
        if same_target and tg != tt:
            continue
        if k in ((k0 + dk) % L, (k0 - dk) % L):
            vals.append(e[0])
    return min(vals) - r['E0'] if vals else np.nan


# ---------------------------------------------------------------- estimators
def fit_cas(Ls, Es, e0=None, quartic=True):
    """Least-squares fit E/L = e0 - (pi/6) A / L^2 [+ b / L^4]. Returns (A, e0, b)."""
    L = np.asarray(Ls, float); y = np.asarray(Es, float) / L
    cols = [] if e0 is not None else [np.ones_like(L)]
    cols.append(-np.pi / (6 * L ** 2))
    if quartic:
        cols.append(1 / L ** 4)
    X = np.vstack(cols).T
    yy = y - (e0 if e0 is not None else 0)
    coef = np.linalg.lstsq(X, yy, rcond=None)[0]
    if e0 is None:
        return coef[1], coef[0], (coef[2] if quartic else 0.0)
    return coef[0], e0, (coef[1] if quartic else 0.0)


def cas_local(Ls, Es):
    """Three-point (consecutive) Casimir estimates A(L_{i}, L_{i+1}, L_{i+2}) to show drift,
    and two-point estimates without the 1/L^4 term."""
    Ls, Es = list(Ls), list(Es)
    three = [fit_cas(Ls[i:i + 3], Es[i:i + 3])[0] for i in range(len(Ls) - 2)]
    two = [fit_cas(Ls[i:i + 2], Es[i:i + 2], quartic=False)[0] for i in range(len(Ls) - 1)]
    return three, two


def c_half(S1, L1, S2, L2):
    """c from the half-ring entropies of two periodic ground states."""
    return 3 * (S2 - S1) / np.log(L2 / L1)
