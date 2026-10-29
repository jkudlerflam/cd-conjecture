"""Scales of the slow-vacancy fluid of the d=12 bond (zero port, K=infinity,
bits zero, which is exact for every bit string by bit blindness).

Open segment  R_1 E [gap_L] R_0^p E [gap_R] R_0 : one cell of p data sites
between two inter-cell gaps.  R_1 E on the left and the lone R_0 on the right
are inert walls (no vacancy can pass or sit beyond them).
  eps_L, eps_R : energy gain of one vacancy in the left/right gap
  X            : cross term E(1,1)-E(1,0)-E(0,1)+E(0,0)  (gap-gap interaction
                 mediated by the cell; NN interaction on the gap lattice)
  U            : E(2,0)-2E(1,0)+E(0,0)  (extra cost of a doubly occupied gap)
Two-cell ring with one vacancy: tunnelling splitting, compared with perturb.py.
"""
import os, sys, json
import numpy as np
from scipy.sparse.linalg import eigsh
import bond12 as b
from closure import closure, nbad_open, nbad

HERE = os.path.dirname(os.path.abspath(__file__))
R0, R1, E, V = ('R', 0), ('R', 1), ('e', 0), b.V


def lowest(H, k=1):
    n = H.shape[0]
    if n <= 3000:
        return np.linalg.eigvalsh(H.toarray())[:k]
    return np.sort(eigsh(H, k=k, sigma=-3., which='LM', tol=1e-14)[0])


VI = b.IX[V]


def left_count(w, k):
    """number of vacancies before the k-th material letter (0-based)."""
    m = -1
    for i, x in enumerate(w):
        if x != VI:
            m += 1
            if m == k:
                return sum(1 for y in w[:i] if y == VI)
    return None


def segment(p, nL, nR, t):
    """Gap occupations are frozen by keeping only states in which exactly nL
    vacancies lie left of the cell's middle letter; this deletes the
    order-t^(p+1) tunnelling paths and nothing at the orders computed here
    when p >= 5 (a vacancy needs about p/2 hops to reach the middle)."""
    seed = [R1, E] + [V] * nL + [R0] * p + [E] + [V] * nR + [R0]
    pr = b.params(max(p, 1), K=0., musel=0.)
    h = b.bond(pr, lam=0., t=t)
    bonds = [h] * (len(seed) - 1) + [None]
    H, pen, mb, st = closure(seed, bonds, legal_only=True, penfun=nbad_open)
    keep = [i for i, w in enumerate(st) if left_count(w, 2 + p // 2) == nL]
    H = H[keep][:, keep]
    return lowest(H)[0], len(keep)


def ring2(p, t):
    seed = [R0] * p + [E] + [V] + [R0] * p + [E]
    pr = b.params(p, K=0., musel=0.)
    H, pen, mb, st = closure(seed, [b.bond(pr, lam=0., t=t)] * len(seed), legal_only=True)
    ev = lowest(H, 2)
    return ev[1] - ev[0], len(st)


def main():
    """ED cross-checks of perturb.py at small p (where t can be made perturbative):
    (i) segment energies: X/t^4 and B_L/t^2 against the Rayleigh-Schroedinger
        coefficients X4, BL2; (ii) one vacancy on a two-cell ring: the
        degenerate gap states form a ring of N_g = 2p+3 states with hopping tau,
        so the lowest splitting is 2 tau (1-cos(2 pi/N_g)) with tau = c_p t^(p+1)."""
    import perturb
    out = dict(segment=[], tunnelling=[])
    for p in (5, 6):
        _, E2, E4, _ = zip(*[perturb.seg_coeffs(p, *cfg, 0.) for cfg in [(0, 0), (1, 0), (0, 1), (1, 1)]])
        X4 = E4[3] - E4[1] - E4[2] + E4[0]; BL2 = E2[0] - E2[1]
        for t in (0.01, 0.005, 0.0025):
            e00, _ = segment(p, 0, 0, t); e10, _ = segment(p, 1, 0, t); e01, _ = segment(p, 0, 1, t)
            e11, _ = segment(p, 1, 1, t)
            rec = dict(p=p, t=t, X_over_t4_ED=float((e11 - e10 - e01 + e00) / t ** 4), X4_PT=float(X4),
                       BL_over_t2_ED=float((e00 - e10) / t ** 2), BL2_PT=float(BL2))
            out['segment'].append(rec); print(json.dumps(rec), flush=True)
    for p in (2, 3):
        c, _, _ = perturb.ring_teff(p, 0.)
        Ng = 2 * p + 3
        for t in (0.02, 0.01):
            dlt, dim = ring2(p, t)
            pred = 2 * c * t ** (p + 1) * (1 - np.cos(2 * np.pi / Ng))
            rec = dict(p=p, t=t, splitting_ED=float(dlt), splitting_PT=float(pred), ratio=float(dlt / pred))
            out['tunnelling'].append(rec); print(json.dumps(rec), flush=True)
    json.dump(out, open(os.path.join(HERE, 'vacancy_scales.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
