"""Item 5: does the frozen selector (nu at the midpoint, cell energy nu - 2cos(pi/(b+1)) < 0)
survive in h10, where the on-site energies are genuine but every two-site material term
(in particular the walker's hop) needs physical contact?

Literal ED of the referee's h10 on a ring: one cell of n data sites (all bits 0, lam = g = t_d = 0,
so only the letters L0, R0, W, V occur) plus Nv vacancies, against the headless ring R0^(n+1) with
the same Nv.  Both are F = 0 sectors.  eps_cell = E0(one cell) - E0(headless) is the energy of
putting a head into the material; the README's selector needs eps_cell = J(nu - 2cos(pi/(n+2))) < 0.
"""
import sys, os, json, itertools
import numpy as np
import scipy.sparse.linalg as spla
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ref_bonds as rb, chain as ch

def nu_sel(b):
    a = lambda x: 2 * (1 - np.cos(np.pi / (x + 1)))
    dl = lambda x: (x + 1) * a(x) - x * a(x + 1)
    return 2 - (dl(b) + dl(b - 1)) / 2

def ground(seeds, par):
    gen = dict(par)
    words, index = ch.closure(seeds, rb.h10(gen), 10, True)
    H, leak = ch.hamiltonian(words, index, rb.h10(par), 10, True)
    assert leak == 0
    if H.shape[0] <= 400:
        return float(np.linalg.eigvalsh(H.toarray())[0]), H.shape[0]
    return float(spla.eigsh(H, k=1, which="SA", tol=1e-11)[0][0]), H.shape[0]

def placements(material, Nv):
    """all words with the material letters in cyclic order starting at site 0 or later (ring)."""
    M = len(material); L = M + Nv
    out = []
    for vac in itertools.combinations(range(L), Nv):
        w = [rb.V] * L
        it = iter(material)
        for i in range(L):
            if i not in vac:
                w[i] = next(it)
        out.append(tuple(w))
    return out

res = []
for n in (9, 4):
    b = n + 1
    nu = nu_sel(b)
    for (J, t) in [(1.0, 1.0), (1.0, 5.0), (1.0, 0.2)]:
        par = dict(nu=nu, Delta=0.5, D=3.0, f=0.0, K=2.0, td=0.0, g=0.0, lam=0.0, phiA=0.1, psiA=np.pi, psiB=np.pi / 2,
                   J=J, t=t, mu=0.0)
        for Nv in range(0, 7 if n == 9 else 9):
            cell = (rb.W,) + (rb.R0,) * n
            dead = (rb.R0,) * b
            EA, dimA = ground(placements(cell, Nv), par)
            EB, dimB = ground(placements(dead, Nv), par)
            d = dict(n=n, b=b, nu=float(nu), J=J, t=t, Nv=Nv, vacancy_density=Nv / (b + Nv), dim_cell=dimA, dim_headless=dimB,
                     E_cell_sector=EA, E_headless_sector=EB, eps_cell=EA - EB,
                     readme_formula=J * (nu - 2 * np.cos(np.pi / (b + 1))))
            res.append(d); print(json.dumps(d), flush=True)
json.dump(res, open(os.path.join(HERE, "item5_mobile_selector.json"), "w"), indent=1)
