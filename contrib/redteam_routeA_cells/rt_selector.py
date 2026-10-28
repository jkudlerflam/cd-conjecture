"""Selector and competing-segment checks for the literal 23-letter bond.

(a) Ring of two cells built directly from the split literal bond (on-site
    halves included) against the cell decomposition plus the port.
(b) Full 23^4 ring (all words, illegal ones included) at n = 1 with mu_sel
    tuned to select b = 2: the global ground state must be the coded,
    selected sector.
(c) For every cell length b' = n'+1 up to nmax: vacuum-walker energy versus
    the lowest excited-token energy of the same legal cell, i.e. whether the
    selector's competitors are what C22's argument assumes.
Usage: python3 rt_selector.py [nmax] [Delta] [g]
"""
import sys, json, time
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from rt_bond import P, bond, LETTERS, IX, D0, legal_pair, is_exc, onsite_matrix
from rt_cell import bond_dict, enum_words, build, cell, port_factors, port_V, code_basis


def ring_two_cells(n, p, lam):
    """Two-cell ring from the literal bond (split on-site), legal words only."""
    L = 2 * (n + 1)
    h = bond(p, lam=lam)
    hd = bond_dict(h)
    words = enum_words(L, {n, 2 * n + 1}, left_fixed=None, cyclic=True)
    bonds = [(k, (k + 1) % L, hd) for k in range(L)]
    H = build(words, bonds, hd, [], np.zeros((D0, D0)))
    return H.toarray(), words


def decomposed(n, p, lam):
    Hd, info = cell(n, p)
    FM, FX = port_factors(n, info)
    D = Hd.shape[0]
    I = np.eye(D)
    # bond (cell1 marker, cell2 site0) and (cell2 marker, cell1 site0)
    H = np.kron(Hd, I) + np.kron(I, Hd)
    for (A, B, c) in port_V(p, FM, FX):
        H = H + lam * c * np.kron(A, B) + lam * c * np.kron(B, A)
    return H


def full_ring(L, p, lam, k=12):
    """All 23^L words, literal bond on every ring bond."""
    h = sp.csr_matrix(bond(p, lam=lam))
    I = lambda m: sp.identity(D0 ** m, format='csr')
    H = None
    for s in range(L - 1):
        term = sp.kron(sp.kron(I(s), h), I(L - s - 2), format='csr')
        H = term if H is None else H + term
    # wrap bond (L-1, 0): permute sites by a cyclic shift
    dims = [D0] * L
    idx = np.arange(D0 ** L).reshape(dims)
    perm = np.moveaxis(idx, L - 1, 0).reshape(-1)       # new order: site L-1 first
    Pm = sp.csr_matrix((np.ones(D0 ** L), (np.arange(D0 ** L), perm)), shape=(D0 ** L, D0 ** L))
    wrap = sp.kron(h, I(L - 2), format='csr')            # acts on (site L-1, site 0) in shifted order
    H = H + Pm.T @ wrap @ Pm
    w = spla.eigsh(H, k=k, which='SA', tol=1e-12)[0]
    return np.sort(w)


def competitors(nmax, p):
    out = []
    for nn in range(0, nmax + 1):
        t0 = time.time()
        if nn == 0:
            # cell = marker only: vacuum M_v vs excited M_s
            Evac = p['musel']
            Eexc = p['musel'] + p['Delta'] + p['DM'] - abs(p['delta'])
            out.append(dict(n=nn, b=1, E_vac=Evac, E_exc_min=Eexc, diff=Eexc - Evac))
            continue
        hk = bond(p, lam=0.0, onsite=False)
        hd = bond_dict(hk)
        words = enum_words(nn + 1, {nn}, left_fixed=IX[('M', None, None)])
        exc = np.array([any(is_exc(LETTERS[c]) for c in w) for w in words])
        bonds = [(-1, 0, hd)] + [(k, k + 1, hd) for k in range(nn)]
        H = build(words, bonds, hd, range(nn + 1), onsite_matrix(p), left_fixed=IX[('M', None, None)])
        ie = np.where(exc)[0]; iv = np.where(~exc)[0]
        assert abs(H[ie][:, iv]).max() == 0 if len(ie) and len(iv) else True
        He = H[ie][:, ie]
        Hv = H[iv][:, iv]
        if He.shape[0] < 3000:
            Eexc = np.linalg.eigvalsh(He.toarray())[0]
            Evac = np.linalg.eigvalsh(Hv.toarray())[0]
        else:
            Eexc = spla.eigsh(He, k=1, which='SA', tol=1e-12)[0][0]
            Evac = spla.eigsh(Hv, k=1, which='SA', tol=1e-12)[0][0]
        Evac_formula = p['musel'] - 2 * p['t'] * np.cos(np.pi / (nn + 2)) + 0.0
        out.append(dict(n=nn, b=nn + 1, dim=len(words), E_vac=float(Evac), E_vac_formula=float(Evac_formula),
                        E_exc_min=float(Eexc), diff=float(Eexc - Evac), sec=time.time() - t0))
        print(json.dumps(out[-1]), flush=True)
    return out


def selection_margin(Ecell, b):
    """gamma_b = min_{k != b} k (e_k - e_b), e_k = E_cell(k)/k, from the TRUE cell
    ground energies (vacuum or excited, whichever is lower)."""
    e = {k: Ecell[k] / k for k in Ecell}
    return min(k * (e[k] - e[b]) for k in Ecell if k != b)


if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'ring':          # (a) literal two-cell ring vs decomposition, n = 1, 2
        res = {}
        for n in (1, 2):
            p = P(td=1e-3)
            Hr, words = ring_two_cells(n, p, 0.1)
            wr = np.linalg.eigvalsh(Hr); wd = np.linalg.eigvalsh(decomposed(n, p, 0.1))
            res['n%d' % n] = dict(dim=len(words), max_eig_diff=float(np.abs(wr - wd).max()))
            print(res['n%d' % n], flush=True)
        json.dump(res, open('ring_check.json', 'w'), indent=1)
    elif mode == 'full':        # (b) all 23^4 words, n = 1, mu_sel selecting b = 2
        p = P(td=1e-3, musel=0.1)
        res = {}
        for lam in (0.0, 0.1):
            t0 = time.time()
            res['lam%g' % lam] = dict(lowest=list(full_ring(4, p, lam, k=12)), sec=time.time() - t0)
            print(lam, res['lam%g' % lam], flush=True)
        res['expected_lam0_ground'] = 2 * (1 - 1.9)
        json.dump(res, open('full_ring_L4.json', 'w'), indent=1, default=float)
    else:                       # (c) competing cell lengths
        nmax, Delta, g = int(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])
        res = competitors(nmax, P(td=1e-3, Delta=Delta, g=g))
        json.dump(res, open('competitors_D%g_g%g.json' % (Delta, g), 'w'), indent=1, default=float)
