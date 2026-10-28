#!/usr/bin/env python3
"""Second-order Schrieffer-Wolff logical chain of the encoded chain.

With H0 = sum_j H_cell(j), V = sum_j sum_a F_a(pR_j) F_a(pL_{j+1}), P the code projector and
R = Q (E0 - H0)^{-1} Q,
    H_eff = P H0 P + lam P V P + lam^2 P V R V P + O(lam^3).
Because H0 is a sum over cells, P V R V P has exactly two kinds of terms:
  T2 (bond j,j+1 twice; intermediate states with cell j and/or j+1 excited)
     T2 = sum_{mn not both code} Y[ik,mn] Y[jl,mn]^* / (2 e_c - E_m - E_n),  Y = sum_a <g|F^R_a|m><g|F^L_a|n>,
  T3 (bonds j,j+1 and j+1,j+2; only the middle cell excited)
     T3 = sum_ab alpha eta_a (x) G^dag (F^R_b R_c F^L_a + F^L_a R_c F^R_b) G (x) alpha gamma_b.
The first-order term is T1 = alpha^2 sum_a eta_a (x) gamma_a.  All on the 2^n-dim logical sites.
Running this file checks T1 + T2 against exact diagonalization of two cells (residual O(lam^3))
and prints the Majorana-monomial content of T2 and T3.

Written by Claude (Opus 5.5) for contrib/encoded_chain_numerics, 2026-10-27.
"""
import os, sys, itertools
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np
from scipy import linalg
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh
import cells as C


def second_order(c):
    H = C.bigsite_H(c)
    E, Vv = linalg.eigh(H)
    N = c['N']
    o = C.cell_spectrum_and_alpha(c)
    G = o['G']
    B = np.hstack([G, Vv[:, N:]])
    Eb = np.concatenate([np.full(N, E[0]), E[N:]])
    ec = E[0]
    FL, FR = C.bigsite_ports(c)
    r = len(FL)
    Lm = [G.conj().T @ f @ B for f in FL]          # (N, D)
    Rm = [G.conj().T @ f @ B for f in FR]
    D = H.shape[0]
    Y = sum(np.einsum('im,kn->ikmn', Rm[a], Lm[a]) for a in range(r))   # (N,N,D,D)
    w = 1.0 / (2 * ec - Eb[:, None] - Eb[None, :] + 0j)
    w[:N, :N] = 0
    T2 = np.einsum('ikmn,mn,jlmn->ikjl', Y, w, Y.conj()).reshape(N * N, N * N)
    Rc = np.zeros(D); Rc[N:] = 1.0 / (ec - E[N:])
    Rop = B[:, N:] @ np.diag(Rc[N:]) @ B[:, N:].conj().T
    al = o['alphaL'][0]
    T3 = np.zeros((N ** 3, N ** 3), complex)
    for a in range(r):
        for b in range(r):
            M = G.conj().T @ (FR[b] @ Rop @ FL[a] + FL[a] @ Rop @ FR[b]) @ G
            T3 += al * al * np.kron(np.kron(c['eta'][a], M), c['gam'][b])
    T1 = al * al * sum(np.kron(e, g) for e, g in zip(c['eta'], c['gam']))
    return dict(T1=T1, T2=(T2 + T2.conj().T) / 2, T3=(T3 + T3.conj().T) / 2, alpha=al, ec=ec, N=N)


def eff_terms(c, lam, ncells, so=None, order=2):
    so = so or second_order(c)
    N = so['N']
    h2 = lam * so['T1'] + (lam ** 2 * so['T2'] if order >= 2 else 0)
    terms = [(j, 2, h2) for j in range(ncells)]
    if order >= 2:
        terms += [(j, 3, lam ** 2 * so['T3']) for j in range(ncells)]
    return [N] * ncells, terms


# ----------------------------------------------------------------------------- checks
def two_cell_exact(c, lam, k):
    H = sp.csr_matrix(C.bigsite_H(c))
    D = H.shape[0]
    FL, FR = C.bigsite_ports(c)
    I = sp.identity(D, format='csr')
    Ht = sp.kron(H, I) + sp.kron(I, H) + lam * sum(sp.kron(sp.csr_matrix(fr), sp.csr_matrix(fl)) for fr, fl in zip(FR, FL))
    if D * D <= 2000:
        return linalg.eigvalsh(Ht.toarray())[:k]
    return np.sort(eigsh(Ht, k=k, which='SA', tol=1e-12)[0])


def jw_majoranas(c, ncell):
    """Global Majoranas c_{a,j} = (prod_{k<j} P_k) gamma_{a,j} on ncell logical sites, a = 1..2n."""
    from majorana_bank_check import majoranas
    gam, P = majoranas(c['n'])
    N = c['N']
    out = []
    for j in range(ncell):
        for a in range(len(gam)):
            ops = [P] * j + [gam[a]] + [np.eye(N)] * (ncell - j - 1)
            M = ops[0]
            for o in ops[1:]:
                M = np.kron(M, o)
            out.append(((a, j), M))
    return out


def monomials(T, c, ncell, tol=1e-9):
    """Expand T in products of JW Majoranas (the full Clifford basis on ncell sites)."""
    mj = jw_majoranas(c, ncell)
    m = len(mj)
    dim = T.shape[0]
    res = []
    for k in range(0, m + 1, 2):
        for idx in itertools.combinations(range(m), k):
            M = np.eye(dim, dtype=complex)
            for i in idx:
                M = M @ mj[i][1]
            coef = np.trace(M.conj().T @ T) / dim
            if abs(coef) > tol:
                res.append((tuple(mj[i][0] for i in idx), coef))
    return res


if __name__ == '__main__':
    np.set_printoptions(precision=5, suppress=True, linewidth=150)
    for (q, qp, n) in [(2, 1, 1), (3, 1, 2)]:
        for g in (0.2, 0.5):
            c = C.cell(q, qp, n, 1.0, g)
            so = second_order(c)
            N = so['N']
            print(f'=== q={q} r={c["r"]} g={g}: alpha={so["alpha"]:.5f} t1/lam=alpha^2={so["alpha"]**2:.5f}')
            for lam in (0.01, 0.02, 0.04):
                ex = two_cell_exact(c, lam, N * N) - 2 * so['ec']
                e1 = linalg.eigvalsh(lam * so['T1'])
                e2 = linalg.eigvalsh(lam * so['T1'] + lam ** 2 * so['T2'])
                print(f'  lam={lam}: max|exact - 1st order| = {np.max(abs(ex - e1)):.2e}, '
                      f'max|exact - 2nd order| = {np.max(abs(ex - e2)):.2e}')
            if q == 2:
                print('  T2 monomials (cells 0,1; (flavour, cell)):')
                for mono, cf in monomials(so['T2'], c, 2):
                    print('    ', mono, np.round(cf, 6))
                print('  T3 monomials (cells 0,1,2):')
                for mono, cf in monomials(so['T3'], c, 3):
                    print('    ', mono, np.round(cf, 6))
            else:
                mon2 = monomials(so['T2'], c, 2)
                quart = [(mm, cf) for mm, cf in mon2 if len(mm) == 4]
                quad = [(mm, cf) for mm, cf in mon2 if len(mm) == 2]
                print(f'  T2: {len(quad)} bilinears, {len(quart)} quartic terms; largest quartic:')
                for mm, cf in sorted(quart, key=lambda x: -abs(x[1]))[:12]:
                    print('    ', mm, np.round(cf, 6))
                print('  T2 bilinears:')
                for mm, cf in quad:
                    print('    ', mm, np.round(cf, 6))


# ----------------------------------------------------------------------------- all-order two-cell (CORE)
def core_bond(c, lam, so=None):
    """All-order nearest-neighbour logical bond from exact diagonalization of two cells:
    with Psi the lowest N^2 two-cell eigenstates and A = (G x G)^dag Psi, the canonical
    (des Cloizeaux) effective Hamiltonian H2 = (A A^dag)^{-1/2} A E A^dag (A A^dag)^{-1/2};
    h_bond = H2 - 2 e_c.  Agrees with lam T1 + lam^2 T2 + O(lam^3).  Also returns the gap from
    the N^2-dim manifold to the next two-cell level and the smallest singular value of A."""
    so = so or second_order(c)
    N = c['N']
    G = C.cell_spectrum_and_alpha(c)['G']
    H = sp.csr_matrix(C.bigsite_H(c))
    D = H.shape[0]
    FL, FR = C.bigsite_ports(c)
    I = sp.identity(D, format='csr')
    Ht = sp.kron(H, I) + sp.kron(I, H) + lam * sum(sp.kron(sp.csr_matrix(fr), sp.csr_matrix(fl)) for fr, fl in zip(FR, FL))
    k = N * N
    if D * D <= 3000:
        E, Psi = linalg.eigh(Ht.toarray())
        gap = E[k] - E[k - 1]
        E, Psi = E[:k], Psi[:, :k]
    else:
        E, Psi = eigsh(Ht.tocsc(), k=k + 4, which='SA', tol=1e-12)
        o = np.argsort(E); E, Psi = E[o], Psi[:, o]
        gap = E[k] - E[k - 1]
        E, Psi = E[:k], Psi[:, :k]
    GG = np.kron(G, G)
    A = GG.conj().T @ Psi
    w, U = linalg.eigh(A @ A.conj().T)
    Sm = U @ np.diag(w ** -0.5) @ U.conj().T
    H2 = Sm @ A @ np.diag(E) @ A.conj().T @ Sm
    return (H2 + H2.conj().T) / 2 - 2 * so['ec'] * np.eye(k), float(gap), float(np.sqrt(w.min()))

