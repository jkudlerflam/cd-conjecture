#!/usr/bin/env python3
"""Marginal current-current couplings of the second-order logical chain (r = 4, so(4) = su(2)+ + su(2)-).

The first-order chain is r Majorana chains i t c_{a,j} c_{a,j+1}, with right movers at k = 0 and
left movers at k = pi, so c_{a,j} -> psi_{aR}(x) + (-1)^j psi_{aL}(x).  Each quartic monomial of
T2 and T3 (written in global Jordan-Wigner Majoranas by effective.monomials) is expanded this way;
the uniform terms with two R and two L fields give
    H_marg = sum_{X,Y} N_{XY} J^X_R J^Y_L,   J^{ab} = i psi_a psi_b,
    J^{+-}_i = (J^{0i} +- (1/2) eps_{ijk} J^{jk}) / 2.
The U(1) charge of the cell is J^{01} + J^{23} = 2 J^+_1, so su(2)+ is the "charge" (eta) algebra
and su(2)- the spin algebra.  Up to one positive normalization (common to all couplings) N is
compared with the half-filled Hubbard chain, U sum (n_up - 1/2)(n_dn - 1/2) = -(U/4) sum c0 c1 c2 c3,
whose physics is known (U > 0: spin coupling marginally irrelevant, charge coupling marginally
relevant, Mott gap), which fixes the sign convention for "relevant".

Written by Claude (Opus 5.5) for contrib/encoded_chain_numerics, 2026-10-27.
"""
import itertools
import numpy as np
import cells as C
import effective as EF


def perm_sign(seq):
    s = 1
    seq = list(seq)
    for i in range(len(seq)):
        for j in range(i + 1, len(seq)):
            if seq[i] > seq[j]:
                s = -s
    return s


def project(monos, r=4):
    """monos: list of (((a, j), ...), coef) quartic monomials (positions j = cell offsets).
    Returns G[a,b,c,d] = coefficient of psi_aR psi_bR psi_cL psi_dL (antisymmetric in ab and cd)."""
    G = np.zeros((r, r, r, r))
    for mono, cf in monos:
        if len(mono) != 4:
            continue
        cf = float(np.real(cf))
        for X in itertools.combinations(range(4), 2):          # positions that take the L field
            sgn = (-1) ** sum(mono[i][1] for i in X)
            Rpos = [i for i in range(4) if i not in X]
            order = Rpos + list(X)                               # reorder to R R L L
            ps = perm_sign(order)
            a, b = mono[Rpos[0]][0], mono[Rpos[1]][0]
            c, d = mono[X[0]][0], mono[X[1]][0]
            if a == b or c == d:
                continue                                          # psi^2 = const: not a four-fermion term
            G[a, b, c, d] += cf * sgn * ps
    # antisymmetrize
    Ga = (G - G.transpose(1, 0, 2, 3) - G.transpose(0, 1, 3, 2) + G.transpose(1, 0, 3, 2)) / 4
    return Ga


def couplings(Ga):
    """N_{XY} with H = sum_abcd Ga_abcd psi_aR psi_bR psi_cL psi_dL = sum_XY N_XY J^X_R J^Y_L.
    J^{+-}_1 = (J^{01} +- J^{23})/2, J^{+-}_2 = (J^{02} -+ J^{13})/2, J^{+-}_3 = (J^{03} +- J^{12})/2."""
    pairs = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    pi = {pq: k for k, pq in enumerate(pairs)}
    A = np.zeros((6, 6))                    # J^X = sum_p A[X, p] J^p
    for X, sgn in ((0, 1), (3, -1)):
        A[X, pi[(0, 1)]] = 0.5; A[X, pi[(2, 3)]] = 0.5 * sgn
        A[X + 1, pi[(0, 2)]] = 0.5; A[X + 1, pi[(1, 3)]] = -0.5 * sgn
        A[X + 2, pi[(0, 3)]] = 0.5; A[X + 2, pi[(1, 2)]] = 0.5 * sgn
    B = np.linalg.inv(A)                    # J^p = sum_X B[p, X] J^X
    Mp = np.zeros((6, 6))
    for p, (a, b) in enumerate(pairs):
        for q_, (c, d) in enumerate(pairs):
            Mp[p, q_] = -4 * Ga[a, b, c, d]  # psi_a psi_b = -i J^{ab}; 4 orderings
    return B.T @ Mp @ B


def describe(N, label):
    Npp, Nmm, Npm = N[:3, :3], N[3:, 3:], N[:3, 3:]
    print(f'{label}:')
    print('   N(+,+) [charge su(2), axis 1 = U(1)] =', np.round(np.diag(Npp), 6), ' offdiag max', np.round(np.abs(Npp - np.diag(np.diag(Npp))).max(), 8))
    print('   N(-,-) [spin su(2)]                  =', np.round(np.diag(Nmm), 6), ' offdiag max', np.round(np.abs(Nmm - np.diag(np.diag(Nmm))).max(), 8))
    print('   N(+,-) max |.|                        =', np.round(np.abs(Npm).max(), 8))


if __name__ == '__main__':
    # calibration: Hubbard U > 0, -(U/4) c0 c1 c2 c3 on one site (U = 1)
    hub = [(((0, 0), (1, 0), (2, 0), (3, 0)), -0.25)]
    describe(couplings(project(hub)), 'Hubbard U = +1 (reference: spin irrelevant, charge relevant)')
    for g in (0.2, 0.5):
        c = C.cell(3, 1, 2, 1.0, g)
        so = EF.second_order(c)
        m2 = EF.monomials(so['T2'], c, 2)
        m3 = EF.monomials(so['T3'], c, 3)
        bil2 = [(m, cf) for m, cf in m2 if len(m) == 2]
        bil3 = [(m, cf) for m, cf in m3 if len(m) == 2]
        print(f'\nq=3 g={g}: alpha^2 = {so["alpha"]**2:.5f}; T2 bilinears {[(m, np.round(cf, 6)) for m, cf in bil2]}')
        print(f'   T3 bilinears {[(m, np.round(cf, 6)) for m, cf in bil3]}; T3 quartic count {sum(len(m) == 4 for m, _ in m3)}, '
              f'higher {sum(len(m) > 4 for m, _ in m3)}')
        N2 = couplings(project(m2)); N3 = couplings(project(m3))
        describe(N2, '   from T2 (per lam^2)')
        describe(N3, '   from T3 (per lam^2)')
        describe(N2 + N3, '   total (per lam^2)')
