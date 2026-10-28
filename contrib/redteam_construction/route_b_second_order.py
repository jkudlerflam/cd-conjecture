"""Second-order effective chain of the route-B "dressed Bell-pair" joint cell.

Red-team check (independent implementation, not a rerun of majorana_bank_check.py).

Cell = (left port x left partner) x (right port x right partner) x logical register of n
qubits with Clifford generators gamma_1..gamma_2n, parity P, eta_a = i gamma_a P.
Port channels F_a, a = 1..r = 2 q+ q-, are the odd (graded) Hermitian matrix units,
tr(F_a F_b) = q delta_ab, and r = 2n so every logical Majorana is used.

  H_cell = Delta(1 - Pi_Omega)_L + Delta(1 - Pi_Omega)_R
           + g sum_a (|Omega_a><Omega| + h.c.)_L x gamma_a
           + g sum_a (|Omega_a><Omega| + h.c.)_R x eta_a

Chain: identical cells, bond lambda sum_a F^R_a(j) F^L_a(j+1).

What is computed
 1. Exact 2^n degeneracy of the cell code and the compressions alpha gamma_b, alpha eta_b.
 2. The complete second-order effective Hamiltonian per unit cell (same-bond two-cell
    terms and adjacent-bond three-cell terms), expanded in global Jordan-Wigner Majoranas.
 3. A check of the whole second-order formula, adjacent-bond term included, against exact
    diagonalization of an open chain of three cells (q = 2 only, dimension 32^3).
    Optional (--ed): the same-bond term against two q = 3 cells (slow).
 4. The quadratic part: the Bloch matrix at k = 0 and k = pi, i.e. whether second order
    splits the 2n Majorana Fermi points (it cannot gap them, by chirality).
 5. The quartic part projected on the continuum current-current couplings
    G_{(ab),(cd)} J_R^{ab} J_L^{cd}, with J^{ab} = i psi_a psi_b and
    c_{j,a} ~ psi_R,a + (-1)^j psi_L,a. Reports the singlet Gross-Neveu coefficient
    g_s (H contains g_s sum_{a<b} eps_a eps_b, eps_a = i psi_R,a psi_L,a; g_s < 0 is the
    attractive, marginally relevant sign in the convention of codex/work_C8 Sec. 11).
 6. One-loop flows dG/dl = B(G,G): at the unshifted point, and after the second-order
    chemical potential has moved the Fermi points, keeping only zero-momentum channels
    (resonant_flow), with and without symmetry-breaking noise of relative size 1e-2.

Run:  OPENBLAS_NUM_THREADS=1 python3 route_b_second_order.py [--big] [--ed]
About 3 minutes on one core; --big adds q=4, q+=2, r=8 (cell dimension 4096, ~10 min).
"""
import sys
import itertools
import numpy as np
from functools import reduce

X = np.array([[0, 1], [1, 0]], complex)
Y = np.array([[0, -1j], [1j, 0]], complex)
Z = np.diag([1.0, -1.0]).astype(complex)
kron = lambda *o: reduce(np.kron, o)


def majoranas(n):
    g = []
    for j in range(n):
        g.append(kron(*([Z] * j + [X] + [np.eye(2)] * (n - j - 1))))
        g.append(kron(*([Z] * j + [Y] + [np.eye(2)] * (n - j - 1))))
    return g, kron(*([Z] * n))


def graded_channels(q, qp):
    F = []
    for j in range(qp):
        for k in range(qp, q):
            A = np.zeros((q, q), complex); A[j, k] = A[k, j] = 1; F.append(A * np.sqrt(q / 2))
            B = np.zeros((q, q), complex); B[j, k] = -1j; B[k, j] = 1j; F.append(B * np.sqrt(q / 2))
    return F


class Cell:
    def __init__(self, q, qp, g=0.2, Delta=1.0):
        self.q, self.qp = q, qp
        F = graded_channels(q, qp)
        self.F, self.r = F, len(F)
        n = self.r // 2
        self.n, self.N = n, 2 ** n
        gam, P = majoranas(n)
        self.gam, self.P = gam, P
        self.eta = [1j * x @ P for x in gam]
        Dp, N = q * q, self.N
        self.Dp, self.D = Dp, Dp * Dp * N
        om = np.eye(q).reshape(Dp) / np.sqrt(q)
        self.om = om
        S = [np.kron(f, np.eye(q)) @ om for f in F]
        Pom = np.outer(om, om.conj())

        def side(ops):
            H = Delta * np.kron(np.eye(Dp) - Pom, np.eye(N))
            for a in range(self.r):
                K = np.outer(S[a], om.conj())
                H = H + g * np.kron(K + K.conj().T, ops[a])
            return H.reshape(Dp, N, Dp, N)
        HL, HR = side(self.gam), side(self.eta)
        H = (np.einsum('aibj,cd->acibdj', HL, np.eye(Dp))
             + np.einsum('cidj,ab->acibdj', HR, np.eye(Dp))).reshape(self.D, self.D)
        self.H = (H + H.conj().T) / 2
        import scipy.linalg as sla
        E, V = sla.eigh(self.H, driver='evr')    # numpy's zheevd fails at D = 4096 here
        self.E, self.V = E, V
        G = V[:, :N]
        self.split = E[N - 1] - E[0]
        self.gap = E[N] - E[0]
        # fix the logical frame: the |Omega>_L|Omega>_R component must be prop. to identity
        blk = np.einsum('a,b,abik->ik', om.conj(), om.conj(), G.reshape(Dp, Dp, N, N))
        U, s, Wh = np.linalg.svd(blk)
        self.G = G @ (U @ Wh).conj().T
        self.eps = E - E[0]

    # port operators applied to a D x k matrix
    def FL(self, a, M):
        q, Dp, N = self.q, self.Dp, self.N
        T = M.reshape(q, q, Dp, N, -1)
        return np.einsum('xy,yabik->xabik', self.F[a], T).reshape(self.D, -1)

    def FR(self, a, M):
        q, Dp, N = self.q, self.Dp, self.N
        T = M.reshape(Dp, q, q, N, -1)
        return np.einsum('xy,ayzik->axzik', self.F[a], T).reshape(self.D, -1)

    def second_order_data(self, tol=1e-7):
        G, V, eps, N, r = self.G, self.V, self.eps, self.N, self.r
        AL = [V.conj().T @ self.FL(a, G) for a in range(r)]
        AR = [V.conj().T @ self.FR(a, G) for a in range(r)]
        CL = [G.conj().T @ self.FL(a, G) for a in range(r)]
        CR = [G.conj().T @ self.FR(a, G) for a in range(r)]
        alpha = np.trace(self.gam[0] @ CL[0]).real / N
        resL = max(np.linalg.norm(CL[a] - alpha * self.gam[a]) for a in range(r))
        resR = max(np.linalg.norm(CR[a] - alpha * self.eta[a]) for a in range(r))
        # group eigenvalues into levels; level 0 = code
        order = np.argsort(eps)
        levels, cur = [], [order[0]]
        for k in order[1:]:
            if abs(eps[k] - eps[cur[-1]]) < tol:
                cur.append(k)
            else:
                levels.append(cur); cur = [k]
        levels.append(cur)
        assert len(levels[0]) == N, 'code level is not exactly N-fold'
        lev_e = np.array([eps[l].mean() for l in levels])
        lev_e[0] = 0.0
        K = len(levels)
        XR = np.zeros((K, r, r, N, N), complex)
        XL = np.zeros((K, r, r, N, N), complex)
        for e, idx in enumerate(levels):
            for a in range(r):
                for b in range(r):
                    XR[e, a, b] = AR[a][idx].conj().T @ AR[b][idx]
                    XL[e, a, b] = AL[a][idx].conj().T @ AL[b][idx]
        Cm = np.zeros((K, K))
        for e in range(K):
            for f in range(K):
                if e or f:
                    Cm[e, f] = 1.0 / (lev_e[e] + lev_e[f])
        # same-bond term on (left member, right member), N^2 x N^2
        S2 = np.zeros((N * N, N * N), complex)
        for a in range(r):
            for b in range(r):
                XLw = np.einsum('ef,fij->eij', Cm, XL[:, a, b])
                S2 -= np.einsum('eij,ekl->ikjl', XR[:, a, b], XLw).reshape(N * N, N * N)
        # adjacent-bond resolvents on the middle cell
        exc = np.concatenate(levels[1:])
        w = 1.0 / eps[exc]
        YRL = np.zeros((r, r, N, N), complex)   # [b,a] = G^ F^R_b R F^L_a G
        YLR = np.zeros((r, r, N, N), complex)   # [b,a] = G^ F^L_b R F^R_a G
        for a in range(r):
            for b in range(r):
                YRL[b, a] = (AR[b][exc].conj().T * w) @ AL[a][exc]
                YLR[b, a] = (AL[b][exc].conj().T * w) @ AR[a][exc]
        self.S2, self.YRL, self.YLR, self.CL, self.CR = S2, YRL, YLR, CL, CR
        self.alpha, self.resL, self.resR, self.nlev = alpha, resL, resR, K
        return self


# ---------------------------------------------------------------- Majorana bookkeeping
def spin_form(cell, mono, ncell=3):
    """Per-cell factors of the global JW monomial c_{x1} c_{x2} ... (x = (j, a))."""
    N = cell.N
    fac = [np.eye(N, dtype=complex) for _ in range(ncell)]
    for (j, a) in mono:
        for k in range(j):
            fac[k] = fac[k] @ cell.P
        fac[j] = fac[j] @ cell.gam[a]
    return fac


def coefficient(cell, mono):
    """Coefficient of the monomial in the per-unit-cell second-order Hamiltonian
    (lambda^2 stripped), using the window logic described in the README."""
    N, r = cell.N, cell.r
    cells = {j for (j, a) in mono}
    f0, f1, f2 = spin_form(cell, mono)
    tr = lambda A, B: np.trace(A.conj().T @ B) / N
    c = 0.0
    if cells <= {1}:          # on-site: both bonds touching cell 1
        c += np.trace(np.kron(f0, f1).conj().T @ cell.S2) / N ** 2 * tr(f2, np.eye(N))
        c += np.trace(np.kron(f1, f2).conj().T @ cell.S2) / N ** 2 * tr(f0, np.eye(N))
        return c
    if cells == {0, 1}:       # one bond
        return np.trace(np.kron(f0, f1).conj().T @ cell.S2) / N ** 2 * tr(f2, np.eye(N))
    if 0 in cells and 2 in cells:   # adjacent-bond term
        t0R = np.array([tr(f0, cell.CR[a]) for a in range(r)])
        t2L = np.array([tr(f2, cell.CL[a]) for a in range(r)])
        for a in range(r):
            for b in range(r):
                c -= t0R[a] * tr(f1, cell.YRL[b, a]) * t2L[b]
                c -= t0R[b] * tr(f1, cell.YLR[b, a]) * t2L[a]
        return c
    return 0.0               # supports {1,2} are translates of {0,1}; avoid double count


def unit_cell_terms(cell, degree):
    """All monomials of given degree in the per-unit-cell set, with coefficients."""
    m = 2 * cell.n
    sites = [(j, a) for j in range(3) for a in range(m)]
    out = {}
    for mono in itertools.combinations(sites, degree):
        cells = {j for (j, a) in mono}
        if not (cells <= {1} or cells == {0, 1} or (0 in cells and 2 in cells)):
            continue
        c = coefficient(cell, mono)
        if abs(c) > 1e-13:
            out[mono] = c
    return out


def bloch_zero_pi(quad, m):
    """h(k) = i sum_delta A(delta) e^{i k delta} at k = 0, pi from quadratic terms
    kappa c_x c_y (kappa imaginary), A_xy = -2 i kappa, antisymmetric."""
    h0 = np.zeros((m, m), complex); hpi = np.zeros((m, m), complex)
    for ((j1, a), (j2, b)), kap in quad.items():
        A = -2j * kap
        d = j2 - j1
        for (u, v, val, dd) in ((a, b, A, d), (b, a, -A, -d)):
            h0[u, v] += 1j * val
            hpi[u, v] += 1j * val * (-1) ** dd
    return h0, hpi


def current_couplings(quart, m):
    """Project quartic Majorana terms onto G_{(ab),(cd)} J_R^{ab} J_L^{cd}."""
    pairs = list(itertools.combinations(range(m), 2))
    idx = {p: i for i, p in enumerate(pairs)}
    Gm = np.zeros((len(pairs), len(pairs)))
    for mono, kap in quart.items():
        assert abs(kap.imag) < 1e-10
        kap = kap.real
        for Lset in itertools.combinations(range(4), 2):
            sgn = (-1) ** sum(mono[i][0] for i in Lset)
            labels = [('L' if i in Lset else 'R', mono[i][1]) for i in range(4)]
            # reorder to R R L L (stable), Grassmann sign from inversions
            perm = sorted(range(4), key=lambda i: (labels[i][0] == 'L', i))
            inv = sum(1 for x in range(4) for y in range(x + 1, 4) if perm[x] > perm[y])
            s = sgn * (-1) ** inv
            (_, u), (_, v), (_, s1), (_, t1) = [labels[i] for i in perm]
            if u == v or s1 == t1:
                continue
            if u > v: u, v = v, u; s = -s
            if s1 > t1: s1, t1 = t1, s1; s = -s
            # psi_Ru psi_Rv psi_Ls psi_Lt = -J_R^{uv} J_L^{st}
            Gm[idx[(u, v)], idx[(s1, t1)]] += -kap * s
    return Gm, pairs


def so_structure(m):
    pairs = list(itertools.combinations(range(m), 2))
    L = []
    for (a, b) in pairs:
        M = np.zeros((m, m), complex); M[a, b] = -1j; M[b, a] = 1j; L.append(M)
    nA = len(L)
    f = np.zeros((nA, nA, nA))
    for A in range(nA):
        for B in range(nA):
            Cm = L[A] @ L[B] - L[B] @ L[A]     # = i f_ABC L_C, tr(L_C L_D) = 2 delta
            for C in range(nA):
                f[A, B, C] = (np.trace(Cm @ L[C]) / (2j)).real
    return f


def one_loop(G0, f, h=0.002, lmax=300.0, maxsteps=200000, proj=None):
    """dG_EF/dl = sum f_ACE f_BDF G_AB G_CD (Cardy, OPE of J_R J_L with itself). Sign check:
    G = +|g| 1 is H = -|g| sum_{a<b} eps_a eps_b (attractive GN) and grows, as it should.
    G0 is normalized to 1, so l is in units of 1/|G0|. Adaptive step dl = h/|G|."""
    def B(G):
        T = np.einsum('ab,bdf,cd->acf', G, f, G, optimize=True)
        return np.einsum('ace,acf->ef', f, T, optimize=True)
    G = G0 / np.linalg.norm(G0)
    ell = 0.0
    for s in range(maxsteps):
        nG = np.linalg.norm(G)
        if nG > 1e3:
            return 'runaway', ell
        if ell > lmax:
            return 'bounded to l=%g (final |G| %.3f)' % (lmax, nG), ell
        dl = min(h / nG, 0.5)
        k1 = B(G); k2 = B(G + 0.5 * dl * k1)
        G = G + dl * k2
        if proj is not None:          # drop round-off in symmetry-forbidden / off-resonant channels
            G = proj(G)
        ell += dl
    return 'undecided', ell


def resonant_flow(Gm, h0, m, seed=1):
    """The second-order on-site term Q = h2(0) = h2(pi) shifts right movers (k ~ 0) by +dk and
    left movers (k ~ pi) by -dk per unit charge. A coupling Jt_R^mu Jt_L^nu between ad_Q
    eigencurrents (charges nu) carries momentum (nu_nu - nu_mu) dk, so only nu_mu = nu_nu
    survives (2x2 check on a charged pair: J1J1 - J2J2 survives, J1J1 + J2J2 does not).
    Three one-loop flows are reported:
      (a) neutral currents only (exact U(1)_Q and exact gadget flavor symmetry),
      (b) (a) plus or minus a random neutral perturbation of relative size 1e-2,
      (c) the full resonant projection plus a random resonant perturbation (no symmetry)."""
    rng = np.random.default_rng(seed)
    pairs = list(itertools.combinations(range(m), 2))
    L = []
    for (a, b) in pairs:
        M = np.zeros((m, m), complex); M[a, b] = -1j; M[b, a] = 1j; L.append(M)
    nA = len(L)
    K = np.array([[(np.trace(L[B] @ (h0 @ L[A] - L[A] @ h0)) / 2j).real
                   for A in range(nA)] for B in range(nA)])
    nu, U = np.linalg.eig(K)
    nu = nu.imag
    scale = max(np.abs(nu).max(), 1e-30)
    Ui = np.linalg.inv(U)

    def project(G, neutral_only=False):
        Ct = U.T @ G @ U
        mask = np.abs(nu[:, None] - nu[None, :]) < 1e-6 * scale
        if neutral_only:
            mask &= (np.abs(nu[:, None]) < 1e-6 * scale)
        return (Ui.T @ (Ct * mask) @ Ui).real
    Gres = project(Gm)
    Gneu = project(Gm, True)
    f = so_structure(m)
    out = dict(nres=np.linalg.norm(Gres) / np.linalg.norm(Gm),
               nneutral=int((np.abs(nu) < 1e-6 * scale).sum()),
               charges=np.round(np.sort(nu) / scale, 4) if scale > 1e-12 else 'no shift')
    # U(1)_Q forward-scattering coefficient: Q direction in the current basis
    qv = np.array([(np.trace(L[A] @ h0) / 2).real for A in range(nA)])
    out['GQQ'] = (qv @ Gm @ qv) / (qv @ qv) if np.linalg.norm(qv) > 1e-12 else float('nan')
    R = rng.normal(size=Gm.shape); R = (R + R.T) / 2
    pn = lambda G: project(G, True)
    pr = lambda G: project(G)
    pa = one_loop(Gneu, f, proj=pn)
    Rn = project(R, True); Rn *= 1e-2 * np.linalg.norm(Gneu) / np.linalg.norm(Rn)
    pb = one_loop(Gneu + Rn, f, proj=pn)
    pb2 = one_loop(Gneu - Rn, f, proj=pn)
    Rr = project(R); Rr *= 1e-2 * np.linalg.norm(Gres) / np.linalg.norm(Rr)
    pc = one_loop(Gres + Rr, f, proj=pr)
    out['flows'] = (pa, pb, pb2, pc)
    return out


def two_cell_ed_check(cell, lams=(0.02, 0.01)):
    """Exact low band of two cells vs lambda alpha^2 T + lambda^2 S2 (same-bond)."""
    import scipy.sparse.linalg as sla
    D, N, r = cell.D, cell.N, cell.r
    Hc = cell.H - cell.E[0] * np.eye(D)
    G2 = np.kron(cell.G, cell.G)
    T = sum(np.kron(cell.CR[a], cell.CL[a]) for a in range(r)) / cell.alpha ** 2
    rows = []
    for lam in lams:
        def mv(v):
            V = v.reshape(D, D)
            out = Hc @ V + V @ Hc.T
            for a in range(r):
                FR = cell.FR(a, V)                         # right port of cell 1 (rows)
                out += lam * cell.FL(a, FR.T).T            # left port of cell 2 (cols)
            return out.reshape(-1)
        op = sla.LinearOperator((D * D, D * D), matvec=mv, dtype=complex)
        w, U = sla.eigsh(op, k=N * N + 2, which='SA', tol=1e-10, v0=G2.sum(1), ncv=4 * N * N + 8)
        o = np.argsort(w); w, U = w[o], U[:, o]
        Ub = U[:, :N * N]
        Bm = G2.conj().T @ Ub
        Ws, s, Vh = np.linalg.svd(Bm)
        W = Ws @ Vh
        Heff = W @ np.diag(w[:N * N]) @ W.conj().T
        pred = lam * cell.alpha ** 2 * T + lam ** 2 * cell.S2
        rows.append((lam, np.linalg.norm(Heff - lam * cell.alpha ** 2 * T),
                     np.linalg.norm(Heff - pred), w[N * N] - w[N * N - 1]))
    return rows


def analyse(q, qp, g=0.2, ed=False):
    cell = Cell(q, qp, g=g).second_order_data()
    m = 2 * cell.n
    print(f'\n=== q={q} q+={qp} r={cell.r} n={cell.n} cell dim={cell.D} g={g} ===')
    print(f'code splitting {cell.split:.1e}, cell gap {cell.gap:.4f}, levels {cell.nlev}, '
          f'alpha={cell.alpha:.6f}, compression residuals L {cell.resL:.1e} R {cell.resR:.1e}')
    if ed:
        print('two-cell ED: lambda, |Heff - first order|, |Heff - (first + second)|, band isolation')
        for row in two_cell_ed_check(cell):
            print('   %.4f  %.3e  %.3e  %.4f' % row)
    quad = unit_cell_terms(cell, 2)
    bad = max([abs(k.real) for k in quad.values()] + [0])
    h0, hpi = bloch_zero_pi(quad, m)
    w0 = np.linalg.eigvalsh(h0); wpi = np.linalg.eigvalsh(hpi)
    a2 = cell.alpha ** 2
    print(f'quadratic second-order terms: {len(quad)} (max real part {bad:.1e}, must be 0)')
    print('  h2(0)  eigenvalues / (4 alpha^2):', np.round(w0 / (4 * a2), 6))
    print('  h2(pi) eigenvalues / (4 alpha^2):', np.round(wpi / (4 * a2), 6))
    print('  [Fermi-point shift of each mode is delta k = omega/(lambda * 4 alpha^2) at order lambda]')
    odd = unit_cell_terms(cell, 3)
    print(f'parity-odd cubic terms: {len(odd)}')
    quart = unit_cell_terms(cell, 4)
    Gm, pairs = current_couplings(quart, m)
    gs = -np.trace(Gm) / len(pairs)
    print(f'quartic terms: {len(quart)}; singlet GN coefficient g_s = {gs:.6f} '
          f'(= {gs / cell.alpha ** 4:.4f} alpha^4); sign {"ATTRACTIVE" if gs < 0 else "repulsive"}')
    sym = np.linalg.norm(Gm - Gm.T) / max(np.linalg.norm(Gm), 1e-30)
    print(f'  |G - G^T|/|G| = {sym:.2e}; eigenvalues of sym(G): {np.round(np.linalg.eigvalsh((Gm + Gm.T) / 2), 6)}')
    if np.linalg.norm(Gm) > 1e-12:
        f = so_structure(m)
        fate, ell = one_loop(Gm, f)
        print(f'  one-loop flow at the symmetric point (all SO({m}) currents resonant): {fate} (l = {ell:.2f})')
        mult = np.round(np.abs(w0) / (4 * a2), 6)
        o = resonant_flow(Gm, h0, m)
        print(f'  |omega| values {sorted(set(mult.tolist()))}; ad_Q charges {o["charges"]}')
        print(f'  neutral currents {o["nneutral"]}; U(1)_Q forward coupling G_QQ = {o["GQQ"]:.6f}; '
              f'resonant fraction |G_res|/|G| = {o["nres"]:.3f}')
        for lab, (fa, la) in zip(('(a) neutral only, symmetric', '(b) neutral + noise 1e-2', "(b') neutral - same noise",
                                  '(c) all resonant + noise 1e-2'), o['flows']):
            print(f'    flow {lab}: {fa} (l = {la:.1f})')
    return cell, quad, quart, Gm


def three_cell_ed_check(q=2, qp=1, g=0.3, lams=(0.04, 0.02, 0.01)):
    """Open chain of three cells (two bonds): exact low band vs first + second order,
    including the adjacent-bond (three-cell) term. Small q only (dimension D^3)."""
    import scipy.sparse.linalg as sla
    cell = Cell(q, qp, g=g).second_order_data()
    D, N, r = cell.D, cell.N, cell.r
    Hc = cell.H - cell.E[0] * np.eye(D)
    I = np.eye(N)
    T01 = sum(np.kron(np.kron(cell.CR[a], cell.CL[a]), I) for a in range(r))
    T12 = sum(np.kron(I, np.kron(cell.CR[a], cell.CL[a])) for a in range(r))
    Xc = np.zeros((N ** 3, N ** 3), complex)
    for a in range(r):
        for b in range(r):
            Xc -= kron(cell.CR[a], cell.YRL[b, a], cell.CL[b])
            Xc -= kron(cell.CR[b], cell.YLR[b, a], cell.CL[a])
    H2 = np.kron(cell.S2, I) + np.kron(I, cell.S2) + Xc
    G3 = kron(cell.G, cell.G, cell.G)
    out = []
    for lam in lams:
        def mv(v):
            V = v.reshape(D, D, D)
            o = (np.einsum('ij,jkl->ikl', Hc, V) + np.einsum('kj,ijl->ikl', Hc, V)
                 + np.einsum('lj,ikj->ikl', Hc, V))
            for a in range(r):
                # bond (1,2): right port of cell 1, left port of cell 2
                W = cell.FR(a, V.reshape(D, -1)).reshape(D, D, D)
                W = np.moveaxis(cell.FL(a, np.moveaxis(W, 1, 0).reshape(D, -1)).reshape(D, D, D), 0, 1)
                o += lam * W
                W = cell.FR(a, np.moveaxis(V, 1, 0).reshape(D, -1)).reshape(D, D, D)
                W = np.moveaxis(W, 0, 1)
                W = np.moveaxis(cell.FL(a, np.moveaxis(W, 2, 0).reshape(D, -1)).reshape(D, D, D), 0, 2)
                o += lam * W
            return o.reshape(-1)
        op = sla.LinearOperator((D ** 3, D ** 3), matvec=mv, dtype=complex)
        w, U = sla.eigsh(op, k=N ** 3 + 2, which='SA', tol=1e-11, v0=G3.sum(1))
        o_ = np.argsort(w); w, U = w[o_], U[:, o_]
        Bm = G3.conj().T @ U[:, :N ** 3]
        Ws, s_, Vh = np.linalg.svd(Bm)
        W_ = Ws @ Vh
        Heff = W_ @ np.diag(w[:N ** 3]) @ W_.conj().T
        H1 = lam * (T01 + T12)
        out.append((lam, np.linalg.norm(Heff - H1), np.linalg.norm(Heff - H1 - lam ** 2 * H2),
                    np.linalg.norm(Heff - H1 - lam ** 2 * (np.kron(cell.S2, I) + np.kron(I, cell.S2)))))
    print(f'three-cell ED check (q={q}, q+={qp}, g={g}): lambda, |Heff-H1|, |Heff-H1-H2|, |Heff-H1-H2 without adjacent term|')
    for row in out:
        print('   %.3f  %.3e  %.3e  %.3e' % row)
    return out


if __name__ == '__main__':
    three_cell_ed_check()
    analyse(3, 1, g=0.2, ed='--ed' in sys.argv)
    analyse(3, 1, g=0.5)
    analyse(4, 1, g=0.2)
    if '--big' in sys.argv:
        analyse(4, 2, g=0.2)
