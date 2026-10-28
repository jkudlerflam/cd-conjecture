#!/usr/bin/env python3
"""A parity-graded route-A cell whose second-order logical bond is n-flavour fermion hopping.

Physical cell (one big site): two qubit ports L, R and an auxiliary register,
    K = C^2_L (x) C^2_R (x) C^A,   A = 4N + extra,   N = 2^n,   D = 4A.
The first 4N auxiliary states are (partner of L) (x) (partner of R) (x) (logical C^N).  Code
    G|i> = sum_{x,y} |x>_L |y>_R |x, y, i>_aux / 2       (both ports in Bell pairs with partners)
so G^dag s^L_a G = G^dag s^R_a G = 0 for every Pauli and the six syndrome families s^{L,R}_a G are
orthonormal (C9, codex/work_C9/RESULT.md Sec. 2.1).

Grading.  U = Z_L Z_R Z_pL Z_pR P_log (and +-1 on the extra states) commutes with H_A and with the
bond lam sum_a s^R_a(j) s^L_a(j+1), and acts on the code as the logical fermion parity P.  X and Y
are odd, Z is even.  Because of this, the second-order bond
    h2 = -(lam^2/Delta) sum_{e,f} sum_{a,b} M^R_e[a,b] (x) M^L_f[a,b] / (e_e + e_f),
    M_e[a,b] = G^dag s_a P_e s_b G,
is odd (x) odd only through the mixed blocks [X,Z], [Y,Z]: hopping must come from excited states
that mix an odd-channel syndrome with the even one.  C9's block-diagonal recipe cannot do this, so
here the weights are
    M^R_e = p_e 1 + rho_R sum_c z^c_e A_c,     M^L_e = p_e 1 + rho_L sum_c z^c_e B_c,
with A_c = sum_{flavours alpha in column c} |ch_alpha><Z| (x) e^{i phi_alpha} a_alpha^dag P + h.c.
(+ |Z><Z| (x) D_mu in column 0) and B_c = sum_alpha |ch_alpha><Z| (x) t_alpha a_alpha + h.c.;
z^c in R^k are orthogonal to (1,...,1) and orthonormal in the metric C_ef = 1/(e_e + e_f).
Flavours sharing a column use different odd channels (X or Y), so no cross-flavour hopping appears.
Then exactly
    h2 = (lam^2/Delta) [ -3 p^T C p - rho_R (z^0.C p) D_mu (x) 1
                         - J0 sum_alpha t_alpha (e^{i phi_alpha} a_alpha^dag P (x) a_alpha + h.c.) ],
    J0 = rho_R rho_L,
i.e. (Jordan-Wigner across cells: c^dag_{j,a} c_{j+1,a} = (a^dag_a P)_j (x) (a_a)_{j+1})
    H_log = sum_j [ -J0 sum_a t_a (e^{i phi_a} c^dag_{j,a} c_{j+1,a} + h.c.) - J0 sum_a mu_a n_{j,a} ].
The weights are realized by an explicit graded Naimark dilation (C9 Sec. 2.3) into the
complement of the code and syndromes, with random (seeded) partner vectors in each parity sector,
so that everything beyond second order is generic up to the exact Z_2.

Conventions match contrib/ir_control/toy_ir.py: target -t (e^{-i th} c^dag_j c_{j+1} + h.c.) - mu n,
i.e. phi = -th, Fermi points th +- k_F with cos k_F = -mu/(2t).

Written by Claude (Opus 5.5) for contrib/routeA_incommensurate, 2026-10-27.
"""
import numpy as np
from scipy import linalg

X = np.array([[0, 1], [1, 0]], complex)
Y = np.array([[0, -1j], [1j, 0]], complex)
Z = np.diag([1.0, -1.0]).astype(complex)
I2 = np.eye(2, dtype=complex)
PAULI = [X, Y, Z]

INC = dict(t=[1.0, 1.0], th=[-2.617, 0.298], kF=[0.962, 0.742])          # ir_control target, n = 2
SYM = dict(t=[1.0, 1.0], th=[np.pi / 2, np.pi / 2], kF=[np.pi / 2, np.pi / 2])   # SO(4) Majorana point
INC3 = dict(t=[1.0, 1.0, 1.0], th=[-2.617, 0.298, 1.571], kF=[0.962, 0.742, 1.25])


def kron(*ops):
    out = np.eye(1, dtype=complex)
    for o in ops:
        out = np.kron(out, o)
    return out


def fermion_ops(n):
    """Annihilators a_k on C^(2^n) (intra-cell Jordan-Wigner, index 0 = empty) and parity P."""
    s = np.array([[0, 1], [0, 0]], complex)
    a = [kron(*([Z] * k + [s] + [I2] * (n - k - 1))) for k in range(n)]
    return a, kron(*([Z] * n))


def psd_sqrt(M):
    w, V = linalg.eigh(0.5 * (M + M.conj().T))
    return (V * np.sqrt(np.clip(w, 0, None))) @ V.conj().T


def target_params(par):
    t = np.asarray(par['t'], float)
    th = np.asarray(par['th'], float)
    kF = np.asarray(par['kF'], float)
    mu = -2 * t * np.cos(kF)
    return t, th, mu, kF


def fermi_data(par):
    t, th, mu, kF = target_params(par)
    pts = np.concatenate([th + kF, th - kF]) % (2 * np.pi)
    v = 2 * t * np.sin(kF)
    # delta_2: smallest distance between any two of {+-p}; delta_4: min |p_a +- p_b +- p_c +- p_d|
    allp = np.sort(np.concatenate([pts, -pts]) % (2 * np.pi))
    gaps = np.diff(np.concatenate([allp, [allp[0] + 2 * np.pi]]))
    return dict(points=pts.tolist(), v=v.tolist(), delta2=float(gaps.min()), filling=(kF / np.pi).tolist())


def target_bond(par, J=1.0):
    """The logical target as a two-cell operator on C^N (x) C^N (on-site terms on the left cell)."""
    t, th, mu, _ = target_params(par)
    n = len(t)
    a, P = fermion_ops(n)
    N = 2 ** n
    h = np.zeros((N * N, N * N), complex)
    for k in range(n):
        hop = np.kron(a[k].conj().T @ P, a[k])
        h += -J * t[k] * (np.exp(-1j * th[k]) * hop + np.exp(1j * th[k]) * hop.conj().T)
        h += -J * mu[k] * np.kron(a[k].conj().T @ a[k], np.eye(N))
    return h


def build(par, E2=5.0, Estar=2.5, frac=0.9, seed=7, columns=None, extra=0):
    """Returns a dict with H0 (D x D), FL, FR (lists of 3 port Paulis), G (D x N), U, and the
    predicted second-order data.  columns[c] = list of (flavour, channel) with channel 0 (X) or 1 (Y)."""
    t, th, mu, _ = target_params(par)
    n = len(t)
    N = 2 ** n
    a, P = fermion_ops(n)
    if columns is None:
        columns = [[(0, 0), (1, 1)]] if n == 2 else [[(0, 0), (1, 1)], [(2, 0)]]
    ncol = len(columns)
    k = ncol + 1
    e = np.array([1.0] + [E2 ** (i / (k - 1)) for i in range(1, k)])  # 1 .. E2, geometric
    C = 1.0 / (e[:, None] + e[None, :])
    p = np.full(k, 1.0 / k)
    # z columns: orthogonal to 1, C-orthonormal; for k = 2 this is z = s (1, -1)
    Bz = linalg.null_space(np.ones((1, k)))[:, :ncol]
    Zc = Bz @ linalg.inv(linalg.sqrtm(Bz.T @ C @ Bz)).real
    if Zc[0, 0] < 0:
        Zc[:, 0] *= -1
    assert np.allclose(Zc.T @ C @ Zc, np.eye(ncol)) and np.allclose(Zc.sum(0), 0)

    # syndrome-space operators (channel (x) logical), channels 0 = X, 1 = Y, 2 = Z
    def E_(i, j):
        m = np.zeros((3, 3), complex); m[i, j] = 1; return m
    Acols, Bcols = [], []
    for col in columns:
        A = np.zeros((3 * N, 3 * N), complex)
        B = np.zeros((3 * N, 3 * N), complex)
        for (f, ch) in col:
            Aop = np.exp(-1j * th[f]) * a[f].conj().T @ P
            Bop = t[f] * a[f]
            A += np.kron(E_(ch, 2), Aop); A += np.kron(E_(ch, 2), Aop).conj().T
            B += np.kron(E_(ch, 2), Bop); B += np.kron(E_(ch, 2), Bop).conj().T
        Acols.append(A); Bcols.append(B)

    def max_rho(ops):
        """largest rho with p_e 1 + rho sum_c z^c_e ops_c >= 0 for all e."""
        r = np.inf
        for ee in range(k):
            Mx = sum(Zc[ee, c] * ops[c] for c in range(ncol))
            w = linalg.eigvalsh(Mx)
            if w[0] < 0:
                r = min(r, p[ee] / (-w[0]))
        return r
    rhoL = frac * max_rho(Bcols)
    # chemical potential through the linear term of column 0 on the even channel
    zCp = float(Zc[:, 0] @ C @ p)
    nops = [a[f].conj().T @ a[f] for f in range(n)]
    Dmu = sum(mu[f] * nops[f] for f in range(n))
    dvals = np.linalg.eigvalsh(Dmu)
    shift = 0.5 * (dvals[0] + dvals[-1])
    Dmu_c = Dmu - shift * np.eye(N)
    Acols[0] = Acols[0] + np.kron(E_(2, 2), rhoL / zCp * Dmu_c)
    rhoR = frac * max_rho(Acols)
    J0 = rhoR * rhoL
    MR = [p[ee] * np.eye(3 * N) + rhoR * sum(Zc[ee, c] * Acols[c] for c in range(ncol)) for ee in range(k)]
    ML = [p[ee] * np.eye(3 * N) + rhoL * sum(Zc[ee, c] * Bcols[c] for c in range(ncol)) for ee in range(k)]
    assert np.allclose(sum(MR), np.eye(3 * N)) and np.allclose(sum(ML), np.eye(3 * N))

    # ---- physical space
    A_aux = 4 * N + extra
    D = 4 * A_aux
    # aux parity: Z_pL Z_pR P on the partner block, alternating on extra states
    Waux = np.ones(A_aux)
    Waux[:4 * N] = np.real(np.diag(kron(Z, Z, P)))
    Waux[4 * N:] = [(-1) ** i for i in range(extra)]
    U = np.kron(np.kron(Z, Z), np.diag(Waux)).astype(complex)
    FL = [np.kron(np.kron(s, I2), np.eye(A_aux)) for s in PAULI]
    FR = [np.kron(np.kron(I2, s), np.eye(A_aux)) for s in PAULI]
    G = np.zeros((D, N), complex)
    for i in range(N):
        for x in range(2):
            for y in range(2):
                aux = (x * 2 + y) * N + i
                G[(x * 2 + y) * A_aux + aux, i] = 0.5
    SL = np.hstack([F @ G for F in FL])           # columns ordered (channel, logical)
    SR = np.hstack([F @ G for F in FR])
    # syndrome grading Gamma = diag(s_ch) (x) P, s = (-1, -1, +1)
    Gam = np.kron(np.diag([-1.0, -1.0, 1.0]), P).real
    assert np.allclose(U @ SR, SR @ Gam) and np.allclose(U @ SL, SL @ Gam) and np.allclose(U @ G, G @ P)

    # complement of span(G, SL, SR), split by U parity, random orthonormal bases (seeded)
    rng = np.random.default_rng(seed)
    Bocc = np.hstack([G, SL, SR])
    Qc = np.eye(D) - Bocc @ Bocc.conj().T
    udiag = np.real(np.diag(U))
    comp = {}
    for sgn in (1, -1):
        Ps = np.diag((udiag == sgn).astype(float))
        M = Ps @ Qc @ Ps
        w, V = linalg.eigh(M)
        Vs = V[:, w > 0.5]
        Rm = rng.normal(size=(Vs.shape[1], Vs.shape[1])) + 1j * rng.normal(size=(Vs.shape[1], Vs.shape[1]))
        Qr, _ = linalg.qr(Rm)
        comp[sgn] = list((Vs @ Qr).T)

    def dilate(Ms, S):
        """Graded Naimark dilation: returns the spectral projectors P_e (D x D) with S^dag P_e S = M_e."""
        m = S.shape[1]
        Wd = np.vstack([psd_sqrt(Mx) for Mx in Ms])               # (k m) x m isometry
        assert np.allclose(Wd.conj().T @ Wd, np.eye(m), atol=1e-10)
        Gam2 = np.kron(np.eye(k), Gam)
        K = np.eye(k * m) - Wd @ Wd.conj().T
        cols_dil, cols_phys = [Wd], [S]
        for sgn in (1, -1):
            Ps = np.diag((np.diag(Gam2) == sgn).astype(float))
            w, V = linalg.eigh(Ps @ K @ Ps)
            Vs = V[:, w > 0.5]
            cols_dil.append(Vs)
            cols_phys.append(np.array([comp[sgn].pop(0) for _ in range(Vs.shape[1])]).T)
        Ud = np.hstack(cols_dil)
        Tp = np.hstack(cols_phys)
        assert np.allclose(Ud.conj().T @ Ud, np.eye(k * m), atol=1e-10)
        J = Tp @ Ud.conj().T                                      # J Wd = S, J isometric
        Pe = []
        for ee in range(k):
            Pi = np.zeros((k * m, k * m)); Pi[ee * m:(ee + 1) * m, ee * m:(ee + 1) * m] = np.eye(m)
            Pe.append(J @ Pi @ J.conj().T)
        return Pe
    PeR = dilate(MR, SR)
    PeL = dilate(ML, SL)
    H0 = np.zeros((D, D), complex)
    used = np.zeros((D, D), complex)
    for ee in range(k):
        H0 += e[ee] * (PeR[ee] + PeL[ee])
        used += PeR[ee] + PeL[ee]
    rest = np.eye(D) - G @ G.conj().T - used
    H0 += Estar * rest
    H0 = 0.5 * (H0 + H0.conj().T)
    pred = dict(const=-3 * float(p @ C @ p), J0=float(J0), rhoR=float(rhoR), rhoL=float(rhoL),
                mu_shift=float(shift), e=e.tolist(), Estar=Estar, zCp=zCp)
    return dict(par=par, n=n, N=N, D=D, H0=H0, FL=FL, FR=FR, G=G, U=U, P=P, a=a, pred=pred,
                MR=MR, ML=ML, Zc=Zc, C=C, p=p)


def second_order_direct(cell):
    """h2 / (lam^2/Delta) from the spectral decomposition of H0 alone (independent of the recipe)."""
    H0, G, FL, FR, N = cell['H0'], cell['G'], cell['FL'], cell['FR'], cell['N']
    w, V = linalg.eigh(H0)
    exc = w > 1e-9
    w, V = w[exc], V[:, exc]
    # matrix elements <i| G^dag F_a |m> for excited eigenvectors m
    R = [G.conj().T @ F @ V for F in FR]          # N x M
    L = [G.conj().T @ F @ V for F in FL]
    den = 1.0 / (w[:, None] + w[None, :])         # M x M
    h = np.zeros((N * N, N * N), complex)
    for aa in range(3):
        for bb in range(3):
            # sum_{m,m'} R_a[i,m] conj(R_b[k,m]) L_a[i',m'] conj(L_b[k',m']) / (w_m + w_m')
            Ra = R[aa]; Rb = R[bb].conj(); La = L[aa]; Lb = L[bb].conj()
            T = np.einsum('im,km,mn,jn,ln->ijkl', Ra, Rb, den, La, Lb, optimize=True)
            h -= T.reshape(N * N, N * N)
    return h


def expected_h2(cell):
    """Predicted h2 / (lam^2/Delta): constant + J0 x target (with the chemical-potential shift)."""
    pr = cell['pred']
    N = cell['N']
    t, th, mu, _ = target_params(cell['par'])
    h = pr['const'] * np.eye(N * N) + target_bond(cell['par'], J=pr['J0'])
    # the chemical potential was centred: D_mu - shift  ->  + J0 shift on the left cell
    h += pr['J0'] * pr['mu_shift'] * np.eye(N * N)
    return h


def checks(cell, verbose=True):
    H0, G, U, FL, FR, N, D = cell['H0'], cell['G'], cell['U'], cell['FL'], cell['FR'], cell['N'], cell['D']
    w = linalg.eigvalsh(H0)
    out = dict(D=D, N=N,
               herm=float(np.abs(H0 - H0.conj().T).max()),
               code=float(np.linalg.norm(H0 @ G)),
               ground_deg=int(np.sum(np.abs(w) < 1e-9)),
               gap=float(w[np.abs(w) > 1e-9].min()),
               top=float(w.max()),
               levels=sorted(set(np.round(w, 9).tolist())),
               sym_U=float(np.linalg.norm(U @ H0 - H0 @ U)),
               # U F U = s F with the same sign s on both ports, so every bond term is U (x) U even
               bond_U=float(max(min(np.linalg.norm(U @ fr @ U - s * fr) + np.linalg.norm(U @ fl @ U - s * fl)
                                    for s in (1, -1)) for fr, fl in zip(FR, FL))),
               UG_eq_GP=float(np.linalg.norm(U @ G - G @ cell['P'])),
               first_order=float(max(np.linalg.norm(G.conj().T @ F @ G) for F in FL + FR)),
               weights_min_eig=float(min(min(linalg.eigvalsh(M)[0] for M in cell['MR']),
                                         min(linalg.eigvalsh(M)[0] for M in cell['ML']))))
    h2 = second_order_direct(cell)
    h2e = expected_h2(cell)
    out['h2_vs_target'] = float(np.linalg.norm(h2 - h2e))
    out['h2_herm'] = float(np.linalg.norm(h2 - h2.conj().T))
    out.update(cell['pred'])
    if verbose:
        for kk, v in out.items():
            print(f'  {kk}: {v}')
    return out, h2


if __name__ == '__main__':
    import json, sys
    for name, par in (('INC', INC), ('SYM', SYM)):
        print(name, fermi_data(par))
        c = build(par)
        o, h2 = checks(c)
    if '--n3' in sys.argv:
        print('INC3', fermi_data(INC3))
        c = build(INC3, extra=6)
        checks(c)
