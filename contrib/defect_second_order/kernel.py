"""Exact tagged contact resolvent kernel chi_N(r) of the free vacancy fluid on a ring.

Gap representation (C94, C96): N material sites in cyclic order, V vacancies, L = N + V,
hopping -t, A_j = [g_j = 0]. chi_N(r) = <0|(A_0 - a) Q (H - E_0)^{-1} Q (A_r - a)|0>.

Method. The material label shift T commutes with H. Its sector m (theta_m = 2 pi m / N) is N free
fermions on the physical L-ring with momenta p = 2 pi (n + xi + m/N)/L (xi = 0 for N odd, 1/2 for
N even) and total momentum 0 mod 2 pi (C96). In the frame anchored at material 0, every sector-m
eigenfunction is sqrt(L/N) times a Slater determinant restricted to configurations containing site 0,
and A_0 = n_0 n_1. Hence, with g_s = U_m diag(exp(i p s - tau eps_p)) U_m^dagger,

  C_m(tau) = <A_0 P_m e^{-tau(H-E_0)} P_m A_0>
           = (L/N)^2 (1/L) sum_{s=0}^{L-1} e^{tau E_0} <FS| n_0 n_1 Gamma(g_s) n_0 n_1 |FS>,

where the s sum projects on total momentum 0. n_0 n_1 |FS> = c |Slater(Psi_R)> with Psi_R =
(delta_0, delta_1, orbitals of the Fermi sea vanishing at sites 0 and 1), |c|^2 = <n_0 n_1>, so each
term is one N x N determinant det(Psi_R^dagger g_s Psi_R). The s dependence is done by FFT and the
determinant is stabilized by factoring out the N lowest modes of grid m. Then

  chi_N(r) = (1/N)... = sum_m e^{i theta_m r} int_0^inf dtau [C_m(tau) - delta_{m0} a^2].

numpy only.
"""
import numpy as np


def contact_mean(N, V):
    L = N + V
    return N / L - np.sin(np.pi * N / L) ** 2 / (N * L * np.sin(np.pi / L) ** 2)


def setup(N, V, t=1.0):
    L = N + V
    xi = 0.0 if N % 2 else 0.5
    x = np.arange(L)
    n = np.arange(L)
    p0 = 2 * np.pi * (n + xi) / L
    e0 = -2 * t * np.cos(p0)
    occ = np.argsort(e0, kind='stable')[:N]
    E0 = e0[occ].sum()
    Phi = np.exp(1j * np.outer(x, p0[occ])) / np.sqrt(L)
    R2 = Phi[:2, :]
    U, S, Vh = np.linalg.svd(R2)
    PV = Phi @ Vh.conj().T
    PsiR = np.zeros((L, N), complex)
    PsiR[0, 0] = 1.0; PsiR[1, 1] = 1.0
    PsiR[:, 2:] = PV[:, 2:]
    pref2 = (S[0] * S[1]) ** 2
    return dict(L=L, xi=xi, E0=E0, PsiR=PsiR, pref2=pref2, x=x, n=n, t=t)


def sector_correlator(st, N, m, taus):
    """C_m(tau) for an array of tau."""
    L = st['L']; x = st['x']; n = st['n']; t = st['t']
    Phim = st['xi'] + m / N
    q = 2 * np.pi * (n + Phim) / L
    eq = -2 * t * np.cos(q)
    Um = np.exp(1j * np.outer(x, q)) / np.sqrt(L)
    O = st['PsiR'].conj().T @ Um                       # N x L
    order = np.argsort(eq, kind='stable')
    F = np.sort(order[:N]); R = np.sort(order[N:])
    OF = O[:, F]
    sign, logdetOF = np.linalg.slogdet(OF)
    absdetOF2 = np.exp(2 * logdetOF.real)
    Ot = np.linalg.solve(OF, O[:, R])                  # N x (L-N)
    mu = eq[F].max()
    EF = eq[F].sum()
    s = np.arange(L)
    phaseF = np.exp(2j * np.pi * s * (n[F].sum() + N * Phim) / L)   # exp(i s sum_F q)
    diagph = np.exp(-2j * np.pi * np.outer(s, n[F]) / L)           # (L, N)
    out = []
    for tau in taus:
        w = np.exp(-tau * (eq[R] - mu))
        Z = np.zeros((L, N, N), complex)
        Ow = (Ot * w).T
        Z[R] = Ow[:, :, None] * Ot.conj().T[:, None, :]
        Y = np.fft.ifft(Z, axis=0) * L                 # Y[s] = sum_n Z[n] e^{2 pi i s n / L}
        dF = np.exp(tau * (eq[F] - mu))                # <= 1
        X = (diagph * dF)[:, :, None] * Y
        X[:, np.arange(N), np.arange(N)] += 1.0
        dets = np.linalg.det(X)
        tot = np.sum(phaseF * dets)
        val = (L / N) ** 2 / L * st['pref2'] * absdetOF2 * np.exp(-tau * (EF - st['E0'])) * tot
        out.append(val)
    return np.array(out)


def _sector_job(args):
    N, V, t, m, taus = args
    st = setup(N, V, t)
    return sector_correlator(st, N, m, taus)


def kernel(N, V, t=1.0, h=0.25, umin=-20.0, umax=None, return_sectors=False, procs=1):
    """Return chi_N(r), r = 0..N-1, and diagnostics."""
    st = setup(N, V, t)
    L = st['L']
    if umax is None:
        umax = np.log(12.0 * L / t)
    u = np.arange(umin, umax + 1e-12, h)
    taus = np.exp(u)
    a = contact_mean(N, V)
    chihat = np.zeros(N, complex); C0sum = 0.0; maximag = 0.0; tailmax = 0.0
    allt = np.concatenate([[0.0], taus])
    jobs = [(N, V, t, m, allt) for m in range(N // 2 + 1)]
    if procs > 1:
        import multiprocessing as mp
        with mp.get_context('fork').Pool(procs) as pool:
            results = pool.map(_sector_job, jobs, chunksize=1)
    else:
        results = [sector_correlator(st, N, m, allt) for m in range(N // 2 + 1)]
    for m, Cm in enumerate(results):
        maximag = max(maximag, np.max(np.abs(Cm.imag)))
        Cm = Cm.real
        c0 = Cm[0]; Cm = Cm[1:]
        if m == 0:
            Cm = Cm - a * a
        mult = 1 if (m == 0 or 2 * m == N) else 2
        C0sum += mult * c0
        f = Cm * taus
        tailmax = max(tailmax, abs(f[-1]))
        val = h * f.sum() + (c0 - (a * a if m == 0 else 0.0)) * np.exp(umin)   # tau < e^umin

        chihat[m] = val
        if m and 2 * m != N:
            chihat[N - m] = val
    r = np.arange(N)
    th = 2 * np.pi * np.outer(r, np.arange(N)) / N
    chi = (np.cos(th) @ chihat.real)
    diag = dict(N=N, V=V, L=L, a=a, sum_C0=C0sum, sum_rule_error=abs(C0sum - a),
                max_imag=maximag, tail_integrand=tailmax, h=h, umin=umin, umax=umax)
    if return_sectors:
        diag['chihat'] = chihat.real.tolist()
    return chi, diag
