#!/usr/bin/env python3
"""Route-A scale bookkeeping (numpy/scipy/stdlib only; seconds).

Recomputes, from the formulas in the cited sources, every scale used in
README.md section 3, and writes scales_results.json next to this file.

  A. C19 Jacobi filters: exact rational r_max^2 and rho^2(n) (C19 eqs. 8-16).
  B. Conditioning of the C19 filters: how errors in pole energies or in the
     relative size of individual spectral weights leak into cross-flavour
     hopping (new here).
  C. Fermi-point margins at n = 2, 3: C19's geometric design (C19 eq. 7),
     the Sidon design of contrib/ir_control (astra_A/margins.py recipe),
     C13's mixed margin, and C20's existence margin 1/(4M(n)).
  D. Length selection: C7 walker (selector_wrapper eq. 6) and the C18
     clock-selected length (C18 section 5), plus the helix equalizer.
  E. A worked hierarchy at n = 2, 3 for an assumed cell length.
"""
import json
import math
import itertools as it
from fractions import Fraction
from pathlib import Path

import numpy as np

OUT = {}


# ---------------------------------------------------------------- A. filters
def jacobi_integer(j):
    """Integer coefficients of Q_j(x) (C19 eq. 8), lowest power first."""
    out = [0] * (j + 1)
    for l in range(j + 1):
        c = math.comb(j + 2, l) * math.comb(j + 1, j - l)
        for r in range(j - l + 1):
            out[l + r] += c * math.comb(j - l, r) * (-1) ** (j - l - r)
    return out


def raw_filters(n):
    """c[e][j] with x(1-x)Q_j(x) = sum_e c[e][j] x^e, e = 1..n+1; h_j."""
    k = n + 1
    c = [[0] * n for _ in range(k)]
    h = []
    for j in range(n):
        q = jacobi_integer(j)
        f = [0] * (n + 3)
        for l, v in enumerate(q):
            f[l + 1] += v
            f[l + 2] -= v
        for e in range(1, k + 1):
            c[e - 1][j] = f[e]
        h.append(Fraction(j + 1, (2 * j + 4) * (j + 3)))
    return c, h


def filter_data(n):
    c, h = raw_filters(n)
    k = n + 1
    # exact checks of C19 eq. 9 in rationals
    for j in range(n):
        assert sum(c[e][j] for e in range(k)) == 0
        for i in range(n):
            g = sum(Fraction(c[e][i] * c[f][j], (e + 1) + (f + 1))
                    for e in range(k) for f in range(k))
            assert g == (h[i] if i == j else 0)
    r2 = [2 * sum(Fraction(c[e][j] ** 2) / h[j] for j in range(n)) for e in range(k)]
    r2max = max(r2)
    rho2_hop = Fraction(1, 16 * k * k) / r2max
    X = np.array([[c[e][j] / math.sqrt(float(h[j])) for j in range(n)] for e in range(k)])
    return dict(k=k, r2max=r2max, rho2_hop=rho2_hop, X=X)


def geometric_design(n):
    """C19 eq. 7."""
    al = math.exp(-1)
    t, th, mu, pts = [], [], [], []
    for a in range(1, n + 1):
        pl, pr = math.pi * al ** (2 * a), math.pi * al ** (2 * a - 1)
        ka = (pr - pl) / 2
        ta = 1 + al ** (2 * a)
        t.append(ta); th.append(-(pr + pl) / 2); mu.append(2 * ta * math.cos(ka))
        pts += [pl, pr]
    return np.array(t), np.array(th), np.array(mu), np.array(pts)


def rho2_with_fields(n, fd, mu):
    k = n + 1
    g = 1 / (2 * (k + 2))
    Fnorm = float(np.sum(np.abs(mu)) / 2)  # ||-sum mu_a Z_a/2||
    return min(float(fd['rho2_hop']), g / (4 * k * Fnorm))


A = []
for n in range(1, 13):
    fd = filter_data(n)
    _, _, mu, _ = geometric_design(n)
    A.append(dict(n=n, k=fd['k'], r2max=float(fd['r2max']),
                  rho2_hopping_only=float(fd['rho2_hop']),
                  rho2_with_fields_geometric=rho2_with_fields(n, fd, mu),
                  cell_dim_reflection=16 * (n + 2) * 2 ** n,
                  cell_dim_separated=16 * (n + 2) ** 2 * 2 ** n))
assert abs(A[0]['rho2_hopping_only'] - 1 / 1536) < 1e-18
assert abs(A[1]['rho2_hopping_only'] - 1 / 172800) < 1e-18
OUT['A_filters'] = A


# ------------------------------------------------------- B. conditioning
def sensitivities(n, rng, samples=200):
    fd = filter_data(n)
    X, k = fd['X'], fd['k']
    E = np.arange(1, k + 1, dtype=float)
    C = 1 / (E[:, None] + E[None, :])
    assert np.allclose(X.T @ C @ X, np.eye(n), atol=1e-6 * max(1, np.abs(X).max() ** 2))
    # (i) relative error eps_e on the whole weight M_e (both p_e and D_e):
    #     hopping matrix -> X^T diag(1+eps) C diag(1+eps) X.  Linear map in eps.
    Lw = [X.T @ (np.outer(np.eye(k)[e], C[e]) + np.outer(C[e], np.eye(k)[e])) @ X for e in range(k)]
    # remove the trace-preserving uniform part: eps_e = const only rescales.
    worst_w = sum(np.linalg.norm(m, 2) for m in Lw)
    typ_w = []
    # (ii) error delta_e in pole energy e (units of Delta):
    #     dC_ef = -(delta_e + delta_f)/(E_e+E_f)^2.
    Le = [X.T @ (-(np.outer(np.eye(k)[e], C[e] ** 2) + np.outer(C[e] ** 2, np.eye(k)[e]))) @ X
          for e in range(k)]
    worst_e = sum(np.linalg.norm(m, 2) for m in Le)
    typ_e = []
    for _ in range(samples):
        eps = rng.uniform(-1, 1, k)
        m = sum(eps[e] * Lw[e] for e in range(k))
        off = m - np.diag(np.diag(m))
        typ_w.append(np.linalg.norm(off, 2))
        d = rng.uniform(-1, 1, k)
        m = sum(d[e] * Le[e] for e in range(k))
        off = m - np.diag(np.diag(m))
        typ_e.append(np.linalg.norm(off, 2))
    cond = np.linalg.cond(C)
    return dict(n=n, cond_C=float(cond),
                weight_error_worst=float(worst_w), weight_error_typical_offdiag=float(np.median(typ_w)),
                pole_error_worst=float(worst_e), pole_error_typical_offdiag=float(np.median(typ_e)))


rng = np.random.default_rng(23)
OUT['B_conditioning'] = [sensitivities(n, rng) for n in range(1, 9)]


# ------------------------------------------------------------- C. margins
def circ(x):
    return np.abs((np.asarray(x) + np.pi) % (2 * np.pi) - np.pi)


def margins(p):
    p = np.asarray(p, float)
    N = len(p)
    bil = [s * p[i] + t * p[j] for i, j in it.combinations(range(N), 2)
           for s, t in it.product([-1, 1], repeat=2)]
    d2 = float(min(np.min(circ(2 * p)), np.min(circ(bil))))
    d4 = None
    if N >= 4:
        d4 = float(min(circ(sum(s * p[i] for i, s in zip(ii, ss)))
                       for ii in it.combinations(range(N), 4)
                       for ss in it.product([-1, 1], repeat=4)))
    return d2, d4, np.array(bil)


def mixed(bil):
    roots = np.unique((bil / (2 * np.pi)) % 1)
    cuts = np.r_[.25, roots[(roots > .25) & (roots < .75)], .75]
    j = int(np.argmax(np.diff(cuts)))
    nu = (cuts[j] + cuts[j + 1]) / 2
    dm = float(min(np.min(circ(bil + 2 * np.pi * nu)), np.min(circ(bil - 2 * np.pi * nu))))
    return nu, dm


def primes(count):
    out, k = [], 2
    while len(out) < count:
        if all(k % q for q in out if q * q <= k):
            out.append(k)
        k += 1
    return out


def sidon(n):
    """contrib/ir_control/astra_A/margins.py designed(), jitter included."""
    size = 2 * n
    q = next(x for x in primes(size + 10) if x >= size)
    i = np.arange(size)
    a = 2 * q * i + i ** 2 % q
    d = int(np.min(np.diff(a)))
    b = (d + 1) // 2
    T = 2 * int(a[-1]) + 4 * b + 2
    pts = 2 * np.pi * (a + b + 1 / 6) / T + 1e-10 * np.sqrt(np.array(primes(size), float))
    k = (pts[n:] - pts[:n]) / 2
    th = (pts[n:] + pts[:n]) / 2
    return pts, 2 * np.sin(k), th


def M_C20(n):
    return 1 + 8 * n + 32 * math.comb(n, 2) + 96 * math.comb(n, 3) + 128 * math.comb(n, 4)


C = []
for n in [1, 2, 3, 4]:
    t, th, mu, pg = geometric_design(n)
    d2g, d4g, bilg = margins(pg)
    kg = (pg[1::2] - pg[0::2]) / 2
    vg = 2 * t * np.sin(kg)
    nug, dmg = mixed(bilg)
    ps, vs, ths = sidon(n)
    d2s, d4s, bils = margins(ps)
    nus, dms = mixed(bils)
    Kc = 4 * math.comb(2 * n, 2)
    C.append(dict(n=n,
                  geometric=dict(delta2=d2g, delta4=d4g, delta_mix=dmg, nu=nug,
                                 v_min=float(vg.min()), pair_momentum_min=float(np.min(circ(2 * th)))),
                  sidon=dict(delta2=d2s, delta4=d4s, delta_mix=dms, nu=nus,
                             v_min=float(vs.min()), pair_momentum_min=float(np.min(circ(2 * ths))),
                             bound_delta2=math.pi / (6 * n), bound_delta4=math.pi / (576 * n * n),
                             bound_delta_mix=math.pi / (4 * (Kc + 1))),
                  C20=dict(M=M_C20(n), delta_normalized=1 / (4 * M_C20(n)),
                           delta_angle=2 * math.pi / (4 * M_C20(n)))))
assert [M_C20(n) for n in range(1, 9)] == [9, 49, 217, 737, 1961, 4369, 8569, 15297]
OUT['C_margins'] = C


# ---------------------------------------------------------- D. selectors
def a_walk(l, t=1.0):
    return 4 * t * math.sin(math.pi / (2 * l)) ** 2   # 2t(1-cos(pi/l)), stable


def c7(b, t=1.0):
    a = lambda l: a_walk(l, t)
    eta = min((a(b - 1) - 2 * a(b) + a(b + 1)) / 2, (a(b) - 2 * a(b + 1) + a(b + 2)) / 3)
    delta_b = (b + 1) * a(b) - b * a(b + 1)
    e_c = a(b + 1) - a(b)
    mu_c = 2 * t - delta_b
    omega = 2 * t * (math.cos(math.pi / b) - math.cos(2 * math.pi / b))
    # selector_wrapper eq. 7: t = 2U[1/eta + b/mu_c + 1/|e_c|] (t=1 units)
    U_over_t = 1 / (2 * (1 / eta + b / mu_c + 1 / abs(e_c)))
    return dict(b=b, a_b=a(b), eta_b=eta, omega_head=omega, U_over_t_max=U_over_t,
                margin_over_U=eta / U_over_t - 1)


def clock(nb, J=1.0):
    """C18 section 5: nb bulk sites, physical period nb+2, 15-letter clock."""
    Cn = lambda m: (2 * m + 4) * 2 ** m + m - 4
    a = lambda m: 4 * J * math.sin(math.pi / (2 * (Cn(m) + 1))) ** 2
    C = Cn(nb)
    gamma_len = (a(nb - 1) - 2 * a(nb) + a(nb + 1)) / 2
    omega = 4 * J * math.sin(3 * math.pi / (2 * (C + 1))) * math.sin(math.pi / (2 * (C + 1)))
    f_min = 2 * math.sin(math.pi / (C + 1)) ** 2 / (C + 1)
    g_in = omega * J * f_min / (omega + J)   # eta = J pin, conservative bound
    pin_lift = J * math.pi ** 2 / C ** 3     # sharper endpoint example
    return dict(bulk_sites=nb, period=nb + 2, C=C, a=a(nb), gamma_len=gamma_len,
                Delta_density=gamma_len / (nb + 3), omega=omega,
                init_gap_lower_bound=g_in, endpoint_pin_lift=pin_lift)


def helix(b, r1=0.5, r2=0.5):
    phi = 1 / b
    beta = math.cos(phi / 2) ** (2 * b)
    gap = 1 / ((b + 1) * (2 * (b + 1) - 1))
    comp = (1 / r1 + 1 / r2) / math.sin(phi / 2) ** 2
    return dict(b=b, beta_helix=beta, helix_gap_bound=gap, compensator_norm_bound=comp)


tab = clock(2)
assert abs(tab['gamma_len'] - 3.94531e-2) < 1e-6 and abs(tab['omega'] - 3.06788e-2) < 1e-6
tab = clock(8)
assert abs(tab['gamma_len'] - 5.89491e-7) < 1e-11 and abs(tab['omega'] - 1.12728e-6) < 1e-10
OUT['D_selector_C7'] = [c7(b) for b in [8, 10, 12, 16, 20, 32]]
OUT['D_clock_C18'] = [clock(m) for m in [2, 4, 6, 8, 10, 12, 16]]
OUT['D_helix'] = [helix(b) for b in [8, 10, 12, 16, 32]]


# ------------------------------------------------------- E. worked hierarchy
def ell_min(n):
    """Register-count floor on the cell length: ports, partners and the
    n + log2(16(n+2)) qubits of the C19 cell, one qubit per site."""
    return n + math.ceil(math.log2(16 * (n + 2))) + 2


try:   # third-order coefficients from third_order.py, if it has been run
    C3 = {r['n']: r['c3_traceless'] for r in
          json.loads(Path(__file__).with_name('third_order_results.json').read_text())}
except FileNotFoundError:
    C3 = {2: 4.465, 3: 6.521}
SAFETY = 0.1   # every "<<" below is implemented as a factor 0.1
E = []
for n in [2, 3]:
    l = ell_min(n)
    cl = clock(l - 2)
    sel = c7(l)
    fd = filter_data(n)
    rho2 = float(fd['rho2_hop'])
    marg = [c for c in C if c['n'] == n][0]
    sens = [s for s in OUT['B_conditioning'] if s['n'] == n][0]
    vdelta = marg['sidon']['v_min'] * marg['sidon']['delta2']
    eps_w = SAFETY * vdelta / sens['weight_error_worst']     # allowed relative weight error
    lam_over_src = SAFETY * vdelta / C3[n]                   # third order below margins
    chain = {}
    for label, g_in in [('pin_lift', cl['endpoint_pin_lift']), ('proved_bound', cl['init_gap_lower_bound'])]:
        # Schur error ~ M/g_in with M ~ Delta_src (C18 section 6)  ->  Delta_src = eps_w * g_in
        d_src = eps_w * g_in
        lam = lam_over_src * d_src
        W = lam ** 2 * rho2 / d_src
        chain[label] = dict(g_in=g_in, Delta_src=d_src, lam=lam, W=W, W_phys_velocity=l * W)
    E.append(dict(n=n, ell_min=l, C=cl['C'], clock_gap=cl['omega'], clock_length_margin=cl['gamma_len'],
                  C7_eta=sel['eta_b'], C7_U_over_t=sel['U_over_t_max'], rho2=rho2,
                  v_delta2_sidon=vdelta, S_worst=sens['weight_error_worst'], c3=C3[n],
                  eps_w=eps_w, lam_over_Delta_src=lam_over_src, frozen_chain=chain,
                  sidon=marg['sidon'], C20=marg['C20']))
OUT['E_worked'] = E

Path(__file__).with_name('scales_results.json').write_text(json.dumps(OUT, indent=1, default=float))

# --------------------------------------------------------------- printout
print('A. rho^2(n), hopping-only and with geometric fields')
for r in A:
    print(f"  n={r['n']:2d} k={r['k']:2d} r2max={r['r2max']:.4e} rho2={r['rho2_hopping_only']:.4e} "
          f"rho2(fields)={r['rho2_with_fields_geometric']:.4e} D_ref={r['cell_dim_reflection']}")
print('B. conditioning (cross-flavour hopping per unit relative weight / pole error)')
for r in OUT['B_conditioning']:
    print(f"  n={r['n']} cond(C)={r['cond_C']:.3e} weight worst={r['weight_error_worst']:.3e} "
          f"typ={r['weight_error_typical_offdiag']:.3e} pole worst={r['pole_error_worst']:.3e} "
          f"typ={r['pole_error_typical_offdiag']:.3e}")
print('C. margins')
for r in C:
    print(' ', r['n'], json.dumps({k: (round(v, 6) if isinstance(v, float) else v)
                                   for k, v in r['geometric'].items()}))
    print('   sidon', json.dumps({k: (round(v, 6) if isinstance(v, float) else v)
                                  for k, v in r['sidon'].items()}))
    print('   C20', r['C20'])
print('D. C7 walker'); [print(' ', r) for r in OUT['D_selector_C7']]
print('D. C18 clock'); [print(' ', r) for r in OUT['D_clock_C18']]
print('D. helix'); [print(' ', r) for r in OUT['D_helix']]
print('E. worked'); [print(' ', r) for r in E]
