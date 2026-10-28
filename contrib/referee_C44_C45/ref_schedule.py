"""Carry out C45 PROOF Sec. 8 (with C44's explicit wrapper bounds) at n=1,2,3.

Inputs are C42's certified data (budgets_n*.json: tau, quartic margin delta,
normalized speed) and C42's witness eta = 1/100, u = 1/64, w = 1/4,
mu = 10^(-8n-10). Everything is in Decimal with unbounded exponents.

Placeholders, stated in the README: the S constants A_n = A_dec,n = 1 (any
larger value only lowers lambda), the core grammar penalty K_c = D, the helix
scale J_h = U = 1, and the root-patch data replaced by the margins themselves
(q_patch = v_* delta, e_det = v_* delta / 32, e_top = 4). The point is to
check that each choice is finite and fixed before the next, and to record the
resulting scales; it is not a certificate.
"""
import json, math
from decimal import Decimal as Dm, getcontext
from pathlib import Path

ctx = getcontext(); ctx.prec = 40; ctx.Emax = 10 ** 15; ctx.Emin = -10 ** 15
C42 = Path(__file__).resolve().parents[2] / 'codex' / 'work_C42'
L10 = lambda x: float(x.log10())


def schedule(n):
    bud = json.loads((C42 / f'budgets_n{n}.json').read_text())
    tau = Dm(10) ** Dm(bud['tau_log10'])
    delta = Dm(10) ** Dm(bud['delta_log10'])
    vstar = min(Dm(1), Dm(10) ** Dm(bud['speed_log10']))
    eta, w, u = Dm(1) / 100, Dm(1) / 4, Dm(1) / 64
    mu = Dm(10) ** (-8 * n - 10)
    q = (w * u) ** 2
    p = 2 * n + 4
    k, b = n + 2, n + 3
    U = Dm(1)
    sigma = Dm(3 * math.sqrt(3 / 5) - 2)
    # Sec. 7: epsilon_n with A_n = 1
    e_det = vstar * delta / 32
    ell = (Dm(4) / e_det).ln()
    eps = min(Dm(1), vstar * delta / 1024, sigma * vstar / (1024 * (1 + ell)))
    # step 1, (27)
    a = 16 * (n + 1)
    Cb = 16 * (n + 1) ** 2 * Dm(3 * a).exp()
    theta = min(tau, eps / (16 * Cb))
    # step 2, (28)
    r, B0 = Dm(1) / 16384, Dm(128)
    D = max(2 / r, 8 * n * B0 * (1 + eta) / (theta * r ** (p + 1)))
    # step 3, (29)
    g2 = min(Dm(1) / 16, theta / (192 * n * (1 + eta) * D ** p),
             eps / (16 * Cb * (2 + 64 * n * n) * 3 * (1 + eta) * D ** p))
    g = g2.sqrt()
    Sq = sum(q ** (2 * j) for j in range(n + 2))
    V0 = q ** (2 * n + 2) / Sq ** 2
    J0 = 2 * V0 * g2 / D ** p
    # step 4, (30)-(31), C44 Sec. 6
    Kc = D
    hmax, amax = Kc + D + 12, Kc + D + 5
    Mc = (k - 1) * hmax + 2 * amax
    Lc, Lbond, ac = 2 * Mc, 2 * (k - 1) * hmax, 1 - g
    x = (3 * (1 + eta) * V0).sqrt()
    Mh0 = Dm(16 / (3 * math.sin(1 / (2 * b)) ** 2))
    y = x * Mh0 ** 2
    deq = J0 * eps / (16 * Cb)
    z0 = min(Dm(1) / 2, (deq * ac ** 2 / (2 * x ** 2 * Lbond)).sqrt())
    rho = z0 / (2 * (k + 1) * (1 + z0))
    dch = {}
    for pp in (1, 2):
        Np, mp = k + pp, math.comb(k + pp, pp)
        dch[pp] = rho * U * (rho / (Np - 1)) ** (2 * k * pp) / (mp * k * pp)
    Jh = Dm(1)
    dfast = min(dch[1], dch[2], Jh / ((b + 1) * (2 * b + 1)))
    A1 = 2 * Lc ** 2 * x ** 2 / ac ** 2 + 4 * Lc * x * y / ac + y ** 2
    A2 = 2 * Lc ** 2 * y ** 2 / ac
    t = min(Dm(1), ac / (2 * Lc ** 2), deq / (2 * (A1 + A2)))
    s = t * dfast
    # step 5: Delta from C44 (3)-(4), filter and selector constants
    kap = Dm(1) / ((k + 2) * (2 * k + 3)); bet = Dm(math.sin(math.pi / (k + 2)) ** 2)
    epsf = kap * bet / ((1 + kap.sqrt()) ** 2 + bet)
    Jf = 2 * U / epsf
    gam_c = 1 - g
    gocc = min(s * gam_c, Jf / (k * (2 * k - 1)))
    Mocc = s * Mc + Jf * (k + 1)
    dmin = 1 - ((rho * (k + 1)) / (1 - rho * (k + 1))) ** 2
    Docc = min(dch[pp] * dmin * gocc / ((dch[pp].sqrt() + Mocc.sqrt()) ** 2 + dmin * gocc) for pp in (1, 2))
    Delta = min(Docc, U * (1 - rho * (k + 1)), Jh / ((b + 1) * (2 * b + 1)), U)
    a1, a2 = rho / (k + 1 - 1), rho / (k + 2 - 1)
    bv = (k + 1) ** 2 * a1 ** (2 * k) * a2 ** (4 * k) / (math.comb(k + 1, 1) * math.comb(k + 2, 2))
    Ab = sum(math.sqrt(2 / b) * math.sin(math.pi * j / b) * math.sqrt(2 / (b + 1)) * math.sin(math.pi * j / (b + 1))
             for j in range(1, b))
    beta = Dm(Ab ** 2 * math.cos(1 / (2 * b)) ** (2 * b)) * bv
    cn = J0 / s
    f_mark, f_read = cn / beta, cn * mu * Sq / q ** 2
    # step 6, (11) and (35), weights e^{a|X|}
    v0 = max(Dm(1), Dm(2 * a).exp() * 3 * (1 + eta.sqrt()) * Mh0 ** 2)
    f0 = Dm(2 * a).exp() * f_mark + Dm(a).exp() * k * f_read
    G = f0 + 48 * v0 ** 2 / Delta
    CSW = 2 ** 16 * (v0 * G / Delta + v0 ** 3 / Delta ** 2)
    lam_caps = dict(sw_v0_over_G=v0 / G, sw_Delta=Delta / (192 * v0),
                    remainder=eps * cn / (16 * CSW),
                    bandwidth=(Delta / (16 * cn * (b + 1) * 4)).sqrt())
    lam = min(lam_caps.values())
    W = lam ** 2 * cn
    vals = dict(eps_n=eps, theta=theta, D=D, g=g, J0=J0, d_eq=deq, rho=rho,
                delta_ch1=dch[1], delta_ch2=dch[2], s=s, Delta=Delta, beta_n=beta,
                c_n=cn, t_m_coefficient=f_mark, v0=v0, f0=f0, C_SW=CSW, lam=lam, W=W,
                slowest_speed=W * vstar, E_det=W * e_det, W_over_Delta=W / Delta,
                remainder_over_W=CSW * lam ** 3 / W)
    res = {k_: L10(v) for k_, v in vals.items()}
    res['binding_lambda_cap'] = min(lam_caps, key=lambda c: lam_caps[c])
    res['ell_star'] = float(ell)
    res['crossover_cells_length_mode'] = L10(1 / e_det)
    # checks: each later quantity positive and fixed; ordering inequalities hold
    res['checks'] = dict(all_positive=all(v > 0 for v in vals.values()),
                         s_below_dfast=bool(s < dfast), W_below_Delta=bool(W < Delta / 16),
                         remainder_below_epsW=bool(CSW * lam ** 3 <= eps * W / 16))
    return res


out = {f'n{n}': schedule(n) for n in (1, 2, 3)}
Path('ref_schedule.json').write_text(json.dumps(out, indent=1))
keys = ['eps_n', 'D', 'g', 'J0', 'rho', 'delta_ch2', 's', 'Delta', 'beta_n', 'c_n', 'C_SW', 'lam', 'W',
        'slowest_speed', 'E_det', 'crossover_cells_length_mode']
print('log10 of', ' '.join(keys))
for n, r in out.items():
    print(n, ' '.join('%s=%.4g' % (k_, r[k_]) for k_ in keys), r['binding_lambda_cap'], r['checks'])
