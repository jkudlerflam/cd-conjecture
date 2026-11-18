"""Independent recomputation of C148 Section 7 and Theorem 4.4 numbers (mpmath intervals where a certificate is claimed)."""
import numpy as np, json
from mpmath import iv, mp, mpf
iv.dps = 40; mp.dps = 40
out = {}
# ---- Theorem 4.4, my own interval version, finer and coarser partitions, plus explicit monotonicity certificates
coth = lambda x: (iv.exp(2 * x) + 1) / (iv.exp(2 * x) - 1)
c = 2 * iv.pi * iv.sin(iv.pi / 8)
def br(l):
    t = iv.exp(-2 * l); return iv.mpf(1) / 2 + 2 * (t / (1 + t) + t ** 2 / (1 + t ** 2) + t ** 3 / (1 + t ** 3))
def cellmax(ncell, lo='1', hi='1.36'):
    worst = mpf(0)
    lo = iv.mpf(lo); hi = iv.mpf(hi)
    for k in range(ncell):
        l0 = lo + (hi - lo) * k / ncell; l1 = lo + (hi - lo) * (k + 1) / ncell
        # monotone enclosure, as in the proof
        num = l1 * (2 * coth(l0) - br(l1)); den = c * (1 - 2 * iv.exp(-iv.pi ** 2 / l1)) ** 2
        assert (2 * coth(l0) - br(l1)).a > 0 and (1 - 2 * iv.exp(-iv.pi ** 2 / l1)).a > 0
        v = num / den
        worst = max(worst, mpf(v.b))
    return worst
out['cell_max_36'] = float(cellmax(36)); out['cell_max_9'] = float(cellmax(9)); out['cell_max_3600'] = float(cellmax(3600))
# direct interval evaluation with the whole cell as an interval (no monotonicity assumed at all)
def cellmax_nomono(ncell):
    worst = mpf(0)
    for k in range(ncell):
        l = iv.mpf([mpf(1) + mpf('0.36') * k / ncell, mpf(1) + mpf('0.36') * (k + 1) / ncell])
        v = l * (2 * coth(l) - br(l)) / (c * (1 - 2 * iv.exp(-iv.pi ** 2 / l)) ** 2)
        worst = max(worst, mpf(v.b))
    return worst
out['cell_max_36_no_monotonicity_assumed'] = float(cellmax_nomono(36))
l = iv.mpf('1.36'); out['direct_margin_1.36'] = float(((iv.exp(l) - 1) / (iv.exp(l) + 1) - coth(l) / 2).a)
# is the direct margin positive below 1.36? root of tanh(l/2) = coth(l)/2
f = lambda x: mp.tanh(x / 2) - mp.coth(x) / 2
out['direct_margin_root'] = float(mp.findroot(f, 1.3))
l = iv.mpf(1); out['f_small_1'] = float(((2 * l * coth(l) - iv.log(2)) / (c * (1 - 2 * iv.exp(-iv.pi ** 2 / l)) ** 2)).b)
out['Phi_limit_0'] = float((2 - mp.log(2)) / (2 * mp.pi * mp.sin(mp.pi / 8)))
# f_small monotone on (0,1]: check on grid
g = [float((2 * x * mp.coth(x) - mp.log(2)) / (2 * mp.pi * mp.sin(mp.pi / 8) * (1 - 2 * mp.exp(-mp.pi ** 2 / x)) ** 2)) for x in np.linspace(1e-3, 1, 2000)]
out['f_small_increasing_on_grid'] = bool(np.all(np.diff(g) > 0))
# the true sup of the Phi majorant with the full series for pi J(0) (not a proof input, shows the margin)
def Phimaj(x):
    t = mp.exp(-2 * x); piJ0 = mpf(1) / 2 + 2 * mp.nsum(lambda j: t ** j / (1 + t ** j), [1, mp.inf])
    return x * (2 * mp.coth(x) - piJ0) / (2 * mp.pi * mp.sin(mp.pi / 8) * (1 - 2 * mp.exp(-mp.pi ** 2 / x)) ** 2)
out['Phi_majorant_exact_piJ0'] = {str(x): float(Phimaj(mpf(x))) for x in ['0.01', '0.5', '1', '1.2', '1.36', '1.5']}
# theta_2 bound and the elliptic identities
for lam in ['0.05', '0.5', '1', '1.36', '3']:
    q = mp.exp(-mpf(lam)); m = mp.mfrom(q=q); Ke = mp.ellipk(m); k = mp.sqrt(m)
    th2 = mp.jtheta(2, 0, q)
    out['kKe_check_' + lam] = [float(k * Ke), float(mp.pi / 2 * th2 ** 2),
                               float(mp.sqrt(mp.pi / mpf(lam)) * (1 - 2 * mp.exp(-mp.pi ** 2 / mpf(lam)))), float(th2)]
# ---- Section 7 majorant
pi = mp.pi; s = mpf(3) / 2; g_ = mpf('13.5') * s
t1 = g_ / pi; t2 = (s / 2) * (pi / 2) * mp.coth(pi) / pi
ends = 2 * s * (2 * g_ + 4 * (t1 + t2))
cp = 2 + 4 / (mp.e ** 2 - 1); x0 = pi ** 2 / 2
kA_term = max(cp * x0 * mp.coth(x0) / (2 * pi), pi / 4 * mp.coth(x0))
ratio = (4 + s / (2 * g_)) * kA_term
kw = 1 + mp.log(2 * cp)
Vf = 4 * cp + 4 * ratio + 8 + 24 * kw
curv = s ** 2 / 4
Ktot = Vf / 8 + ends + curv
out['sec7'] = dict(D_factor=float(t1 + t2), ends=float(ends), cp=float(cp), kappaA_D1_over_bD=float(ratio), kappa_w=float(kw),
                   Vf=float(Vf), core=float(Vf / 8), curv=float(curv), K_over_ab=float(Ktot),
                   K_over_M0=float(Ktot * 2 * pi ** 2), stated_218p4_times_2pi2=float(mpf('218.4') * 2 * pi ** 2),
                   doubling=float(mpf(5) / 2 * Ktot * 2 * pi ** 2), final=float(mpf(5) / 2 * Ktot * 2 * pi ** 2 + pi ** 2 / 4),
                   final_from_stated=float(mpf(5) / 2 * mpf('218.4') * 2 * pi ** 2 + pi ** 2 / 4))
# ---- sup over lambda of the pieces the majorant bounds, on a grid (checks each "<=" of Section 7)
def pieces(lam, s):
    a, b = np.pi * np.sinh(lam), np.pi / lam; g = 13.5 * s
    K0 = 1 / (np.pi * np.tanh(lam)); Kp1 = 4 / (np.pi * np.sinh(2 * lam))
    xx = np.linspace(0, np.pi / 2, 200001)
    Kp = np.abs(-2 * np.sinh(2 * lam) * np.sin(2 * xx) / (np.pi * (np.cosh(2 * lam) - np.cos(2 * xx)) ** 2)); M1 = Kp.max()
    D = g * Kp1 + K0 * s / 2; D1 = M1 * (4 * g + s / 2)
    cp = 2 + 4 / np.expm1(np.pi ** 2 / lam); kA = min(np.pi ** 2 / (2 * lam), cp); kw = min(np.pi ** 2 / (2 * lam), 1 + np.log(2 * cp))
    Vf = 4 * kA * (1 + D1 / (b * D)) + 8 + 24 * kw
    ends = 2 * s * (2 * g + 4 * min(lam / 2, np.pi / 2) * D)
    K = Vf / 8 + ends + s * s / 4
    abM0 = 4 * np.pi ** 2 * np.tanh(lam / 2) / lam
    return dict(K_over_ab=K, K_over_M0=K * abM0, kA_D1_bD=kA * D1 / (b * D), minD=min(lam / 2, np.pi / 2) * D, M1_ratio=M1 / Kp1 / (0.5 / np.tanh(lam)))
lams = np.r_[np.geomspace(1e-3, 0.05, 20), np.linspace(0.05, 12, 800), np.geomspace(12, 300, 30)]
P = [pieces(l, 1.5) for l in lams]
out['grid_s_3_2'] = {k: float(max(p[k] for p in P)) for k in P[0]}
P2 = [pieces(l, 0.5) for l in lams if l >= 1.36]
out['grid_s_1_2_lam_ge_1.36'] = {k: float(max(p[k] for p in P2)) for k in P2[0]}
json.dump(out, open('check_constants.json', 'w'), indent=1)
print(json.dumps(out, indent=1))
