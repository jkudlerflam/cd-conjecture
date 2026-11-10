"""C191: low-energy (Luttinger + g1) parameters of the two-valley gas of model (3.2) from the exact
two-body S-matrix at the Fermi-surface kinematics.  |w| = 1, phi = pi/2, rho = n/2 per valley.

Frame: valley fields mutually bosonic, fermionic within a valley (the frame in which the unlike pair's
transmission is close to its free value +1; see PROOF section 3).  Landau parameters (first order in the
phases, f = -(v_j - v_l) tau for an approaching pair with transmission phase tau relative to free):
  f_par = (aR,aL) = -(vR + vL) * 2 delta_s          (like pair, fermion reference S = -1)
  f_pp  = (aR,bL) = -2 vR * arg t,   f_mm = (aL,bR) = -2 vL * arg t'
  f_4   = (aR,bR) = (aL,bL): variants, see below
g1 = 2 sqrt(vR vL) |r|, y1 = g1 / (pi vbar).
"""
import sys
import numpy as np
from scat191 import Model, unlike_amplitudes, like_amplitude
from gauss191 import metric, dim, OPS, W_from


def threshold_g4(m, ka, kb, pL, pR):
    """K-matrix (scattering-length) pseudopotential for the (aR, bR) pair: total K = ka+kb+pR-pL,
    threshold a_+ from the even eigenphase, g4 = -(eps''_1 + eps''_2) / a_+."""
    K = ka + kb + pR - pL
    ps = np.linspace(-np.pi, np.pi, 400001)
    Et = m.eps(ps) + m.eps(K - ps)
    i = np.argmin(Et)
    p0 = ps[i]
    for _ in range(40):
        g1 = m.deps(p0) - m.deps(K - p0); g2 = m.d2eps(p0) + m.d2eps(K - p0)
        p0 -= g1 / g2
    Emin = m.eps(p0) + m.eps(K - p0)
    inv_mu = m.d2eps(p0) + m.d2eps(K - p0)
    out = []
    for dE in (1e-7, 1e-6):
        S, pin, pout, v = m.smatrix(K, Emin + dE)
        q = np.sqrt(dE / (inv_mu / 2))          # E_rel = q^2 (eps''1+eps''2)/2
        ev = np.linalg.eigvals(-S)
        a = np.sort(-np.angle(ev) / (2 * q))
        out.append(a)
    a_m, a_p = out[0]
    return -inv_mu / a_p, a_p, a_m, inv_mu


def analyse(n, m=None, verbose=True):
    m = m or Model(w=1j)
    ka, kb = m.valleys()[:2]
    rho = n / 2
    pL, pR, eF = m.fermi_points(ka, rho)
    vR, vL = m.deps(ka + pR), -m.deps(ka + pL)
    Lk = like_amplitude(m, ka, pL, pR)
    U = unlike_amplitudes(m, ka, kb, pL, pR)
    ds = Lk['delta']
    tau, taup = np.angle(U['t']), np.angle(U['tp'])
    r = abs(U['r'])
    vbar = 0.5 * (vR + vL)
    g1 = 2 * np.sqrt(vR * vL) * r
    y1 = g1 / (np.pi * vbar)
    g4_KRR, ap_RR, am_RR, invmu = threshold_g4(m, ka, kb, pL, pR)
    g4_Kpi = -2 * m.d2eps(ka) / 20.7111  # K = pi threshold, a_+ = 20.711 (closed form, C175)
    res = dict(n=n, pL=pL, pR=pR, vR=vR, vL=vL, eF=eF - m.eps(ka), delta_s=ds, tau=tau, taup=taup, r=r,
               y1=y1, aplus_RR=ap_RR, g4_KRR=g4_KRR, g4_Kpi=g4_Kpi)
    rows = []
    for like in ('delta', 'tan'):
        dsv = ds if like == 'delta' else np.arctan(np.tan(ds))  # same angle; tan variant below
        f_par = -(vR + vL) * 2 * (ds if like == 'delta' else np.tan(ds))
        f_pp = -2 * vR * (tau if like == 'delta' else np.tan(tau))
        f_mm = -2 * vL * (taup if like == 'delta' else np.tan(taup))
        for lab4, f4 in (('g4=0', 0.0), ('g4=K(pi)', g4_Kpi), ('g4=K(RR)', g4_KRR)):
            W = W_from(vR, vL, f_par, f_pp, f_mm, f4)
            G, lam, sig = metric(W)
            x1 = dim(G, OPS['O1 = cos 2(phi_a - phi_b)'])
            xp = dim(G, OPS['pair Psi_a Psi_b = e^{i(th_a+th_b)}'])
            xb = dim(G, OPS['boson Psi_a = e^{i th_a}'])
            xf1 = dim(G, OPS['fermion e^{i th_a + i(phi_a - phi_b)}'])
            xf2 = dim(G, OPS['fermion e^{i th_a + i(phi_a + phi_b)}'])
            Ks, Kc = x1 / 2, 1 / (2 * xp)
            rows.append(dict(variant=f"{like},{lab4}", Ks=Ks, Kc=Kc, exp_pair=2 * xp, exp_boson=2 * xb,
                             exp_fermion=2 * min(xf1, xf2), margin=(x1 - 2) - abs(y1), vel=sorted(np.abs(lam))))
    res['rows'] = rows
    if verbose:
        print(f"n={n:.5f} pL={pL:+.4f} pR={pR:+.4f} vR={vR:.4f} vL={vL:.4f} eF-emin={res['eF']:.4f}")
        print(f"   delta_s={ds:+.4f}  arg t={tau:+.4f} arg t'={taup:+.4f} |r|={r:.4f}  y1={y1:.4f}  "
              f"a+(K_RR)={ap_RR:.3f} a-(K_RR)={am_RR:.3f} g4_KRR={g4_KRR:+.4f} g4_Kpi={g4_Kpi:+.4f}")
        for row in rows:
            print(f"   {row['variant']:16s} K_s={row['Ks']:.3f} K_c={row['Kc']:.3f}  exps: pair={row['exp_pair']:.3f} "
                  f"boson={row['exp_boson']:.3f} fermion={row['exp_fermion']:.3f}  KT margin 2(Ks-1)-|y1|={row['margin']:+.3f}"
                  f"  velocities={np.round(row['vel'],3)}")
    return res


if __name__ == "__main__":
    ns = [float(eval(a)) for a in sys.argv[1:]] or [1/16, 3/32, 1/8, 5/32, 3/16]
    for n in ns:
        analyse(n)
