"""Free-fermion (XX) checks of the modular-route identities.

XX chain h = SxSx + SySy = (1/2)(c^dag_x c_{x+1} + h.c.): e_0 = -1/pi per bond, v = 1, c = 1,
A = c v = 1, M = 1/2, B = c/v = 1, C = c = 1.

Interval I of ell sites (0..ell-1), bonds x = 0..ell-2, quadratic weight
    w_x = (x+1)(ell-1-x)/ell          (vanishes at the virtual bonds x=-1 and x=ell-1),
K_w = sum_x w_x h_x, W = sum_x w_x.

Checked quantities (all Gaussian, so exact at every ell):
 (a) defect  Delta(ell) = e_0 W - E_0(K_w)            -> slope in log ell should be A/(12 pi)
 (b) defect of the top Schmidt vector chi of rho_I    -> same slope (CFT claim)
 (c) energy-density profile of the K_w ground state   -> delta_x = (A/24 pi) ell^2/(x^2 (ell-x)^2)
 (d) relative entropy S(rho_I || tau_beta), tau_beta = exp(-beta K_w)/Z
        -> slope in log ell  beta A/(12 pi) + pi B/(3 beta) - C/3
 (e) same for two decoupled XX chains with t = 1/2 and t = 1/4 (velocities 1 and 1/2, d = 4)
 (f) layer-cake inequality  E_0(K_w) >= sum_k Delta_k E_0^open(n_k)  and the boundary energy f_b
"""
import numpy as np

def weights(ell):
    x = np.arange(ell - 1)
    return (x + 1) * (ell - 1 - x) / ell

def kmat(wts, t=0.5):
    n = len(wts) + 1
    k = np.zeros((n, n))
    k[np.arange(n - 1), np.arange(1, n)] = t * wts
    k[np.arange(1, n), np.arange(n - 1)] = t * wts
    return k

def corr_xx(ell):
    """<c^dag_x c_y> of the half-filled infinite XX chain restricted to ell sites."""
    x = np.arange(ell)
    D = x[:, None] - x[None, :]
    with np.errstate(divide='ignore', invalid='ignore'):
        C = (-1.0) ** D * np.sin(np.pi * D / 2) / (np.pi * D)   # occupied |k|>pi/2 (ground state of +t hopping)
    C[np.arange(ell), np.arange(ell)] = 0.5
    return C

def gs_projector(k):
    e, U = np.linalg.eigh(k)
    occ = U[:, e < 0]
    return occ @ occ.T, e

def chain_data(ell, t=0.5, beta_list=()):
    """All quantities for one XX chain with hopping t (velocity 2t, e_0 = -2t/pi per bond)."""
    wts = weights(ell); W = wts.sum()
    k = kmat(wts, t)
    G, e = gs_projector(k)
    E0 = e[e < 0].sum()
    e0 = -2 * t / np.pi
    C = corr_xx(ell)
    nu, V = np.linalg.eigh(C)
    nu = np.clip(nu, 1e-16, 1 - 1e-16)
    S_rho = -(nu * np.log(nu) + (1 - nu) * np.log(1 - nu)).sum()
    P = V[:, nu > 0.5] @ V[:, nu > 0.5].T          # top Schmidt vector (Gaussian)
    out = dict(W=W, defect=e0 * W - E0, defect_chi=e0 * W - np.trace(k @ P), S_rho=S_rho,
               E_rho=np.trace(k @ C), E0=E0, e0=e0)
    # profile of the K_w ground state
    bond = np.real(G[np.arange(ell - 1), np.arange(1, ell)]) * 2 * t   # <h_x> = 2t Re G_{x,x+1}
    out['delta_x'] = e0 - bond
    rel = {}
    for beta in beta_list:
        logZ = np.logaddexp(0, -beta * e).sum()
        rel[beta] = -S_rho + beta * out['E_rho'] + logZ
    out['rel'] = rel
    return out

def open_energy(n_bonds, t=0.5):
    k = kmat(np.ones(n_bonds), t)
    e = np.linalg.eigvalsh(k)
    return e[e < 0].sum()

if __name__ == '__main__':
    betas = [np.pi, 2 * np.pi, 4 * np.pi, 8 * np.pi]
    ells = [64, 128, 256, 512, 1024, 2048]
    print("== (a),(b),(d): single XX chain, A=B=C=1, M=1/2.  Predicted slopes: defect A/(12pi)=%.6f" % (1 / (12 * np.pi)))
    print("   S(rho||tau) slope beta/(12pi)+pi/(3beta)-1/3:", {f"{b/np.pi:.0f}pi": round(b / (12 * np.pi) + np.pi / (3 * b) - 1 / 3, 6) for b in betas})
    prev = None
    rows = []
    for ell in ells:
        d = chain_data(ell, 0.5, betas)
        rows.append((ell, d))
        if prev is not None:
            lp, dp = prev
            dl = np.log(ell / lp)
            s1 = (d['defect'] - dp['defect']) / dl
            s2 = (d['defect_chi'] - dp['defect_chi']) / dl
            sS = (d['S_rho'] - dp['S_rho']) / dl
            sr = {f"{b/np.pi:.0f}pi": round((d['rel'][b] - dp['rel'][b]) / dl, 6) for b in betas}
            print(f"ell={ell:5d}  defect={d['defect']:.6f} slope={s1:.6f}   chi-defect slope={s2:.6f}   S slope={sS:.5f} (c/3={1/3:.5f})   S(rho||tau) slopes={sr}")
        prev = (ell, d)
    print()
    print("== (c): profile x^2 (ell-x)^2 delta_x / ell^2  vs  A/(24 pi) = %.6f   (ell = 2048; x measured from the virtual bond, x -> x+1)" % (1 / (24 * np.pi)))
    ell, d = rows[-1]
    x = np.arange(1, ell)          # bond x sits at position x in units where w = x(ell-x)/ell
    dl = d['delta_x']
    sm = (dl[:-2] + 2 * dl[1:-1] + dl[2:]) / 4        # kills the alternating (Friedel) part and its first difference
    prof = x[1:-1]**2 * (ell - x[1:-1])**2 * sm / ell**2
    print("   (smoothed over three bonds; raw delta_x alternates in sign with amplitude ~ 0.03/x)")
    for xx in [2, 4, 8, 16, 64, 256, 512, 1024]:
        print(f"   x={xx:5d}  {prof[xx-2]:.6f}")
    print()
    print("== (e): two decoupled XX chains t=1/2 and t=1/4 (v=1 and 1/2), d=4, M=3/4: A=1.5, B=3, C=2")
    print("   predicted S(rho||tau) slopes:", {f"{b/np.pi:.0f}pi": round(1.5 * b / (12 * np.pi) + 3 * np.pi / (3 * b) - 2 / 3, 6) for b in betas}, "  defect slope A/(12pi)=%.6f" % (1.5 / (12 * np.pi)))
    prev = None
    for ell in ells:
        d1 = chain_data(ell, 0.5, betas); d2 = chain_data(ell, 0.25, betas)
        tot = dict(defect=d1['defect'] + d2['defect'], rel={b: d1['rel'][b] + d2['rel'][b] for b in betas})
        if prev is not None:
            lp, dp = prev; dl = np.log(ell / lp)
            s1 = (tot['defect'] - dp['defect']) / dl
            sr = {f"{b/np.pi:.0f}pi": round((tot['rel'][b] - dp['rel'][b]) / dl, 6) for b in betas}
            print(f"ell={ell:5d}  defect slope={s1:.6f}   S(rho||tau) slopes={sr}")
        prev = (ell, tot)
    print()
    print("== (f): layer cake and boundary energy, single XX chain")
    e0 = -1 / np.pi
    for n in [9, 17, 33, 65, 129, 257]:   # odd bond numbers: even site numbers, no zero mode
        fb = open_energy(n) - n * e0
        print(f"   n_bonds={n:4d}  f_b(n)=E_open-n e0={fb:.6f}   exact limit 1/2-2/pi={0.5-2/np.pi:.6f}   (n+2)*(f_b(n)-f_b)={(n+2)*(fb-(0.5-2/np.pi)):.5f}  (CFT: -pi A/24 = {-np.pi/24:.5f})")
    for ell in [32, 64, 128, 256]:
        wts = weights(ell); k = kmat(wts); e = np.linalg.eigvalsh(k); E0 = e[e < 0].sum()
        levels = np.unique(wts)
        lc = 0.0; prev_t = 0.0
        for tlev in levels:
            n_k = int((wts >= tlev).sum())
            lc += (tlev - prev_t) * open_energy(n_k); prev_t = tlev
        print(f"   ell={ell:4d}  E_0(K_w)={E0:.5f}  layer-cake lower bound={lc:.5f}  e0 W={e0*wts.sum():.5f}  defect={e0*wts.sum()-E0:.5f}  layer-cake allowance={e0*wts.sum()-lc:.5f}")
