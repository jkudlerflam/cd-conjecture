"""C167 (Claude Opus 5.5, 2026-11-08): U(1)_m charge decomposition of the harmonics of a bond.

H = sum_j [s_j^T J s_{j+1} + f . s_j].  Polar J = S O, O = R_m(theta) (tk.spiral_parent).  In the spiral frame the
bond is tau_j^T R_m(j theta) S R_m(-j theta) tau_{j+1} (C163 PROOF 3.1 with S in place of -1), so the harmonics are
the parts of S with U(1)_m charge q = 1 (S_{m,perp}) and q = 2 (traceless in-plane part), and of f with q = 1 (f_perp).

Output per bond: parent Delta (AF convention), theta, and the norms
  h1x = |S_{m perp}| (Frobenius norm of the charge-1 symmetric part) / |S0|,
  h1f = |f_perp| / |S0|,
  h2  = |charge-2 part| / |S0|,
with |S0| the Frobenius norm of the axial parent part.  Also checks the identity
  "h1x = 0 and h1f = 0  <=>  the bond commutes with the pi rotation about m" (parity class), Lemma R0 of PROOF.md.
usage: python3 charges.py   (self test + the C159 attractive points)
"""
import json, numpy as np
import tk, est as E


def decompose(J, f=np.zeros(3)):
    p = tk.spiral_parent(J, f)
    S, m = p['S'], p['m']
    P = np.eye(3) - np.outer(m, m)
    S0 = p['s_m'] * np.outer(m, m) + p['s_perp'] * P
    Sq1 = P @ S @ np.outer(m, m) + np.outer(m, m) @ S @ P          # charge +-1 (m-perp block)
    Sq2 = P @ S @ P - p['s_perp'] * P                              # charge +-2 (traceless in-plane)
    fperp = P @ f
    n0 = np.linalg.norm(S0)
    res = np.linalg.norm(S - S0 - Sq1 - Sq2)
    return dict(Delta=p['Delta'], theta=p['theta'], m=m, h1x=np.linalg.norm(Sq1) / n0, h1f=np.linalg.norm(fperp) / n0,
                h2=np.linalg.norm(Sq2) / n0, hm=p['hm'] / n0, harm=p['harm'], resid=res, s_perp=p['s_perp'], s_m=p['s_m'])


def pi_rot(m):
    m = np.asarray(m, float) / np.linalg.norm(m)
    return 2 * np.outer(m, m) - np.eye(3)


if __name__ == '__main__':
    rng = np.random.default_rng(167)
    # Lemma R0 check: parity iff no charge-1 parts
    worst = 0.0
    for _ in range(200):
        J = rng.normal(size=(3, 3)); f = rng.normal(size=3)
        d = decompose(J, f)
        R = pi_rot(d['m'])
        Jp = R @ (d['m'] @ d['m'] * 0 + J) @ R    # R^T J R with R = R^T
        # remove charge-1 parts and check the pi rotation becomes a symmetry
        p = tk.spiral_parent(J, f)
        P = np.eye(3) - np.outer(d['m'], d['m'])
        mm = np.outer(d['m'], d['m'])
        S = p['S']; Sn = S - (P @ S @ mm + mm @ S @ P)
        Jn = Sn @ p['O']; fn = mm @ f
        worst = max(worst, np.max(np.abs(R @ Jn @ R - Jn)), np.max(np.abs(R @ fn - fn)))
        worst = max(worst, d['resid'])
    print('R0: removing the charge-1 parts gives a pi_m-symmetric bond, max violation', f'{worst:.1e}')
    # conversely: a parity bond XYZ + Dz has no charge-1 part
    w = 0.0
    for _ in range(100):
        a = rng.normal(size=4)
        J = np.diag(a[:3]); J[0, 1] += a[3]; J[1, 0] -= a[3]
        d = decompose(J)
        if abs(abs(d['m'][2]) - 1) < 1e-8:
            w = max(w, d['h1x'])
    print('R0 converse: XYZ + D_z bonds (spiral axis z) have charge-1 part', f'{w:.1e}')
    names = ['hi_TFM_59', 'hi_TFM_60', 'hi_TFM_62', 'hi_TFM_63', 'hi_TFM_64', 'hi_TFM_65', 'hi_TFM_66', 'dist_TFM_77',
             'dist_TFM_78', 'hi_T_34', 'hi_T_36', 'hi_T_37', 'hi_T_39', 'ising_T_43', 'ising_T_46', 'ising_T_41', 'hi_K_1']
    for n in names:
        d0 = json.load(open(f'c159_refine/{n}.json'))
        sc = d0['scan'] if isinstance(d0['scan'], dict) else eval(d0['scan'])
        v = np.asarray(sc['vec'], float); v /= np.linalg.norm(v)
        J, f = E.J_f(v)
        d = decompose(J, f)
        print(f"{n:12s} Delta={d['Delta']:+.2f} theta={d['theta']:.2f} harm={d['harm']:.2f}  h1(exch)={d['h1x']:.2f} "
              f"h1(field)={d['h1f']:.2f} h2={d['h2']:.2f}")
