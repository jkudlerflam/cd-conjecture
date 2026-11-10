"""C191 predictions. For each density: Landau (delta, tan) and hybrid (ABA intra-valley) estimates with four
variants of the R-R inter-valley parameter g4: 0, K-matrix at K = pi, K-matrix at the R-R kinematics, and +0.4
(an in-medium 'background' value with the shallow-pole part Pauli-blocked).  Outputs K_s, K_c, exponents,
u_c/K_c, u_s/K_s, L*Delta_P = (pi/2)(u_c/K_c + u_s/K_s), dp/ddelta = 1/(pi n u_s/K_s), y1, KT margin."""
import sys, json
import numpy as np
from scat191 import Model
from gauss191 import metric, dim, OPS, W_from
from lowE import analyse
from aba_valley import ABA

m = Model(w=1j)
ka = m.valleys()[0]
aba = ABA(m, ka, 0.75, npts=121)
vc = np.array([.5, .5, .5, .5]); vs = np.array([.5, .5, -.5, -.5])


def evalW(W, n, y1):
    G, lam, sig = metric(W)
    x1 = dim(G, OPS['O1 = cos 2(phi_a - phi_b)']); xp = dim(G, OPS['pair Psi_a Psi_b = e^{i(th_a+th_b)}'])
    xb = dim(G, OPS['boson Psi_a = e^{i th_a}'])
    xf = min(dim(G, OPS['fermion e^{i th_a + i(phi_a - phi_b)}']), dim(G, OPS['fermion e^{i th_a + i(phi_a + phi_b)}']))
    ucK, usK = vc @ W @ vc, vs @ W @ vs
    return dict(Ks=x1 / 2, Kc=1 / (2 * xp), e_pair=2 * xp, e_bos=2 * xb, e_ferm=2 * xf, ucoK=ucK, usoK=usK,
                LDP=np.pi / 2 * (ucK + usK), dpdd=1 / (np.pi * n * usK), margin=x1 - 2 - abs(y1),
                vel=sorted(np.abs(lam))[::2])


def run(n, Lr=512):
    res = analyse(n, m, verbose=False)
    vR, vL, y1 = res['vR'], res['vL'], res['y1']
    f_pp = -2 * vR * res['tau']; f_mm = -2 * vL * res['taup']
    N = int(round(n / 2 * Lr)); Kp, up = aba.KU(Lr, N)
    out = dict(n=n, y1=y1, r=res['r'], delta_s=res['delta_s'], tau=res['tau'], taup=res['taup'], Kpar=Kp, upar=up, rows=[])
    g4s = (('0', 0.0), ('K(pi)', res['g4_Kpi']), ('K(RR)', res['g4_KRR']), ('+0.4', 0.4))
    for meth in ('landau-delta', 'landau-tan', 'hybrid-ABA'):
        for lab, f4 in g4s:
            if meth == 'hybrid-ABA':
                d = 0.5 * (vR - vL); uoK, uK = up / Kp, up * Kp
                wRL = 0.5 * (uoK - uK); ws = uoK + uK
                W = W_from(0.5 * (ws + 2 * d), 0.5 * (ws - 2 * d), 2 * np.pi * wRL, f_pp, f_mm, f4)
            else:
                tn = (lambda x: x) if meth == 'landau-delta' else np.tan
                W = W_from(vR, vL, -(vR + vL) * 2 * tn(res['delta_s']), -2 * vR * tn(res['tau']), -2 * vL * tn(res['taup']), f4)
            row = evalW(W, n, y1); row.update(method=meth, g4=lab, g4val=f4)
            out['rows'].append(row)
    return out


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else 'main'
    if mode == 'scan':
        print("# n  |r|  y1  delta_s  tau  tau'  | margin range over 12 variants | K_s range")
        for n in np.arange(0.05, 0.2501, 0.0078125):
            o = run(n)
            ms = [r['margin'] for r in o['rows']]; ks = [r['Ks'] for r in o['rows']]
            print(f"{n:.4f} {o['r']:.4f} {o['y1']:.4f} {o['delta_s']:+.3f} {o['tau']:+.3f} {o['taup']:+.3f} | "
                  f"{min(ms):+.3f} .. {max(ms):+.3f} | {min(ks):.3f} .. {max(ks):.3f}")
    else:
        allout = []
        for n in (1/8, 5/32, 3/16):
            o = run(n); allout.append(o)
            print(f"\n=== n = {n:.5f}: |r|={o['r']:.4f} y1={o['y1']:.4f} delta_s={o['delta_s']:+.4f} arg t={o['tau']:+.4f} arg t'={o['taup']:+.4f} K_par(ABA)={o['Kpar']:.3f} u_par={o['upar']:.3f}")
            for r in o['rows']:
                print(f"  {r['method']:13s} g4={r['g4']:6s} K_s={r['Ks']:.3f} K_c={r['Kc']:.3f} | pair {r['e_pair']:.3f} boson {r['e_bos']:.3f} fermion {r['e_ferm']:.3f} |"
                      f" uc/Kc={r['ucoK']:.3f} us/Ks={r['usoK']:.3f} L*DP={r['LDP']:.2f} dp/dd={r['dpdd']:.2f} | margin {r['margin']:+.3f} | u={np.round(r['vel'],2)}")
            ks = [r['Ks'] for r in o['rows']]; kc = [r['Kc'] for r in o['rows']]
            print(f"  RANGE K_s {min(ks):.3f}..{max(ks):.3f}  median {np.median(ks):.3f};  K_c {min(kc):.3f}..{max(kc):.3f} median {np.median(kc):.3f}")
        json.dump(allout, open('predict.json', 'w'), indent=1, default=float)
