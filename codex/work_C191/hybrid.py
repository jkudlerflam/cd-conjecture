"""C191: hybrid estimate. Intra-valley block from the one-valley ABA (u_par, K_par; drift from the bare
Fermi velocities), inter-valley Landau parameters at first order in the exact transmission phases.
Variants of the R-R inter-valley parameter g4 as in lowE.py.  Prints K_s, K_c, exponents, y1, KT margin."""
import sys
import numpy as np
from scat191 import Model
from gauss191 import metric, dim, OPS
from lowE import analyse
from aba_valley import ABA

m = Model(w=1j)
ka = m.valleys()[0]
aba = ABA(m, ka, 0.75, npts=121)
ns = [float(eval(a)) for a in sys.argv[1:]] or [1/16, 3/32, 1/8, 5/32, 3/16]
print("# n | variant | K_par u_par | K_s K_c | exponents pair boson fermion | y1 | margin 2(K_s-1)-|y1| | u_s u_c")
for n in ns:
    res = analyse(n, m, verbose=False)
    vR, vL = res['vR'], res['vL']
    L = 512; N = int(round(n / 2 * L))
    Kp, up = aba.KU(L, N)
    d = 0.5 * (vR - vL)
    uoK, uK = up / Kp, up * Kp
    wRL = 0.5 * (uoK - uK)
    wsum = uoK + uK
    wRR, wLL = 0.5 * (wsum + 2 * d), 0.5 * (wsum - 2 * d)
    f_pp = -2 * vR * res['tau']
    f_mm = -2 * vL * res['taup']
    for lab4, f4 in (('g4=0', 0.0), ('g4=K(pi)', res['g4_Kpi']), ('g4=K(RR)', res['g4_KRR'])):
        W = np.zeros((4, 4))
        W[0, 0], W[1, 1], W[2, 2], W[3, 3] = wRR, wLL, wLL, wRR      # bR = image of aL, bL = image of aR
        W[0, 1] = W[1, 0] = wRL; W[2, 3] = W[3, 2] = wRL
        W[0, 3] = W[3, 0] = f_pp / (2 * np.pi); W[1, 2] = W[2, 1] = f_mm / (2 * np.pi)
        W[0, 2] = W[2, 0] = f4 / (2 * np.pi); W[1, 3] = W[3, 1] = f4 / (2 * np.pi)
        G, lam, sig = metric(W)
        x1 = dim(G, OPS['O1 = cos 2(phi_a - phi_b)']); xp = dim(G, OPS['pair Psi_a Psi_b = e^{i(th_a+th_b)}'])
        xb = dim(G, OPS['boson Psi_a = e^{i th_a}'])
        xf = min(dim(G, OPS['fermion e^{i th_a + i(phi_a - phi_b)}']), dim(G, OPS['fermion e^{i th_a + i(phi_a + phi_b)}']))
        Ks, Kc = x1 / 2, 1 / (2 * xp)
        # spin/charge velocities (symmetric-part estimate)
        vel = sorted(np.abs(lam))
        print(f"{n:.5f} | {lab4:9s} | {Kp:.3f} {up:.3f} | K_s={Ks:.3f} K_c={Kc:.3f} | {2*xp:.3f} {2*xb:.3f} {2*xf:.3f} | "
              f"{res['y1']:.4f} | {x1-2-abs(res['y1']):+.3f} | vel={np.round(vel,3)}")
