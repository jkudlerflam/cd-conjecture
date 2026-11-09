"""C181: atom-dimer threshold data from the three-body box (gap coordinates, every gap <= R, nbody_box.py of C175).
At total momentum K = K_pair + k_a the lowest levels are the aab trimer (if bound) and atom-dimer scattering states.
For a scattering level E_i(R), the hard-wall relation E_i - E_th = pi^2/(2 mu (R - a_i)^2), 1/mu = 1/m_a + 1/M_pair,
gives an effective length a_i (eigenvalues of the 2x2 atom-dimer length matrix on the channels atom-left, atom-right,
up to a common offset of order the dimer size from the wall convention). Usage: python3 ad_scat.py wa R1,R2 phi1,phi2"""
import numpy as np, sys
from scat import valleys
from fewbody import pair_min, epp
import nbody_box
wa = float(sys.argv[1]); Rs = [int(x) for x in sys.argv[2].split(',')]; phis = [float(x) for x in sys.argv[3].split(',')]
print("# |w| phi R E_th(pair+a) levels-E_th (4 lowest) | a_i from the scattering levels", flush=True)
for phi in phis:
    w = wa * np.exp(1j * phi)
    (ka, ea), (kb, eb) = valleys(w)[:2]
    Kp, Ep, Mp = pair_min(w, ka + kb)
    mu = 1.0 / (epp(ka, w) + Mp)
    Eth = Ep + ea
    for R in Rs:
        e = nbody_box.lowest(3, Kp + ka, w, R, k=4) - Eth
        sc = e[e > 0]
        a = [R - np.pi / np.sqrt(2 * mu * x) for x in sc[:2]]
        print(f"{wa:.1f} {phi:.4f} {R} {Eth:.6f} " + " ".join(f"{x:+.6f}" for x in e) + " | " + " ".join(f"{x:.2f}" for x in a), flush=True)
