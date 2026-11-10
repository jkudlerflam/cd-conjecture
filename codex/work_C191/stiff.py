"""C191: spin and charge stiffness u K from the hybrid W (central variant), and the open-chain parity gap under
the two boundary conditions: phi_s pinned (intra-valley wall reflection) and theta_s pinned (valley-swapping wall)."""
import numpy as np
from scat191 import Model
from gauss191 import W_from
from lowE import analyse
from aba_valley import ABA
m = Model(w=1j); ka = m.valleys()[0]; aba = ABA(m, ka, 0.75, npts=121)
vc = np.array([.5,.5,.5,.5]); vs = np.array([.5,.5,-.5,-.5]); jc = np.array([.5,-.5,.5,-.5]); js = np.array([.5,-.5,-.5,.5])
for n in (1/8, 5/32, 3/16):
    res = analyse(n, m, verbose=False); vR, vL = res['vR'], res['vL']
    Kp, up = aba.KU(512, int(round(n/2*512)))
    d = 0.5*(vR-vL); uoK, uK = up/Kp, up*Kp; wRL = 0.5*(uoK-uK); ws = uoK+uK
    W = W_from(0.5*(ws+2*d), 0.5*(ws-2*d), 2*np.pi*wRL, -2*vR*res['tau'], -2*vL*res['taup'], res['g4_Kpi'])
    a = {k: v@W@v for k, v in (('uc/Kc', vc), ('us/Ks', vs), ('ucKc', jc), ('usKs', js))}
    LDP_phi = np.pi/2*(a['uc/Kc'] + a['us/Ks']); LDP_th = np.pi/2*(a['uc/Kc'] + a['usKs'])
    print(f"n={n:.5f} " + " ".join(f"{k}={v:.3f}" for k, v in a.items()) +
          f" | K_s(sym)={np.sqrt(a['usKs']/a['us/Ks']):.3f} K_c(sym)={np.sqrt(a['ucKc']/a['uc/Kc']):.3f}"
          f" | L*Delta_P: phi_s pinned {LDP_phi:.2f}, theta_s pinned {LDP_th:.2f}")
