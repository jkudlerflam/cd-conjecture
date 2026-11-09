"""C182: pair bands at phi = pi/2 against |w|: minima of the two-body lower edge, curvature at K = pi,
binding of ab (K = pi) and like (K = 2k_a) pairs."""
import numpy as np
from pairband import analyse, lowest, econt, valleys
for wa in [0.3,0.5,0.75,1.0,1.5,2.0,2.5,3.0,3.25,3.5,3.6,3.7,3.8,3.9]:
    R = 4000 if wa < 0.6 else 2000
    res = analyse(wa, np.pi/2, R=R, nK=361)
    w = 1j*wa; v = valleys(w); ka = v[0][0]; ea=v[0][1]
    # like pair at K=2ka
    Kl = 2*ka
    el = lowest(Kl, w, R, k=1)[0]; ecl = econt(Kl, w)
    # pair curvature at pi
    h=2e-3; f=lambda K: lowest(K,w,R,k=1)[0]
    curv = (f(np.pi+h)-2*f(np.pi)+f(np.pi-h))/h**2
    print(f"|w|={wa:4.2f} valleys k={v[0][0]:+.4f},{v[1][0]:+.4f} eps_min={ea:.5f}  E_b(ab,K=pi)={econt(np.pi,w)-f(np.pi):.6f}  1/m_pair(pi)={curv:.5f} free={(8-wa**2/2)/2:.4f}  like(K=2ka): E-Econt={el-ecl:+.2e}", flush=True)
    for m in res['mins']:
        print(f"    min K={m['K']:+.4f} E-2eps={m['E_minus_2ea']:+.6f} depth={m['depth']:+.2e} dR={m['E_2R']-m['E']:+.1e} curv={m['curv']:.4f}", flush=True)
