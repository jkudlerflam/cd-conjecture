import numpy as np
from scat191 import Model
m = Model(w=1j)
ka, kb = m.valleys()[:2]
def eig_phases(S):
    ev = np.linalg.eigvals(S)
    return ev
for n in (1/8, 5/32, 3/16):
    rho = n/2
    pL, pR, eF = m.fermi_points(ka, rho)
    # (aR, bR): a at ka+pR, b at kb-pL  (bR is U-image of aL)
    pa, pb = ka+pR, kb-pL
    K = pa+pb; E = m.eps(pa)+m.eps(pb)
    S, pin, pout, v = m.smatrix(K, E)
    print(f"n={n:.5f} K-(ka+kb)={K-(ka+kb):.4f} v_aR={m.deps(pa):.4f} v_bR={m.deps(pb):.4f} open={len(pin)} pin={np.round(pin,4)} pout={np.round(pout,4)}")
    print("  |S|=\n", np.round(np.abs(S),4), "\n  argS=\n", np.round(np.angle(S),4))
    # threshold scattering lengths at this K: eigenphases of -S at small relative momentum
    # find continuum minimum at this K
    ps = np.linspace(-np.pi, np.pi, 200001)
    Etot = m.eps(ps)+m.eps(K-ps); i = np.argmin(Etot); Emin = Etot[i]
    for dE in (1e-6, 1e-4, 1e-3):
        S2, pin2, pout2, v2 = m.smatrix(K, Emin+dE)
        if len(pin2) == 2:
            ev = np.linalg.eigvals(-S2)
            # relative momentum estimate
            q = np.sqrt(dE/ (m.d2eps(ka)))  # E_rel = q^2/m*... rough: 2 * eps''/2 q^2
            print(f"   dE={dE:g} eig(-S) phases {np.round(np.angle(ev),5)}  -> a ~ {np.round(-np.angle(ev)/(2*q),3)} (q~{q:.4g})")
