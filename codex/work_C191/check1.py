import numpy as np
from scat191 import Model, unlike_amplitudes, like_amplitude
m = Model(w=1j)
ka, kb = m.valleys()[:2]
print("valleys", ka, kb, ka+kb, m.eps(ka), m.eps(kb), "eps''", m.d2eps(ka))
# small-q check against C175 eigenphase table: a_pm(q) from eigenphases of -S at E = 2 eps(kb+q)
for q in (0.01, 0.07, 0.13, 0.19, 0.25, 0.31):
    S, pin, pout, v = m.smatrix(ka+kb, 2*m.eps(kb+q))
    ev = np.linalg.eigvals(-S)
    print(q, "unitarity", np.linalg.norm(S.conj().T@S-np.eye(2)), "a(q) =", np.sort(-np.angle(ev)/(2*q)))
# hard rods: NN hop only, no V: phase -q*a with a=1
h = Model(w=1.0, nnn=False, V2=0, V3=0)
kh = h.valleys()[0]
for rho in (0.05, 0.1, 0.2):
    pL, pR, eF = h.fermi_points(kh, rho)
    L = like_amplitude(h, kh, pL, pR)
    q = (pR-pL)/2
    print("hard rods rho", rho, "delta", L['delta'], "-q", -q, "nopen", L['nopen'])
