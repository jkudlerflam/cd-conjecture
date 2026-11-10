"""C191 check: the construction applied to the attractive Gaudin-Yang gas (contact g < 0, equal populations).
Inputs exactly as for the cusp model: like phase delta_s = 0 (contact is invisible to identical fermions),
unlike channels delta_e(q) = arctan(1/(q a)), a = -2/(m g) > 0, delta_o = 0, at q = k_F = pi n / 2.
Fermion frame here (spin-1/2 fermions): transmission relative to the free value -1 is t = (e^{2i de} + 1)/2.
Expected: K_s - 1 = -|y1|/2 + O(g^2) (the SU(2) separatrix, flowing to the gap), K_c - 1 = +|y1|/2,
and the KT time l* = 1/|y1| gives Delta ~ E_F exp(-pi^2/(2|gamma|)), gamma = m g / n (exact BA asymptotics)."""
import numpy as np
from gauss191 import metric, dim, OPS, W_from
m = 1.0; n = 1.0
kF = np.pi * n / 2; v = kF / m
print("# gamma   y1     2(K_s-1)   2(K_c-1)   margin(=2(Ks-1)-|y1|)   l*=1/|y1|   pi^2/(2|gamma|)")
for gamma in (-0.05, -0.1, -0.2, -0.4):
    g = gamma * n / m
    a = -2 / (m * g)
    de = np.arctan(1 / (kF * a))
    t = (np.exp(2j * de) + 1) / 2; r = (np.exp(2j * de) - 1) / 2
    tau = np.angle(t)
    f_pp = f_mm = -2 * v * tau
    f4 = g                       # K-matrix value at q -> 0: -(2q/m) tan(de) -> -2/(m a) = g
    W = W_from(v, v, 0.0, f_pp, f_mm, f4)
    G, lam, sig = metric(W)
    Ks = dim(G, OPS['O1 = cos 2(phi_a - phi_b)']) / 2
    Kc = 1 / (2 * dim(G, OPS['pair Psi_a Psi_b = e^{i(th_a+th_b)}']))
    y1 = 2 * v * abs(r) / (np.pi * v)
    print(f"{gamma:+.3f}  {y1:.5f}  {2*(Ks-1):+.5f}  {2*(Kc-1):+.5f}  {2*(Ks-1)-abs(y1):+.6f}  {1/abs(y1):8.3f}  {np.pi**2/(2*abs(gamma)):8.3f}")
