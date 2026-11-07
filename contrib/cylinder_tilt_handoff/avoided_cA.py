# d = 4 quasi-free chain: band eps1 = -2t cos k crossing a flat level eps2 = 0 at k_c = pi/2, on-site hybridization D.
# Fermi level at the kink, mu = -D. IR: one Dirac fermion with v = t (A = t). Ring: APBC momenta, L = 0 mod 4
# (Fermi points midway between grid points, as for XX). c_A(L) = (6L/pi)(e0 L - E0(L)).
import numpy as np
t = 1.0
def bands(k, D):
    e1 = -2*t*np.cos(k); s = np.sqrt(e1**2/4 + D**2)
    return e1/2 - s, e1/2 + s
def E0(k, D):
    lo, hi = bands(k, D); mu = -D
    return np.minimum(lo - mu, 0) + np.minimum(hi - mu, 0)
for D in [1e-2, 1e-3]:
    kk = np.linspace(-np.pi, np.pi, 2_000_001)[:-1]; e0 = E0(kk, D).mean()
    out = []
    for L in [8, 16, 32, 64, 128, 256, 512, 1024, 4096, 16384]:
        k = 2*np.pi*(np.arange(L) + 0.5)/L
        out.append((L, 6*L*(e0*L - E0(k, D).sum())/np.pi))
    print(f"D={D}: crossover t/D={t/D:.0f};  c_A(L)/t:", "  ".join(f"{L}:{c:.3f}" for L, c in out))
