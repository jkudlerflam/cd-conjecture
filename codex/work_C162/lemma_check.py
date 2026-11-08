"""C162. Numerical check of the rearrangement lemma (PROOF section 3) on a K bond whose field is OFF the in-plane
principal axes, so that the extra classical symmetry x -> -x, n_y -> -n_y of C161's family is absent and the
tilted soliton is left-right asymmetric.  W = -h(cos psi n_x + sin psi n_y) + (kappa/2) n_z^2 + (mu/2) n_y^2,
kappa = 0.05, mu = 2.3 (C161), psi = argv[1].
Checks: (1) tilt instability at the planar onset; (2) tilted soliton: n_z single-signed, both tail coefficients
B_L, B_R > 0 and unequal; (3) two-soliton ring energies U_same(T) < U_alt(T) at every T (lemma), and U_same < 0
somewhere (like-valley binding, first-order onset)."""
import sys
import numpy as np
from scipy.optimize import brentq
from kfun import KFun, planar, tilt_lowest, ring_min

kappa, mu = 0.05, 2.3
psi = float(sys.argv[1]) if len(sys.argv) > 1 else 0.3
def mk(h):
    return KFun(-kappa, mu - kappa, h * np.cos(psi), h * np.sin(psi), 1.0)

def Epl(h):
    F = mk(h); pb, _ = F.background()
    return planar(F, pb, X=15.0, dx=0.08)[0] - 2 * np.pi
hc = brentq(Epl, 0.005, 0.2, xtol=1e-5)
F = mk(hc); pb, glob = F.background(); mu2, mw2 = F.masses(pb)
E, xs, p = planar(F, pb, X=15.0, dx=0.08)
lam0, _ = tilt_lowest(F, pb, xs, p)
print(f'psi={psi}: planar onset h_c={hc:.5f}, background phi_bg={pb:.4f} (global {glob}), m_u^2={mu2:.4f}, m_w^2={mw2:.4f}, lam0={lam0:+.4f}')

def Etilt(h, T=60.0):
    F = mk(h); pb, _ = F.background()
    E, n, dx = ring_min(F, pb, T, [0.6], dx=0.1)
    return E, n, dx, F, pb
for h in (hc, hc + 0.004, hc + 0.008):
    E1, n1, dx1, F1, pb1 = Etilt(h)
    print(f'h={h:.5f}: isolated tilted soliton E={E1:+.6f}, max n_z={n1[:,2].max():.3f}, min n_z={n1[:,2].min():+.2e}')
# choose h with small positive soliton energy
h = hc + 0.004
E1, n1, dx1, F1, pb1 = Etilt(h)
# tails of n_z on both sides of the soliton centre (ring of length 60, soliton near the middle)
x = np.arange(len(n1)) * dx1
ic = np.argmax(n1[:, 2]); m_w = np.sqrt(F1.masses(pb1)[1])
w = n1[:, 2]
def tailfit(idx, sgn):
    xx = x[idx] - x[ic]; ww = w[idx]
    return np.exp(np.polyfit(xx, np.log(ww), 1)), np.polyfit(xx, np.log(ww), 1)[0]
iR = np.arange(ic + int(8 / dx1), ic + int(14 / dx1)); iL = np.arange(ic - int(14 / dx1), ic - int(8 / dx1))
(cR, kR), (cL, kL) = tailfit(iR, 1), tailfit(iL, -1)
print(f'h={h:.5f}: E_sol={E1:+.6f}; n_z tails: right slope {kR:+.4f}, left slope {kL:+.4f} (m_w = {m_w:.4f}); '
      f'B_R={cR[1]:.4f}, B_L={cL[1]:.4f}; min n_z/max n_z = {w.min()/w.max():+.2e}')
print('T      U_same        U_alt        U_alt-U_same')
for T in (8.0, 10.0, 12.0, 15.0, 18.0, 22.0):
    Es, ns, _ = ring_min(F1, pb1, T, [0.6, 0.6], dx=0.1)
    Ea, na, _ = ring_min(F1, pb1, T, [0.6, -0.6], dx=0.1)
    alt_ok = (na[:, 2].max() > 0.05) and (na[:, 2].min() < -0.05)
    print(f'{T:5.1f}  {Es/2-E1:+.6e}  {Ea/2-E1:+.6e}  {(Ea-Es)/2:+.6e}  {"" if alt_ok else "(alternating start relaxed to same tilt)"}')
    sys.stdout.flush()
