"""Closed form for D_2 on the XXZ critical line, from the gradient identity
    L^2 (rho_2^inf - rho_2^L) -> (pi/6) grad_h A   (A = c v, fixed charge sector),
and the Bethe-ansatz velocity v(Delta) = pi sin(gamma) / (2 gamma), Delta = cos(gamma).

Bond h = S^x S^x + S^y S^y + Delta S^z S^z, d = 2, c = 1, A = v.
By U(1) x Z_2 x reflection symmetry the two-site marginal is
    rho_2 = 1/4 + 4 g_xx (S^xS^x + S^yS^y) + 4 g_zz S^zS^z,
so rho_2^inf - rho_2^L is fixed by two numbers.  Hellmann-Feynman gives
    L^2 (g_zz^inf - g_zz^L) -> (pi/6) A'(Delta)                  (pairing with S^zS^z)
    L^2 (2 g_xx^inf - 2 g_xx^L) -> (pi/6) (A - Delta A')        (pairing with the xy part, Euler)
and the eigenvalues of X = rho_2^inf - rho_2^L are a, a, b - a, -b - a (times 1/L^2) with
    a = (pi/6) A',   b = (pi/6)(A - Delta A'),
so  D_2(Delta) = 2|a| + |b - a| + |b + a| = 2|a| + 2 max(|a|, |b|).

Also printed: the Yang-Yang ground-state energy density e_0 and its Delta-derivative
(the exact infinite-chain <S^zS^z>), used by xxz_ed_check.py, and the bond half-spread.
"""
import numpy as np
from scipy.integrate import quad

def v_of(Delta):
    g = np.arccos(Delta)
    return np.pi*np.sin(g)/(2*g)

def vprime_of(Delta):
    # dv/dDelta = (dv/dgamma)/(dDelta/dgamma),  dDelta/dgamma = -sin(gamma)
    g = np.arccos(Delta)
    dv_dg = (np.pi/2)*(g*np.cos(g) - np.sin(g))/g**2
    return dv_dg/(-np.sin(g))

def e0_of(Delta):
    """Yang-Yang: e_0 = Delta/4 - (sin g/2) int dx sinh((pi-g)x)/(sinh(pi x) cosh(g x)), Delta=cos g."""
    g = np.arccos(Delta)
    a = np.pi - g
    # same integrand written with decaying exponentials: 2 e^{-2gx}(1-e^{-2ax}) / ((1-e^{-2 pi x})(1+e^{-2gx}))
    def f(x):
        if x < 1e-8: return a/np.pi
        return 2*np.exp(-2*g*x)*(1-np.exp(-2*a*x))/((1-np.exp(-2*np.pi*x))*(1+np.exp(-2*g*x)))
    val, err = quad(f, 0, np.inf, limit=400)
    return Delta/4 - (np.sin(g)/2)*2*val

def gzz_inf(Delta, h=1e-5):
    return (e0_of(Delta+h) - e0_of(Delta-h))/(2*h)

def gxx_inf(Delta):
    # e_0 = 2 g_xx + Delta g_zz
    return (e0_of(Delta) - Delta*gzz_inf(Delta))/2

def D2_closed(Delta):
    A = v_of(Delta); Ap = vprime_of(Delta)
    a = np.pi/6*Ap; b = np.pi/6*(A - Delta*Ap)
    return 2*abs(a) + 2*max(abs(a), abs(b)), a, b

def M_half_spread(Delta):
    ev = np.array([Delta/4, Delta/4, 0.5 - Delta/4, -0.5 - Delta/4])
    return (ev.max() - ev.min())/2

if __name__ == "__main__":
    print("check e_0(0) = -1/pi:", e0_of(0.0), -1/np.pi)
    print("check e_0(1/2) = -3/8 (Razumov-Stroganov):", e0_of(0.5))
    print("check e_0(1) -> 1/4 - ln 2:", e0_of(0.999999), 0.25 - np.log(2))
    print()
    print(" Delta      v=A       A'(Delta)    e_0        g_zz^inf    g_xx^inf     M     D_2 closed   (pi/6)A/M   (pi^2/6)log2 d   (6/pi) M D_2")
    for Delta in [0.999, 0.9, 0.5, 0.0, -0.5, -0.8, -0.9, -0.95, -0.99, -0.999, -0.9999]:
        A = v_of(Delta); Ap = vprime_of(Delta); D2, a, b = D2_closed(Delta); M = M_half_spread(Delta)
        print(f"{Delta:8.4f} {A:9.6f} {Ap:11.6f} {e0_of(Delta):10.6f} {gzz_inf(Delta):10.6f} {gxx_inf(Delta):10.6f} {M:6.3f} {D2:11.5f} {np.pi/6*A/M:10.5f} {np.pi**2/6:12.5f} {6/np.pi*M*D2:12.4f}")
    print()
    print("Asymptotics near the ferromagnetic point, eps = 1 + Delta:")
    for eps in [1e-2, 1e-3, 1e-4, 1e-5, 1e-6]:
        Delta = -1 + eps
        D2, a, b = D2_closed(Delta)
        print(f"  eps={eps:8.1e}  A={v_of(Delta):10.3e}  A'={vprime_of(Delta):10.4f}  D_2={D2:10.4f}   D_2*sqrt(eps)={D2*np.sqrt(eps):8.5f}  (limit pi sqrt(2)/6 = {np.pi*np.sqrt(2)/6:.5f})")
    print()
    print("Where the closed form crosses the Heisenberg value pi^2/6 = 1.644934 (D_2 is minimal at Delta=1):")
    Ds = np.linspace(-0.999, 0.999, 4000)
    vals = np.array([D2_closed(D)[0] for D in Ds])
    print("  min over the line:", vals.min(), "at Delta =", Ds[vals.argmin()])
    print("  D_2 at Delta=0 (XX):", D2_closed(0.0)[0], " expected (pi+2)/3 =", (np.pi+2)/3)
    print("  D_2 at Delta->1 (Heisenberg):", D2_closed(0.999999)[0], " expected pi^2/6 =", np.pi**2/6)
